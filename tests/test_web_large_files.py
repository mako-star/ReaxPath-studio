from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import httpx


ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = ROOT / "web"
if str(WEB_ROOT) not in sys.path:
    sys.path.insert(0, str(WEB_ROOT))

from app.main import create_app  # noqa: E402
from app.services import analyzer, job_manager, storage  # noqa: E402


class ChunkedStorageTests(unittest.TestCase):
    def test_chunks_are_written_in_place_and_can_arrive_out_of_order(self) -> None:
        payload = b"0123456789abcdefghijklmnopqrstuvwxyz"
        with tempfile.TemporaryDirectory() as temporary:
            upload_root = Path(temporary) / "uploads"
            upload_root.mkdir()
            with patch.object(storage, "UPLOAD_ROOT", upload_root):
                meta = storage.create_upload(
                    "upload-test-001",
                    "trajectory.bonds",
                    len(payload),
                    chunk_size=10,
                )
                self.assertEqual(meta.total_chunks, 4)
                storage.receive_chunk("upload-test-001", 2, payload[20:30])
                storage.receive_chunk("upload-test-001", 0, payload[0:10])
                repeated = storage.receive_chunk("upload-test-001", 0, payload[0:10])
                self.assertEqual(repeated.received_chunks, [0, 2])
                storage.receive_chunk("upload-test-001", 3, payload[30:])
                storage.receive_chunk("upload-test-001", 1, payload[10:20])
                completed = storage.finalize_upload("upload-test-001")

                merged = Path(completed.merged_path or "")
                self.assertEqual(merged.read_bytes(), payload)
                self.assertFalse((merged.parent / "merged.part").exists())
                self.assertFalse((merged.parent / "chunks").exists())


class WebApiTests(unittest.IsolatedAsyncioTestCase):
    async def asyncSetUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.runtime_root = Path(self.temporary.name)
        self.job_root = self.runtime_root / "jobs"
        self.upload_root = self.runtime_root / "uploads"
        self.job_root.mkdir()
        self.upload_root.mkdir()

        self.patches = [
            patch.object(job_manager, "JOB_ROOT", self.job_root),
            patch.object(analyzer, "JOB_ROOT", self.job_root),
            patch.object(storage, "UPLOAD_ROOT", self.upload_root),
        ]
        for active_patch in self.patches:
            active_patch.start()
        job_manager._singleton = None

        app = create_app()
        transport = httpx.ASGITransport(app=app, client=("127.0.0.1", 43210))
        self.client = httpx.AsyncClient(
            transport=transport,
            base_url="http://127.0.0.1",
        )

    async def asyncTearDown(self) -> None:
        await self.client.aclose()
        job_manager._singleton = None
        for active_patch in reversed(self.patches):
            active_patch.stop()
        self.temporary.cleanup()

    async def test_removed_connection_options_are_rejected(self) -> None:
        for route, source in (
            ('/api/jobs/local', {'bonds_path': str(ROOT / 'examples' / 'ch4_dissociation.bonds')}),
            ('/api/jobs', {'upload_id': 'unused'}),
        ):
            response = await self.client.post(route, params={
                **source, 'element_mapping': '1:C,2:H', 'bond_cutoff': 0.5,
            })
            self.assertEqual(response.status_code, 422)
            self.assertIn('removed', response.text)
        with self.assertRaisesRegex(ValueError, 'removed'):
            analyzer.build_analyzer_config({'bond_cutoff': 0.5})

    async def test_local_path_job_runs_streaming_analysis_and_previews_results(self) -> None:
        source = ROOT / "examples" / "ch4_dissociation.bonds"
        response = await self.client.post(
            "/api/jobs/local",
            params={
                "bonds_path": str(source),
                "element_mapping": "1:C,2:H",
                "timestep_scale": 0.001,
                "recrossing_window": 0.0,
                "generate_identifiers": "true",
            },
        )
        self.assertEqual(response.status_code, 201, response.text)
        job = response.json()
        snapshot = (await self.client.get(f"/api/jobs/{job['job_id']}")).json()

        self.assertEqual(snapshot["status"], "completed")
        self.assertEqual(snapshot["source_mode"], "local")
        self.assertTrue(snapshot["summary"]["streaming_analysis"])
        self.assertEqual(snapshot["summary"]["frames"], 2)
        self.assertEqual(snapshot["summary"]["events_accepted"], 1)

        registry = await self.client.get(
            f"/api/jobs/{job['job_id']}/preview/molecule_registry",
            params={"limit": 1},
        )
        self.assertEqual(registry.status_code, 200, registry.text)
        self.assertEqual(len(registry.json()["rows"]), 1)
        self.assertTrue(registry.json()["possibly_truncated"])

        network = await self.client.get(
            f"/api/jobs/{job['job_id']}/outputs/reaction_network"
        )
        self.assertEqual(network.status_code, 200, network.text)
        self.assertEqual(network.json()["metadata"]["frame_count"], 2)

    async def test_raw_chunk_upload_round_trip(self) -> None:
        payload = (ROOT / "examples" / "ch4_dissociation.bonds").read_bytes()
        chunk_size = max(1, len(payload) // 3)
        initialized = await self.client.post(
            "/api/uploads",
            params={
                "filename": "small.bonds",
                "total_size": len(payload),
                "chunk_size": chunk_size,
            },
        )
        self.assertEqual(initialized.status_code, 201, initialized.text)
        upload = initialized.json()

        for index in reversed(range(upload["total_chunks"])):
            start = index * chunk_size
            end = min(start + chunk_size, len(payload))
            response = await self.client.put(
                f"/api/uploads/{upload['upload_id']}/chunks/{index}",
                content=payload[start:end],
                headers={"Content-Type": "application/octet-stream"},
            )
            self.assertEqual(response.status_code, 200, response.text)

        completed = await self.client.post(
            f"/api/uploads/{upload['upload_id']}/complete"
        )
        self.assertEqual(completed.status_code, 200, completed.text)
        merged = Path(completed.json()["merged_path"])
        self.assertEqual(merged.read_bytes(), payload)


if __name__ == "__main__":
    unittest.main()
