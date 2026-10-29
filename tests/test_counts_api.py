from __future__ import annotations

import asyncio
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import httpx


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "web") not in sys.path:
    sys.path.insert(0, str(ROOT / "web"))

from app.main import create_app
from app.services import analyzer, count_analysis, job_manager, storage


BONDS = ("C-C", "C-H", "C-O", "O-H", "H-H", "O-O")
GROUPS = tuple(count_analysis.GROUP_LABELS)


class CountApiTests(unittest.IsolatedAsyncioTestCase):
    async def asyncSetUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.patches = [patch.object(job_manager, "JOB_ROOT", self.root / "jobs"),
                        patch.object(analyzer, "JOB_ROOT", self.root / "jobs"),
                        patch.object(storage, "UPLOAD_ROOT", self.root / "uploads")]
        for active in self.patches:
            active.start()
        job_manager._singleton = None
        app = create_app()
        self.manager = job_manager.get_job_manager()
        self.client = httpx.AsyncClient(transport=httpx.ASGITransport(app=app,
                                        client=("127.0.0.1", 1234)), base_url="http://127.0.0.1")

    async def asyncTearDown(self):
        await self.client.aclose()
        job_manager._singleton = None
        for active in reversed(self.patches):
            active.stop()
        self.temporary.cleanup()

    async def make_job(self, *, files=None, frames=101, config=None, job_id="count-test-001"):
        job = await self.manager.create(job_id, "unused", "source.bonds", {"_raw": "1:H,2:C,3:O"},
                                        config or {"timestep_scale": 0.0002, "expected_interval": 1000})
        await self.manager.update(job.job_id, status=job_manager.JobStatus.COMPLETED,
                                  summary={"frames": frames}, result_files=files or {})
        return job.job_id

    def make_outputs(self, frames=101):
        target = self.root / "counts"
        target.mkdir()
        files = {}
        series = {}
        for name in ("bond_counts", "functional_groups", "functional_group_molecules"):
            categories = BONDS if name == "bond_counts" else GROUPS
            rows = []
            for i in range(frames):
                values = {key: 0 for key in categories}
                if name == "bond_counts":
                    values.update({"C-C": 5000 if i == 49 else i,
                                   "C-H": 9000 if i == 77 else 2})
                else:
                    values["hydroxyl"] = 2 if name == "functional_groups" else 1
                rows.append({"timestep": i * 1000, "time_ps": i / 5, **values})
            series[name] = rows
            path = target / (name + ".csv")
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["timestep", "time_ps", *categories])
                writer.writeheader(); writer.writerows(rows)
            files[name] = str(path)
        changes = target / "bond_changes.csv"
        columns = [key + suffix for key in BONDS for suffix in ("_formed", "_broken")]
        with changes.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["timestep", "time_ps", *columns])
            writer.writeheader()
            for i in range(frames):
                writer.writerow({"timestep": i * 1000, "time_ps": i / 5,
                                 **{key: None if i == 0 else 1 for key in columns}})
        files["bond_changes"] = str(changes)
        metadata = {"status": "completed", "frame_count": frames,
                    "first_timestep": 0, "last_timestep": (frames - 1) * 1000,
                    "time_start_ps": 0.0, "time_end_ps": (frames - 1) / 5,
                    "timestep_scale_ps": 0.0002, "rule_version": "test-topology",
                    "bond_categories": [{"key": key, "label": key.replace("-", "–")} for key in BONDS],
                    "functional_group_categories": [{"key": key, "label": count_analysis.GROUP_LABELS[key]} for key in GROUPS]}
        path = target / "count_analysis.json"
        path.write_text(json.dumps(metadata), encoding="utf-8")
        files["count_metadata"] = str(path)
        return files

    async def test_all_frames_modes_labels_and_initial_changes_are_not_zero_filled(self):
        files = self.make_outputs(3)
        job_id = await self.make_job(files=files, frames=3)
        base = f"/api/jobs/{job_id}/counts"
        response = await self.client.get(base)
        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["frame_count"], 3)
        self.assertEqual(payload["returned_points"], 3)
        self.assertFalse(payload["downsampled"])
        self.assertEqual(payload["rows"][-1]["time_ps"], 0.4)
        self.assertEqual(payload["categories"][0], {"key": "C-C", "label": "C–C"})
        self.assertTrue(payload["csv_url"].endswith("/outputs/bond_counts"))
        for metric in ("formed", "broken"):
            changed = (await self.client.get(base, params={"metric": metric})).json()
            self.assertIsNone(changed["rows"][0]["C-C"])
            self.assertEqual(changed["rows"][1]["C-C"], 1)
        group_count = (await self.client.get(base, params={"kind": "functional"})).json()
        carriers = (await self.client.get(base, params={"kind": "functional", "metric": "molecules"})).json()
        self.assertEqual(group_count["rows"][0]["hydroxyl"], 2)
        self.assertEqual(carriers["rows"][0]["hydroxyl"], 1)

    async def test_sampling_covers_whole_range_and_preserves_each_curve_extrema(self):
        files = self.make_outputs()
        job_id = await self.make_job(files=files)
        response = await self.client.get(f"/api/jobs/{job_id}/counts", params={"max_points": 50})
        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertTrue(payload["downsampled"])
        self.assertLessEqual(len(payload["rows"]), 50)
        steps = {row["timestep"] for row in payload["rows"]}
        self.assertTrue({0, 49000, 77000, 100000}.issubset(steps))
        self.assertEqual(max(row["C-C"] for row in payload["rows"]), 5000)
        self.assertEqual(max(row["C-H"] for row in payload["rows"]), 9000)

    async def test_missing_running_and_failed_are_explicit_without_empty_zero_series(self):
        job_id = await self.make_job(frames=2)
        url = f"/api/jobs/{job_id}/counts"
        missing = (await self.client.get(url)).json()
        self.assertEqual(missing["status"], "missing")
        self.assertEqual(missing["rows"], [])
        for status in ("running", "failed"):
            await self.manager.update(job_id, count_analysis={"status": status, "error": "evidence"})
            result = (await self.client.get(url)).json()
            self.assertEqual(result["status"], status)
            self.assertEqual(result["rows"], [])
            self.assertEqual(self.manager.snapshot(job_id)["status"], "completed")

    async def test_invalid_combinations_and_corrupt_series_are_rejected(self):
        files = self.make_outputs(3)
        job_id = await self.make_job(files=files, frames=3)
        url = f"/api/jobs/{job_id}/counts"
        for parameters in ({"kind": "functional", "metric": "formed"},
                           {"kind": "bonds", "metric": "molecules"},
                           {"kind": "other"}, {"max_points": 0}):
            self.assertEqual((await self.client.get(url, params=parameters)).status_code, 422)
        path = Path(files["bond_counts"])
        path.write_text(path.read_text().replace("0,0.0,0,2", "0,0.0,-1,2"), encoding="utf-8")
        self.assertEqual((await self.client.get(url)).status_code, 409)
        self.assertEqual((await self.client.get("/api/jobs/count-missing-001/counts")).status_code, 404)

    async def test_claim_is_idempotent_under_concurrent_requests_and_no_unknown_scope(self):
        source = self.root / "input.bonds"
        source.write_text("# source", encoding="utf-8")
        job_id = await self.make_job(frames=2, config={"source_path": str(source), "timestep_scale": 0.0002})
        left, right = await asyncio.gather(count_analysis.prepare_backfill(self.manager, job_id),
                                           count_analysis.prepare_backfill(self.manager, job_id))
        self.assertEqual(sum(parameters is not None for _, parameters in (left, right)), 1)
        claimed = next(parameters for _, parameters in (left, right) if parameters is not None)
        self.assertEqual(claimed[4], 2)
        await self.manager.update(job_id, count_analysis={})
        seen = []
        async def no_scan(manager, key, parameters):
            seen.append(parameters[4])
        with patch.object(count_analysis, "run_backfill", no_scan):
            response = await self.client.post(f"/api/jobs/{job_id}/counts")
            again = await self.client.post(f"/api/jobs/{job_id}/counts")
        self.assertEqual(response.status_code, 202)
        self.assertEqual(again.status_code, 202)
        self.assertEqual(seen, [2])
        unknown = await self.make_job(frames=0, config={"source_path": str(source)}, job_id="count-unknown-001")
        self.assertEqual((await self.client.post(f"/api/jobs/{unknown}/counts")).status_code, 409)

    async def test_corrupt_monotone_time_axis_and_metadata_range_are_rejected(self):
        files = self.make_outputs(3)
        job_id = await self.make_job(files=files, frames=3)
        url = f"/api/jobs/{job_id}/counts"
        csv_path = Path(files["bond_counts"])
        original_csv = csv_path.read_text(encoding="utf-8")
        csv_path.write_text(original_csv.replace("1000,0.2,", "1000,0.21,"), encoding="utf-8")
        self.assertEqual((await self.client.get(url)).status_code, 409)
        csv_path.write_text(original_csv, encoding="utf-8")
        metadata_path = Path(files["count_metadata"])
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        metadata["time_end_ps"] = 0.41
        metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
        self.assertEqual((await self.client.get(url)).status_code, 409)

    async def test_backfill_caps_existing_job_scope_and_preserves_reaction_outputs(self):
        source = self.root / "input.bonds"
        source.write_text("\n".join(f"# Timestep {i * 1000}\n# Number of particles 1\n1 1 0 0 0 0 0" for i in range(3)), encoding="utf-8")
        stat = source.stat()
        marker = self.root / "events.csv"
        marker.write_text("original reaction evidence", encoding="utf-8")
        job_id = await self.make_job(frames=2, files={"reaction_events": str(marker)}, config={
            "source_path": str(source), "source_size_bytes": stat.st_size,
            "source_mtime_ns": stat.st_mtime_ns, "timestep_scale": 0.0002, "expected_interval": 1000})
        response = await self.client.post(f"/api/jobs/{job_id}/counts")
        self.assertEqual(response.status_code, 202, response.text)
        completed = self.manager.snapshot(job_id)
        self.assertEqual(completed["status"], "completed")
        self.assertEqual(completed["count_analysis"]["status"], "completed", completed)
        self.assertEqual(marker.read_text(), "original reaction evidence")
        counts = (await self.client.get(f"/api/jobs/{job_id}/counts")).json()
        self.assertEqual([row["timestep"] for row in counts["rows"]], [0, 1000])
        self.assertEqual(counts["frame_count"], 2)
        with patch.object(count_analysis, "run_backfill", side_effect=AssertionError("must not rerun")):
            duplicate = await self.client.post(f"/api/jobs/{job_id}/counts")
        self.assertEqual(duplicate.json()["status"], "completed")

    async def test_source_drift_failed_worker_and_restart_keep_primary_result(self):
        source = self.root / "input.bonds"
        source.write_text("invalid source", encoding="utf-8")
        stat = source.stat()
        job_id = await self.make_job(frames=2, config={"source_path": str(source), "source_size_bytes": stat.st_size})
        response = await self.client.post(f"/api/jobs/{job_id}/counts")
        self.assertEqual(response.status_code, 202)
        self.assertEqual(self.manager.snapshot(job_id)["count_analysis"]["status"], "failed")
        self.assertEqual(self.manager.snapshot(job_id)["status"], "completed")
        source.write_text("changed", encoding="utf-8")
        self.assertEqual((await self.client.post(f"/api/jobs/{job_id}/counts")).status_code, 409)
        await self.manager.update(job_id, count_analysis={"status": "running"})
        restored = job_manager.JobManager().snapshot(job_id)
        self.assertEqual(restored["status"], "completed")
        self.assertEqual(restored["count_analysis"]["status"], "failed")


if __name__ == "__main__":
    unittest.main()
