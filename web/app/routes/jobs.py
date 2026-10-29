"""Job lifecycle routes: create from a completed upload, poll status, download."""

from __future__ import annotations

import csv
import json
import secrets
from pathlib import Path as FilePath
from typing import Annotated

from fastapi import APIRouter, BackgroundTasks, HTTPException, Path as ApiPath, Query, Request
from fastapi.responses import FileResponse, Response

from ..services import event_queries, storage
from ..services.analyzer import cleanup_job_outputs, run_analysis
from ..services.job_manager import JobStatus, get_job_manager

router = APIRouter(prefix="/api/jobs", tags=["jobs"])


def reject_connection_filters(request: Request) -> None:
    removed = {"bond_cutoff", "valence_overrides", "coordination_overrides"} & request.query_params.keys()
    if removed:
        raise HTTPException(422, f"Connection filtering options were removed: {sorted(removed)}")


@router.post("", status_code=201)
async def create_job(
    request: Request,
    background_tasks: BackgroundTasks,
    upload_id: Annotated[str, "ID returned by the upload completion endpoint"],
    element_mapping: Annotated[str, 'LAMMPS atom-type mapping, e.g. "1:C,2:H,3:O"'],
    timestep_scale: Annotated[float, Query(gt=0, allow_inf_nan=False)] = 1.0,
    recrossing_window: Annotated[float, Query(ge=0, allow_inf_nan=False)] = 0.0,
    expected_interval: Annotated[int | None, Query(gt=0)] = None,
    generate_identifiers: Annotated[bool, "Generate SMILES/InChI (requires RDKit)"] = True,
) -> dict:
    reject_connection_filters(request)
    try:
        upload_meta = storage.get_upload(upload_id)
    except storage.UploadError as exc:
        raise HTTPException(404, f"Upload not found: {exc}") from exc
    if upload_meta.completed_at is None or upload_meta.merged_path is None:
        raise HTTPException(409, "Upload has not been finalised yet")
    if not FilePath(upload_meta.merged_path).exists():
        raise HTTPException(409, "Merged upload file is missing on disk")

    job_id = secrets.token_urlsafe(12)
    config_payload = {
        "timestep_scale": timestep_scale,
        "recrossing_window": recrossing_window,
        "generate_identifiers": generate_identifiers,
        "expected_interval": expected_interval,
    }
    manager = get_job_manager()
    job = await manager.create(
        job_id=job_id,
        upload_id=upload_id,
        filename=upload_meta.filename,
        element_mapping={},  # filled from string by the analyzer via parse_element_mapping
        config=config_payload,
        source_mode="upload",
    )
    # Stash the raw mapping string where the analyzer can find it.
    job.element_mapping = {"_raw": element_mapping}
    manager._persist(job)  # noqa: SLF001 - intentional, the raw mapping lives in the meta

    background_tasks.add_task(
        run_analysis,
        manager,
        job_id,
        FilePath(upload_meta.merged_path),
        element_mapping,
        config_payload,
    )

    return job.to_dict()


@router.post("/local", status_code=201)
async def create_local_job(
    request: Request,
    background_tasks: BackgroundTasks,
    bonds_path: Annotated[str, "Absolute path to a bonds file on the API host"],
    element_mapping: Annotated[str, 'LAMMPS atom-type mapping, e.g. "1:C,2:H,3:O"'],
    timestep_scale: Annotated[float, Query(gt=0, allow_inf_nan=False)] = 1.0,
    recrossing_window: Annotated[float, Query(ge=0, allow_inf_nan=False)] = 0.0,
    expected_interval: Annotated[int | None, Query(gt=0)] = None,
    generate_identifiers: Annotated[bool, "Generate SMILES/InChI"] = True,
) -> dict:
    reject_connection_filters(request)
    client_host = request.client.host if request.client else ""
    if client_host not in {"127.0.0.1", "::1", "localhost", "testclient"}:
        raise HTTPException(403, "Local-path analysis is only available from this machine")

    source = FilePath(bonds_path.strip().strip('"')).expanduser().resolve()
    if not source.is_file():
        raise HTTPException(404, f"Bonds file not found: {source}")

    config_payload = {
        "timestep_scale": timestep_scale,
        "recrossing_window": recrossing_window,
        "generate_identifiers": generate_identifiers,
        "expected_interval": expected_interval,
    }
    job_id = secrets.token_urlsafe(12)
    manager = get_job_manager()
    job = await manager.create(
        job_id=job_id,
        upload_id="local-path",
        filename=source.name,
        element_mapping={},
        config=config_payload,
        source_mode="local",
    )
    job.element_mapping = {"_raw": element_mapping}
    manager._persist(job)  # noqa: SLF001
    background_tasks.add_task(
        run_analysis,
        manager,
        job_id,
        source,
        element_mapping,
        config_payload,
    )
    return job.to_dict()


@router.get("/{job_id}")
async def get_job(
    job_id: Annotated[str, ApiPath(min_length=8, max_length=48)],
) -> dict:
    manager = get_job_manager()
    snapshot = manager.snapshot(job_id)
    if snapshot is None:
        raise HTTPException(404, f"Job {job_id!r} not found")
    return snapshot


@router.get("")
async def list_jobs() -> dict:
    manager = get_job_manager()
    return {"jobs": [manager.snapshot(jid) for jid in manager._jobs.keys()]}  # noqa: SLF001


@router.get("/{job_id}/outputs/{name}")
async def download_output(
    job_id: Annotated[str, ApiPath(min_length=8, max_length=48)],
    name: Annotated[str, ApiPath(pattern=r"^[A-Za-z0-9_.-]+$")],
) -> FileResponse:
    manager = get_job_manager()
    snapshot = manager.snapshot(job_id)
    if snapshot is None:
        raise HTTPException(404, f"Job {job_id!r} not found")
    if snapshot.get("status") != JobStatus.COMPLETED.value:
        raise HTTPException(409, "Job is not completed yet")
    files: dict[str, str] = snapshot.get("result_files", {})
    if name not in files:
        raise HTTPException(404, f"Output {name!r} not produced for this job")
    path = FilePath(files[name])
    if not path.is_file():
        raise HTTPException(410, "Output file vanished from disk")
    media_type = {
        '.json': 'application/json', '.jsonl': 'application/x-ndjson',
        '.sqlite': 'application/vnd.sqlite3', '.sqlite3': 'application/vnd.sqlite3',
        '.csv': 'text/csv',
    }.get(path.suffix.lower(), 'application/octet-stream')
    return FileResponse(path, filename=path.name, media_type=media_type)


@router.get("/{job_id}/preview/{name}")
async def preview_csv(
    job_id: Annotated[str, ApiPath(min_length=8, max_length=48)],
    name: Annotated[str, ApiPath(pattern=r"^(molecule_registry|reaction_events|species_timeseries)$")],
    limit: Annotated[int, Query(ge=1, le=2000)] = 500,
    offset: Annotated[int, Query(ge=0)] = 0,
) -> dict:
    manager = get_job_manager()
    snapshot = manager.snapshot(job_id)
    if snapshot is None:
        raise HTTPException(404, f"Job {job_id!r} not found")
    if snapshot.get("status") != JobStatus.COMPLETED.value:
        raise HTTPException(409, "Job is not completed yet")
    files: dict[str, str] = snapshot.get("result_files", {})
    path_text = files.get(name)
    if not path_text:
        raise HTTPException(404, f"Output {name!r} not produced for this job")
    path = FilePath(path_text)
    if not path.is_file():
        raise HTTPException(410, "Output file vanished from disk")
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for index, row in enumerate(reader):
            if index < offset:
                continue
            if len(rows) >= limit + 1:
                break
            rows.append(row)
    has_more = len(rows) > limit
    return {"rows": rows[:limit], "limit": limit, "offset": offset,
            "has_more": has_more, "possibly_truncated": has_more}


def _event_index_path(job_id: str) -> FilePath:
    snapshot = get_job_manager().snapshot(job_id)
    if snapshot is None:
        raise HTTPException(404, f"Job {job_id!r} not found")
    if snapshot.get('status') != JobStatus.COMPLETED.value:
        raise HTTPException(409, 'Job is not completed yet')
    filename = snapshot.get('result_files', {}).get('event_index')
    if not filename:
        raise HTTPException(409, '此旧任务没有事件详情索引，请使用当前版本重新分析。')
    return FilePath(filename)


@router.get('/{job_id}/events')
def query_events(
    job_id: Annotated[str, ApiPath(min_length=8, max_length=48)],
    offset: Annotated[int, Query(ge=0)] = 0,
    limit: Annotated[int, Query(ge=1, le=1000)] = 50,
    q: Annotated[str, Query(max_length=256)] = '',
    formula: Annotated[str, Query(max_length=128)] = '',
    species_key: Annotated[str, Query(max_length=256)] = '',
    start_step: Annotated[int | None, Query(ge=0)] = None,
    end_step: Annotated[int | None, Query(ge=0)] = None,
    status: Annotated[str, Query(pattern=r'^(all|accepted|raw|recrossing)$')] = 'all',
) -> dict:
    path = _event_index_path(job_id)
    try:
        return event_queries.list_events(path, offset=offset, limit=limit, q=q, formula=formula,
                                         species_key=species_key, start_step=start_step,
                                         end_step=end_step, status=status)
    except event_queries.EventIndexUnavailable as exc:
        raise HTTPException(409, str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(422, str(exc)) from exc


def _read_event(job_id: str, event_id: str) -> dict:
    try:
        payload = event_queries.get_event(_event_index_path(job_id), event_id)
    except event_queries.EventIndexUnavailable as exc:
        raise HTTPException(409, str(exc)) from exc
    if payload is None:
        raise HTTPException(404, f'Event {event_id!r} not found')
    return payload


@router.get('/{job_id}/events/{event_id}')
def event_detail(
    job_id: Annotated[str, ApiPath(min_length=8, max_length=48)],
    event_id: Annotated[str, ApiPath(pattern=r'^[A-Za-z0-9_.-]{1,80}$')],
) -> dict:
    return _read_event(job_id, event_id)


@router.get('/{job_id}/events/{event_id}/export')
def export_event(
    job_id: Annotated[str, ApiPath(min_length=8, max_length=48)],
    event_id: Annotated[str, ApiPath(pattern=r'^[A-Za-z0-9_.-]{1,80}$')],
) -> Response:
    payload = _read_event(job_id, event_id)
    return Response(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False),
                    media_type='application/json',
                    headers={'Content-Disposition': f'attachment; filename="{event_id}.json"'})


@router.delete("/{job_id}", status_code=204)
async def delete_job(
    job_id: Annotated[str, ApiPath(min_length=8, max_length=48)],
) -> None:
    manager = get_job_manager()
    if manager.snapshot(job_id) is None:
        raise HTTPException(404, f"Job {job_id!r} not found")
    cleanup_job_outputs(job_id)
    manager._jobs.pop(job_id, None)  # noqa: SLF001
    return None


@router.get("/{job_id}/outputs")
async def list_outputs(
    job_id: Annotated[str, ApiPath(min_length=8, max_length=48)],
) -> dict:
    manager = get_job_manager()
    snapshot = manager.snapshot(job_id)
    if snapshot is None:
        raise HTTPException(404, f"Job {job_id!r} not found")
    return {
        "files": [
            {"name": name, "url": f"/api/jobs/{job_id}/outputs/{name}"}
            for name in snapshot.get("result_files", {}).keys()
        ]
    }


@router.get("/{job_id}/counts")
def count_series(
    job_id: Annotated[str, ApiPath(min_length=8, max_length=48)],
    kind: Annotated[str, Query(pattern=r"^(bonds|functional)$")] = "bonds",
    metric: Annotated[str, Query(pattern=r"^(count|formed|broken|molecules)$")] = "count",
    max_points: Annotated[int, Query(ge=50, le=20000)] = 5000,
) -> dict:
    from ..services.count_analysis import CountOutputError, read_series

    snapshot = get_job_manager().snapshot(job_id)
    if snapshot is None:
        raise HTTPException(404, f"Job {job_id!r} not found")
    if snapshot.get("status") != JobStatus.COMPLETED.value:
        raise HTTPException(409, "原分析尚未完成")
    try:
        return read_series(snapshot, kind=kind, metric=metric, max_points=max_points)
    except CountOutputError as exc:
        raise HTTPException(409, str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(422, str(exc)) from exc


@router.post("/{job_id}/counts", status_code=202)
async def backfill_counts(
    job_id: Annotated[str, ApiPath(min_length=8, max_length=48)],
    background_tasks: BackgroundTasks,
) -> dict:
    from ..services.count_analysis import CountOutputError, prepare_backfill, run_backfill

    manager = get_job_manager()
    try:
        state, parameters = await prepare_backfill(manager, job_id)
    except KeyError as exc:
        raise HTTPException(404, f"Job {job_id!r} not found") from exc
    except (ValueError, CountOutputError, OSError) as exc:
        raise HTTPException(409, str(exc)) from exc
    if parameters is not None:
        background_tasks.add_task(run_backfill, manager, job_id, parameters)
    return {"job_id": job_id, **state}
