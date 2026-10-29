"""Background, bounded-memory ReaxFF analysis service."""

from __future__ import annotations

import asyncio
import logging
import shutil
import time
from pathlib import Path
from typing import Any

from reaxnet.analysis import AnalyzerConfig
from reaxnet.parser import BondsParseError, parse_element_mapping
from reaxnet.streaming import analyze_file_streaming

from ..config import JOB_ROOT
from .job_manager import JobManager, JobStatus

logger = logging.getLogger(__name__)


def build_analyzer_config(payload: dict[str, Any]) -> AnalyzerConfig:
    removed = {"bond_cutoff", "valence_overrides", "coordination_overrides"} & payload.keys()
    if removed:
        raise ValueError(f"Connection filtering options were removed: {sorted(removed)}")
    return AnalyzerConfig(
        timestep_scale=float(payload.get("timestep_scale", 1.0)),
        recrossing_window=float(payload.get("recrossing_window", 0.0)),
        expected_interval=payload.get("expected_interval"),
        generate_identifiers=bool(payload.get("generate_identifiers", True)),
    )


async def run_analysis(
    job_manager: JobManager,
    job_id: str,
    bonds_path: Path,
    element_mapping_payload: str,
    config_payload: dict[str, Any],
) -> None:
    """Stream one trajectory into CSV/JSON outputs without retaining all frames."""

    try:
        element_mapping = parse_element_mapping(element_mapping_payload)
        config = build_analyzer_config(config_payload)
    except (TypeError, ValueError) as exc:
        await job_manager.update(
            job_id,
            status=JobStatus.FAILED,
            error=f"Invalid configuration: {exc}",
            finished_at=time.time(),
            message="Aborted: invalid configuration",
        )
        return

    source_stat = bonds_path.stat()
    source_size = source_stat.st_size
    saved_config = dict((job_manager.snapshot(job_id) or {}).get("config", config_payload))
    saved_config.update(source_path=str(bonds_path.resolve()),
                        source_size_bytes=source_size, source_mtime_ns=source_stat.st_mtime_ns)
    await job_manager.update(
        job_id,
        status=JobStatus.RUNNING,
        stage="analysing",
        progress=1.0,
        started_at=time.time(),
        config=saved_config,
        summary={
            "frames_processed": 0,
            "bytes_processed": 0,
            "bytes_total": source_size,
            "streaming_analysis": True,
        },
        message="Streaming trajectory frames",
    )

    output_dir = JOB_ROOT / job_id / "outputs"
    output_dir.mkdir(parents=True, exist_ok=True)
    loop = asyncio.get_running_loop()
    last_progress_future = None

    def report_progress(progress_payload: dict[str, int | float]) -> None:
        nonlocal last_progress_future
        bytes_total = int(progress_payload.get("bytes_total", source_size)) or 1
        bytes_processed = int(progress_payload.get("bytes_processed", 0))
        percent = 1.0 + min(1.0, bytes_processed / bytes_total) * 96.0
        last_progress_future = asyncio.run_coroutine_threadsafe(
            job_manager.update(
                job_id,
                progress=percent,
                stage="analysing",
                summary={
                    "frames_processed": int(progress_payload.get("frames_processed", 0)),
                    "bytes_processed": bytes_processed,
                    "bytes_total": bytes_total,
                    "streaming_analysis": True,
                },
                message=(
                    f"Processed {int(progress_payload.get('frames_processed', 0)):,} frames"
                ),
            ),
            loop,
        )

    try:
        result = await asyncio.to_thread(
            analyze_file_streaming,
            bonds_path,
            element_mapping,
            output_dir,
            config,
            progress_callback=report_progress,
        )
        if last_progress_future is not None:
            await asyncio.wrap_future(last_progress_future)
    except (BondsParseError, ValueError, OSError) as exc:
        logger.exception("Job %s: streaming analysis failed", job_id)
        await job_manager.update(
            job_id,
            status=JobStatus.FAILED,
            error=f"{type(exc).__name__}: {exc}",
            finished_at=time.time(),
            stage="analysing",
            message="Analysis failed",
        )
        return
    except Exception as exc:  # pragma: no cover - defensive
        logger.exception("Job %s: unexpected failure", job_id)
        await job_manager.update(
            job_id,
            status=JobStatus.FAILED,
            error=f"{type(exc).__name__}: {exc}",
            finished_at=time.time(),
            stage="analysing",
            message="Analysis failed",
        )
        return

    summary = {
        "frames": result.summary["frame_count"],
        "species": result.summary["species_count"],
        "events_total": result.summary["detected_event_count"],
        "events_accepted": result.summary["accepted_event_count"],
        "events_filtered": result.summary["filtered_event_count"],
        "reactions": result.summary["reaction_count"],
        "source_size_bytes": result.summary["source_size_bytes"],
        "peak_buffered_states": result.summary["peak_buffered_states"],
        "peak_pending_events": result.summary["peak_pending_events"],
        "streaming_analysis": True,
    }
    result_paths = {name: str(path) for name, path in result.paths.items()}
    await job_manager.update(
        job_id,
        status=JobStatus.COMPLETED,
        progress=100.0,
        stage="done",
        finished_at=time.time(),
        result_files=result_paths,
        summary=summary,
        count_analysis=({"status": "completed", "frame_count": summary["frames"],
                         "progress": 100.0, "finished_at": time.time()}
                        if "count_metadata" in result_paths else {}),
        message=(
            f"Done: {summary['frames']:,} frames, {summary['species']} species, "
            f"{summary['events_accepted']}/{summary['events_total']} events accepted"
        ),
    )
    logger.info("Job %s completed: %s", job_id, summary)


def cleanup_job_outputs(job_id: str) -> bool:
    base = JOB_ROOT / job_id
    if not base.exists():
        return False
    shutil.rmtree(base, ignore_errors=True)
    return True
