"""In-memory + on-disk job tracking for reaxnet analysis tasks.

A job is the unit of work that runs after a chunked upload completes. The
state is mirrored to ``jobs/<job_id>/meta.json`` so a service restart can pick
up where it left off. For the MVP the analysis itself runs inside a FastAPI
``BackgroundTask``; the file mirror keeps things observable and makes it easy
to swap in a real queue (ARQ / Celery) later without changing the public API.
"""

from __future__ import annotations

import asyncio
import json
import os
import time
from dataclasses import asdict, dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

from ..config import JOB_ROOT, JOB_TTL_SECONDS


class JobStatus(str, Enum):
    PENDING = "pending"
    UPLOADING = "uploading"
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


@dataclass
class JobInfo:
    job_id: str
    upload_id: str
    filename: str
    element_mapping: dict[int, str]
    config: dict[str, Any]
    source_mode: str = "upload"
    status: JobStatus = JobStatus.PENDING
    progress: float = 0.0
    stage: str = ""
    message: str = ""
    created_at: float = field(default_factory=time.time)
    started_at: float | None = None
    finished_at: float | None = None
    result_files: dict[str, str] = field(default_factory=dict)
    summary: dict[str, Any] = field(default_factory=dict)
    error: str | None = None
    count_analysis: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict:
        data = asdict(self)
        data["status"] = self.status.value
        return data

    @classmethod
    def from_dict(cls, data: dict) -> "JobInfo":
        data = dict(data)
        data["status"] = JobStatus(data.get("status", JobStatus.PENDING.value))
        return cls(**data)


class JobManager:
    """Process-wide singleton that tracks all jobs."""

    def __init__(self) -> None:
        self._jobs: dict[str, JobInfo] = {}
        self._lock = asyncio.Lock()
        JOB_ROOT.mkdir(parents=True, exist_ok=True)
        self._restore_from_disk()

    # ------------------------------------------------------------------
    # CRUD
    # ------------------------------------------------------------------

    async def create(
        self,
        job_id: str,
        upload_id: str,
        filename: str,
        element_mapping: dict[int, str],
        config: dict[str, Any],
        source_mode: str = "upload",
    ) -> JobInfo:
        async with self._lock:
            if job_id in self._jobs:
                raise ValueError(f"Job {job_id!r} already exists")
            job = JobInfo(
                job_id=job_id,
                upload_id=upload_id,
                filename=filename,
                element_mapping=dict(element_mapping),
                config=dict(config),
                source_mode=source_mode,
                status=JobStatus.QUEUED,
                stage="waiting",
                message="Queued for analysis",
            )
            self._jobs[job_id] = job
            self._persist(job)
            (JOB_ROOT / job_id).mkdir(parents=True, exist_ok=True)
            return job

    async def get(self, job_id: str) -> JobInfo | None:
        return self._jobs.get(job_id)

    def snapshot(self, job_id: str) -> dict | None:
        job = self._jobs.get(job_id)
        return job.to_dict() if job else None

    async def update(self, job_id: str, **fields: Any) -> JobInfo | None:
        async with self._lock:
            job = self._jobs.get(job_id)
            if job is None:
                return None
            for key, value in fields.items():
                if not hasattr(job, key):
                    raise AttributeError(f"JobInfo has no field {key!r}")
                setattr(job, key, value)
            self._persist(job)
            return job

    # ------------------------------------------------------------------
    # Garbage collection
    # ------------------------------------------------------------------

    def collect_garbage(self, now: float | None = None) -> int:
        now = now or time.time()
        removed = 0
        for entry in JOB_ROOT.iterdir():
            if not entry.is_dir():
                continue
            meta_path = entry / "meta.json"
            if not meta_path.is_file():
                continue
            try:
                data = json.loads(meta_path.read_text(encoding="utf-8"))
                anchor = data.get("finished_at") or data.get("created_at") or 0
            except (OSError, ValueError, json.JSONDecodeError):
                anchor = 0
            if now - anchor > JOB_TTL_SECONDS:
                import shutil
                shutil.rmtree(entry, ignore_errors=True)
                self._jobs.pop(entry.name, None)
                removed += 1
        return removed

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _persist(self, job: JobInfo) -> None:
        base = JOB_ROOT / job.job_id
        base.mkdir(parents=True, exist_ok=True)
        meta_path = base / "meta.json"
        tmp = meta_path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(job.to_dict(), indent=2), encoding="utf-8")
        os.replace(tmp, meta_path)

    def _restore_from_disk(self) -> None:
        if not JOB_ROOT.exists():
            return
        for entry in JOB_ROOT.iterdir():
            if not entry.is_dir():
                continue
            meta_path = entry / "meta.json"
            if not meta_path.is_file():
                continue
            try:
                data = json.loads(meta_path.read_text(encoding="utf-8"))
            except (OSError, ValueError, json.JSONDecodeError):
                continue
            try:
                job = JobInfo.from_dict(data)
            except (TypeError, ValueError):
                continue
            # Anything mid-flight before a restart cannot be resumed safely
            # in the BackgroundTask MVP; mark it as failed with a hint.
            if job.status in (JobStatus.RUNNING, JobStatus.QUEUED):
                job.status = JobStatus.FAILED
                job.error = (
                    "Service restarted while job was running. "
                    "Please resubmit the upload."
                )
                job.finished_at = time.time()
                self._persist(job)
            if job.count_analysis.get("status") == "running":
                job.count_analysis = {
                    **job.count_analysis,
                    "status": "failed",
                    "error": "Service restarted during count analysis; request counts again.",
                    "finished_at": time.time(),
                }
                self._persist(job)
            self._jobs[job.job_id] = job


_singleton: JobManager | None = None


def get_job_manager() -> JobManager:
    global _singleton
    if _singleton is None:
        _singleton = JobManager()
    return _singleton
