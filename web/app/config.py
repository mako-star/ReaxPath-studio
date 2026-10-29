"""Runtime configuration for the web service."""

from __future__ import annotations

from pathlib import Path
from typing import Final

# Project layout: web/ holds pyproject.toml, the FastAPI package lives in app/.
WEB_ROOT: Final[Path] = Path(__file__).resolve().parent.parent
UPLOAD_ROOT: Final[Path] = WEB_ROOT / "uploads"
JOB_ROOT: Final[Path] = WEB_ROOT / "jobs"
STATIC_ROOT: Final[Path] = Path(__file__).resolve().parent / "static"

# Resumable upload defaults. 8 MiB matches the default chunk size used by most
# browser chunked-upload libraries; tune via the HTML frontend if needed.
DEFAULT_CHUNK_SIZE: Final[int] = 8 * 1024 * 1024
MAX_FILE_SIZE: Final[int] = 16 * 1024 * 1024 * 1024  # 16 GiB hard cap
MAX_CHUNK_SIZE: Final[int] = 32 * 1024 * 1024        # 32 MiB per chunk

# How long an unfinished upload session lingers before cleanup sweeps it.
UPLOAD_TTL_SECONDS: Final[int] = 24 * 3600

# How long finished job outputs are kept on disk before cleanup.
JOB_TTL_SECONDS: Final[int] = 7 * 24 * 3600


def ensure_directories() -> None:
    """Make sure the upload and job directories exist."""
    UPLOAD_ROOT.mkdir(parents=True, exist_ok=True)
    JOB_ROOT.mkdir(parents=True, exist_ok=True)