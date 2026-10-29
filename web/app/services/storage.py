"""Disk-backed chunked upload storage.

Each upload session lives in ``uploads/<upload_id>/`` and contains metadata plus
one ``merged.part`` file. Chunks are written directly at their final byte offset,
so a 4 GiB upload does not require a second 4 GiB merge copy.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

from ..config import (
    DEFAULT_CHUNK_SIZE,
    MAX_CHUNK_SIZE,
    MAX_FILE_SIZE,
    UPLOAD_ROOT,
    UPLOAD_TTL_SECONDS,
)


@dataclass
class UploadMeta:
    upload_id: str
    filename: str
    total_size: int
    chunk_size: int
    total_chunks: int
    received_chunks: list[int] = field(default_factory=list)
    content_sha256: str = ""
    created_at: float = field(default_factory=time.time)
    completed_at: float | None = None
    merged_path: str | None = None

    def to_dict(self) -> dict:
        return asdict(self)


class UploadError(Exception):
    """Raised on any user-facing upload problem (400/404/409)."""


class ChunkValidationError(UploadError):
    """Raised when a chunk payload is malformed or oversized."""


def _upload_dir(upload_id: str) -> Path:
    safe = "".join(c for c in upload_id if c.isalnum() or c in "-_")
    if not safe or safe != upload_id:
        raise UploadError(f"Invalid upload id: {upload_id!r}")
    path = UPLOAD_ROOT / safe
    if not path.resolve().is_relative_to(UPLOAD_ROOT.resolve()):
        raise UploadError(f"Invalid upload id: {upload_id!r}")
    return path


def create_upload(
    upload_id: str,
    filename: str,
    total_size: int,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> UploadMeta:
    if total_size <= 0:
        raise UploadError("total_size must be positive")
    if total_size > MAX_FILE_SIZE:
        raise UploadError(
            f"File too large: {total_size} bytes > {MAX_FILE_SIZE} bytes"
        )
    if chunk_size <= 0 or chunk_size > MAX_CHUNK_SIZE:
        raise UploadError(
            f"chunk_size must be in (0, {MAX_CHUNK_SIZE}], got {chunk_size}"
        )

    base = _upload_dir(upload_id)
    if base.exists():
        raise UploadError(f"Upload id {upload_id!r} already exists")

    base.mkdir(parents=False)

    total_chunks = (total_size + chunk_size - 1) // chunk_size
    meta = UploadMeta(
        upload_id=upload_id,
        filename=filename,
        total_size=total_size,
        chunk_size=chunk_size,
        total_chunks=total_chunks,
    )
    _write_meta(base, meta)
    return meta


def get_upload(upload_id: str) -> UploadMeta:
    base = _upload_dir(upload_id)
    meta_path = base / "meta.json"
    if not meta_path.is_file():
        raise UploadError(f"Upload {upload_id!r} not found")
    return _read_meta(base)


def receive_chunk(
    upload_id: str,
    index: int,
    payload: bytes,
) -> UploadMeta:
    base = _upload_dir(upload_id)
    if not base.exists():
        raise UploadError(f"Upload {upload_id!r} not found")
    meta = _read_meta(base)

    if not (0 <= index < meta.total_chunks):
        raise ChunkValidationError(
            f"chunk index {index} out of range [0, {meta.total_chunks})"
        )
    expected = _expected_chunk_size(meta, index)
    actual = len(payload)
    if actual > MAX_CHUNK_SIZE:
        raise ChunkValidationError(
            f"chunk {index} too large: {actual} > {MAX_CHUNK_SIZE}"
        )
    if actual != expected:
        raise ChunkValidationError(
            f"chunk {index} size mismatch: got {actual}, expected {expected}"
        )

    partial_path = base / "merged.part"
    if index in meta.received_chunks:
        if not partial_path.is_file():
            raise ChunkValidationError(f"chunk {index} metadata exists but data is missing")
        with partial_path.open("rb") as handle:
            handle.seek(index * meta.chunk_size)
            existing = handle.read(expected)
        if existing != payload:
            raise ChunkValidationError(f"chunk {index} was already received with different data")
        return meta

    mode = "r+b" if partial_path.exists() else "w+b"
    with partial_path.open(mode) as handle:
        handle.seek(index * meta.chunk_size)
        handle.write(payload)

    meta.received_chunks.append(index)
    meta.received_chunks.sort()
    _write_meta(base, meta)
    return meta


def _expected_chunk_size(meta: UploadMeta, index: int) -> int:
    if index < meta.total_chunks - 1:
        return meta.chunk_size
    return meta.total_size - meta.chunk_size * (meta.total_chunks - 1)


def finalize_upload(upload_id: str, expected_sha256: str = "") -> UploadMeta:
    base = _upload_dir(upload_id)
    meta = _read_meta(base)
    if meta.completed_at is not None:
        return meta

    missing = sorted(set(range(meta.total_chunks)) - set(meta.received_chunks))
    if missing:
        raise UploadError(
            f"Upload {upload_id!r} incomplete, missing chunks: {missing[:10]}"
            f"{'...' if len(missing) > 10 else ''}"
        )

    partial_path = base / "merged.part"
    if not partial_path.is_file() or partial_path.stat().st_size != meta.total_size:
        raise UploadError(
            f"Upload data size mismatch: expected {meta.total_size} bytes"
        )

    hasher = hashlib.sha256() if expected_sha256 else None
    if hasher is not None:
        with partial_path.open("rb") as src:
            while buf := src.read(8 * 1024 * 1024):
                hasher.update(buf)

    if hasher is not None:
        digest = hasher.hexdigest()
        if expected_sha256 and digest.lower() != expected_sha256.lower():
            raise UploadError(
                f"SHA256 mismatch: declared {expected_sha256}, got {digest}"
            )
        meta.content_sha256 = digest

    merged_path = base / "merged"
    os.replace(partial_path, merged_path)
    meta.completed_at = time.time()
    meta.merged_path = str(merged_path)
    _write_meta(base, meta)
    (base / "complete").touch()
    return meta


def discard_upload(upload_id: str) -> None:
    base = _upload_dir(upload_id)
    if base.exists():
        shutil.rmtree(base, ignore_errors=True)


def collect_garbage(now: float | None = None) -> int:
    """Remove upload sessions older than ``UPLOAD_TTL_SECONDS``.

    Returns the number of sessions removed.
    """
    if not UPLOAD_ROOT.exists():
        return 0
    now = now or time.time()
    removed = 0
    for entry in UPLOAD_ROOT.iterdir():
        if not entry.is_dir():
            continue
        meta_path = entry / "meta.json"
        if not meta_path.is_file():
            continue
        try:
            meta = _read_meta(entry)
        except (OSError, ValueError, json.JSONDecodeError):
            shutil.rmtree(entry, ignore_errors=True)
            removed += 1
            continue
        anchor = meta.completed_at or meta.created_at
        if now - anchor > UPLOAD_TTL_SECONDS:
            shutil.rmtree(entry, ignore_errors=True)
            removed += 1
    return removed


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _write_meta(base: Path, meta: UploadMeta) -> None:
    tmp = base / "meta.json.tmp"
    with tmp.open("w", encoding="utf-8") as handle:
        json.dump(meta.to_dict(), handle, indent=2)
    os.replace(tmp, base / "meta.json")


def _read_meta(base: Path) -> UploadMeta:
    meta_path = base / "meta.json"
    data = json.loads(meta_path.read_text(encoding="utf-8"))
    return UploadMeta(**data)
