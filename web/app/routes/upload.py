"""Resumable upload routes for multi-gigabyte bonds trajectories.

The browser sends raw fixed-size chunks. Each chunk is written directly to its
final offset in one partial file, avoiding both whole-file request buffering
and a second full-size merge copy.
"""

from __future__ import annotations

import secrets
from typing import Annotated

from fastapi import APIRouter, HTTPException, Path, Request, status
from fastapi.responses import JSONResponse

from ..services import storage
from ..services.storage import UploadError, UploadMeta

router = APIRouter(prefix="/api/uploads", tags=["uploads"])


@router.post("", status_code=status.HTTP_201_CREATED)
async def init_upload(
    filename: Annotated[str, "Original filename"],
    total_size: Annotated[int, "Total file size in bytes"],
    chunk_size: Annotated[int, "Chunk size in bytes (must be positive)"],
    content_sha256: Annotated[str, ""] = "",
) -> JSONResponse:
    if total_size <= 0:
        raise HTTPException(400, "total_size must be positive")
    if chunk_size <= 0:
        raise HTTPException(400, "chunk_size must be positive")
    upload_id = secrets.token_urlsafe(16)
    try:
        meta = storage.create_upload(
            upload_id=upload_id,
            filename=filename,
            total_size=total_size,
            chunk_size=chunk_size,
        )
    except UploadError as exc:
        raise HTTPException(400, str(exc)) from exc
    payload = meta.to_dict()
    payload["content_sha256_expected"] = content_sha256 or None
    return JSONResponse(payload, status_code=201)


@router.put("/{upload_id}/chunks/{index}")
async def upload_chunk(
    request: Request,
    upload_id: Annotated[str, Path(min_length=8, max_length=64)],
    index: Annotated[int, Path(ge=0)],
) -> dict:
    # The browser sends one raw 8 MiB slice. This avoids multipart temporary
    # files while keeping request memory bounded by MAX_CHUNK_SIZE.
    payload = await request.body()
    try:
        meta = storage.receive_chunk(upload_id, index, payload)
    except storage.UploadError as exc:
        status_code = 404 if "not found" in str(exc) else 400
        raise HTTPException(status_code, str(exc)) from exc
    return _progress_payload(meta)


@router.post("/{upload_id}/complete")
async def complete_upload(
    upload_id: Annotated[str, Path(min_length=8, max_length=64)],
    content_sha256: Annotated[str, ""] = "",
) -> dict:
    try:
        meta = storage.finalize_upload(upload_id, expected_sha256=content_sha256)
    except UploadError as exc:
        raise HTTPException(400, str(exc)) from exc
    return _progress_payload(meta)


@router.delete("/{upload_id}", status_code=status.HTTP_204_NO_CONTENT)
async def discard_upload(
    upload_id: Annotated[str, Path(min_length=8, max_length=64)],
) -> None:
    storage.discard_upload(upload_id)
    return None


@router.get("/{upload_id}")
async def get_upload_status(
    upload_id: Annotated[str, Path(min_length=8, max_length=64)],
) -> dict:
    try:
        meta = storage.get_upload(upload_id)
    except UploadError as exc:
        raise HTTPException(404, str(exc)) from exc
    return _progress_payload(meta)


def _progress_payload(meta: UploadMeta) -> dict:
    received = len(meta.received_chunks)
    total = meta.total_chunks
    return {
        "upload_id": meta.upload_id,
        "filename": meta.filename,
        "total_size": meta.total_size,
        "chunk_size": meta.chunk_size,
        "total_chunks": total,
        "received_chunks": received,
        "missing_chunks": sorted(set(range(total)) - set(meta.received_chunks)),
        "progress": (received / total * 100.0) if total else 0.0,
        "completed": meta.completed_at is not None,
        "content_sha256": meta.content_sha256,
        "merged_path": meta.merged_path,
    }
