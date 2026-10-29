"""Static-file serving for the HTML frontend."""

from __future__ import annotations

from fastapi import APIRouter
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from ..config import STATIC_ROOT

router = APIRouter(tags=["frontend"])


@router.get("/", include_in_schema=False)
async def index() -> FileResponse:
    return FileResponse(STATIC_ROOT / "index.html", media_type="text/html")


def mount_static(app) -> None:
    """Mount the static dir at /static for any non-index assets we add later."""
    app.mount("/static", StaticFiles(directory=str(STATIC_ROOT)), name="static")