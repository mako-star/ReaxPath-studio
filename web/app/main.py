"""FastAPI application entry point."""

from __future__ import annotations

import logging

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .config import DEFAULT_CHUNK_SIZE, MAX_FILE_SIZE, ensure_directories
from .routes import frontend as frontend_route
from .routes import jobs as jobs_route
from .routes import upload as upload_route
from .services.job_manager import get_job_manager


def create_app() -> FastAPI:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    ensure_directories()
    # Eagerly construct the singleton so any on-disk job state is restored.
    get_job_manager()

    app = FastAPI(
        title="smiletosvg web",
        description="Upload a ReaxFF bonds trajectory and run reaxnet analysis in the background.",
        version="0.1.0",
    )

    # Permissive CORS so a separate dev frontend (Vite, etc.) can talk to us.
    # Tighten before deploying.
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_methods=["*"],
        allow_headers=["*"],
        allow_credentials=False,
    )

    app.include_router(upload_route.router)
    app.include_router(jobs_route.router)
    frontend_route.mount_static(app)
    app.include_router(frontend_route.router)

    @app.get("/healthz", tags=["meta"])
    async def healthz() -> dict:
        return {
            "status": "ok",
            "max_file_size": MAX_FILE_SIZE,
            "default_chunk_size": DEFAULT_CHUNK_SIZE,
            "streaming_analysis": True,
            "local_path_enabled": True,
        }

    return app


app = create_app()


def main() -> None:
    """Console-script entry: ``smiletosvg-web``."""
    import uvicorn

    uvicorn.run(
        "app.main:app",
        host="127.0.0.1",
        port=8000,
        reload=False,
        log_level="info",
    )


if __name__ == "__main__":
    main()
