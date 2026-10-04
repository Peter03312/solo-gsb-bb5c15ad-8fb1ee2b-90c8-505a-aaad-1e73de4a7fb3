"""FastAPI 入口：上传、预览、下载。

结果只在校验完全通过后保存在进程内存中（BytesIO 等价物），
任何解析错误都直接 400 返回，不创建、更不提供任何可下载文件。
"""

from __future__ import annotations

import io
import os
import secrets
import time
from dataclasses import dataclass

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles

from .bed import parse_bed
from .chain import parse_chain
from .engine import LiftOverEngine
from .report import build_tsv, rows_to_preview

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
STATIC_DIR = os.path.join(BASE_DIR, "..", "static")

MAX_UPLOAD_BYTES = int(os.getenv("LIFTOVER_MAX_UPLOAD_MB", "50")) * 1024 * 1024
JOB_TTL_SECONDS = int(os.getenv("LIFTOVER_JOB_TTL_SECONDS", "3600"))
MAX_JOBS = int(os.getenv("LIFTOVER_MAX_JOBS", "100"))

app = FastAPI(title="离线基因组区间迁移 (BED liftover)", version="1.0.0")


@dataclass
class _Job:
    tsv_bytes: bytes
    filename: str
    created_at: float


_jobs: dict[str, _Job] = {}


def _prune_jobs() -> None:
    now = time.time()
    expired = [k for k, j in _jobs.items() if now - j.created_at > JOB_TTL_SECONDS]
    for k in expired:
        _jobs.pop(k, None)
    # 兜底：超过容量时淘汰最早的任务
    if len(_jobs) > MAX_JOBS:
        for k in sorted(_jobs, key=lambda x: _jobs[x].created_at)[: len(_jobs) - MAX_JOBS]:
            _jobs.pop(k, None)


def _read_upload(f: UploadFile) -> str:
    data = f.file.read()
    if not data:
        raise HTTPException(status_code=400, detail={"errors": [f"文件 {f.filename} 为空"]})
    if len(data) > MAX_UPLOAD_BYTES:
        raise HTTPException(
            status_code=413,
            detail={"errors": [f"文件 {f.filename} 超过大小上限 "
                               f"{MAX_UPLOAD_BYTES // (1024 * 1024)} MB"]},
        )
    try:
        # BOM 仅允许出现在文件最前面
        return data.decode("utf-8-sig")
    except UnicodeDecodeError:
        try:
            return data.decode("latin-1")
        except UnicodeDecodeError:
            raise HTTPException(
                status_code=400,
                detail={"errors": [f"文件 {f.filename} 不是有效的 UTF-8 文本"]},
            )


@app.exception_handler(HTTPException)
async def _http_exc_handler(request: Request, exc: HTTPException) -> JSONResponse:
    return JSONResponse(status_code=exc.status_code, content=exc.detail
                        if isinstance(exc.detail, dict) else {"errors": [str(exc.detail)]})


@app.get("/healthz")
async def healthz() -> dict[str, str]:
    return {"status": "ok"}


@app.post("/api/liftover")
async def liftover(
    bed: UploadFile = File(...),
    chain: UploadFile = File(...),
    side: str = Form(...),
):
    if side not in ("target", "query"):
        raise HTTPException(status_code=400,
                            detail={"errors": ["side 必须是 target 或 query"]})

    bed_text = _read_upload(bed)
    chain_text = _read_upload(chain)

    bed_result = parse_bed(bed_text)
    chain_result = parse_chain(chain_text)

    errors = list(bed_result.errors) + list(chain_result.errors)
    if errors:
        # 校验失败：只返回错误，绝不写入 _jobs —— 不存在可下载文件
        raise HTTPException(status_code=400, detail={"errors": errors})

    engine = LiftOverEngine(chain_result.chains, side)
    rows = engine.convert_all(bed_result.records)
    preview = rows_to_preview(rows)

    tsv = build_tsv(rows)
    # Excel 友好的 UTF-8 BOM；HTTP 下载本质上仍是内存字节，不落盘
    payload = b"\xef\xbb\xbf" + tsv.encode("utf-8")
    token = secrets.token_urlsafe(24)
    stamp = time.strftime("%Y%m%d-%H%M%S", time.gmtime())
    filename = f"liftover_{side}_to_{'query' if side == 'target' else 'target'}_{stamp}.tsv"
    _prune_jobs()
    _jobs[token] = _Job(tsv_bytes=payload, filename=filename, created_at=time.time())

    return {"download_token": token, "download_name": filename, "preview": preview}


@app.get("/api/download/{token}")
async def download(token: str) -> Response:
    _prune_jobs()
    job = _jobs.get(token)
    if job is None:
        raise HTTPException(status_code=404,
                            detail={"errors": ["结果不存在或已过期，请重新提交转换"]})
    return Response(
        content=job.tsv_bytes,
        media_type="text/tab-separated-values; charset=utf-8",
        headers={"Content-Disposition": f'attachment; filename="{job.filename}"'},
    )


# 静态资源（纯本地，无任何外部 CDN 依赖）
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")


@app.get("/", response_class=HTMLResponse)
async def index() -> FileResponse:
    return FileResponse(os.path.join(STATIC_DIR, "index.html"))
