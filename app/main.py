"""FastAPI 入口：离线基因组区间迁移服务。

接口：
  GET  /                 前端页面
  GET  /api/health       健康检查
  POST /api/convert      上传 BED + chain，返回 JSON 预览
  POST /api/export       上传 BED + chain，返回 TSV 下载（失败时 422 JSON，不产生文件）
"""

from __future__ import annotations

import os
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, Response
from pydantic import BaseModel

from .liftover import MAX_BYTES, SIDES, LiftError, build_tsv, run_liftover

PREVIEW_LIMIT = int(os.environ.get("PREVIEW_LIMIT", "1000"))

app = FastAPI(title="离线基因组区间迁移", version="1.0.0")

STATIC_DIR = Path(__file__).resolve().parent / "static"


async def _read_upload(f: UploadFile, label: str) -> bytes:
    if not f or not f.filename:
        raise HTTPException(status_code=422, detail=[{"error": f"未选择{label}"}])
    data = await f.read()
    if not data:
        raise HTTPException(status_code=422, detail=[{"error": f"{label}为空文件"}])
    if len(data) > MAX_BYTES:
        raise HTTPException(
            status_code=422,
            detail=[{"error": f"{label}超过 {MAX_BYTES // (1024 * 1024)} MiB 上限"}],
        )
    return data


def _do_convert(bed_bytes: bytes, chain_bytes: bytes, side: str) -> dict:
    if side not in SIDES:
        raise HTTPException(
            status_code=422,
            detail=[{"error": f"源侧必须为 'target' 或 'query'，得到 {side!r}"}],
        )
    try:
        return run_liftover(bed_bytes, chain_bytes, side)
    except LiftError as e:
        # 逐条错误拆分，前端可分行展示
        errors = [line for line in str(e).splitlines() if line.strip()]
        raise HTTPException(status_code=422, detail=[{"error": x} for x in errors])


class ConvertResponse(BaseModel):
    side: str
    bed_rows: int
    result_rows: int
    stats: dict
    rows: list
    truncated: bool
    preview_limit: int


@app.get("/api/health")
def health() -> dict:
    return {"status": "ok"}


@app.post("/api/convert")
async def convert(
    bed: UploadFile = File(...),
    chain: UploadFile = File(...),
    side: str = Form(...),
) -> ConvertResponse:
    bed_bytes = await _read_upload(bed, "BED 文件")
    chain_bytes = await _read_upload(chain, "chain 文件")
    payload = _do_convert(bed_bytes, chain_bytes, side)
    truncated = payload["result_rows"] > PREVIEW_LIMIT
    return ConvertResponse(
        side=payload["side"],
        bed_rows=payload["bed_rows"],
        result_rows=payload["result_rows"],
        stats=payload["stats"],
        rows=payload["rows"][:PREVIEW_LIMIT],
        truncated=truncated,
        preview_limit=PREVIEW_LIMIT,
    )


@app.post("/api/export")
async def export(
    bed: UploadFile = File(...),
    chain: UploadFile = File(...),
    side: str = Form(...),
):
    bed_bytes = await _read_upload(bed, "BED 文件")
    chain_bytes = await _read_upload(chain, "chain 文件")
    # 校验或转换失败：返回 422 JSON，浏览器拿到的是错误响应而非下载文件
    payload = _do_convert(bed_bytes, chain_bytes, side)
    tsv = build_tsv(payload)
    return Response(
        # 加 UTF-8 BOM 便于 Excel 直接打开
        content=b"\xef\xbb\xbf" + tsv.encode("utf-8"),
        media_type="text/tab-separated-values; charset=utf-8",
        headers={"Content-Disposition": 'attachment; filename="liftover_result.tsv"'},
    )


@app.get("/", response_class=HTMLResponse)
def index() -> FileResponse:
    return FileResponse(STATIC_DIR / "index.html")
