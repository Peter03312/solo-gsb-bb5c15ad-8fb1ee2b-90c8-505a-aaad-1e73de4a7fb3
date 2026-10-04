"""结果导出：TSV 文本与预览 JSON。

TSV 列（英文表头，避免下载文件里混入列名歧义）::

    source_line, status, source_chrom, source_start, source_end,
    candidate_index, chain_id, dest_chrom, dest_strand,
    fragment_index, src_start, src_end, dest_start, dest_end, reason

每个源行号可能对应多行（多个片段 / 多个歧义候选）。
unmapped 行没有候选，相关列写 '.'。
"""

from __future__ import annotations

import csv
import io
from dataclasses import asdict
from typing import Any

from .engine import RowResult

STATUS_LABELS = {
    "mapped": "mapped",
    "partial": "partial",
    "ambiguous": "ambiguous",
    "unmapped": "unmapped",
}

TSV_HEADER = [
    "source_line",
    "status",
    "source_chrom",
    "source_start",
    "source_end",
    "candidate_index",
    "chain_id",
    "dest_chrom",
    "dest_strand",
    "fragment_index",
    "src_start",
    "src_end",
    "dest_start",
    "dest_end",
    "reason",
]


def build_tsv(rows: list[RowResult]) -> str:
    buf = io.StringIO(newline="")
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n",
                        quoting=csv.QUOTE_MINIMAL)
    writer.writerow(TSV_HEADER)
    for row in rows:
        if not row.candidates:
            writer.writerow([
                row.line_no, row.status, row.src_chrom, row.src_start, row.src_end,
                ".", ".", ".", ".", ".", ".", ".", ".", ".", row.reason,
            ])
            continue
        for ci, cand in enumerate(row.candidates, start=1):
            for fi, frag in enumerate(cand.fragments, start=1):
                writer.writerow([
                    row.line_no, row.status, row.src_chrom, row.src_start, row.src_end,
                    ci, cand.chain_id, cand.dest_chrom, cand.dest_strand,
                    fi, frag.src_start, frag.src_end, frag.dest_start, frag.dest_end,
                    row.reason,
                ])
    return buf.getvalue()


def row_to_dict(row: RowResult) -> dict[str, Any]:
    return asdict(row)


def rows_to_preview(rows: list[RowResult]) -> dict[str, Any]:
    counts = {"mapped": 0, "partial": 0, "ambiguous": 0, "unmapped": 0}
    for r in rows:
        counts[r.status] += 1
    return {
        "total": len(rows),
        "counts": counts,
        "rows": [row_to_dict(r) for r in rows],
    }
