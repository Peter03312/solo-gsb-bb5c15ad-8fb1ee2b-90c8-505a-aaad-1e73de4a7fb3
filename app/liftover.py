"""离线基因组区间坐标迁移引擎（BED + UCSC chain）。

坐标约定：BED 与 chain 均按零基左闭右开 [start, end) 解释，
转换后的结果同样保持零基左闭右开（负链区间已在算法内部反转）。
"""

from __future__ import annotations

from bisect import bisect_left
from dataclasses import dataclass, field
from typing import Optional

TARGET = "target"
QUERY = "query"
SIDES = (TARGET, QUERY)

MAX_ISSUES = 50
MAX_BYTES = 50 * 1024 * 1024  # 与路由层保持一致


class LiftError(Exception):
    """输入数据存在致命错误（坏行/格式错误）。"""


@dataclass(frozen=True)
class Block:
    """一段无缺口的对齐块，两侧坐标均为零基左闭右开。"""

    t_start: int
    t_end: int
    q_start: int  # 始终为数值较小的坐标（query 负链时对应原链末端）
    q_end: int
    same: bool  # 两侧相对方向是否一致（qstrand == '+'）


@dataclass
class Chain:
    chain_id: str
    t_name: str
    t_size: int
    t_strand: str
    t_start: int
    t_end: int
    q_name: str
    q_size: int
    q_strand: str
    q_start: int
    q_end: int
    blocks: list[Block] = field(default_factory=list)

    def side_chrom(self, side: str) -> str:
        return self.t_name if side == TARGET else self.q_name


# --------------------------------------------------------------------------- #
# chain 解析
# --------------------------------------------------------------------------- #


def _parse_int(tok: str, lineno: int, what: str) -> int:
    try:
        return int(tok)
    except ValueError:
        raise LiftError(f"chain 第 {lineno} 行：{what} 不是整数（得到 {tok!r}）")


class _ChainParser:
    def __init__(self, text: str):
        self.lines = text.splitlines()
        self.issues: list[str] = []
        self.chains: list[Chain] = []
        self._cur: Optional[Chain] = None
        self._cur_lines: list[Block] = []
        # 头部记录的已消耗长度（含块与 gap），用于最终一致性校验
        self._t_consumed = 0
        self._q_consumed = 0
        self._header_line = 0
        self._terminated = False

    def fail(self, msg: str) -> None:
        self.issues.append(msg)

    def run(self) -> list[Chain]:
        for i, raw in enumerate(self.lines, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.split(maxsplit=1)[0] == "chain":
                self._finish_chain(i)
                self._start_chain(i, line)
            elif self._cur is None:
                self.fail(f"chain 第 {i} 行：孤立数据行 {line!r}，缺少 chain 头部")
            else:
                self._data_line(i, line)
        self._finish_chain(None)
        if self.issues:
            head = self.issues[:MAX_ISSUES]
            if len(self.issues) > MAX_ISSUES:
                head.append(f"……另有 {len(self.issues) - MAX_ISSUES} 个错误未展开")
            raise LiftError("\n".join(head))
        return self.chains

    def _start_chain(self, lineno: int, line: str) -> None:
        toks = line.split()
        if len(toks) != 13:
            self.fail(f"chain 第 {lineno} 行：头部应有 13 个字段，实际 {len(toks)} 个")
            self._cur = None
            return
        if toks[0] != "chain":
            self.fail(f"chain 第 {lineno} 行：行首应为 'chain'，得到 {toks[0]!r}")
            self._cur = None
            return

        # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
        def pi(idx: int, what: str) -> Optional[int]:
            try:
                return int(toks[idx])
            except ValueError:
                self.fail(
                    f"chain 第 {lineno} 行：{what} 不是整数（得到 {toks[idx]!r}）"
                )
                return None

        score = pi(1, "score")
        t_size = pi(3, "target size")
        t_start = pi(5, "target start")
        t_end = pi(6, "target end")
        q_size = pi(8, "query size")
        q_start = pi(10, "query start")
        q_end = pi(11, "query end")

        t_name, q_name, chain_id = toks[2], toks[7], toks[12]
        t_strand, q_strand = toks[4], toks[9]
        if t_strand != "+":
            self.fail(
                f"chain 第 {lineno} 行：本工具仅支持 target 正链（target strand='+'），"
                f"得到 {t_strand!r}"
            )
        if q_strand not in ("+", "-"):
            self.fail(f"chain 第 {lineno} 行：query strand 只能为 + 或 -，得到 {q_strand!r}")
        if score is not None and int(score) < 0:
            self.fail(f"chain 第 {lineno} 行：score 不能为负")
        vals = {
            "target size": t_size,
            "target start": t_start,
            "target end": t_end,
            "query size": q_size,
            "query start": q_start,
            "query end": q_end,
        }
        bad = False
        for name, v in vals.items():
            if v is not None and v < 0:
                self.fail(f"chain 第 {lineno} 行：{name} 不能为负（{v}）")
                bad = True
        # 边界合法性
        pairs_t = (("target", t_size, t_start, t_end),)
        pairs_q = (("query", q_size, q_start, q_end),)
        for label, size, start, end in pairs_t + pairs_q:
            if size is None or start is None or end is None or bad:
                continue
            if not (0 <= start < end <= size):
                self.fail(
                    f"chain 第 {lineno} 行：{label} 坐标不满足 0 <= start < end <= size"
                    f"（start={start}, end={end}, size={size}）"
                )
        if self.issues:
            self._cur = None
            return

        self._cur = Chain(
            chain_id=chain_id,
            t_name=t_name,
            t_size=t_size,
            t_strand=t_strand,
            t_start=t_start,
            t_end=t_end,
            q_name=q_name,
            q_size=q_size,
            q_strand=q_strand,
            q_start=q_start,
            q_end=q_end,
        )
        self._cur_lines = []
        self._t_consumed = 0
        self._q_consumed = 0
        self._header_line = lineno
        self._terminated = False

    def _data_line(self, lineno: int, line: str) -> None:
        assert self._cur is not None
        if self._terminated:
            self.fail(
                f"chain 第 {lineno} 行（chain {self._cur.chain_id}）："
                f"终止行（仅 size）之后不允许再出现数据行"
            )
            return
        toks = line.split()
        if len(toks) not in (1, 3):
            self.fail(
                f"chain 第 {lineno} 行（chain {self._cur.chain_id}）："
                f"数据行应为 size 或 size dt dq，实际 {len(toks)} 个字段"
            )
            return
        size = _safe_int(toks[0], lineno, "size", self.issues)
        dt = _safe_int(toks[1], lineno, "dt", self.issues) if len(toks) == 3 else 0
        dq = _safe_int(toks[2], lineno, "dq", self.issues) if len(toks) == 3 else 0
        if self.issues:
            return
        if size <= 0:
            self.fail(
                f"chain 第 {lineno} 行（chain {self._cur.chain_id}）：对齐块大小须 >= 1"
            )
            return
        if dt < 0 or dq < 0:
            self.fail(
                f"chain 第 {lineno} 行（chain {self._cur.chain_id}）：gap 不能为负"
            )
            return

        c = self._cur
        t_off = c.t_start + self._t_consumed
        q_off = c.q_start + self._q_consumed
        # 本块之前的 gap 不能把游标推过头部边界（否则块起点已越界）
        if t_off > c.t_end or (c.q_strand == "+" and q_off > c.q_end):
            self.fail(
                f"chain 第 {lineno} 行（chain {c.chain_id}）："
                f"累计 gap 使对齐块起点超出头部声明的范围"
            )
            return
        if c.q_strand == "-" and c.q_end - self._q_consumed < c.q_start:
            self.fail(
                f"chain 第 {lineno} 行（chain {c.chain_id}）："
                f"累计 gap 使对齐块起点超出头部声明的范围（query 负链）"
            )
            return
        if c.q_strand == "+":
            qs, qe = q_off, q_off + size
        else:
            # query 负链：从原链末端向起点消耗
            qs, qe = c.q_end - self._q_consumed - size, c.q_end - self._q_consumed
        block = Block(
            t_start=t_off,
            t_end=t_off + size,
            q_start=qs,
            q_end=qe,
            same=(c.q_strand == "+"),
        )
        # 块必须落在头部声明的范围内
        t_ok = c.t_start <= block.t_start and block.t_end <= c.t_end
        q_ok = c.q_start <= block.q_start and block.q_end <= c.q_end
        if not t_ok or not q_ok:
            self.fail(
                f"chain 第 {lineno} 行（chain {c.chain_id}）："
                f"累计块超出头部声明的范围，请检查 size/gap"
            )
            return
        self._cur_lines.append(block)
        self._t_consumed += size + dt
        self._q_consumed += size + dq
        if len(toks) == 1:
            self._terminated = True

    def _finish_chain(self, next_header_lineno: Optional[int]) -> None:
        c = self._cur
        if c is None:
            return
        where = f"；下一头部在第 {next_header_lineno} 行" if next_header_lineno else "（文件结束）"
        if not self._cur_lines:
            self.fail(f"chain 第 {self._header_line} 行（chain {c.chain_id}）：缺少任何对齐块")
        else:
            last = self._cur_lines[-1]
            first = self._cur_lines[0]
            # 块序列必须恰好覆盖头部声明的范围（末行 size 后不允许再有 gap）：
            # target 恒为正链，故首块起于 tStart、末块止于 tEnd；
            # query 正链同理，负链时首块（先消耗的）终止于 qEnd、末块起于 qStart。
            t_ok = first.t_start == c.t_start and last.t_end == c.t_end
            if c.q_strand == "+":
                q_ok = first.q_start == c.q_start and last.q_end == c.q_end
            else:
                q_ok = first.q_end == c.q_end and last.q_start == c.q_start
            if not t_ok or not q_ok:
                self.fail(
                    f"chain 第 {self._header_line} 行（chain {c.chain_id}）："
                    f"块/gap 累计范围与头部不符"
                    f"（实际 target {first.t_start}..{last.t_end}，"
                    f"query {last.q_start}..{first.q_end}；"
                    f"头部 target {c.t_start}..{c.t_end}，query {c.q_start}..{c.q_end}）{where}"
                )
            # 同一条 chain 内块不得在任一侧相互重叠
            t_sorted = sorted(self._cur_lines, key=lambda b: (b.t_start, b.t_end))
            q_sorted = sorted(self._cur_lines, key=lambda b: (b.q_start, b.q_end))
            for seq, label in ((t_sorted, "target"), (q_sorted, "query")):
                for a, b in zip(seq, seq[1:]):
                    a_end = a.t_end if label == "target" else a.q_end
                    b_start = b.t_start if label == "target" else b.q_start
                    if a_end > b_start:
                        self.fail(
                            f"chain 第 {self._header_line} 行（chain {c.chain_id}）："
                            f"{label} 侧对齐块之间存在重叠（{a_end} > {b_start}）"
                        )
                        break
            c.blocks = self._cur_lines
            self.chains.append(c)
        self._cur = None
        self._cur_lines = []


def _safe_int(tok: str, lineno: int, what: str, issues: list[str]) -> int:
    try:
        return int(tok)
    except ValueError:
        issues.append(f"chain 第 {lineno} 行：{what} 不是整数（得到 {tok!r}）")
        return 0


def parse_chain(text: str) -> list[Chain]:
    return _ChainParser(text).run()


# --------------------------------------------------------------------------- #
# BED 解析
# --------------------------------------------------------------------------- #


@dataclass
class BedRow:
    lineno: int  # 物理行号（从 1 开始）
    chrom: str
    start: int
    end: int
    name: str = "."
    strand: str = "+"
    raw: str = ""


def parse_bed(text: str) -> list[BedRow]:
    rows: list[BedRow] = []
    issues: list[str] = []
    for i, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        first = line.split(maxsplit=1)[0]
        if first in ("track", "browser") or line.startswith("#"):
            continue
        toks = line.split()
        if len(toks) < 3:
            issues.append(f"BED 第 {i} 行：至少需要 chrom/start/end 三个字段，实际 {len(toks)} 个")
            continue
        chrom = toks[0]
        try:
            start = int(toks[1])
            end = int(toks[2])
        except ValueError:
            issues.append(
                f"BED 第 {i} 行：start/end 不是整数（得到 {toks[1]!r}/{toks[2]!r}）"
            )
            continue
        if start < 0 or end < 0:
            issues.append(f"BED 第 {i} 行：坐标不能为负（start={start}, end={end}）")
            continue
        if start >= end:
            issues.append(f"BED 第 {i} 行：空区间或反向区间（start={start}, end={end}），要求 start < end")
            continue
        name = toks[3] if len(toks) >= 4 else "."
        strand = "+"
        if len(toks) >= 6:
            strand = toks[5]
            if strand not in ("+", "-", "."):
                issues.append(f"BED 第 {i} 行：strand 只能为 +、- 或 .，得到 {strand!r}")
                continue
        rows.append(BedRow(i, chrom, start, end, name, strand, line))
    if issues:
        head = issues[:MAX_ISSUES]
        if len(issues) > MAX_ISSUES:
            head.append(f"……另有 {len(issues) - MAX_ISSUES} 个错误未展开")
        raise LiftError("\n".join(head))
    return rows


# --------------------------------------------------------------------------- #
# 迁移
# --------------------------------------------------------------------------- #


@dataclass
class Fragment:
    src_start: int
    src_end: int
    dst_chrom: str
    dst_start: int
    dst_end: int
    dst_strand: str


class _Index:
    """按染色体聚合对齐块，供区间相交查询。"""

    def __init__(self, chains: list[Chain], side: str):
        self.entries: dict[str, list[tuple[Block, Chain]]] = {}
        for c in chains:
            self.entries.setdefault(c.side_chrom(side), []).extend((b, c) for b in c.blocks)
        for key in self.entries:
            self.entries[key].sort(key=lambda e: (e[0].t_start if side == TARGET else e[0].q_start))
        self.side = side

    def candidates(self, chrom: str, start: int, end: int) -> list[tuple[Block, Chain]]:
        out = []
        for b, c in self.entries.get(chrom, ()):
            bs = b.t_start if self.side == TARGET else b.q_start
            be = b.t_end if self.side == TARGET else b.q_end
            if be <= start:
                continue
            if bs >= end:
                break
            out.append((b, c))
        return out


def _map_block(side: str, b: Block, lo: int, hi: int) -> tuple[int, int]:
    """块内 [lo, hi)（源侧坐标）映射到对侧坐标（零基左闭右开）。"""
    if side == TARGET:
        off = lo - b.t_start
        span = hi - lo
        if b.same:
            return b.q_start + off, b.q_start + off + span
        return b.q_end - off - span, b.q_end - off
    off = lo - b.q_start
    span = hi - lo
    return b.t_start + off, b.t_start + off + span


def _out_strand(src_strand: str, same: bool) -> str:
    if src_strand == ".":
        return "."
    return src_strand if same else ("-" if src_strand == "+" else "+")


@dataclass
class ResultRow:
    source_line: int
    source_chrom: str
    source_start: int
    source_end: int
    source_strand: str
    status: str  # mapped / partial / ambiguous / unmapped
    frag_index: int
    frag_count: int
    dest_chrom: str
    dest_start: str
    dest_end: str
    dest_strand: str
    note: str

    def as_dict(self) -> dict:
        return {
            "source_line": self.source_line,
            "source_chrom": self.source_chrom,
            "source_start": self.source_start,
            "source_end": self.source_end,
            "source_strand": self.source_strand,
            "status": self.status,
            "frag_index": self.frag_index,
            "frag_count": self.frag_count,
            "dest_chrom": self.dest_chrom,
            "dest_start": self.dest_start,
            "dest_end": self.dest_end,
            "dest_strand": self.dest_strand,
            "note": self.note,
        }


def run_liftover(
    bed_bytes: bytes,
    chain_bytes: bytes,
    side: str,
) -> dict:
    if side not in SIDES:
        raise LiftError(f"源侧必须为 {TARGET} 或 {QUERY}")
    if len(bed_bytes) > MAX_BYTES:
        raise LiftError("BED 文件超过 50 MiB 上限")
    if len(chain_bytes) > MAX_BYTES:
        raise LiftError("chain 文件超过 50 MiB 上限")

    def decode(data: bytes, label: str) -> str:
        if data.startswith(b"\xef\xbb\xbf"):
            data = data[3:]
        try:
            return data.decode("utf-8")
        except UnicodeDecodeError as e:
            raise LiftError(f"{label} 不是合法 UTF-8：第 {e.start} 字节附近解码失败")

    bed_text = decode(bed_bytes, "BED 文件")
    chain_text = decode(chain_bytes, "chain 文件")

    bed_rows = parse_bed(bed_text)
    chains = parse_chain(chain_text)  # 有错误时直接抛出，阻止整次导出
    if not chains:
        raise LiftError("chain 文件不含任何 chain 记录")

    index = _Index(chains, side)

    results: list[ResultRow] = []
    stats = {"mapped": 0, "partial": 0, "ambiguous": 0, "unmapped": 0}

    for row in bed_rows:
        cands = index.candidates(row.chrom, row.start, row.end)
        if not cands:
            chroms = sorted(index.entries.keys())
            if row.chrom not in chroms:
                reason = f"无映射：chain 中不存在源染色体 {row.chrom}"
            else:
                reason = "无映射：区间与任何对齐块都不相交（可能全部落在缺口上）"
            stats["unmapped"] += 1
            results.append(
                ResultRow(
                    row.lineno, row.chrom, row.start, row.end, row.strand,
                    "unmapped", 0, 0, "", "", "", "", reason,
                )
            )
            continue

        # 按 chain 归组，一个 chain 贡献一组有序片段
        groups: dict[str, list[Fragment]] = {}
        chain_of: dict[str, Chain] = {}
        mapped_bp: dict[str, int] = {}
        for b, c in cands:
            lo = max(row.start, b.t_start if side == TARGET else b.q_start)
            hi = min(row.end, b.t_end if side == TARGET else b.q_end)
            if lo >= hi:
                continue
            d_lo, d_hi = _map_block(side, b, lo, hi)
            groups.setdefault(c.chain_id, []).append(
                Fragment(
                    lo, hi,
                    c.side_chrom(QUERY if side == TARGET else TARGET),
                    d_lo, d_hi,
                    _out_strand(row.strand, b.same),
                )
            )
            chain_of[c.chain_id] = c
            mapped_bp[c.chain_id] = mapped_bp.get(c.chain_id, 0) + (hi - lo)

        # 块按源坐标本就有序；query 负链时目标坐标严格递减，这里保持源序输出
        for fl in groups.values():
            fl.sort(key=lambda f: (f.src_start, f.src_end))

        total = row.end - row.start
        if len(groups) == 1:
            cid = next(iter(groups))
            frags = groups[cid]
            bp = mapped_bp[cid]
            if bp < total:
                status = "partial"
                stats["partial"] += 1
                note = f"部分映射：{bp}/{total} bp 命中，其余跨缺口未映射"
            else:
                status = "mapped"
                stats["mapped"] += 1
                note = "完整映射" if len(frags) == 1 else f"完整映射：跨缺口拆分为 {len(frags)} 个有序片段"
            for k, f in enumerate(frags, start=1):
                results.append(
                    ResultRow(
                        row.lineno, row.chrom, row.start, row.end, row.strand,
                        status, k, len(frags),
                        f.dst_chrom, str(f.dst_start), str(f.dst_end), f.dst_strand,
                        note,
                    )
                )
        else:
            stats["ambiguous"] += 1
            total_frags = sum(len(v) for v in groups.values())
            detail = "；".join(
                f"chain {cid} 命中 {mapped_bp[cid]}/{total} bp"
                for cid in sorted(groups)
            )
            for cid in sorted(groups):
                for f in groups[cid]:
                    results.append(
                        ResultRow(
                            row.lineno, row.chrom, row.start, row.end, row.strand,
                            "ambiguous", 0, total_frags,
                            f.dst_chrom, str(f.dst_start), str(f.dst_end), f.dst_strand,
                            f"歧义：命中 {len(groups)} 个目标位置（{detail}），未自动挑选",
                        )
                    )

    return {
        "side": side,
        "bed_rows": len(bed_rows),
        "result_rows": len(results),
        "stats": stats,
        "rows": [r.as_dict() for r in results],
    }


# --------------------------------------------------------------------------- #
# 结果导出（TSV）
# --------------------------------------------------------------------------- #

TSV_COLUMNS = [
    ("source_line", "source_line"),
    ("source_chrom", "source_chrom"),
    ("source_start", "source_start"),
    ("source_end", "source_end"),
    ("source_strand", "source_strand"),
    ("status", "status"),
    ("frag_index", "frag_index"),
    ("frag_count", "frag_count"),
    ("dest_chrom", "dest_chrom"),
    ("dest_start", "dest_start"),
    ("dest_end", "dest_end"),
    ("dest_strand", "dest_strand"),
    ("note", "note"),
]


def _tsv_escape(value: str) -> str:
    return value.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n").replace("\r", "\\r")


def build_tsv(payload: dict) -> str:
    header = "\t".join(h for _, h in TSV_COLUMNS)
    lines = [header]
    for r in payload["rows"]:
        lines.append("\t".join(_tsv_escape(str(r[k])) for k, _ in TSV_COLUMNS))
    return "\n".join(lines) + "\n"
