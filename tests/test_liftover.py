"""引擎正确性测试：pytest tests/test_liftover.py，或 python tests/test_liftover.py"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.liftover import LiftError, build_tsv, parse_chain, run_liftover  # noqa: E402

# 正链 chain：chr1(target) -> chrA(query)
# 数据行: size [dt dq]（dt/dq 为本块之后的缺口，末行只有 size）
#   行1: size30, 其后 t gap 10 / q gap 10
#   行2: size30, 其后无缺口
#   行3: size30（末块）
# block1 t0-30   <-> q10-40
# block2 t40-70  <-> q50-80
# block3 t70-100 <-> q80-110（target 侧 block2/3 相邻）
CHAIN_PLUS = """chain 5000 chr1 100 + 0 100 chrA 120 + 10 110 1
30	10	10
30	0	0
30

"""

# 负链 chain：chr2(target) -> chrB(query, -)
#   行1: size30, 其后 t gap 10 / q gap 10
#   行2: size40（末块）
# block1 t10-40 <-> q 原坐标 80..110（升序存储，same=False）
# block2 t50-90 <-> q 原坐标 30..70
CHAIN_MINUS = """chain 4000 chr2 100 + 10 90 chrB 120 - 30 110 2
30	10	10
40
"""

CHAIN_DUP = CHAIN_PLUS + "\n" + CHAIN_PLUS.replace(" 1\n", " 3\n")


def lift(bed_text, chain_text=CHAIN_PLUS, side="target"):
    return run_liftover(bed_text.encode(), chain_text.encode(), side)


def statuses(payload):
    return [r["status"] for r in payload["rows"]]


def test_chain_structure_plus():
    chains = parse_chain(CHAIN_PLUS)
    assert len(chains) == 1
    c = chains[0]
    assert c.chain_id == "1"
    blocks = c.blocks
    assert (blocks[0].t_start, blocks[0].t_end) == (0, 30)
    assert (blocks[0].q_start, blocks[0].q_end) == (10, 40)
    assert (blocks[1].t_start, blocks[1].t_end) == (40, 70)
    assert (blocks[1].q_start, blocks[1].q_end) == (50, 80)
    assert (blocks[2].t_start, blocks[2].t_end) == (70, 100)
    assert (blocks[2].q_start, blocks[2].q_end) == (80, 110)
    assert all(b.same for b in blocks)


def test_chain_structure_minus():
    c = parse_chain(CHAIN_MINUS)[0]
    b0, b1 = c.blocks
    assert (b0.t_start, b0.t_end) == (10, 40)
    # 负链第一块对应 query 原坐标 80..110（存为升序），same=False
    assert (b0.q_start, b0.q_end) == (80, 110)
    assert b0.same is False
    assert (b1.t_start, b1.t_end) == (50, 90)
    assert (b1.q_start, b1.q_end) == (30, 70)
    assert b1.same is False


def test_simple_plus_mapping():
    p = lift("chr1\t5\t15\n")
    r = p["rows"][0]
    assert r["status"] == "mapped"
    assert (r["dest_chrom"], r["dest_start"], r["dest_end"], r["dest_strand"]) == (
        "chrA", "15", "25", "+",
    )
    assert p["stats"]["mapped"] == 1


def test_split_across_gap_is_complete_ordered():
    # 0..50：块1 0-30 映射 10-40；30-40 落缺口；40-50 在块2 映射 50-60
    p = lift("chr1\t0\t50\n")
    rows = p["rows"]
    assert statuses(p) == ["partial", "partial"]
    assert rows[0]["frag_index"] == 1 and rows[1]["frag_index"] == 2
    assert rows[0]["frag_count"] == 2
    assert (rows[0]["dest_start"], rows[0]["dest_end"]) == ("10", "40")
    assert (rows[1]["dest_start"], rows[1]["dest_end"]) == ("50", "60")
    assert p["stats"]["partial"] == 1  # 源区间只计一次

def test_all_in_gap_unmapped():
    p = lift("chr1\t31\t39\n")
    assert statuses(p) == ["unmapped"]
    assert "缺口" in p["rows"][0]["note"]
    assert p["rows"][0]["dest_chrom"] == ""


def test_unknown_chrom_reason_and_line_preserved():
    p = lift("chr9\t0\t10\n")
    r = p["rows"][0]
    assert r["status"] == "unmapped"
    assert "chr9" in r["note"]
    assert r["source_line"] == 1
    assert r["source_chrom"] == "chr9"


def test_ambiguous_two_chains_not_picked():
    # 两条 chain 都能映射 chr1:0-10
    p = lift("chr1\t0\t10\n", CHAIN_DUP)
    assert set(statuses(p)) == {"ambiguous"}
    assert len(p["rows"]) == 2
    assert all("歧义" in r["note"] and "未自动挑选" in r["note"] for r in p["rows"])
    assert p["stats"]["ambiguous"] == 1


def test_minus_chain_target_side():
    # chr2:10..20 -> query 原坐标：qEnd=110，off=0 => 100..110，方向 -
    p = run_liftover(b"chr2\t10\t20\tgene\t0\t+\n", CHAIN_MINUS.encode(), "target")
    r = p["rows"][0]
    assert r["status"] == "mapped"
    assert (r["dest_start"], r["dest_end"], r["dest_strand"]) == ("100", "110", "-")
    # 源 BED 本身是负链，穿过负链 chain 后方向翻回 +
    p2 = run_liftover(b"chr2\t10\t20\tg\t0\t-\n", CHAIN_MINUS.encode(), "target")
    assert p2["rows"][0]["dest_strand"] == "+"


def test_minus_chain_query_side_half_open():
    # query 侧源区间 chrB:90..100（升序坐标），落第一块 80..110 内：
    # off=10, span=10 -> target 20..30，方向 -
    p = run_liftover(b"chrB\t90\t100\n", CHAIN_MINUS.encode(), "query")
    r = p["rows"][0]
    assert r["status"] == "mapped"
    assert (r["dest_chrom"], r["dest_start"], r["dest_end"], r["dest_strand"]) == (
        "chr2", "20", "30", "-",
    )
    # 左闭右开边界：第一块紧贴边界 80..110 -> target 10..40
    p2 = run_liftover(b"chrB\t80\t110\n", CHAIN_MINUS.encode(), "query")
    assert (p2["rows"][0]["dest_start"], p2["rows"][0]["dest_end"]) == ("10", "40")


def test_strand_dot_preserved():
    p = lift("chr1\t5\t15\t.\t0\t.\n")
    assert p["rows"][0]["dest_strand"] == "."


def test_skip_directives_count_physical_lines():
    bed = "# comment\ntrack name=x\n\nchr1\t5\t15\n"
    p = lift(bed)
    assert p["rows"][0]["source_line"] == 4


def test_bad_bed_reports_line():
    for bad, needle in [
        ("chr1\t10\n", "BED 第 1 行"),
        ("chr1\tabc\t15\n", "第 1 行"),
        ("chr1\t-1\t10\n", "第 1 行"),
        ("chr1\t10\t10\n", "第 1 行"),
        ("chr1\t10\t20\tx\t0\t?\n", "strand"),
    ]:
        try:
            lift(bad)
        except LiftError as e:
            assert needle in str(e), (bad, str(e))
        else:
            raise AssertionError(f"应当报错: {bad}")


def test_bad_bad_line_does_not_block_other_lines_collected():
    try:
        lift("chr1\t0\t10\nchr1\t5\nchr1\t20\t30\n")
    except LiftError as e:
        msg = str(e)
        assert "BED 第 2 行" in msg
    else:
        raise AssertionError("应阻止整次导出")


def test_bad_chain_reports_positions():
    cases = [
        ("10\n", "孤立数据行"),
        ("chain 1 a 100 + 0 100 b 100 + 0 100\n", "13 个字段"),
        ("chain 1 a 100 + 0 100 b 100 + 0 100 9\n0\n", "对齐块大小"),
        ("chain 1 a 100 + 0 100 b 100 + 0 100 9\n60\t0\t0\n60\n", "超出头部"),
        ("chain 1 a 100 - 0 100 b 100 + 0 100 9\n100\n", "target 正链"),
        ("chain 1 a 100 + 0 100 b 100 x 0 100 9\n100\n", "query strand"),
        ("chain 1 a abc + 0 100 b 100 + 0 100 9\n100\n", "target size"),
    ]
    for text, needle in cases:
        try:
            parse_chain(text)
        except LiftError as e:
            assert needle in str(e), (text, needle, str(e))
        else:
            raise AssertionError(f"应当报错: {text!r}")


def test_chain_size_mismatch_detected():
    # 头部声明 t 0..100，但块只覆盖 50
    text = "chain 1 a 100 + 0 100 b 100 + 0 50 9\n50\n"
    try:
        parse_chain(text)
    except LiftError as e:
        assert "累计" in str(e)
    else:
        raise AssertionError("应检测出累计长度不符")


def test_chain_data_after_terminator_rejected():
    text = (
        "chain 1 a 100 + 0 100 b 100 + 0 100 9\n"
        "50\t0\t0\n50\n"
        "40\n"  # 终止行之后又来数据
    )
    try:
        parse_chain(text)
    except LiftError as e:
        assert "终止行" in str(e)
    else:
        raise AssertionError("应拒绝终止行后的数据行")


def test_chain_dangling_gap_rejected():
    # 最后一个数据行带 gap 却没有后续末块
    text = "chain 1 a 100 + 0 100 b 100 + 0 100 9\n60\t0\t0\n"
    try:
        parse_chain(text)
    except LiftError as e:
        assert "累计" in str(e)
    else:
        raise AssertionError("应拒绝悬空 gap")


def test_tsv_roundtrip_shape():
    p = lift("chr1\t0\t50\nchr9\t0\t3\n")
    tsv = build_tsv(p)
    lines = tsv.splitlines()
    assert lines[0].startswith("source_line\tsource_chrom")
    # 2 条源区间：第一条 2 片段 + 第二条 1 行
    assert len(lines) == 1 + 3
    assert "无映射" in lines[-1]


def test_chain_block_partial_boundary():
    # 65..90：块2 65..70 (5bp) 与块3 70..90 (20bp)，两侧相邻无 gap => 完整但拆分
    p = lift("chr1\t65\t90\n")
    assert p["stats"]["mapped"] == 1
    rows = p["rows"]
    assert len(rows) == 2
    assert (rows[0]["dest_start"], rows[0]["dest_end"]) == ("75", "80")
    assert (rows[1]["dest_start"], rows[1]["dest_end"]) == ("80", "100")


def test_split_adjacent_fragments_half_open_touch():
    # 相邻片段的目标坐标必须首尾相接，保证零基左闭右开无重叠无空洞
    p = lift("chr1\t65\t90\n")
    rows = p["rows"]
    assert rows[0]["dest_end"] == rows[1]["dest_start"]


def test_full_split_across_gap_marked_partial():
    # 25..55：块1 25-30 + 块2 40-55，中间 30-40 落缺口
    p = lift("chr1\t25\t55\n")
    assert p["stats"]["partial"] == 1
    rows = p["rows"]
    assert len(rows) == 2
    assert (rows[0]["dest_start"], rows[0]["dest_end"]) == ("35", "40")
    assert (rows[1]["dest_start"], rows[1]["dest_end"]) == ("50", "65")


if __name__ == "__main__":
    fns = [(n, f) for n, f in sorted(globals().items()) if n.startswith("test_")]
    failed = 0
    for n, f in fns:
        try:
            f()
            print(f"PASS {n}")
        except Exception as e:  # noqa: BLE001
            failed += 1
            print(f"FAIL {n}: {e}")
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    sys.exit(1 if failed else 0)
