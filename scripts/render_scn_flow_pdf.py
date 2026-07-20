#!/usr/bin/env python3
"""Render the SCN architecture as a single-page vector PDF without dependencies."""

from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "docs" / "scn_workflow.pdf"
W, H = 842, 595  # A4 landscape, points


def pdf_string(value: str) -> str:
    return "<FEFF" + value.encode("utf-16-be").hex().upper() + ">"


parts: list[str] = []


def line(x1, y1, x2, y2, width=0.8, color="0.30 0.35 0.42"):
    parts.append(f"{color} RG {width} w {x1} {y1} m {x2} {y2} l S")


def arrow(x1, y1, x2, y2):
    line(x1, y1, x2, y2)
    if x2 >= x1:
        line(x2, y2, x2 - 5, y2 + 3)
        line(x2, y2, x2 - 5, y2 - 3)
    else:
        line(x2, y2, x2 + 5, y2 + 3)
        line(x2, y2, x2 + 5, y2 - 3)


def text(x, y, value, size=9, color="0.10 0.12 0.16"):
    parts.append(f"{color} rg BT /F1 {size} Tf {x} {y} Td {pdf_string(value)} Tj ET")


def box(x, y, width, height, title, body, fill="0.94 0.97 1.00"):
    parts.append(f"{fill} rg {x} {y} {width} {height} re f")
    parts.append(f"0.32 0.43 0.56 RG 0.8 w {x} {y} {width} {height} re S")
    text(x + 7, y + height - 14, title, 9.6, "0.05 0.18 0.32")
    for i, row in enumerate(body):
        text(x + 7, y + height - 27 - i * 10, row, 7.2, "0.18 0.22 0.27")


text(34, 560, "GeneMiner2-UCE  SCN：单拷贝核基因恢复与鉴定架构", 16, "0.04 0.16 0.30")
text(34, 542, "原则：样本内只报告 copy-state；跨样本与基因树才判定 pseudoSCO 或 resolved_1to1。", 8.5, "0.27 0.31 0.38")

# Main recovery row
box(34, 430, 132, 72, "1 参考包（离线）", ["多物种 CDS / 蛋白", "BUSCO / OrthoDB 先验", "UniProt 仅作功能注释"], "0.92 0.96 1.00")
box(198, 430, 132, 72, "2 Reads 招募", ["MainFilter：多 bait", "不修改 MainFilter", "按 orthogroup 归属"], "0.92 0.98 0.96")
box(362, 430, 132, 72, "3 精炼与组装", ["refilter：paired reads", "original-rust", "保留有限候选 contigs"], "1.00 0.97 0.90")
box(526, 430, 132, 72, "4 样本 copy-state", ["coverage + read-pair", "翻译一致性 + alternatives", "变异结构 / phased alleles"], "1.00 0.94 0.91")
box(690, 430, 118, 72, "5 样本结果", ["single_candidate", "allelic / multicopy", "collapsed / ambiguous"], "0.98 0.93 0.98")
for x in (166, 330, 494, 658):
    arrow(x + 3, 466, x + 29, 466)

# Cohort row
box(78, 285, 156, 72, "6 跨样本基础 QC", ["orthogroup 确认", "嵌合检查、codon-aware MSA", "局部异常位点 masking"], "0.92 0.96 1.00")
box(290, 285, 156, 72, "7 pseudoSCO", ["single_candidate + 可信 allele", "按 occupancy 选择", "不强制所有样本交集"], "0.92 0.98 0.96")
box(502, 285, 156, 72, "8 问题基因树解析", ["只处理 multicopy / ambiguous", "species-overlap 树拆分", "OrthoSNAP / UPhO 思路"], "1.00 0.97 0.90")
box(690, 285, 118, 72, "9 严格结果", ["resolved_1to1", "multicopy families", "保留审计证据"], "0.98 0.93 0.98")
arrow(749, 430, 749, 358)
arrow(234, 321, 290, 321)
arrow(446, 321, 502, 321)
arrow(658, 321, 690, 321)

# Outputs
text(34, 242, "输出面板", 11, "0.04 0.16 0.30")
box(34, 155, 178, 64, "样本审计", ["copy_state.tsv", "候选 contigs 与 reads 支持", "不把未检出 paralog 当作 SCO"], "0.95 0.96 0.98")
box(238, 155, 178, 64, "快速系统发育", ["core_pseudoSCO", "occupancy_pseudoSCO", "保留缺失与判定记录"], "0.92 0.98 0.96")
box(442, 155, 178, 64, "严格系统发育", ["resolved_1to1", "树支持与拆分记录", "适合严格物种树分析"], "1.00 0.97 0.90")
box(646, 155, 162, 64, "多拷贝路线", ["multicopy_families", "DISCO / ASTRAL-Pro", "不浪费复杂家族"], "0.98 0.93 0.98")
arrow(368, 285, 327, 219)
arrow(560, 285, 531, 219)
arrow(749, 285, 727, 219)

text(34, 105, "判定层级：", 9, "0.04 0.16 0.30")
text(103, 105, "reference_SCO  →  single_candidate  →  pseudoSCO  →  resolved_1to1", 9, "0.13 0.18 0.24")
text(34, 82, "默认模式：仅对有冲突的 locus 建 gene tree；复杂阈值由 conservative / balanced / sensitive preset 管理。", 8, "0.27 0.31 0.38")
text(34, 58, "模块边界：CLI 调度与缓存 ｜ MainFilter/refilter 保持现状 ｜ original-rust 候选组装 ｜ Rust scn_workflow 判定与队列分析。", 8, "0.27 0.31 0.38")

stream = "\n".join(parts) + "\n"
objects = [
    "<< /Type /Catalog /Pages 2 0 R >>",
    "<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
    f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {W} {H}] /Resources << /Font << /F1 5 0 R >> >> /Contents 4 0 R >>",
    f"<< /Length {len(stream.encode('ascii'))} >>\nstream\n{stream}endstream",
    "<< /Type /Font /Subtype /Type0 /BaseFont /STSong-Light /Encoding /UniGB-UCS2-H /DescendantFonts [<< /Type /Font /Subtype /CIDFontType0 /BaseFont /STSong-Light /CIDSystemInfo << /Registry (Adobe) /Ordering (GB1) /Supplement 2 >> >>] >>",
]

pdf = "%PDF-1.4\n%\xe2\xe3\xcf\xd3\n"
offsets = [0]
for number, obj in enumerate(objects, 1):
    offsets.append(len(pdf.encode("latin-1")))
    pdf += f"{number} 0 obj\n{obj}\nendobj\n"
xref = len(pdf.encode("latin-1"))
pdf += f"xref\n0 {len(objects) + 1}\n0000000000 65535 f \n"
for offset in offsets[1:]:
    pdf += f"{offset:010d} 00000 n \n"
pdf += f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\nstartxref\n{xref}\n%%EOF\n"

OUT.write_bytes(pdf.encode("latin-1"))
print(OUT)
