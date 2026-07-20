#!/usr/bin/env python3
"""Render the SCN architecture as a UTF-8 PNG with a bundled system CJK font."""

from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "scn_workflow.png"
FONT = "/usr/share/fonts/opentype/noto/NotoSansCJK-Medium.ttc"
BOLD = "/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc"
W, H = 2400, 1700

image = Image.new("RGB", (W, H), "#FAFCFF")
draw = ImageDraw.Draw(image)


def font(size, bold=False):
    return ImageFont.truetype(BOLD if bold else FONT, size, index=2)


def centered(x, y, value, fnt, fill="#182432"):
    box = draw.textbbox((0, 0), value, font=fnt)
    draw.text((x - (box[2] - box[0]) / 2, y), value, font=fnt, fill=fill)


def wrapped_lines(value, fnt, width):
    rows, row = [], ""
    for char in value:
        candidate = row + char
        if draw.textlength(candidate, font=fnt) > width and row:
            rows.append(row)
            row = char
        else:
            row = candidate
    if row:
        rows.append(row)
    return rows


def box(x, y, width, height, title, rows, fill):
    draw.rounded_rectangle((x, y, x + width, y + height), radius=20, fill=fill, outline="#4A6480", width=3)
    draw.text((x + 24, y + 20), title, font=font(30, True), fill="#103451")
    y_pos = y + 70
    body = font(21)
    for value in rows:
        for row in wrapped_lines(value, body, width - 48):
            draw.text((x + 24, y_pos), row, font=body, fill="#263440")
            y_pos += 28
        y_pos += 4


def arrow(x1, y1, x2, y2):
    draw.line((x1, y1, x2, y2), fill="#53677B", width=4)
    if abs(x2 - x1) >= abs(y2 - y1):
        if x2 >= x1:
            points = [(x2, y2), (x2 - 18, y2 - 10), (x2 - 18, y2 + 10)]
        else:
            points = [(x2, y2), (x2 + 18, y2 - 10), (x2 + 18, y2 + 10)]
    elif y2 >= y1:
        points = [(x2, y2), (x2 - 10, y2 - 18), (x2 + 10, y2 - 18)]
    else:
        points = [(x2, y2), (x2 - 10, y2 + 18), (x2 + 10, y2 + 18)]
    draw.polygon(points, fill="#53677B")


draw.text((85, 55), "GeneMiner2-UCE  SCN：单拷贝核基因恢复与鉴定架构", font=font(50, True), fill="#0B2945")
draw.text((88, 122), "原则：样本内只报告 copy-state；跨样本与基因树才判定 pseudoSCO 或 resolved_1to1。", font=font(27), fill="#4A5564")

main_y, main_h, main_w = 280, 220, 390
main_x = [70, 550, 1030, 1510, 1990]
main = [
    ("1 参考包（离线）", ["多物种 CDS / 蛋白", "BUSCO / OrthoDB 先验", "UniProt 仅作功能注释"], "#EAF3FF"),
    ("2 Reads 招募", ["MainFilter：多 bait", "不修改 MainFilter", "按 orthogroup 归属"], "#E9FAF0"),
    ("3 精炼与组装", ["refilter：paired reads", "original-rust", "保留有限候选 contigs"], "#FFF4DF"),
    ("4 样本 copy-state", ["coverage + read-pair", "翻译一致性 + alternatives", "变异结构 / phased alleles"], "#FDECE8"),
    ("5 样本结果", ["single_candidate", "allelic / multicopy", "collapsed / ambiguous"], "#F7ECFA"),
]
for x, (title, rows, fill) in zip(main_x, main):
    box(x, main_y, main_w, main_h, title, rows, fill)
for i in range(4):
    arrow(main_x[i] + main_w + 8, main_y + main_h / 2, main_x[i + 1] - 10, main_y + main_h / 2)

cohort_y, cohort_h, cohort_w = 740, 220, 450
cohort_x = [170, 700, 1230, 1760]
cohort = [
    ("6 跨样本基础 QC", ["orthogroup 确认", "嵌合检查、codon-aware MSA", "局部异常位点 masking"], "#EAF3FF"),
    ("7 pseudoSCO", ["single_candidate + 可信 allele", "按 occupancy 选择", "不强制所有样本交集"], "#E9FAF0"),
    ("8 问题基因树解析", ["仅 multicopy / ambiguous", "species-overlap 树拆分", "OrthoSNAP / UPhO 思路"], "#FFF4DF"),
    ("9 严格结果", ["resolved_1to1", "multicopy families", "保留完整审计证据"], "#F7ECFA"),
]
for x, (title, rows, fill) in zip(cohort_x, cohort):
    box(x, cohort_y, cohort_w, cohort_h, title, rows, fill)
arrow(main_x[-1] + main_w / 2, main_y + main_h + 8, cohort_x[-1] + cohort_w / 2, cohort_y - 8)
for i in range(3):
    arrow(cohort_x[i] + cohort_w + 10, cohort_y + cohort_h / 2, cohort_x[i + 1] - 10, cohort_y + cohort_h / 2)

draw.text((85, 1040), "输出面板", font=font(34, True), fill="#0B2945")
out_y, out_h, out_w = 1120, 205, 500
out_x = [70, 650, 1230, 1810]
outputs = [
    ("样本审计", ["copy_state.tsv", "候选 contigs 与 reads 支持", "不把未检出 paralog 当 SCO"], "#F2F5F9"),
    ("快速系统发育", ["core_pseudoSCO", "occupancy_pseudoSCO", "保留缺失与判定记录"], "#E9FAF0"),
    ("严格系统发育", ["resolved_1to1", "树支持与拆分记录", "适合严格物种树分析"], "#FFF4DF"),
    ("多拷贝路线", ["multicopy_families", "DISCO / ASTRAL-Pro", "不浪费复杂家族"], "#F7ECFA"),
]
for x, (title, rows, fill) in zip(out_x, outputs):
    box(x, out_y, out_w, out_h, title, rows, fill)
arrow(cohort_x[1] + cohort_w / 2, cohort_y + cohort_h + 8, out_x[1] + out_w / 2, out_y - 8)
arrow(cohort_x[2] + cohort_w / 2, cohort_y + cohort_h + 8, out_x[2] + out_w / 2, out_y - 8)
arrow(cohort_x[3] + cohort_w / 2, cohort_y + cohort_h + 8, out_x[3] + out_w / 2, out_y - 8)

draw.line((85, 1420, 2315, 1420), fill="#C6D1DD", width=2)
draw.text((85, 1450), "判定层级：reference_SCO  →  single_candidate  →  pseudoSCO  →  resolved_1to1", font=font(28, True), fill="#16334D")
draw.text((85, 1515), "默认：仅对有冲突的 locus 建 gene tree；复杂阈值由 conservative / balanced / sensitive preset 管理。", font=font(23), fill="#4A5564")
draw.text((85, 1560), "模块边界：CLI 调度与缓存 ｜ MainFilter/refilter 保持现状 ｜ original-rust 候选组装 ｜ Rust scn_workflow 判定与队列分析。", font=font(23), fill="#4A5564")

image.save(OUT, "PNG", optimize=True)
print(OUT)
