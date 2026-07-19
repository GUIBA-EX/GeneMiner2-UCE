# 2. Assembler

[English version](assembler_EN.md)

Assembler 将招募后的逐 locus reads 组装为 contig。应按生物学目标选择流程，而不是按测序平台或文库名称选择。

| 工作流 | 适用目标 | 默认后端 | 结果 |
|---|---|---|---|
| `--assembly-mode original` | exon、SCO、核 marker、线粒体 marker | `original-rust` | 参考引导 contig，默认按参考裁切 |
| `--assembly-mode uce` | genome skimming 或 target capture 中的 UCE | `uce-rust` | UCE core 及有 read 支持的 flank |

## 后端

`original-rust` 是 `original` 的确定性 Rust 兼容后端。只有需要与固定的上游 GeneMiner2 Python 实现严格对照时才使用 `--assembler-implementation original`。`uce-rust` 是 `uce` 的唯一后端，不会回退到 Python。

所有后端均使用参考位置明确的 seed、read k-mer 支持、双向延伸和 read-slice 验证。

## Original 工作流

原版算法优先延伸权重最高的边，将其他分支压栈以便回溯。每侧最多保留三条候选，再组合和评分。该行为适用于 exon 与细胞器 marker 等较短的参考引导目标。

```bash
cli/geneminer2 -f samples.tsv -r references -o output -p 8 \
  --assembly-mode original
```

`original-rust` 可用 `--reuse-reference-cache` 复用带版本的二进制参考 k-mer cache。

## UCE 工作流与 QC

UCE 模式在 refilter 中保留 linked mates，且默认跳过裁切。`backbone` 策略在每个气泡作一次有限前瞻决策，提交胜出分支而不回溯。`--uce-rescue-reads` 以已接受的首轮 contig 再做一次 raw-read 招募和组装；rescue 失败或密度降低时会回退。

应检查 `uce_assembly_summary.csv`、`uce_rescue_summary.csv`、`contigs_all_low/` 和下游比对。主要护栏包括 unique-read density、支持 breadth、最大无支持缺口、k-mer depth CV 与 max/median depth ratio。`--uce-path-strategy search` 只应用于与旧分支枚举的敏感性比较。

文件字段见[输出说明](../manual/ZH_CN/output.md)，参数见[命令行指南](../manual/ZH_CN/command_line.md)。
