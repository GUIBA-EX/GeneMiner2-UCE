# GeneMiner2-UCE assembler 说明

本文只说明当前 CLI 支持的组装路径：`reference` 与 `uce`。marker profiling 是独立的读段定量流程，不使用 assembler，见命令行手册。

## 实现与选择

| 使用场景 | 默认后端 | 可选后端 | 适用目的 |
|---|---|---|---|
| `--assembly-mode original` | `original-rust` | `original`、`uce-rust` | exon、SCO 及核/线粒体 marker 的参考引导恢复 |
| `--assembly-mode uce` | `uce-rust` | 无 | 从 genome skimming 或 target capture 恢复 UCE core 和有 read 支持的 flank |

- `original-rust` 是 `original` 默认后端：单线程 Rust 实现，保留上游 GeneMiner2 的参考引导 seed、候选组合与回溯逻辑。
- `original` 是固定版本的上游 Python 后端，只用于严格结果对照。
- `uce-rust` 是 UCE 定向 Rust 后端，也是 `uce` 模式唯一可用的后端；程序不可用或失败时会直接报错。

## 共用原则

三个后端都以参考序列定位 locus，以过滤后的 reads 支持或否定候选 contig：

1. 参考和 reads 被切分为 k-mer，低频 read k-mer 可被过滤；
2. 参考位置与方向一致的 k-mer 优先作为 seed；
3. 从 seed 双向延伸、合并为候选 contig；
4. 使用 reads 的高质量内部 slice 验证候选，而不只依赖图深度。

因此，参考负责界定起点，reads 决定结果是否可信。

## Original：保留 Python 原版作为对照

`original` 后端固定在上游 GeneMiner2 提交 `36e06feeb99654bdb87f45d4cde225d8c3e311d0` 的 Python 原版。它仍随本项目发布，并可用 `--assembler-implementation original` 显式运行；这是结果对照的基准，不是废弃代码。

Python 原版以 read depth 与参考权重给图节点排序。每次延伸先走当前权重最高的邻接 k-mer，其余分叉压入栈；路径走入死路后回退到最近分叉，继续尝试栈中候选，直至达到 `--search-depth`。每个方向按 read-slice 支持保留最多三个候选，再组合、评分并要求足够的 read 支持覆盖 contig。

这种回溯枚举对 exon、SCO 及核/线粒体 marker 等较短的参考引导目标很稳妥，但在 UCE flank、低复杂度区域或多气泡图中会增加搜索和候选组合。因此 `original + auto` 使用逻辑兼容的 `original-rust`，而不是移除 Python 原版。若出现结果差异、需要复现实验基线或准备与上游比较，应优先运行 Python 原版。

`original-rust` 支持 `--reuse-reference-cache`。它按参考、k 和实现版本验证二进制 k-mer cache；缓存不匹配、损坏或截断时会自动重建。

## UCE：受控延伸优先

UCE 模式默认的 `--uce-path-strategy backbone` 在无分叉区域连续延伸；遇到气泡时，对每条出边做最多 `--uce-backbone-lookahead` 步的有限前瞻，依次比较：

1. 可继续延伸的长度；
2. 前瞻路径累计的 k-mer 与参考权重；
3. 起始边权重。

选中的分支进入主干，其他分支随即丢弃；算法不保存回溯栈，也不重复访问 k-mer。这牺牲了枚举全部变体路径的能力，换取可控的运行时间，并降低由重复 flank、杂合或噪声气泡组合出异常长 contig 的风险。

`--uce-path-strategy search` 保留较旧的分支枚举行为，仅建议用于与 backbone 做 A/B 比较。UCE 输出 `uce_assembly_summary.csv`，提供 read-supported span、breadth、最大缺口、unique-read density、k-mer depth CV 与 max/median depth 等质量指标；可通过 `--assembler-graph-format gfa|dot|both` 导出诊断图。

## Rust 实现要点

Rust 后端以 2-bit `u128` 精确表示 k-mer（支持 `k ≤ 63`），不会因概率哈希碰撞合并不同序列。reads 按 `--assembler-read-chunk-size` 分块处理，先在块内排序和计数，再汇总至 locus 图，以限制内存并减少随机字典更新。

这些工程实现借鉴了短读长 de Bruijn graph assembler 的常见做法；本项目不调用 MaSuRCA、SPAdes 或 Sparrowhawk 的代码，也未实现它们的全基因组、多 k、scaffolding 或概率哈希工作流。

## 实际使用与验证

- 常规基因恢复：使用默认 `reference`；需要严格比较时改用 `--assembler-implementation original`。
- UCE：使用 `--assembly-mode uce`，并检查 `uce_assembly_summary.csv` 与下游比对；放宽边界或 k-mer 参数后尤其应检查长 contig 的 read density、depth CV 和 flank balance。
- 评估 backbone：在同一数据上与 `--uce-path-strategy search` 比较恢复率、contig 长度、read-supported breadth、异常长序列和下游 SNP 缺失率。不能只凭运行时间或 contig 更长判断更好。

## 代码位置

- `original`（Python）：`scripts/main_assembler_original.py`
- `original-rust`：`rust/main_assembler_original/src/main.rs`
- `uce-rust`：`rust/main_assembler/src/main.rs`
- UCE 路径和候选验收：`rust/main_assembler/src/pipeline.rs`
- 图与 I/O：`rust/main_assembler/src/unitig.rs`、`rust/main_assembler/src/io_utils.rs`

## 参考

- [MaSuRCA](https://pmc.ncbi.nlm.nih.gov/articles/PMC3799473/)
- [SPAdes](https://pmc.ncbi.nlm.nih.gov/articles/PMC3342519/)
- [Sparrowhawk](https://github.com/bacpop/sparrowhawk)
- [GeneMiner2 upstream assembler](https://github.com/sculab/GeneMiner2/blob/36e06feeb99654bdb87f45d4cde225d8c3e311d0/scripts/main_assembler.py)
