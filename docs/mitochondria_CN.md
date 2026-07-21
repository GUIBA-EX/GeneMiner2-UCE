# 线粒体工作流：常规环形动物线粒体

`mito` 使用现有 MainFilter、refilter 与 Rust UCE assembler，从带注释的 GenBank 参考招募 reads，并只在样本 reads 支持时报告闭环。

| 从什么开始 | 得到什么 | 主要边界 |
| --- | --- | --- |
| 带注释 GenBank + 样本 reads | 有 read 支持的闭环或 partial 线粒体序列 | 仅适用于常规单条环形动物线粒体 |

## 适用范围

本流程只面向**常规的单条环形动物线粒体基因组**：基因顺序与参考大体可比，不存在明显多分子结构、大片段重排或需要分相报告的强异质性。它不适合植物或真菌复杂线粒体、分裂或多染色体线粒体、严重重排，或需要专门处理异质性与 NUMT 的研究。

## 最小命令

```bash
cli/geneminer2 mito \
  -f samples.tsv \
  -o mito_output \
  -p 8 \
  --mito-genbank mitochondrial_reference.gb
```

`--mito-max-reads 320` 限制自适应流程的最大输入；`--no-mito-adaptive-stop` 改为一次性过滤。

## 流程

```text
GenBank gene/rRNA/tRNA + genome/tile baits
→ MainFilter 招募 paired reads
→ 单一 mitochondrial read pool
→ refilter → Rust UCE 联合图组装
→ 非冗余、信息性 contigs_all 作为样本特异性 seed，从原始 paired reads 再招募一次并联合重组装
→ 严格 overlap、GFA 唯一路径与 mate-link 连接
→ 原 k 失败时，局部尝试 k−10 / k / k+10 的唯一图路径
→ junction-spanning reads 验证闭环
```

所有 bait 写入一个线粒体 locus。参考仅用于招募与 seed；最终序列不按参考坐标拼接或补洞。mate link 只决定候选邻接和方向，断点碱基必须由同一 filtered read pool 的唯一路径恢复；不能恢复时保持断裂，不插入 `N`。多 k 仅在原 k 不存在可接受路径时运行，且不能绕过唯一性或 junction-read 条件。

`mito` 默认启用无限延伸与 GFA 图输出；若参考较远或覆盖较低，建议显式使用较敏感的 `-kf 17–25 -s 1`。每个自适应深度仍受 `--mito-max-reads` 限制；需要扫描完整文库时应将其设为不小于输入量。
无达阈值 mate-link 时不会构建 read graph。

## 成功标准

闭环结果必须同时满足：单一 component、无 `N`、每个连接都来自严格 overlap 或唯一 GFA/read-graph bridge、首尾闭合，并有至少 `--mito-min-junction-support` 条真实 reads 跨越最终 junction。条件不足时保留 partial 输出，但命令返回失败。

## 常用与专家参数

- `--mito-genbank`：必需；带注释的线粒体 GenBank 参考。
- `--mito-max-reads 320`：自适应流程最大约 1.05M paired-read blocks；连续两个阶段产生严格一致的闭环（忽略切点和链方向）时提前停止。
- `--no-mito-adaptive-stop`：关闭分阶段早停，改用普通 `--max-reads` 的一次性过滤。

以下为隐藏的专家覆盖参数，仅在明确诊断招募、图连接或闭环失败时调整：`--mito-min-overlap`、`--mito-min-overlap-identity`、`--mito-terminal-window`、`--mito-link-kmer`、`--mito-min-link-hits`、`--mito-min-pair-support`、`--mito-bridge-kmer`、`--mito-bridge-min-depth`、`--mito-max-bridge`、`--mito-min-junction-support`。

## 输出

- `<sample>/mito/mitochondrial_assembly.fasta`：闭环序列或 partial components。
- `<sample>/mito/mitochondrial_assembly_summary.tsv`：兼容的粗粒度状态、`resolution_reason`、连接与 junction 支持。
- `<sample>/mito/mitochondrial_evidence.json`：结构、图、mate-link 与 junction 的机器可读证据。
- `<sample>/mito/mitochondrial_feature_evidence.tsv`：参考 feature 的 canonical 21-mer 相似度；仅说明与参考的精确锚点共享，不能据此判断远缘样本的基因存在或缺失；`translation_status=not_checked`，不是 CDS 注释或翻译判定。
- `<sample>/mito/mitochondrial_mate_links.tsv`：已接受的 read-supported links。
- `.gm2_mito_reference/metadata/mitochondrial_genes.tsv`：bait 元数据；坐标为 0-based 半开区间，`segments_0_half_open` 保留跨 origin 或 `join(...)` feature 的全部片段。
- <sample>/uce_rescue_round_1/assembly_refs/mito_rescue_seeds.tsv：救援 seed 审计；对含歧义碱基 contig 仅使用其中连续 ACGT 区段；再剔除正反向重复或无信息低复杂度 seed，不按参考相似度筛除远缘序列。
