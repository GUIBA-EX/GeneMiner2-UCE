# 线粒体基因组工作流

`mito` 完全复用 GeneMiner2-UCE 现有技术栈：MainFilter、refilter 和 Rust UCE assembler。没有迭代 baiting、外部组装器、参考坐标拼接或参考补洞。

    cli/geneminer2 mito \
      -f samples.tsv \
      -o mito_output \
      -p 8 \
      --mito-genbank mitochondrial_reference.gb

## 主流程

1. Rust 从 GenBank 提取 gene、rRNA、tRNA、flank、环状 tiles 和全基因组 bait。
2. 所有 bait 写在同一个 `.gm2_mito_reference/mitochondrion.fasta` 中，因此 MainFilter 只产生一个线粒体 locus。
3. MainFilter 招募 reads；任一 mate 命中时，refilter 的 `--keep-linked-mates` 保留整对 reads。
4. 所有线粒体 paired reads 位于 `<sample>/filtered/mitochondrion.fq`。
5. 现有 `main_assembler-rust` 以 UCE backbone 模式组装该单一 read pool，并把全部候选写入 `<sample>/contigs_all/mitochondrion.fasta`。
6. Rust finalizer 先按双向序列 overlap 合并 contigs。
7. 对剩余断点，仅使用 contig 两端的唯一 k-mer 建立 mate-link。至少多个独立 read pairs 支持相同方向后才接受连接关系。
8. mate-link 只决定邻接和方向。断点碱基必须由同一 filtered read pool 的唯一局部 de Bruijn path 恢复；找不到唯一路径时保持断裂，不插入 `N`。
9. 首尾也按同样规则解析，并用真实 junction-spanning reads 验证闭环。

GB 参考只参与 MainFilter 和 UCE seed，不复制到最终序列。没有近缘参考时，只要 gene/tile bait 仍能招募到一部分同源 reads，后续组装和 mate-link 均只使用样本数据。

## 闭环成功条件

`mito` 默认要求：

- 最终只有一个 component；
- 所有连接均为直接 overlap 或唯一的 read-graph bridge；
- 序列不含 `N`；
- 首尾闭合；
- 至少 `--mito-min-junction-support` 条 reads 跨越最终 junction。

条件不足时保留 partial FASTA 和统计文件，但命令返回失败，不能把仅由 mate-link 排序的 scaffold 标记为完整线粒体。

## 常用参数

- `--mito-genbank`：必需；带注释的线粒体 GenBank 参考。
- `--mito-max-reads 320`：自适应流程允许处理的最大约 1.05M paired-read block 数；连续两个阶段得到一致闭环时会提前停止。
- `--no-mito-adaptive-stop`：关闭分阶段早停；此时沿用普通 `--max-reads` 的一次性过滤行为。

## 专家覆盖参数

以下参数为兼容已有流程而保留，但不出现在默认命令帮助中；仅在诊断明确的招募、图连接或闭环失败时使用。

- `--mito-min-overlap 41`：直接合并所需的最短 contig overlap。
- `--mito-min-overlap-identity 0.98`：overlap 最低一致性。
- `--mito-terminal-window 500`：用于 mate 映射的 contig 末端窗口。
- `--mito-link-kmer 31`：mate 末端定位 k-mer。
- `--mito-min-link-hits 2`：read 唯一定位所需的不同 k-mer 数。
- `--mito-min-pair-support 3`：接受邻接关系所需的独立 read pairs。
- `--mito-bridge-kmer 31`：解析断点局部 read graph 的 k-mer。
- `--mito-bridge-min-depth 2`：局部路径节点的最小物理 read 深度。
- `--mito-max-bridge 1000`：允许解析的最大断点长度。
- `--mito-min-junction-support 3`：闭环 junction 的最小 read 支持。

## 输出

`metadata/mitochondrial_genes.tsv` 坐标采用 0-based 半开区间：`start_0_inclusive` 包含，`end_0_exclusive` 不包含；`segments_0_half_open` 保留跨 origin 或 `join(...)` feature 的全部坐标段。

- `<sample>/filtered/mitochondrion.fq`：MainFilter/refilter 得到的 paired mitochondrial read pool。
- `<sample>/contigs_all/mitochondrion.fasta`：UCE assembler 的全部候选 contigs。
- `<sample>/mito/mitochondrial_assembly.fasta`：闭环序列或 partial components。
- `<sample>/mito/mitochondrial_assembly_summary.tsv`：组件、mate bridge、闭环方式和 junction 支持。
- `<sample>/mito/mitochondrial_mate_links.tsv`：被真实 read-graph path 解析的 mate links。
