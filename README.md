# GeneMiner2-UCE

**[English README](README_EN.md)**

GeneMiner2-UCE 是 GeneMiner2 的命令行分支，面向 target-enrichment 和 ultraconserved element（UCE）数据。它保留原有的参考引导 reads 捕获和组装框架，并针对短探针、侧翼序列恢复、结果质控及系统发育分析进行了调整。本仓库仅包含 CLI 源码、构建文件和命令行文档，不包含原 GUI、截图或演示数据。

## 使用声明

使用当前版本时必须引用 [GeneMiner2-UCE GitHub 仓库](https://github.com/GUIBA-EX/GeneMiner2-UCE)。介绍 GeneMiner2-UCE 的预印本正在准备中，发布后将在此补充。如需修改代码，请联系 [xf@g.ecc.u-tokyo.ac.jp](mailto:xf@g.ecc.u-tokyo.ac.jp)。

## 功能概览

- 从二代测序 reads 中参考引导恢复目标分子标记。
- 通过 `--assembly-mode uce` 保留有 reads 支持的 UCE 侧翼序列。
- 通过 `--uce-rescue-reads` 执行一轮受控的 raw-read rescue。
- 根据唯一定位 reads、连续覆盖和深度指标筛除弱支持 contig。
- 导出 PHYLUCE 兼容的 UCE contig，并生成样本和 locus 层面的恢复统计。
- 通过 Rust `population` 流程构建公共 UCE 参考、统一 mapping 和联合 VCF，并输出 PCA 与 ADMIXTURE 遗传成分分析所需面板。
- 支持 MAFFT、MUSCLE、Clustal Omega、trimAl、AliFilter 及多种系统发育树程序。

![GeneMiner2-UCE 流程](docs/images/summary_ZH.png)

## 构建

GeneMiner2-UCE 目前不是可直接通过 `pip install` 安装的 Python 包。安装完整构建依赖后，在仓库根目录运行：

```bash
make
```

构建后的入口为：

```bash
cli/geneminer2
```

每次执行 `git pull` 更新源码后都应重新运行 `make`。完整构建方法和运行时依赖见[命令行说明](manual/ZH_CN/command_line.md)。

## 快速开始

准备 tab 分隔的样本列表；每行依次为样本名、R1 和可选的 R2：

```text
Sample_A	/path/to/Sample_A_R1.fq.gz	/path/to/Sample_A_R2.fq.gz
Sample_B	/path/to/Sample_B_R1.fq.gz	/path/to/Sample_B_R2.fq.gz
```

将每个目标 locus 保存为单独的 FASTA 文件：

```text
references/
  uce-0001.fasta
  uce-0002.fasta
```

下面的命令运行 UCE 流程并启用一轮 raw-read rescue：

```bash
cli/geneminer2 \
  -f samples.tsv \
  -r references \
  -o output \
  --assembly-mode uce \
  --uce-rescue-reads
```

## UCE 流程设计

### 组装模式与候选筛选

`--assembly-mode uce` 降低短参考或探针边界对组装的限制，优先选择延伸更长且仍有 reads 支持的候选 contig。未显式指定子命令时，UCE 模式默认跳过参考引导的 `trim` 步骤，避免将新恢复的侧翼序列再次裁回探针区域；如需参考切齐，可显式加入 `trim`。

针对短 UCE bait 或与参考序列存在一定分化的样本，可考虑：

```bash
--assembly-mode uce \
  -sb unlimited \
  -ka 0 \
  --min-ka 17 \
  --max-ka 31 \
  -e 1
```

这些设置会放宽边界控制、在较低 assembly k-mer 范围内自动选 k，并降低 k-mer 计数阈值，因此也可能增加低支持候选。使用后应检查组装摘要、rescue 摘要和下游比对。

UCE 模式默认不限制 contig 长度。长度至少为 1000 bp 的 contig 仍必须满足 `uniquely_placed_read_count / contig_length >= 0.003`。重复定位 reads 会记录在摘要中，但不提供唯一位置支持。若需要更保守地排除异常延长 contig，可显式设置长度上限：

```bash
--uce-max-contig-length 5000 \
  --uce-min-read-density 0.003 \
  --uce-density-check-min-length 1000
```

以下深度均匀性限制默认关闭；只有在确实需要排除深度极不均匀或疑似重复序列的候选时才建议设为正值：

```bash
--uce-max-depth-cv 0 \
  --uce-max-depth-ratio 0
```

### 参考序列缓存

重复使用同一参考目录时，可通过 `--reuse-reference-cache` 复用 reference k-mer index。缓存默认位于 `output/.gm2_reference_cache`；也可用 `--reference-cache-dir` 指定共享项目目录或临时存储位置。

缓存指纹包含参考文件名、大小、修改时间、`-kf` 和 `-s`。它只减少重复建索引的时间，不改变 contig 选择或组装质量。由于 UCE rescue reference 包含样本特异的第一轮 contig，因此 rescue index 会按样本重新构建。

### 成对 reads 保留与 raw-read rescue

在 UCE 模式的 re-filtering 阶段，只要任一 mate 通过 locus 过滤，整对 paired-end reads 都会保留。这样既能保留落在保守 UCE core 上的 read，也能利用另一端提供的侧翼信息。

`--uce-rescue-reads` 在第一轮组装后执行一次额外招募：

1. 合并原始 locus reference 与第一轮 contig，构建临时 rescue reference。
2. 从 raw reads 中重新捕获匹配 reads。
3. 使用 rescue reference 重新执行 re-filtering 和 assembly。
4. 比较第一轮与 rescue 结果，决定保留或回退。

rescue 阶段最多同时处理 4 个样本，每个样本最多使用 4 个线程，并根据 `-p` 和样本数自动缩减并行规模。

若第一轮 locus 已被接受，而 rescue 结果缺失或未通过接受标准，流程会恢复第一轮 contig，并标记为 `reverted_failed_rescue`。若两轮结果都被接受，但 rescue 后唯一 reads density 明显下降，则根据以下比值判断是否回退：

```text
before_density = before_unique_read_count / before_contig_length
rescue_density = rescue_unique_read_count / rescue_contig_length
density_ratio = rescue_density / before_density
```

默认在 `density_ratio < 0.5` 时恢复第一轮 contig，并标记为 `reverted_density_drop`。阈值可通过 `--uce-rescue-min-density-ratio` 调整。`uce_rescue_summary.csv` 中的 `after_*` 字段描述 rescue 尝试结果；发生回退时，最终序列仍是第一轮 contig。

### 输出与统计

UCE 模式额外生成：

- `uce_assembly_summary.csv`：每个样本和 locus 的接受状态、拒绝原因、contig 长度、唯一位置支持范围、支持碱基并集、支持广度、最大无支持间隔、总 reads、唯一定位 reads、多重定位 reads、density、k-mer depth、候选数及低质量标记。
- `uce_rescue_summary.csv`：rescue 前后比较、density ratio、回退状态和错误信息。
- `uce_contigs/`：按样本导出的 PHYLUCE 兼容 contig FASTA。
- `contigs_all_low/`：保留未通过标准的候选供检查；这些序列不会进入主结果、rescue reference、合并矩阵或 PHYLUCE 导出。

运行结束后可汇总恢复情况：

```bash
cli/geneminer2 stats \
  -f samples.tsv \
  -r references \
  -o output \
  --stats-no-heatmap
```

该命令输出 `uce_stats.tsv`、`uce_locus_stats.tsv`、`uce_seq_lengths.tsv`、`uce_read_counts.tsv` 和 `uce_filtered_read_counts.tsv`。如果环境中安装了 `pandas`、`seaborn` 和 `matplotlib`，且未使用 `--stats-no-heatmap`，还会生成恢复率和 read-count heatmap。

## Population 模式

`population` 面向二倍体 UCE 重测序或 target-enrichment 数据中的 PCA、ADMIXTURE 和物种界定。它不要求单倍型定相，而是从每个样本的已接受 UCE contig 和原始 reads 构建一致的二倍体基因型矩阵：

1. 按 locus 汇总已接受 contig，默认采用类似 SqCL `make_PRG.py` 的策略选择最长合格 contig，构建公共参考；`supported` 策略可改为优先 reads 支持。验证或复用既有参考时，可用 `--population-reference-fasta` 指定固定 FASTA。
2. 使用 minibwa 将所有样本的原始 reads 统一映射到同一公共参考，并用 samtools 排序、去重复和统计 mapping 质量。
3. 使用 bcftools 对全部样本联合检测变异，在基因型层面应用 DP/GQ，在位点层面应用 QUAL、call rate 和 MAC 过滤。
4. 输出全部 SNP、每个 UCE 一个 SNP 和 LD-pruned 三种面板；PLINK 对三种面板分别执行 PCA，ADMIXTURE 默认使用每个 UCE 一个 SNP 的主面板。

已有 UCE 组装结果后运行：

```bash
cli/geneminer2 population \
  -f samples.tsv \
  -r references \
  -o output \
  --assembly-mode uce \
  -p 8 \
  --population-admixture-k-min 2 \
  --population-admixture-k-max 6
```

运行时需要 `minibwa`、samtools、bcftools、PLINK 1.9 和 ADMIXTURE 位于 `PATH`。ADMIXTURE 缺失时，公共参考、VCF、PLINK 面板和 PCA 仍会完成，并在 `population/structure/admixture/status.tsv` 中记录 `unavailable`。`--population-start-at mapping`、`calling` 或 `selection` 只会复用经过检查的既有公共参考、BAM 或过滤 VCF，适合在不重复 mapping 的情况下调整下游阈值。正式推断前应检查 `mapping/mapping_qc.tsv`、`variants/variant_qc.tsv` 和不同 SNP 面板的 PCA 是否一致；仅在内部构建公共参考时，才需结合 `reference/reference_contribution.tsv` 判断是否由少数样本贡献过多而产生参考偏倚。

## 实现与下游工具

主 reads 过滤器现在提供 Rust 实现 `rust/main_filter_new/`。其命令行参数、六种输出模式和 paired-end 保留逻辑与原 `scripts/filter/MainFilterNew.hx` 保持兼容，并使用精确的 k-mer 到 locus 映射和带版本标识的缓存格式。`make` 在 Cargo 可用时默认构建 Rust 版本；原 Haxe 源码完整保留，运行 `make haxe-filter` 可另行生成 `cli/bin/MainFilterNew-haxe`。Rust 可以读取旧 Haxe 缓存；切回 Haxe 版本时应重建缓存，因为 Haxe 不能读取新的 Rust 缓存格式。

二次 reads 过滤器提供 Rust 实现 `rust/main_refilter_new/`，并保持与 `scripts/main_refilter_new.py` 相同的命令行参数和输出结构。构建时如检测到 Cargo，会优先编译 Rust 版本；否则回退到 Python/PyInstaller 实现。

assembler 默认使用 `rust/main_assembler/` 的 Rust 实现。`--assembler-implementation auto`（默认）会先运行 `cli/bin/main_assembler-rust`；若二进制不可用、执行失败或没有生成结果，则清理不完整输出并自动重试由 Git 基线保存的未修改原始实现 `cli/bin/main_assembler-original`。`rust` 可用于严格只跑 Rust，`original` 可直接跳过 Rust。当前修改版 `scripts/main_assembler.py` 不参与 fallback。

population 主流程位于 `rust/main_population/`。公共参考构建、流程编排、SNP 选择、面板汇总和 ADMIXTURE 交叉验证解析均由 Rust 实现；minibwa、samtools、bcftools、PLINK 和 ADMIXTURE 作为经过验证的外部程序调用。

combine 阶段可通过 `--msa-threads` 和 `--filter-processes` 控制并行。`--alignment-filter alifilter` 可调用 AliFilter 替代 trimAl；AliFilter 不随本仓库分发，必须另行安装并确保 `AliFilter` 位于 `PATH`。省略 `--alifilter-model` 或设为 `default` 时使用内置模型，只有自定义模型才需要提供真实的 `model.json` 路径。

## 文档

- [中文命令行说明](manual/ZH_CN/command_line.md)
- [中文输出文件说明](manual/ZH_CN/output.md)
- [English command-line usage](manual/EN_US/command_line.md)
- [English output files](manual/EN_US/output.md)
- [版本更新记录](CHANGELOG.md)

## 相关文献

GeneMiner2 主要引用：

Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111. [https://doi.org/10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

相关前期工具：

Zhang Z, Xie PL, Guo YL, Zhou WB, Liu EY, Yu Y. 2022. **Easy353**: A tool to get Angiosperms353 genes for phylogenomic research. *Molecular Biology and Evolution* 39(12): msac261. [https://doi.org/10.1093/molbev/msac261](https://doi.org/10.1093/molbev/msac261)

Xie PL, Guo YL, Teng Y, Zhou WB, Yu Y. 2024. **GeneMiner**: A tool for extracting phylogenetic markers from next-generation sequencing data. *Molecular Ecology Resources* 24(3): e13924. [https://doi.org/10.1111/1755-0998.13924](https://doi.org/10.1111/1755-0998.13924)

如果使用 `--alignment-filter alifilter`，还应引用：

Bianchini G, Zhu R, Cicconardi F, Moody ERR. 2026. **AliFilter: a machine learning approach to alignment filtering.** *Molecular Biology and Evolution* 43(4): msag097. [https://doi.org/10.1093/molbev/msag097](https://doi.org/10.1093/molbev/msag097)
