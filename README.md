# GeneMiner2-UCE

**[English README](README_EN.md)**

GeneMiner2-UCE 是 GeneMiner2 专门给 UCE 扩出来的版本，主要收拾 target-enrichment、ultraconserved element（UCE）和相关二代测序数据。原来参考引导招募 reads 那套框架还留着，另外又加上了 UCE 侧翼恢复、marker profiling 和 UCE 群体遗传分析，基本上一套就能整明白。主打的就是一个带派。

用这个软件出了结果，别忘了引用 [GeneMiner2-UCE GitHub 仓库](https://github.com/GUIBA-EX/GeneMiner2-UCE)。正式论文发出来以后，这旮沓会马上补上。

![GeneMiner2-UCE 流程](docs/images/summary_ZH.png)

## 这玩意儿能干啥

- 从 genome skimming 或 target-capture reads 里把目标分子标记捞出来。
- 跑 UCE 模式时，把有 reads 撑腰的 core 和侧翼序列都留下。
- `profiling` 从 WGS 或其他 shotgun reads 招募 marker 相关 reads，再直接作 k-mer 伪比对并输出参考序列支持。
- 拿多个 UCE 样本整公共伪参考、联合 VCF、PCA 和 ADMIXTURE 输入。
- 导出 PHYLUCE 能直接接上的 contig，再把样本和 locus 的恢复质量归拢成表。

## 先瞅瞅该用哪个模式

| 模式 | 啥数据适合 | 最后能整出啥 |
| --- | --- | --- |
| `--assembly-mode original` | exon、SCO 及核/线粒体 marker | 参考引导 contig；默认流程会按参考裁切 |
| `--assembly-mode uce` | genome skimming 或 target capture 中的 UCE | UCE core 和有 reads 支持的 flanking sequence |
| `profiling` 子命令 | WGS / metagenome 中任意扩增子 marker | 每条参考序列的命中、共享支持和单例支持 |
| `mito` 子命令 | 常规环形动物线粒体与带注释 GenBank 参考 | 仅由样本 reads 支持的闭环或 partial 线粒体组装 |
| `population` 子命令 | 多个已完成 UCE 组装的样本 | 公共伪参考、联合 VCF、PCA 和 ADMIXTURE 输入 |

## 咋整进你那系统里？

GeneMiner2-UCE 现在得从源码构建。先激活包含 Cython、PyInstaller 和 Rust 的构建环境，再进仓库根目录：

```bash
conda activate geneminer2uce
make build
```

螚完以后，入口搁这儿：

```bash
cli/geneminer2
```

以后要是拉了带源码改动的新版本，记着再跑一遍 `make build`，别拿旧程序硬对付。完整构建依赖和外部工具要求都写在[中文命令行指南](manual/ZH_CN/command_line.md)里了。

## 麻溜儿跑一遍

先准备一个用 tab 分隔的样本表。每行先写样本名，再写 R1，后头的 R2 有就写、没有就拉倒：

```text
Sample_A\t/path/to/Sample_A_R1.fq.gz\t/path/to/Sample_A_R2.fq.gz
Sample_B\t/path/to/Sample_B_R1.fq.gz
```

参考目录里，一个 locus 单独搁一个 FASTA 文件（FASTA内可以有复数参考），文件主名就当 locus 名，别整串了：

```text
references/
  uce-0001.fasta
  uce-0002.fasta
```

下面这条命令跑 UCE 恢复，再额外整一轮受控的 raw-read rescue：

```bash
cli/geneminer2 \
  -f samples.tsv \
  -r references \
  -o output \
  -p 8 \
  --assembly-mode uce \
  --uce-rescue-reads
```

跑完先瞅这仨地方，基本就知道整得咋样了：

- `output/uce_assembly_summary.csv`：每个样本、每个 locus 组没组上，reads 支持够不够；
- `output/uce_rescue_summary.csv`：rescue 前后差多少，最后有没有退回第一轮；
- `output/uce_contigs/`：给 PHYLUCE 预备好的样本 contig。

剩下的目录和文件都在[输出文件说明](manual/ZH_CN/output.md)里，找不着就上那儿翻。

## Original 模式咋回事

`original` 是默认组装模式，面向 exon、SCO 以及核或线粒体 marker 的参考引导恢复；不写 `--assembly-mode` 时即使用它。默认后端是 `original-rust`。如需与固定的上游 GeneMiner2 Python 原版逐项对照，显式加上 `--assembler-implementation original`。`uce-rust` 是 UCE 定向后端；常规 UCE 恢复应使用 `--assembly-mode uce`。

```bash
cli/geneminer2 \
  -f samples.tsv -r references -o output -p 8 \
  --assembly-mode original
```

`original` 是工作流模式；`--assembler-implementation original` 是该模式下可选的 Python 后端，两者不要混淆。

## UCE 模式咋回事

UCE 模式会松开短 probe 边界对组装的限制，默认跳过参考引导的 `trim`，优先留下延伸更长、同时还有 reads 支持的候选。refilter 的时候，只要一对儿 mate 里有一个通过 locus 过滤，整对 paired-end reads 都留下，不能把有用的侧翼信息半道扔了。

默认 Rust assembler 走不反复 backtrack 的 backbone 策略，遇上气泡不会来回来去磨叽。`--uce-rescue-reads` 会拿第一轮 contig 加原始参考再招一遍 reads；要是 rescue 以后质量掉了，就麻溜儿退回第一轮结果。组装策略、质量护栏、cache 和 rescue 规则见 [Assembler 章节](docs/assembler_ZH.md)。

## 线粒体模式咋回事

`mito` 专门给**常规单条环形动物线粒体基因组**用：从带注释的 GenBank 参考招募 reads，组装后只在 overlap、局部 read graph 和 junction reads 都支持时报告闭环。它不拿参考坐标补洞，也不适合多分子、严重重排、强异质性或复杂植物和真菌线粒体。

```bash
cli/geneminer2 mito \
  -f samples.tsv \
  -o mito_output \
  -p 8 \
  --mito-genbank mitochondrial_reference.gb
```

完整流程、成功标准和专家参数见 [线粒体章节](docs/5.mito.md)。

## Profiling 模式咋回事

`profiling` 是**读段证据模式，不组装**：先做一次 k-mer 招募，再用 Themisto 把招募 reads 同 marker 参考库比较。默认主结果是每条参考序列的支持证据；不会运行 `refilter`、`assemble`、`combine` 或 `tree`。

最常用的参考序列级模式就是默认：

```bash
cli/geneminer2 profiling \
  -f samples.tsv -r marker_reference.fasta \
  -o output -p 8
```

主结果 `<output>/<sample>/marker_profile/marker_reference_support.tsv` 会给每条命中的参考序列写出总命中、分数化共享支持与单例支持。一条 read 若有 N 个候选，每个候选只获得 `1/N`，不能把同一条 read 重复算给所有候选。它表示与参考相容的证据，不是唯一鉴定或生物量丰度。

`--profile-group-map` 可选：用于在结果的 `group` 列给参考序列附加物种、clade 或整理过的 DIV 标签，不再计算 group 丰度。

完整参数、输出和 QC 见 [Profiling 章节](docs/profiling_ZH.md)。

## Population 模式咋回事

`population` 拿多个已完成 UCE 组装的样本及其原始 reads，整出公共伪参考、联合 VCF、PCA 和 ADMIXTURE 输入。它先生成或者复用公共伪参考，再把所有样本统一 mapping、联合检测变异，最后给出 all-SNP、每个 UCE 一个 SNP 和 LD-pruned 三套面板。

UCE 已经组装完了，就这么跑：

```bash
cli/geneminer2 population \
  -f samples.tsv \
  -r baits_by_locus \
  -o output \
  -p 8 \
  --assembly-mode uce \
  --engine panref \
  --population-admixture-k-min 2 \
  --population-admixture-k-max 6
```

`panref` 以冻结核心图和受限侧翼救援构建公共参考；默认只将通过QC的locus写入参考。运行时需要 minibwa、samtools、bcftools 和 PLINK 1.9；ADMIXTURE 可选。参数、阶段和QC见 [Population 章节](docs/population_ZH.md)。

## 文档往哪儿找

首页只管选模式和跑起第一条命令；参数、质量控制和输出字段都归到下面的专题文档。`original` 默认用 `original-rust`，需要与上游 Python 对照时才显式加 `--assembler-implementation original`；`uce` 只使用 `uce-rust`。

- [中文命令行指南](manual/ZH_CN/command_line.md)
- [中文输出文件说明](manual/ZH_CN/output.md)
- [Filter 章节](docs/filter_ZH.md)
- [Assembler 章节](docs/assembler_ZH.md)
- [5. Mito：线粒体基因组](docs/5.mito.md)
- [Profiling 章节](docs/profiling_ZH.md)
- [Population 章节](docs/population_ZH.md)
- [版本更新记录](CHANGELOG.md)

## 引用和联系

GeneMiner2 主要引用搁这儿，可别落下：

Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111. [https://doi.org/10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

相关工具还有 [Easy353](https://doi.org/10.1093/molbev/msac261) 和 [GeneMiner](https://doi.org/10.1111/1755-0998.13924)。要是用了 `--alignment-filter alifilter`，也别忘了引用 [AliFilter](https://doi.org/10.1093/molbev/msag097)。

要联系项目维护者，就往 [xf@g.ecc.u-tokyo.ac.jp](mailto:xf@g.ecc.u-tokyo.ac.jp) 发邮件，有啥事儿说明白就行。
