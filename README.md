# GeneMiner2-UCE

**[English README](README_EN.md)**

GeneMiner2-UCE 是 GeneMiner2 专门给 UCE 扩出来的版本，主要收拾 target-enrichment、ultraconserved element（UCE）和相关二代测序数据。原来参考引导招募 reads 那套框架还留着，另外又加上了 UCE 侧翼恢复、ITS2 多候选组装和 UCE 群体遗传分析，基本上一套就能整明白。主打的就是一个带派。

用这个软件出了结果，别忘了引用 [GeneMiner2-UCE GitHub 仓库](https://github.com/GUIBA-EX/GeneMiner2-UCE)。正式论文发出来以后，这旮沓会马上补上。

![GeneMiner2-UCE 流程](docs/images/summary_ZH.png)

## 这玩意儿能干啥

- 从 genome skimming 或 target-capture reads 里把目标分子标记捞出来。
- 跑 UCE 模式时，把有 reads 撑腰的 core 和侧翼序列都留下。
- 跑 ITS2 模式时，不硬挑一条，多个候选都留着，还给出 paired、diagnostic 和 EM 丰度证据。
- 拿多个 UCE 样本整公共伪参考、联合 VCF、PCA 和 ADMIXTURE 输入。
- 导出 PHYLUCE 能直接接上的 contig，再把样本和 locus 的恢复质量归拢成表。

## 先瞅瞅该用哪个模式

| 模式 | 啥数据适合 | 最后能整出啥 |
| --- | --- | --- |
| `--assembly-mode reference` | genome skimming、常规基因恢复 | 参考引导 contig；默认流程会按参考裁切 |
| `--assembly-mode uce` | UCE target capture | UCE core 和有 reads 支持的 flanking sequence |
| `--assembly-mode its2` | 同一个 locus 里可能有多个真 ITS2 变体 | 多条候选序列和 fragment-level 支持 |
| `population` 子命令 | 二倍体 UCE 群体重测序或 target capture | 公共伪参考、联合 SNP、PCA 和 ADMIXTURE 面板 |

## 咋整进你那系统里？

GeneMiner2-UCE 现在得从源码构建。进到仓库根目录，bash螚（make）一下就成了：

```bash
make
```

螚完以后，入口搁这儿：

```bash
cli/geneminer2
```

以后要是拉了带源码改动的新版本，记着再跑一遍 `make`，别拿旧程序硬对付。完整构建依赖和外部工具要求都写在[中文命令行指南](manual/ZH_CN/command_line.md)里了。

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

## UCE 模式咋回事

UCE 模式会松开短 probe 边界对组装的限制，默认跳过参考引导的 `trim`，优先留下延伸更长、同时还有 reads 支持的候选。refilter 的时候，只要一对儿 mate 里有一个通过 locus 过滤，整对 paired-end reads 都留下，不能把有用的侧翼信息半道扔了。

默认 Rust assembler 走不反复 backtrack 的 backbone 策略，遇上气泡不会来回来去磨叽。`--uce-rescue-reads` 会拿第一轮 contig 加原始参考再招一遍 reads；要是 rescue 以后质量掉了，就麻溜儿退回第一轮结果。参数、质量护栏、reference cache 和回退规则都归拢在 [UCE 流程说明](docs/uce-workflow_ZH.md)里。

## ITS2 模式咋回事

ITS2 模式固定用 Rust assembler 和 21-mer。一个 locus 里要是真有多个候选，就都留着，再报告 fragment、paired-fragment、diagnostic-fragment 和 EM abundance。reads 分不清的候选会归到 equivalence group 里，不瞎拍脑袋硬选一个。

常规 `combine` 是按每个 locus 一条序列设计的，所以 ITS2 得明确告诉程序只跑这三步：

```bash
cli/geneminer2 filter refilter assemble \
  -f samples.tsv \
  -r references \
  -o output \
  -p 8 \
  --assembly-mode its2
```

ITS2 不认原始 Python assembler，Rust assembler 要是失败了也不往 Python 回退。结果搁在各样本的 `results/`、`its2_assembly_summary.csv` 和 `<locus>.its2_support.tsv`，上这几个地方找准没错。

## Population 模式咋回事

`population` 会拿已经接受的 UCE contig 和原始 reads，整出一套口径一致的未定相二倍体 SNP 矩阵。它先生成或者复用公共伪参考，再把所有样本统一 mapping、联合检测变异，最后给出 all-SNP、每个 UCE 一个 SNP 和 LD-pruned 三套面板。

UCE 已经组装完了，就这么跑：

```bash
cli/geneminer2 population \
  -f samples.tsv \
  -r references \
  -o output \
  -p 8 \
  --assembly-mode uce \
  --population-admixture-k-min 2 \
  --population-admixture-k-max 6
```

运行时得把 minibwa、samtools、bcftools 和 PLINK 1.9 预备好；ADMIXTURE 没有也能接着跑，只是不出那部分结果。伪参考咋选、断点咋接、三套 SNP 面板咋看、哪些 QC 必须瞅，都写在 [Population 流程说明](docs/population_ZH.md)里了。

## 实现和文档都搁哪儿

默认会构建 Rust MainFilter、Refilter、Assembler、Population 和 Rust 辅助工具。`--assembler-implementation auto` 先试 Rust，真整不动了再退回没改过的原始 Python assembler；ITS2 是例外，它只认 Rust。主 CLI 编排器和 consensus 程序还继续用 Python，这块儿没硬改。

- [中文命令行指南](manual/ZH_CN/command_line.md)
- [中文输出文件说明](manual/ZH_CN/output.md)
- [UCE 流程说明](docs/uce-workflow_ZH.md)
- [Population 流程说明](docs/population_ZH.md)
- [MainFilter 性能与兼容性说明](docs/mainfilter-performance.md)
- [版本更新记录](CHANGELOG.md)

## 引用和联系

GeneMiner2 主要引用搁这儿，可别落下：

Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111. [https://doi.org/10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

相关工具还有 [Easy353](https://doi.org/10.1093/molbev/msac261) 和 [GeneMiner](https://doi.org/10.1111/1755-0998.13924)。要是用了 `--alignment-filter alifilter`，也别忘了引用 [AliFilter](https://doi.org/10.1093/molbev/msag097)。

要联系项目维护者，就往 [xf@g.ecc.u-tokyo.ac.jp](mailto:xf@g.ecc.u-tokyo.ac.jp) 发邮件，有啥事儿说明白就行。
