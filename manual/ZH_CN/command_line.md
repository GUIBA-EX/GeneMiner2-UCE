# GeneMiner2-UCE 命令行指南

GeneMiner2-UCE 是面向 UCE 和其他系统发育标记的命令行流程。本仓库不包含 GUI、内置演示数据或旧版图形界面文档。

## 1. 从源码构建

### 1.1 构建依赖

完整构建需要：

- C/C++ 编译器和 [zlib](https://zlib.net/)；
- [Rust 和 Cargo](https://www.rust-lang.org/tools/install/)；
- Python 3.11、Cython、PyInstaller、Biopython、NumPy、SciPy、pandas、matplotlib 和 seaborn。

Rust/Cargo 是当前完整构建的必需依赖，用于编译主 reads 过滤器、二次过滤器和 population 流程。Haxe 源码仅作为主过滤器的兼容实现保留，不能替代完整 Rust 构建。

Ubuntu 可先安装系统依赖：

```bash
sudo apt install build-essential zlib1g zlib1g-dev
```

然后创建 Python 环境：

```bash
conda create -c conda-forge -n geneminer \
  python=3.11 numpy=2.1.3 biopython cython matplotlib \
  pandas seaborn pyinstaller scipy setuptools wheel
conda activate geneminer
```

### 1.2 运行时依赖

系统发育流程按所选步骤调用部分外部程序：

```bash
conda install -c bioconda \
  aster blast clustalo fasttree iqtree mafft magicblast minimap2 \
  muscle raxml-ng trimal veryfasttree
```

可选工具：

- `--alignment-filter alifilter` 需要单独安装 AliFilter，并确保 `AliFilter` 位于 `PATH`。
- `population` 需要 `minibwa`、samtools、bcftools 和 PLINK 1.9。
- ADMIXTURE 用于遗传成分分析。若未安装，population 仍会完成公共参考、VCF、PLINK 面板和 PCA，并将状态记录为 `unavailable`。

population 外部程序默认从 `PATH` 查找，也可通过 `--population-minibwa`、`--population-samtools`、`--population-bcftools`、`--population-plink` 和 `--population-admixture` 分别指定。

### 1.3 下载与构建

```bash
git clone --depth 1 https://github.com/GUIBA-EX/GeneMiner2-UCE.git
cd GeneMiner2-UCE
make
```

构建后的统一入口为：

```bash
cli/geneminer2 -h
```

源码更新后应重新运行 `make`。

## 2. 输入文件

所有命令均需要以下三个参数：

- `-f FILE`：tab 分隔的样本表；
- `-r DIR`：参考序列目录；
- `-o DIR`：输出目录。

### 2.1 样本表

每行表示一个样本。单端数据使用两列，双端数据使用三列：

```text
Sample_A	/data/reads/Sample_A_R1.fq.gz	/data/reads/Sample_A_R2.fq.gz
Sample_B	/data/reads/Sample_B_R1.fq.gz	/data/reads/Sample_B_R2.fq.gz
Sample_C	/data/reads/Sample_C_R1.fq.gz	/data/reads/Sample_C_R2.fq.gz
```

建议使用绝对路径。样本名应保持唯一；population 会另行生成内部样本名和 VCF 样本名的对应表。

### 2.2 参考序列目录

每个 locus 使用一个 FASTA 文件，文件名即 locus 名；同一文件可包含一条或多条参考序列：

```text
references/
  uce-0001.fasta
  uce-0002.fasta
```

## 3. 子命令与默认流程

可在同一条命令中按执行顺序列出一个或多个子命令：

| 子命令 | 功能 |
| --- | --- |
| `filter` | 根据参考 k-mer 从原始数据招募 reads |
| `refilter` | 进一步分配和过滤每个 locus 的 reads |
| `assemble` | 使用 wDBG 组装目标序列 |
| `population` | 构建公共 UCE 参考并生成群体 SNP、PCA 和 ADMIXTURE 结果 |
| `consensus` | 在杂合位点生成一致性序列 |
| `trim` | 按参考序列裁切侧翼 |
| `combine` | 合并样本、比对、清理和过滤 alignment |
| `tree` | 构建溯祖树或串联树 |
| `stats` | 汇总 UCE 恢复统计并可选生成热图 |

不显式指定子命令时：

- `--assembly-mode reference`（默认）运行 `filter refilter assemble trim combine tree`；
- `--assembly-mode uce` 运行 `filter refilter assemble combine tree`，跳过 `trim`，避免新恢复的 UCE 侧翼再次被裁回参考范围。

默认参考模式示例：

```bash
cli/geneminer2 \
  -f /home/user/project/samples.tsv \
  -r /home/user/project/references \
  -o /home/user/project/output \
  -p 8
```

## 4. UCE 组装

UCE 模式放宽短 probe 的参考边界限制，优先保留更长且仍有 reads 支持的侧翼序列。在 refilter 阶段，只要任一 mate 通过 locus 过滤，整对 paired-end reads 都会保留。

基础运行示例：

```bash
cli/geneminer2 \
  -f samples.tsv \
  -r references \
  -o output \
  -p 8 \
  --assembly-mode uce \
  --uce-rescue-reads
```

`--uce-rescue-reads` 在第一轮组装后，以初步 contig 和原始参考再次招募 raw reads，并只执行一轮 re-filtering 和 assembly。rescue 最多并行处理 4 个样本、每个样本最多使用 4 个线程，并受 `-p` 总量约束。

建议在放宽 `-sb`、`-e` 或 assembly k-mer 范围后检查 `uce_assembly_summary.csv`、`uce_rescue_summary.csv` 和下游 alignment。详细输出见[输出文件说明](output.md)。

## 5. Population 群体遗传分析

### 5.1 适用范围

`population` 从多个样本的 UCE 组装结果和原始 reads 构建未定相的二倍体 SNP 矩阵，主要用于 PCA、ADMIXTURE、遗传成分比较和物种界定。

运行前，每个样本必须已经完成 UCE 组装，并保留：

- `uce_assembly_summary.csv`；
- `results/` 中已接受的 UCE contig；
- 样本表所列的原始 reads。

该模式输出未定相基因型，而不是两条完整单倍型。它不替代需要单倍型序列、重组信息或单 locus 基因树的 phasing 流程。

### 5.2 分析流程

1. 按 locus 汇总已接受 contig，构建公共 UCE 参考；
2. 使用 minibwa 将全部样本统一 mapping 到同一参考；
3. 使用 bcftools 联合检测变异并执行基因型和位点过滤；
4. 输出 all-SNP、每个 UCE 一个 SNP 和 LD-pruned 三种面板。

PLINK 对三种面板分别执行 PCA。ADMIXTURE 默认使用每个 UCE 一个 SNP 的主面板，并在指定 K 范围内计算交叉验证误差。

### 5.3 运行示例

已有 UCE 组装结果后，可单独运行：

```bash
cli/geneminer2 population \
  -f /home/user/project/samples.tsv \
  -r /home/user/project/references \
  -o /home/user/project/output \
  -p 8 \
  --assembly-mode uce \
  --population-admixture-k-min 2 \
  --population-admixture-k-max 6
```

正式解释遗传成分前，应检查：

- `population/mapping/mapping_qc.tsv` 中的 mapping rate、coverage breadth 和 depth；
- `population/reference/reference_contribution.tsv` 中公共参考来源是否过度集中；
- 样本和位点缺失率；
- 三种 SNP 面板的 PCA 是否给出一致的主要结构。

## 6. 分步运行示例

只基于已有结果重建串联树：

```bash
cli/geneminer2 tree \
  -f samples.tsv -r references -o output \
  -m concatenation --phylo-program iqtree
```

生成 consensus，按参考裁切并合并：

```bash
cli/geneminer2 consensus trim combine \
  -f samples.tsv -r references -o output \
  -c 0.75 -ts consensus -tm all -tr 0.5 \
  -cs trimmed --msa-program muscle
```

从已有 UCE 结果生成统计表但不绘图：

```bash
cli/geneminer2 stats \
  -f samples.tsv -r references -o output \
  --stats-no-heatmap
```

## 7. 参数参考

以下列出主要公开参数及当前默认值。运行 `cli/geneminer2 -h` 可查看同版本的完整帮助。

### 7.1 通用输入与并行

| 参数 | 说明 |
| --- | --- |
| `-f FILE` | 样本表，必需 |
| `-r DIR` | 参考序列目录，必需 |
| `-o DIR` | 输出目录，必需 |
| `-p INT` | 并行进程总数，默认 `1` |

### 7.2 Reads 过滤与 re-filtering

| 参数 | 说明 |
| --- | --- |
| `-kf INT` | 过滤 k-mer 大小，默认 `31` |
| `-s, --step-size INT` | reads 扫描步长，默认 `4` |
| `--max-reads INT` | 每个文件最多处理的 reads 数，单位为百万；`0` 表示不限 |
| `--reuse-reference-cache` | 复用带输入指纹的 reference k-mer index |
| `--reference-cache-dir DIR` | reference cache 目录；默认 `output/.gm2_reference_cache`，且必须与上一参数同时使用 |
| `--depth-low-water-mark INT` | 低于该深度时尝试放宽条件招募更多 reads，默认 `50` |
| `--depth-limit INT` | re-filtering 处理的最高深度，默认 `768` |
| `--file-size-limit INT` | re-filtering 文件大小限制，默认 `6` |

### 7.3 组装与 UCE

| 参数 | 说明 |
| --- | --- |
| `-ka INT` | assembly k-mer 大小；默认 `0`，表示自动估算 |
| `--min-ka INT` / `--max-ka INT` | 自动估算范围，默认 `21` / `51` |
| `-e, --error-threshold INT` | k-mer 错误阈值，默认 `2` |
| `-sb, --soft-boundary VALUE` | 软边界：整数、`auto` 或 `unlimited`；默认 `auto` |
| `-i, --search-depth INT` | 搜索深度，默认 `4096` |
| `--min-coverage INT` | contig 最低 read depth，默认 `0` |
| `--assembly-mode MODE` | `reference` 或 `uce`；默认 `reference` |
| `--uce-side-candidates INT` | 每侧参与组合的分支候选数，默认 `8`，最小为 `3` |
| `--uce-max-contig-length INT` | 进入评分的 UCE contig 最大长度，默认 `5000`；`0` 关闭 |
| `--uce-min-read-density FLOAT` | 长 contig 的最低唯一 reads/长度，默认 `0.003` |
| `--uce-density-check-min-length INT` | 启用 density 阈值的最短 contig，默认 `1000` |
| `--uce-max-depth-cv FLOAT` | k-mer depth CV 上限，默认 `0`，表示关闭 |
| `--uce-max-depth-ratio FLOAT` | 最大/中位 k-mer depth 上限，默认 `0`，表示关闭 |
| `--uce-rescue-reads` | 执行一轮 UCE raw-read rescue |
| `--uce-rescue-min-contig-length INT` | rescue reference 的最短 contig，默认 `60`，且不会小于 `-kf` |
| `--uce-rescue-min-density-ratio FLOAT` | 保留 rescue 结果的最低 rescue/首轮 density，默认 `0.5` |

### 7.4 Population

| 参数 | 说明 |
| --- | --- |
| `--population-reference-strategy MODE` | `sqcl-longest`（默认）或 `supported` |
| `--population-min-mapq INT` / `--population-min-baseq INT` | 联合检测最低 MAPQ / base quality，默认 `20` / `20` |
| `--population-min-dp INT` / `--population-min-gq INT` | 低于阈值的基因型设为缺失，默认 `5` / `20` |
| `--population-min-qual FLOAT` | 最低位点 QUAL，默认 `20` |
| `--population-min-call-rate FLOAT` | 最低非缺失基因型比例，默认 `0.8` |
| `--population-min-mac INT` | 最低 minor allele count，默认 `2` |
| `--population-ld-window INT` / `--population-ld-step INT` | LD pruning 窗口和步长，默认 `50` / `5` SNP |
| `--population-ld-r2 FLOAT` | LD pruning 的 r² 阈值，默认 `0.2` |
| `--population-admixture-k-min INT` / `--population-admixture-k-max INT` | ADMIXTURE K 范围，默认 `2` / `6`；最大 K 不超过样本数 |
| `--population-admixture-cv INT` | ADMIXTURE CV 折数，默认 `10`，不超过样本数 |
| `--population-stop-after STAGE` | 在 `reference`、`mapping`、`calling` 或 `selection` 后停止；默认 `selection` |
| `--population-skip-mark-duplicates` | 跳过 samtools duplicate marking |
| `--population-skip-plink` | 不生成 PLINK、PCA、LD-pruned 或 ADMIXTURE 结果 |
| `--population-skip-admixture` | 生成 PLINK 和 PCA，但不运行 ADMIXTURE |
| `--population-minibwa PATH` | minibwa 可执行文件，默认 `minibwa` |
| `--population-samtools PATH` | samtools 可执行文件，默认 `samtools` |
| `--population-bcftools PATH` | bcftools 可执行文件，默认 `bcftools` |
| `--population-plink PATH` | PLINK 1.9 可执行文件，默认 `plink` |
| `--population-admixture PATH` | ADMIXTURE 可执行文件，默认 `admixture` |

### 7.5 Consensus、裁切与合并

| 参数 | 说明 |
| --- | --- |
| `-c, --consensus-threshold FLOAT` | consensus 阈值，默认 `0.75` |
| `-ts, --trim-source SOURCE` | `assembly` 或 `consensus`；默认使用上一步结果 |
| `-tm, --trim-mode MODE` | `all`、`longest`、`terminal` 或 `isoform`；默认 `terminal` |
| `-tr, --trim-retention FLOAT` | 参考裁切的保留长度比例，默认 `0` |
| `-cs, --combine-source SOURCE` | `assembly`、`consensus` 或 `trimmed`；默认使用上一步结果 |
| `-cd, --clean-difference FLOAT` | alignment 最大可接受成对差异，默认 `1.0` |
| `-cn, --clean-sequences INT` | alignment 最少序列数，默认 `0` |
| `--msa-program PROGRAM` | `clustalo`、`mafft` 或 `muscle`；默认 `mafft` |
| `--msa-threads INT` | 每个 MSA 任务的线程数，默认 `1`，且不能大于 `-p` |
| `--alignment-filter PROGRAM` | `trimal`、`alifilter` 或 `none`；默认 `trimal` |
| `--filter-processes INT` | 并行 alignment-filter 任务上限，默认等于 `-p` |
| `--alifilter-model MODEL` | AliFilter 模型名或 `model.json` 路径 |
| `--strict-combine-errors` | 任一 locus 的 alignment、清理或过滤失败时立即停止 |
| `--no-alignment` | 跳过多序列比对 |
| `--no-trimal` | 已弃用；等同于 `--alignment-filter none` |

### 7.6 建树与统计

| 参数 | 说明 |
| --- | --- |
| `-m, --tree-method METHOD` | `coalescent` 或 `concatenation`；默认 `coalescent` |
| `-b, --bootstrap INT` | bootstrap 次数，默认 `1000` |
| `--phylo-program PROGRAM` | `raxmlng`、`iqtree`、`fasttree` 或 `veryfasttree`；默认 `fasttree` |
| `--stats-no-heatmap` | 不生成 UCE 统计热图 |
| `--stats-count-input-reads` | 统计输入 FASTQ reads 以填写 `InputReads` 和 `PctFiltered`；大数据集上较慢 |

比例类参数使用 `0–1` 小数。对正式数据修改过滤阈值前，建议先保留默认值完成小规模测试，并同时检查中间 VCF、mapping QC 和 locus 占有率。
