# 命令行版本

本仓库只保留 GeneMiner2 命令行流程。GUI 工程、图形界面文档、截图和内置演示数据已经移除。

## 从源码构建

从源码构建独立 CLI 包需要 Haxe、C++ 编译器、Python 和 PyInstaller。Cargo/Rust 是可选依赖：如果环境中有 Cargo，`make` 会构建 Rust 版 `main_refilter_new`；如果没有 Cargo，则自动回退到 Python 版。真实分析应先运行 `make`，再使用打包生成的 `cli/geneminer2`。GeneMiner2 的使用方法请跳转到文档后半部分的[用法](#用法)章节。

第一步，安装系统的默认C++编译器和[zlib](https://zlib.net/)。在Ubuntu 24.04 或更新版本上，可以运行下面的命令：

```
sudo apt install build-essential zlib1g zlib1g-dev
```

高级用户可以选择安装[zlib-ng](https://github.com/zlib-ng/zlib-ng)或[cloudflare-zlib](https://github.com/cloudflare/zlib)代替zlib以提高性能。

第二步，安装[Haxe 4](https://haxe.org/)并[配置C++目标](https://haxe.org/manual/target-cpp-getting-started.html)。在Ubuntu 24.04 或更新版本上，运行下面的命令：

```
sudo apt install haxe
haxelib setup ~/haxelib
haxelib install hxcpp
```

第三步，安装Python和依赖项。如果可以使用conda，执行下面的命令即可；否则，需要用系统的软件包管理器手动安装这些包。

```
conda create -c conda-forge -n geneminer python=3.11 numpy=2.1.3 biopython cython matplotlib pandas seaborn pyinstaller scipy setuptools wheel
conda activate geneminer
```

最后，将本分支的源代码下载到GeneMiner2文件夹中，构建可执行文件。

```
git clone --depth 1 https://github.com/GUIBA-EX/GeneMiner2.git
cd GeneMiner2
make
```

编译完成的可执行文件放在`cli`目录下，运行`cli/geneminer2 -h`查看帮助。软件在运行时仍然需要几个独立的生物信息学工具。用conda安装的方法如下：

```
conda install -c bioconda aster blast clustalo fasttree iqtree mafft magicblast minimap2 muscle raxml-ng trimal veryfasttree
```

AliFilter可以作为可选的多序列比对过滤程序使用。上面的conda命令不会安装AliFilter；如需使用`--alignment-filter alifilter`，请单独下载AliFilter可执行文件，并确保命令行中可以直接调用`AliFilter`。

## 用法

GeneMiner2需要tsv格式的样本列表和FASTA格式的参考序列。样本列表的具体格式是`<物种名><Tab><数据文件1>`（单端）或者`<物种名><Tab><数据文件1><Tab><数据文件2>`（双端），每一行代表一个样本。其中的数据文件建议采用绝对路径。

对于双端 reads，样本列表可以写成：

```
Sample_A	/data/reads/Sample_A_R1.fq.gz	/data/reads/Sample_A_R2.fq.gz
Sample_B	/data/reads/Sample_B_R1.fq.gz	/data/reads/Sample_B_R2.fq.gz
Sample_C	/data/reads/Sample_C_R1.fq.gz	/data/reads/Sample_C_R2.fq.gz
```

参考序列需要放在一个独立的文件夹下，每个基因可以有一条或多条参考序列。对于每个基因，将对应的所有参考序列存在`<基因名>.fasta`下。假设需要提取matK和psbA两个基因，则需要在一个空白文件夹下创建`matK.fasta`和`psbA.fasta`两个文件，分别保存对应的参考序列。

假设样本列表保存在`/home/user/project/samples.tsv`，参考序列保存在`/home/user/project/references`，期望的输出文件夹为`/home/user/project/output`，用默认设置运行 GeneMiner2 的命令如下：

```
cli/geneminer2 -f /home/user/project/samples.tsv -r /home/user/project/references -o /home/user/project/output
```

执行这行命令后，GeneMiner2 会构建一棵溯祖树，路径为`/home/user/project/output/Coalescent.tree`。

如果使用`--assembly-mode uce`且不指定子命令，默认流程会跳过基于参考序列的`trim`步骤，运行`filter refilter assemble combine tree`。这样可以避免刚组装出的UCE侧翼序列在后续参考切齐步骤中被再次裁剪。如需仍然执行参考切齐，可以显式加入`trim`子命令。

命令行参数的说明如下：

- `-f`: 样本列表tsv文件
- `-r`: 参考序列文件夹
- `-o`: 输出文件夹
- `-p`: 同时运行的进程数
- `--max-reads`: 每文件限制的最大读长数量（单位为百万），默认禁用
- `-kf`: 过滤kmer大小
- `-s`: 过滤步长
- `--reuse-reference-cache`: 跨运行复用带指纹的 reference k-mer index，避免每次重新构建。
- `--reference-cache-dir`: `--reuse-reference-cache` 使用的缓存目录。默认位于 `output/.gm2_reference_cache`。
- `-ka`: 组装kmer大小，默认为自动
- `--min-ka`: 自动估算k值的下限
- `--max-ka`: 自动估算k值的上限
- `-e`: 错误阈值
- `-sb`: 软边界，支持`0`、`auto`和`unlimited`三个选项
- `-i`: 搜索深度
- `--assembly-mode`: 组装模式，可选`reference`或`uce`。`reference`为默认模式，保持原有基于参考序列的组装和边界控制；`uce`模式会放宽参考边界裁剪，并优先保留更长且有reads支持的UCE侧翼序列。
- UCE组装模式下，refilter阶段会在任一端通过locus过滤时保留整对paired-end reads。这有助于保留与短UCE探针配对的侧翼reads。
- `--uce-side-candidates`: UCE组装中每侧参与组合的分支候选数。数值越大，越可能保留较长的低支持侧翼候选，但会增加运行时间和候选路径数量。
- `--uce-max-contig-length`: UCE候选contig进入评分前允许保留的最大长度。默认值为`5000`；设为`0`可关闭该限制。
- `--uce-min-read-density`: 长UCE候选contig进入评分前所需的最低`read_count / contig_length`。默认值为`0.003`。
- `--uce-density-check-min-length`: 启用UCE read-density保护阈值的最短contig长度。默认值为`1000`。
- `--uce-max-depth-cv`: 可选的UCE候选contig k-mer深度变异系数上限。默认值为`0`，表示关闭该保护阈值。
- `--uce-max-depth-ratio`: 可选的UCE候选contig最大/中位k-mer深度比值上限。默认值为`0`，表示关闭该保护阈值。
- `--uce-rescue-reads`: 仅用于UCE模式。初次组装后，用初步contig和原始参考序列再招募一次raw reads，然后重新进一步过滤并重新组装。
- UCE raw-read rescue 使用受控的样本级并行：最多同时 rescue 4 个样本，每个样本最多 4 个线程；当 `-p` 较小时会自动降低并行度。
- `--uce-rescue-min-contig-length`: 参与UCE raw-read rescue的初步contig最短长度。
- `--uce-rescue-min-density-ratio`: UCE raw-read rescue 后保留结果所需的最低 rescue/第一轮 read-density 比值。默认值为`0.5`；低于该比值的locus会恢复为第一轮contig。
- `-c`: 一致性阈值（介于0-1的小数形式）
- `-ts`: 基于参考切齐的来源序列，可以是`assembly`或`consensus`
- `-tm`: 基于参考切齐的模式，可以是`all`、`longest`、`terminal`或`isoform`
- `-tr`: 基于参考切齐的保留长度阈值（介于0-1的小数形式）
- `-cs`: 合并结果的来源序列，可以是`assembly`、`consensus`或`trimmed`
- `--msa-program`: 多序列比对的软件，支持`mafft`、`muscle`和`clustalo`
- `--msa-threads`: 每个多序列比对任务使用的线程数，默认为1。GeneMiner2会限制同时运行的比对任务数量，使比对阶段请求的总线程数不超过`-p`
- `--alignment-filter`: 建树前的比对列过滤程序，支持`trimal`、`alifilter`和`none`，默认为`trimal`
- `--filter-processes`: 同时运行的trimAl或AliFilter过滤任务数上限，默认等于`-p`
- `--alifilter-model`: 使用`--alignment-filter alifilter`时可选指定AliFilter的`model.json`路径；省略该参数或设为`default`时使用AliFilter内置默认模型
- `--strict-combine-errors`: 如果任一locus在多序列比对、alignment清理或比对列过滤阶段失败，则立即停止`combine`。默认行为是输出警告并跳过失败的locus
- `--no-trimal`: 已弃用，等同于`--alignment-filter none`
- `stats`: 可选子命令，用于在UCE流程结束后汇总恢复结果。该命令会输出样本级统计、locus级统计、长度矩阵、read-count矩阵和过滤read-count矩阵。
- `--stats-no-heatmap`: 不生成统计热图。
- `--stats-count-input-reads`: 统计输入FASTQ reads数量，用于填写`InputReads`和`PctFiltered`；大数据集上可能较慢。
- `-cd`: 清理序列的最大差异
- `-cn`: 清理序列的最小序列数量
- `-m`: 建树的方法，支持`coalescent`和`concatenation`
- `-b`: 支持度计算的重复次数
- `--phylo-program`: 建树的软件，支持`fasttree`、`veryfasttree`、`iqtree`和`raxmlng`

例如，在上述命令运行结束后，使用`tree`子命令和`-m concatenation`参数，可以要求GeneMiner2基于之前的结果构建串联树：

```
cli/geneminer2 tree -f /home/user/project/samples.tsv -r /home/user/project/references -o /home/user/project/output -m concatenation
```

类似，可以通过传递子命令和参数，只运行特定的分析步骤。假设给出这些参数：

| 参数                        | 值                |
| -------------------------- | ----------------- |
| Source Sequence            | Consensus Results |
| 保留长度阈值                 | 50%               |
| Trim Method                | All Fragments     |
| 比对程序                    | muscle            |
| 最大差异大于                 | 0.2               |
| 序列数量小于                 | 5                 |

可以用下面的命令只执行**基于参考切齐**和**合并结果**：

```
cli/geneminer2 trim combine -f /home/user/project/samples.tsv -r /home/user/project/references -o /home/user/project/output -ts consensus -tm all -tr 0.5 -cd 0.2 -cn 5 --msa-program muscle
```

命令行版本在需要填写百分数的位置使用 0-1 的小数。一部分内部参数（例如`--min-coverage`）也可以修改，允许高级用户更灵活地控制软件行为。
