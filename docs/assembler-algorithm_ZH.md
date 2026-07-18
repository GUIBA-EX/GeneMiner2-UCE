# GeneMiner2-UCE 新旧 Assembler 算法说明

本文解释 GeneMiner2 上游原版与 GeneMiner2-UCE Rust 新版之间的算法差异，并说明新版从 MaSuRCA、SPAdes 与 Sparrowhawk 分别借鉴了什么。

最重要的结论是：Rust 新版不是三个外部 assembler 的拼装版，也没有调用它们的代码。它保留了 GeneMiner2 的参考引导种子、k-mer 深度权重和 read-slice 验证，只吸收适合 UCE 小型逐 locus 图的算法思想与工程实现方式。

## 1. 当前实际存在的三个 GeneMiner2 assembler

| 实现 | 可执行文件 | 主要用途 | 分支策略 | 当前默认规则 |
|---|---|---|---|---|
| 上游原版 Python | `main_assembler-original` | genome skimming、普通基因恢复 | 深度优先枚举分支，遇到死路后 backtrack | `reference + auto` 直接使用 |
| Rust 新版 | `main_assembler-rust` | UCE 主组装、ITS2 多候选组装；也可显式测试普通模式 | UCE 默认 bounded-lookahead backbone；另保留 `search` 兼容策略 | `uce + auto` 优先使用；ITS2 强制使用 |
| Rust 原版兼容版 | `main_assembler-original-rust` | reference 模式的性能 A/B 与兼容性验证 | 保留原版的 seed、权重、候选组合和 backtrack；单线程 | 仅 `reference + original-rust` 显式使用，默认不变 |

这里的“原版”是固定到上游 GeneMiner2 提交 `36e06feeb99654bdb87f45d4cde225d8c3e311d0` 的逐字节副本，不含 GeneMiner2-UCE 后来加入的参数或判断。

## 2. 三个实现共有的 GeneMiner2 基础
这三个实现并不是完全不同的 assembler。它们共享同一条核心思路：
1. 将每个 locus 的参考序列切成 k-mer，记录参考位置和方向。
2. 从已经过滤到该 locus 的 reads 建立 k-mer 图，并去掉低频错误 k-mer；参考中出现的 k-mer可以保留并获得参考权重。
3. 优先选择参考中位置明确、方向一致的 k-mer 作为 seed。
4. 从 seed 向左右延伸，再把两侧与 seed 合并为 contig。
5. 用 reads 中间的高质量 slice 检查候选 contig，而不是只相信图上的 k-mer 深度。

所以真正没有改变的是“参考负责指出从哪儿开始，reads 负责证明组装结果靠不靠谱”。新版改变的主要是图如何存、分叉如何走、候选如何评分，以及怎样限制内存和运行时间。

## 3. 上游原版如何组装

原版把每个 k-mer 编码成整数，节点权重由 read depth 与参考权重共同构成。从 seed 向一个方向延伸时，它总是先走当前权重最高的邻接节点；如果存在第二条或更多分支，就把其余分支压入栈中。当前路径走到死路后，算法退回最近的分叉点，再尝试栈中的其他分支，直到搜索深度耗尽。

可以把它想成在迷宫里走路：先走看起来最宽的一条路，但把其他路口记在本子上；走不通时翻本子回来重走。这个办法对普通目标基因很稳妥，因为短区域中的真实变异可以作为多个候选保留下来，但在 UCE target-capture 数据中，旁边的低复杂度 flank、杂合位点和深度噪声可能制造许多小气泡，导致同一段图被反复搜索。

原版每侧按 read-slice 支持度保留前三个候选，最多组合成 `3 × 3 = 9` 条完整 contig。最终主要根据参考种子数、图权重和 read-slice 数排序，并要求 read 支持量足以覆盖 contig。

## 4. Rust 新版改变了什么

### 4.1 数据表示与计数

Rust 版用 2-bit 表示 A、C、G、T，并将 k-mer 精确存入 `u128`，当前支持 `k ≤ 63`。这不是概率哈希：相同整数就是相同 k-mer，不存在仅由哈希碰撞造成的序列合并。

reads 按 `--assembler-read-chunk-size` 分块读取。每条 read 内同一个正向或反向 k-mer 只计数一次，避免重复片段虚抬深度；每个工作块先收集、排序并合并相同 k-mer，最后才汇总到 locus 图中。这减少了大量随机 HashMap 更新，并限制一次驻留内存中的 reads 数量。

### 4.2 UCE backbone：有限前瞻后一次决定

UCE 默认的 `backbone` 模式从 seed 两侧各建立一条主干。没有分叉时连续向前消费整段线性链；遇到气泡时，对每条出边最多向前查看 `--uce-backbone-lookahead` 步，然后依次比较：

1. 能继续延伸的长度；
2. 前瞻路径的累计 k-mer 与参考权重；
3. 第一条边本身的权重。

胜出的分支被写入 backbone，兄弟分支立即标记为 discarded。算法不保存待回溯分支栈，也不会重新进入已经访问的 k-mer，因此气泡只做一次有限决策，环路也会停止。

这相当于在路口先站高一点看有限距离，再选一条路一直走；选完就把其他路从本次主干任务中划掉，而不是走到很远以后再回来重走。它牺牲了穷举所有组合的能力，换取可控运行时间，并降低把多个短气泡组合成超长垃圾 backbone 的机会。

### 4.3 评分、诊断和其他模式

Rust 版在 GeneMiner2 的 read-slice 验证基础上增加唯一定位支持、覆盖碱基数、breadth、最大缺口、unique-read density、k-mer depth CV 和 max/median ratio 等 UCE 护栏，并输出 `uce_assembly_summary.csv`。可选 GFA/DOT 图会把入度和出度均为 1 的线性链压缩为 unitig，方便观察真正的分叉点；当前主干遍历本身通过连续消费线性节点加速，并不是先把整个计算图替换成 unitig 图。

`--uce-path-strategy search` 仍保留分支枚举，便于与旧策略做 A/B 测试。ITS2 模式则保留多条候选，并统计 fragment support、paired-fragment support、diagnostic support、equivalence group 和 EM abundance；这部分不是 UCE 单 backbone 的逻辑。

## 5. 新版参考的三个 assembler 分别贡献了什么

| 参考项目 | 原项目的核心思想 | GeneMiner2-UCE 实际借鉴 | 对 UCE 的具体贡献 | 明确没有照搬的部分 |
|---|---|---|---|---|
| MaSuRCA | 沿唯一 k-mer 邻接关系延伸 read，形成更长的 super-read，以减少后续问题规模 | 无分叉区域连续延伸；把决策集中在真正分叉处 | 大部分 UCE core 和可信 flank 是线性的，不必每个 k-mer 都重复做完整分支决策 | 没有构建 MaSuRCA super-reads，没有 CABOG/OLC 阶段，也没有做全基因组 scaffolding |
| SPAdes | multisized de Bruijn graph、tip/bulge 简化，以及在删除 bulge 时保存映射信息 | 把气泡当成局部竞争路径；比较有限长度与累计支持后选择一条 backbone | 避免只凭路口第一步深度决定，能跨过短噪声气泡，同时把搜索范围限制在局部 | 没有实现 multi-k 图、BayesHammer、bulge corremoval 投影、paired-distance graph 或 SPAdes 全套图简化流程 |
| Sparrowhawk | Rust 短读长 assembler；紧凑 k-mer 表示、批量提取与 sort/count、图压缩和图导出 | Rust 工程结构、分块读取、排序后计数、紧凑整数 k-mer，以及 GFA/DOT 诊断图的设计启发 | 降低 Python 对象和随机字典更新开销；让逐 locus 图的内存上限更容易控制和诊断 | 没有复制 Sparrowhawk 代码；没有采用 ntHash 作为 k-mer 身份、Bloom filter、petgraph 图或自动 k-mer spectrum 拟合 |

### 5.1 MaSuRCA 的贡献：直道就直接走

[MaSuRCA 论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3799473/)的关键做法是：只要一个 k-mer 只有唯一的合法后继，就继续向前延伸，把大量 reads 压缩为数量更少、长度更长的 super-reads。GeneMiner2-UCE 没有复制它的完整流水线，但接受了一个非常适合 UCE 的原则：**线性区域不值得反复搜索，计算应集中在分叉点。**

在本项目中，这表现为 backbone 内部连续消费唯一后继，以及诊断输出中的 unitig 压缩。贡献主要是减少无意义决策，而不是让 GeneMiner2-UCE 变成 super-read/OLC assembler。

### 5.2 SPAdes 的贡献：气泡不能只看第一步

[SPAdes 论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3342519/)把短而相似的替代路径视作 bulge，并结合覆盖度与图结构进行修正和移除。GeneMiner2-UCE 借鉴的是局部 bulge 判断，而不是 SPAdes 的完整算法：在每条候选出边上向前看有限步数，先比较能走多远，再比较累计支持。

这对 UCE 很重要，因为捕获探针附近的真实 SNP、测序错误或 flank 重复都可能形成短气泡。只看第一条边的深度容易选错；枚举所有气泡组合又可能爆炸。有限前瞻位于两者之间：它比单步贪婪多看一些上下文，但不会无限 backtrack。

### 5.3 Sparrowhawk 的贡献：把算法写成跑得稳的 Rust 工程

[Sparrowhawk](https://github.com/bacpop/sparrowhawk)是 Rust 编写的轻量短读长 de Bruijn graph assembler。对本项目最有价值的不是替换 GeneMiner2 的生物学判定，而是说明 Rust 中可以怎样组织 k-mer 热路径：紧凑整数、滚动提取、先批量排序再计数、图压缩和标准图格式输出。

GeneMiner2-UCE 因此采用了精确 `u128` 2-bit k-mer、分块读取和局部 sort/count。不过它没有使用 Sparrowhawk 的 ntHash 作为节点身份。这样做是有意的：UCE 是先按 locus 过滤后再组装，每张图通常比全基因组图小，优先保证精确可复现比用概率哈希再处理碰撞更合适。

## 6. 新旧算法总表

| 比较项 | 上游原版 Python | Rust backbone 新版 |
|---|---|---|
| 主要目标 | 普通参考引导基因恢复 | 快速、受控地恢复 UCE core 与可信 flank |
| k-mer 表示 | Python 大整数与字典 | 精确 2-bit `u128`，`k ≤ 63` |
| reads 读取 | locus 文件整体处理 | 分块流式读取 |
| k-mer 汇总 | 逐条更新字典 | 块内排序计数，再汇总图 |
| 分叉处理 | 压栈并 backtrack | 有限前瞻，提交一条分支，不 backtrack |
| 每侧候选 | 最多 3 | backbone 通常 1；`search` 模式可多候选 |
| read 支持 | 中间 slice 总数和跨度 | 唯一定位、覆盖碱基、breadth、gap、density，并写入摘要 |
| 重复区护栏 | 深度上限与最小覆盖 | 增加 depth CV、max/median ratio 等护栏 |
| 图诊断 | 无标准图输出 | 可选压缩 GFA/DOT |
| 最适合 | genome skimming、普通目标基因 | UCE 主流程和 ITS2 专用流程 |
| 主要风险 | 气泡多时搜索慢；分支组合增加 | 单 backbone 可能丢掉低支持的真实等位路径 |

## 7. 为什么 reference 默认仍使用原版

普通基因恢复与 UCE 的目标不同。普通模式通常希望忠实保留 GeneMiner2 已验证的参考引导搜索行为；目标区域较短时，适度枚举分支的成本可接受，而且能保留不同候选。UCE 更关心数百或上千 loci 的稳定吞吐、core 两侧 flank 的连续 read 支持，以及避免气泡组合产生异常长 contig。

Rust 原版兼容版在启用 `--reuse-reference-cache` 时，会为每个参考 locus 和 k 写入独立的二进制 k-mer cache。header 记录 magic、格式版本、实现版本、k、参考规范路径、大小和修改时间；条目保存 u128 k-mer、depth、position 与 strand。任何不匹配、截断或解析错误都会触发重建，写入采用临时文件加原子替换。

因此当前策略是：

- `reference + auto`：直接使用上游原版，保证普通模式行为稳定；
- `reference + rust`：仅在用户显式要求 A/B 测试时使用 Rust；
- `reference + original-rust`：可选的单线程 Rust 原版兼容后端；保留原版回溯逻辑，用于速度与结果 A/B。它不属于 UCE backbone，不会取代默认原版；在并列 k-mer 的极端重复序列中，容器遍历顺序可能选到同权重的另一 seed，因此应先在目标数据上比较。
- `uce + auto`：只使用 Rust backbone；Rust 不可用或失败时直接报错；
- `its2`：强制 Rust，因为多候选 fragment 统计没有 Python 等价实现。

### 7.1 已完成的 Python 与 `original-rust` 对照

为避免把“跑得快”误当成“结果一样”，我们用同一批已过滤的 UCE reads 做过一次单线程 A/B：40 个 locus、相同参考、`k=31`、`k_min=21`、`k_max=51`、`limit_count=2`、`iteration=4096`、`cov_min=0`。这是一组开发验证数据，不应代替正式数据集上的准确性评估。

| 指标 | 上游 Python 原版 | `original-rust` | 结论 |
|---|---:|---:|---|
| 状态：成功 / low quality / 过滤不足 | 30 / 9 / 1 | 30 / 9 / 1 | 分类完全一致 |
| 写出的最佳 contig | 39 | 39 | 完全一致 |
| 最佳 contig 碱基序列 | — | 38/39 与 Python 完全相同 | `v1__uce-1200` 同为 1119 bp，但序列和 read-support 不同 |
| `contigs_all` 候选序列流 | — | 35/39 与 Python 相同 | 4 个 locus 的候选路径集合不同；其中 3 个的最佳序列仍相同 |
| 墙钟时间 | 0.48 s | 0.05 s | 小型数据下约快 9.6 倍 |
| 峰值常驻内存 | 48,200 KB | 5,432 KB | 小型数据下约低 8.9 倍 |

`result_dict.txt` 的 locus 输出顺序不同，是 Python 依文件遍历顺序写入、Rust 按 locus 名排序的实现差异；按 locus 名排序后，只有 `v1__uce-1200` 的 support 值不同（Python 为 45，Rust 为 50）。39 个最佳 FASTA 中另有 30 个只在表头的候选坐标或统计字段上不同，碱基序列相同。

同一轮还验证了 `original-rust` 的参考 k-mer 二进制 cache：无 cache、首次建 cache、热 cache 的 39 个最佳 FASTA 均逐字节一致。40 个 locus 的小参考集中，首次建 cache 为 0.07 s、热 cache 为 0.05 s；热 cache 与无 cache 都约 0.05 s，说明这个规模不足以证明稳定的性能收益。更大的参考集或重复样本运行时，仍应重新计时。

这次验证支持将 `original-rust` 保留为高效的 A/B 后端，但**不支持**把它替代 `reference + auto` 的上游 Python 原版。若要进一步追求逐字节兼容，应优先追踪 `v1__uce-1200` 的 seed 选择、分叉候选排序和 support 累计。

## 8. 应怎样验证这些改动

新版是否更好不能只看运行时间或 contig 长度。建议同一批 UCE 数据同时运行 `backbone` 与 `search`，至少比较：

1. 每个 locus 是否恢复成功以及缺失率；
2. contig 长度、UCE core 覆盖和两侧 flank balance；
3. unique-read density、supported-bases breadth 与最大无支持缺口；
4. k-mer depth CV 和 max/median ratio；
5. 多样本比对中的异常长分支、疑似 paralog 和 SNP 缺失率；
6. population 伪参考建成后，各样本 mapping rate、coverage breadth 和 depth 是否更均衡。

只有当 backbone 在明显缩短运行时间的同时，没有降低这些 read-supported 与群体一致性指标，才能说它对该数据集真正有帮助。

## 9. 代码对应位置

- 上游原版：`scripts/main_assembler_original.py`
- Rust 原版兼容版：`rust/main_assembler_original/src/main.rs`
- Rust locus 流程与候选验收：`rust/main_assembler/src/pipeline.rs`
- Rust unitig 与 GFA/DOT：`rust/main_assembler/src/unitig.rs`
- Rust 分块读取与参考 cache：`rust/main_assembler/src/io_utils.rs`
- 上游来源和哈希：`scripts/main_assembler_original.UPSTREAM.md`

## 参考资料

- Zimin et al. 2013. [The MaSuRCA genome assembler](https://pmc.ncbi.nlm.nih.gov/articles/PMC3799473/).
- Bankevich et al. 2012. [SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing](https://pmc.ncbi.nlm.nih.gov/articles/PMC3342519/).
- bacpop. [Sparrowhawk: a lightweight short-read assembler](https://github.com/bacpop/sparrowhawk).
- sculab. [GeneMiner2 upstream `main_assembler.py`](https://github.com/sculab/GeneMiner2/blob/36e06feeb99654bdb87f45d4cde225d8c3e311d0/scripts/main_assembler.py).
