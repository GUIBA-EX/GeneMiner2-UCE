# MainFilter 性能优化与兼容性说明

本文记录 Rust `MainFilterNew` 已采用的确定性优化、实测边界和未采用方案。所有已采用优化都以保持命令行参数、k-mer 判定和筛选输出不变为前提；字典格式升级时必须安全失效并重建；性能收益必须经真实 reads 的 A/B 比较确认。

## 已采用的确定性优化

| 修改 | 位置 | 原理 | 兼容性边界 |
| --- | --- | --- | --- |
| DNA 碱基查表 | `base_code` | 用固定的 256 项表把 `A/a/C/c/G/g/T/t/U/u` 转成 2-bit 编码，避免每个碱基做大小写转换和多分支匹配。 | 非 A/C/G/T/U 的字符仍会中断该 k-mer；与原规则一致。 |
| 无取模采样扫描 | `KmerIndex::collect_hits`（`k <= 32`） | 用单调递增的 `next_probe` 记录下一个 `0, step, 2*step...` 采样起点，替代热循环里的 `start % step == 0`。 | 仍检查相同的全局起点，并始终补查 read 尾端 k-mer；`N` 或其他非标准碱基不会重置采样坐标。 |
| `AHashMap` k-mer 索引 | `KmerStore` | 仅将内存中的 k-mer 到 locus 命中表从标准 `HashMap` 改为 `AHashMap`，减少哈希查询成本。 | 哈希表迭代顺序不参与筛选结果；缓存二进制排列不保证稳定。 |
| 按输出模式保留 FASTQ 文本行 | `SequenceReader`、`Record` | 默认 GM2（`-m 5`）和只扫描（`-m 3`）只保留序列与质量值；跳过后续不会使用的 header、`+` 行和规范化文本副本。 | 文本输出模式 `-m 0/1/4` 仍完整保留并写出所有 FASTQ/FASTA 行；筛选判定和 GM2 编码不变。 |
| Canonical k-mer 索引 | `KmerIndex` | 每个 k-mer 与其反向互补统一为同一 key，单次查询即可覆盖双链。 | `-gr` 保留为兼容参数；旧策略与 canonical 策略在 A/C/G/T/U 窗口上的招募集合相同。 |
| Packed locus postings | `ReferenceHits`、`packed_hits` | 多 locus 命中在建索引后压为连续 ID 切片，map value 只保留 offset/length。 | ID 顺序与输出顺序不变；范围超过 `u32` 时明确报错。 |
| v3 内容验证 cache | `reference_content_hash`、dictionary loader | cache 写入 canonical 标志、`k` 和参考内容 SHA-256；首次建索引时哈希与索引共用一次参考遍历。 | 旧、损坏或参考不匹配 cache 自动重建，不会静默复用。 |

默认 UCE 设置通常为 `k=31`。因此 DNA 查表和无取模扫描覆盖常用短 k-mer 路径；按输出模式保留文本行适用于所有 k-mer 长度。

## 实测结果

基准使用 Li et al. UCE 数据中的 DK40（`SRR29729138`）target-capture 样本。输入为 `fastp` 质控后的双端 reads，参考为 Bivalve 2k UCE probe loci，参数为 `k=31`、`step=4`；每次运行最多处理前 1,000,000 个 read pairs。

### 仅扫描

以下为仅扫描模式（`-m 3`）的三次运行平均 wall time：

| 实现 | 平均时间 | 相对原始 Rust MainFilter |
| --- | ---: | ---: |
| 原始 `std::HashMap` | 5.15 s | 基线 |
| 无取模采样 | 5.00 s | 约快 2.9% |
| 加 DNA 查表 | 3.60 s | 约快 30.2% |
| 加 `AHashMap`（未含低分配读取） | 3.07 s | 约快 40.5% |

### 默认 GM2 输出

在默认输出模式（`-m 5`）下，低分配读取将过滤阶段从 3.65–3.73 s 降至 3.39 s，约快 7–9%。旧版与候选版生成的 4,466 个 GM2 文件及 `ref_reads_count_dict.txt` 均逐字节一致。

这证明该改动在此输入、参数和输出模式下保持结果不变；它不构成对所有测序策略、参考库或参数组合的生物学等价性声明。

### 解释边界

上述数据来自 **target capture**，不能外推为 genome-skimming 的性能承诺。对 genome-skimming，应选取有代表性的样本，在相同参考、`k` 和 `step` 下进行至少三次计时，并逐字节比较每个输出文件。

## 明确未采用的方案

| 方案 | 结论 | 原因 |
| --- | --- | --- |
| LRU 持久输出文件句柄 | 不采用 | 在上述实际输出运行中，`openat`/`close` 合计约占筛选时间 0.5%，不足以抵消句柄淘汰、写入顺序和错误处理的复杂度。 |
| MainFilter 内部多线程或共享跨样本索引 | 不采用 | 增加并发行为和 I/O 争用；当前优先保留容易逐字节验证的单线程实现。 |
| `FxHashMap` | 不采用 | 单独 A/B 测试仅比 `AHashMap` 快约 1.6%，未达到 5% 的采用阈值。 |
| 编码缓冲区复用 | 不采用 | 实际输出模式的 A/B 测试没有稳定收益。 |
| 改大 `step` | 不属于性能优化 | 会改变采样密度和潜在检出率，必须用 UCE recovery/准确性验证决定，而不能只看运行时间。 |
| SSHash/GGCAT 重构 | 暂不采用 | 这不是简单替换哈希表：需要静态 k-mer ID 到 locus posting 的新索引和缓存格式，应在确认索引查询仍是主瓶颈后独立设计和验证。 |

## 修改后的验证要求

修改 `MainFilterNew` 后，至少执行：

```
cargo test --manifest-path rust/main_filter_new/Cargo.toml
cargo clippy --manifest-path rust/main_filter_new/Cargo.toml -- -D warnings
cargo build --release --manifest-path rust/main_filter_new/Cargo.toml
```

若改动涉及筛选逻辑、k-mer 编码、缓存或输出路径，还必须用固定的真实双端样本分别运行旧二进制与候选二进制，并比较每个 locus 输出及 `ref_reads_count_dict.txt`：

```
diff -qr old_run/filtered_pe new_run/filtered_pe
cmp old_run/ref_reads_count_dict.txt new_run/ref_reads_count_dict.txt
```

只有在输出一致、测试通过且真实数据计时有稳定收益时，才应将候选二进制安装到 `cli/bin/MainFilterNew` 并发布。
