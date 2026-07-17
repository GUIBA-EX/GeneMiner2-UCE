# MainFilter 性能优化与兼容性说明

本说明记录当前 Rust `MainFilterNew` 的确定性优化、验证范围和暂不采用的方案。目标是降低 UCE 读段筛选的 CPU 开销，**不改变**命令行参数、字典格式、筛选规则或输出格式。

## 当前保留的修改

| 修改 | 位置 | 原理 | 兼容性边界 |
| --- | --- | --- | --- |
| DNA 碱基查表 | `base_code` | 用固定的 256 项表把 `A/a/C/c/G/g/T/t/U/u` 转成 2-bit 编码，避免每个碱基做大小写转换和多分支匹配。 | 非 A/C/G/T/U 的字符仍会中断该 k-mer；与原规则一致。 |
| 无取模采样扫描 | `KmerIndex::collect_hits`（`k <= 32`） | 用单调递增的 `next_probe` 记录下一个 `0, step, 2*step...` 采样起点，替代热循环里的 `start % step == 0`。 | 仍检查相同的全局起点，并始终补查 read 尾端 k-mer；`N` 或其他非标准碱基不会重置采样坐标。 |
| `AHashMap` k-mer 索引 | `KmerStore` | 仅将内存中的 k-mer 到 locus 命中表从标准 `HashMap` 改为 `AHashMap`，减少哈希查询成本。 | 现有 v2 缓存仍可读取，缓存结构未改变；重新生成时条目排列可能不同，因此不承诺缓存二进制逐字节一致。哈希表的迭代顺序不参与筛选结果。 |

长 k-mer 路径（`k > 32`）没有改动。默认 UCE 设置通常为 `k=31`，因此上述前两项覆盖常用路径。

## 实测结果

基准使用 Li et al. UCE 数据中的 DK40（`SRR29729138`）target-capture 样本。输入为 `fastp` 质控后的双端 reads；参考为 Bivalve 2k UCE probe loci；参数为 `k=31`、`step=4`。为使不同版本可重复比较，运行限制为前 1,000,000 个 read pairs。

仅扫描（`-m 3`）的三次运行平均 wall time：

| 实现 | 平均时间 | 相对原始 Rust MainFilter |
| --- | ---: | ---: |
| 原始 `std::HashMap` | 5.15 s | 基线 |
| 无取模采样 | 5.00 s | 约快 2.9% |
| 加 DNA 查表 | 3.60 s | 约快 30.2% |
| 加 `AHashMap`（当前） | 3.07 s | 约快 40.5% |

在实际筛选输出模式（`-m 5`）下，当前实现与改动前实现产生的 4,466 个 GM2 文件以及 `ref_reads_count_dict.txt` 均逐字节一致。该比较证明本次改动没有改变此数据、此参数组合下的输出；它不是对所有测序策略或所有参数的生物学等价性证明。

该基准是 **target capture**，不能外推成 genome-skimming 的性能承诺。后者应另选有代表性的样本，以相同参考、k 和 step 进行输出逐字节比较和至少三次计时。

## 明确未采用的方案

| 方案 | 结论 | 原因 |
| --- | --- | --- |
| LRU 持久输出文件句柄 | 不采用 | 在上述实际输出运行中，`openat`/`close` 合计约占筛选时间 0.5%，不足以抵消句柄淘汰、写入顺序和错误处理的复杂度。 |
| MainFilter 内部多线程或共享跨样本索引 | 不采用 | 增加并发行为和 I/O 争用；当前优先保留容易逐字节验证的单线程实现。 |
| `FxHashMap` | 不采用 | 单独 A/B 测试仅比 `AHashMap` 快约 1.6%，未达到 5% 的采用阈值。 |
| 编码缓冲区复用 | 不采用 | 实际输出模式的 A/B 测试没有稳定收益。 |
| 改大 `step` | 不属于性能优化 | 会改变采样密度和潜在检出率，必须用 UCE recovery/准确性验证决定，而不能只看运行时间。 |
| SSHash/GGCAT 重构 | 暂不采用 | 这不是简单替换哈希表：需要静态 k-mer ID 到 locus posting 的新索引和缓存格式，应在确认索引查询仍是主瓶颈后独立设计和验证。 |

## 维护与复核

修改 `MainFilterNew` 后，至少执行：

```
cargo test --manifest-path rust/main_filter_new/Cargo.toml
cargo clippy --manifest-path rust/main_filter_new/Cargo.toml -- -D warnings
cargo build --release --manifest-path rust/main_filter_new/Cargo.toml
```

涉及筛选逻辑、k-mer 编码、缓存或输出路径时，还应选择固定的真实双端样本，分别运行旧二进制和候选二进制，并比较每个 locus 输出及 `ref_reads_count_dict.txt`：

```
diff -qr old_run/filtered_pe new_run/filtered_pe
cmp old_run/ref_reads_count_dict.txt new_run/ref_reads_count_dict.txt
```

只有在输出一致、测试通过且真实数据计时有稳定收益时，才应将候选二进制安装到 `cli/bin/MainFilterNew` 并发布。
