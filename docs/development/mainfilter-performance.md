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
| 长 k-mer scratch buffers | `long_kmer_into`（`k > 32`） | 复用正反链缓冲；索引中已有 key 时不复制序列；采样扫描直接迭代起点而不创建索引 `Vec`。 | canonical key 与尾端补查规则不变；新 key 仍只复制一次以进入索引。 |
| 定向输出刷新 | `OutputManager` | 达到总缓冲预算时只刷新当前最大的 buffer，保留 64 MiB 内存上限，避免一次刷新所有冷文件。 | 每个文件内容与记录顺序不变；仅改变写入时机。 |
| 编码与日志缓冲复用 | `encode_*_into`、`Logger` | GM2/text 编码复用 scratch buffer；日志采用 64 KiB 缓冲并按时间或退出时刷新。 | 编码格式和日志内容不变。 |
| 常驻输出文件句柄 | `OutputManager`（`handles` 字段） | 每个输出路径只在首次 flush 时 `open` 一次，此后复用同一个 `File` 句柄，不再每次 flush 都重新 open/close；进程启动时额外尝试把 `RLIMIT_NOFILE` 抬到硬上限，避免 locus 数一多、句柄常驻同时打开时撞到系统上限。 | 文件内容与写入顺序不变；只改变 flush 涉及的系统调用次数，此前「明确未采用」一节给出的结论已被更大规模的真实数据推翻，见下方实测结果。 |
| 字节级行读取 | `SequenceReader::read_line`、`Record.lines` | 用 `read_until(b'\n', …)` 直接读进 `Vec<u8>`，手动去掉尾部 `\r\n`，不再经过 `String`（省掉一次 UTF-8 校验和 `trim().to_string()` 的额外拷贝）。 | 只精确剥掉行尾换行符，不再像 `.trim()` 那样顺带清掉行首尾其他空白字符；正常 FASTA/FASTQ 不受影响，含异常内嵌空白的畸形输入容错会略弱于旧实现。 |
| 放大输入端缓冲 | `open_input`、`GzipReader::open`（`gzbuffer`） | 把外层 `BufReader` 容量和 zlib 内部压缩端缓冲都从默认 8 KiB 提到 1 MiB（`READ_BUFFER_SIZE`），减少 `gzread`/`read` 调用频次。 | 不改变读到的字节内容，只改变每次系统调用/FFI 调用搬运的数据量。 |
| zlib-ng 运行时探测 | `detect_zlib_ng`、`ZlibBackend` | 启动时用 `dlopen`/`dlsym` 查找 `libz-ng.so{.2,.1,}` 里的 `zng_gzopen/zng_gzread/zng_gzclose/zng_gzbuffer` 四个符号，全部找到就用 zlib-ng 的 SIMD 加速实现，否则退回静态链接的系统 zlib；结果用 `OnceLock` 缓存，全程只探测一次。 | 两种后端解压出的字节逐位相同（已用真实数据 A/B `diff` 验证）；环境里没有 zlib-ng 时行为、产物与之前完全一致，不引入新的构建期依赖。 |

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

### v1.3.1 合成压力测试

在 release 构建、相同机器的合成 FASTQ A/B 测试中，`k=33`、50 万 150 bp reads、仅扫描模式从 **1.11 s** 降至 **1.01 s**（约快 **9%**）。在 128 loci、约 534 MB GM2 分散输出的压力场景中，耗时从 **0.28 s** 降至 **0.27 s**；该约 4% 差异接近文件系统噪声，但两版输出总字节数一致。该测试证明长 k-mer 路径的改动有效；输出优化主要用于限制 I/O 抖动和避免全局 flush，而非承诺固定加速比例。

### 常驻输出句柄 + 字节级解析 + 输入缓冲 + zlib-ng

基准同样使用 DK40（`SRR29729138`）双端 reads，参考换成 Bivalve 2k UCE probe loci 全集（3,579 个 locus，比此前测试的参考规模大得多，更能放大按 locus 逐文件输出的开销），`k=31`、默认 `step`，每次运行处理前 3,000,000 个 read pairs，默认 GM2/文本混合的 `-m 0` 全流程输出。

先用 `strace -c` 核对常驻句柄的效果（`-m_reads 1`）：

| 版本 | `openat` 总次数 | 打开输出文件的次数 | 平均每个输出文件被 open 几次 |
| --- | ---: | ---: | ---: |
| 每次 flush 重新 open/close（旧策略） | 11,938 | 8,347（对应 2,237 个不同文件） | 约 3.7 次 |
| 常驻句柄（本次） | 5,826 | 2,235（对应 2,235 个不同文件） | 恰好 1 次 |

两版本都还各含约 3,579 次参考 FASTA 读取（用于校验字典缓存的内容哈希），与本次改动无关，已从上表中的输出文件计数里分离。

再看 `-m 0` 全流程 wall time（每组 2 次重复）：

| 版本 | Filtering 耗时 | 相对上一版 |
| --- | ---: | ---: |
| 原始版本（本轮改动前） | ~10.36 s | 基线 |
| + 常驻输出句柄 | ~10.27–10.41 s | 本机热缓存 ext4 下收益在噪声范围内；句柄常驻的价值主要在高延迟文件系统（NFS/云盘）或句柄反复淘汰的场景，此处以 `openat` 计数为准 |
| + 字节级行读取、放大输入缓冲 | ~9.13–9.30 s | 约快 10–12% |
| + zlib-ng（检测到时，本机通过 conda 环境的 `libz-ng.so.2` 验证） | ~7.39–7.47 s | 再快约 20% |

四版本对同一输入跑 `-m 0/1/4/5`，输出目录 `diff -rq` 逐字节一致（仅 `log.txt` 因计时数字不同）；system-zlib 与 zlib-ng 两个后端各自跑一遍同样逐字节一致。

### 解释边界

上述数据来自 **target capture**，不能外推为 genome-skimming 的性能承诺。对 genome-skimming，应选取有代表性的样本，在相同参考、`k` 和 `step` 下进行至少三次计时，并逐字节比较每个输出文件。zlib-ng 相关收益仅在检测到该库时生效，且收益幅度取决于具体的 zlib-ng 版本与是否启用 SIMD 加速构建；未检测到时行为与纯 system zlib 版本一致。

## 明确未采用的方案

| 方案 | 结论 | 原因 |
| --- | --- | --- |
| MainFilter 内部多线程或共享跨样本索引 | 不采用 | 增加并发行为和 I/O 争用；当前优先保留容易逐字节验证的单线程实现。已确认扫描阶段在 32 核机器上完全单线程，是后续最有希望的优化方向，但需要重新设计 `HitCollector` 的多线程复用和 `OutputManager` 的并发写入合并策略，工作量和验证成本明显更高，暂缓。 |
| `FxHashMap` | 不采用 | 单独 A/B 测试仅比 `AHashMap` 快约 1.6%，未达到 5% 的采用阈值。 |
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
