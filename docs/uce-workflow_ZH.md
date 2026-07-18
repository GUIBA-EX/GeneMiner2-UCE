# UCE 流程说明

[English version](uce-workflow_EN.md)

本文解释 `--assembly-mode uce` 的组装策略、reads 支持、raw-read rescue、质量护栏和参考缓存。参数全集见[中文命令行指南](../manual/ZH_CN/command_line.md)，输出字段见[输出文件说明](../manual/ZH_CN/output.md)。

## 默认流程

未显式指定子命令时，UCE 模式运行：

```text
filter → refilter → assemble → combine → tree
```

它默认跳过 `trim`，避免刚恢复的侧翼序列再次被裁回 probe 范围。需要按参考切齐时可显式加入 `trim`。

refilter 阶段只要任一 mate 通过 locus 过滤，整对 paired-end reads 都会保留。这样一端可以命中保守 UCE core，另一端继续提供侧翼信息。

## Backbone 组装

默认 `--uce-path-strategy backbone` 在 seed 两侧各构建一条主干。遇到气泡时，assembler 对每条出边执行受 `--uce-backbone-lookahead` 限制的前瞻，优先选择能继续延伸更远的分支，再以累计 k-mer 支持度打破平局。

选定分支后，同一气泡的其他出边会被丢弃；算法不保存分支栈，也不反复 backtrack。已经访问的 k-mer 不会再次进入，因此环路会终止。无分叉链会作为 unitig 连续消费，减少逐节点决策开销。

如需与旧分支枚举算法比较，可使用：

```bash
--uce-path-strategy search \
  --uce-side-candidates 8
```

`--uce-side-candidates` 只影响 `search` 策略。

## 候选支持和质量护栏

候选 contig 不只按长度排序，还记录：

- 总 reads、唯一定位 reads 和重复定位 reads；
- reads 支持的碱基并集和覆盖广度；
- 最大无支持间隔；
- 唯一 reads density；
- k-mer depth 的中位数、变异系数和 max/median ratio。

UCE 模式默认不限制 contig 长度。长度至少为 1000 bp 的候选默认需要满足：

```text
unique_read_density = uniquely_placed_read_count / contig_length
unique_read_density >= 0.003
```

相关参数为：

```bash
--uce-max-contig-length 0 \
  --uce-min-read-density 0.003 \
  --uce-density-check-min-length 1000 \
  --uce-max-depth-cv 0 \
  --uce-max-depth-ratio 0
```

长度上限、depth CV 和 depth ratio 默认关闭。只有在数据中确实存在异常延长、重复峰或深度极不均匀的候选时才建议启用。

## 放宽短 probe 参数

对于短 bait 或样本与参考存在一定分化的情况，可测试：

```bash
-sb unlimited \
  -ka 0 \
  --min-ka 17 \
  --max-ka 31 \
  -e 1
```

这些设置会放宽边界并允许较低 assembly k-mer，也可能增加弱支持或错误候选。调整后应检查 `uce_assembly_summary.csv`、下游 alignment 和不同参数间的结果稳定性。

## Raw-read rescue

`--uce-rescue-reads` 只执行一轮额外招募：

1. 合并原始 locus reference 与第一轮合格 contig；
2. 从 raw reads 中重新捕获匹配 reads；
3. 重新执行 refilter 和 assemble；
4. 比较第一轮与 rescue 结果并决定保留或回退。

长度小于 `--uce-rescue-min-contig-length`（默认 60）的第一轮 contig 不作为 rescue reference。rescue 最多同时处理 4 个样本，每个样本最多使用 4 个线程，并受全局 `-p` 约束。

第一轮已接受而 rescue 缺失或不合格时，恢复第一轮结果并标记 `reverted_failed_rescue`。两轮都合格时比较唯一 reads density：

```text
before_density = before_unique_read_count / before_contig_length
rescue_density = rescue_unique_read_count / rescue_contig_length
density_ratio = rescue_density / before_density
```

默认 `density_ratio < 0.5` 时恢复第一轮结果并标记 `reverted_density_drop`。阈值由 `--uce-rescue-min-density-ratio` 控制。`uce_rescue_summary.csv` 的 `after_*` 字段始终描述 rescue 尝试；发生回退时最终 FASTA 仍保存第一轮 contig。

## Reference cache

重复使用相同参考时可加入 `--reuse-reference-cache`。缓存默认位于 `output/.gm2_reference_cache`，也可通过 `--reference-cache-dir` 指定共享目录。

缓存指纹包含参考文件名、大小、修改时间、filter k-mer 和扫描步长。缓存只减少重复建索引时间，不改变筛选或组装规则。rescue reference 含样本特异的第一轮 contig，因此 rescue index 仍按样本构建。

## 推荐检查

正式使用 UCE 结果前至少检查：

- `uce_assembly_summary.csv` 中的接受状态、拒绝原因和 reads 支持；
- 启用 rescue 时的 `uce_rescue_summary.csv`；
- `contigs_all_low/` 中被拒绝候选是否呈现异常延长或低支持；
- 合并 alignment 中的缺失、异常长分支和疑似 paralog；
- `stats` 输出的样本与 locus 恢复率。
