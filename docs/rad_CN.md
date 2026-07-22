# RAD：用 WGS 补充已有 ipyrad loci

`rad` 用 WGS reads 补充 ipyrad `.loci` 中缺失的新样本；它不从缺失的 RAD reads 伪造等位基因。

| 起点 | 产物 | 边界 |
| --- | --- | --- |
| 已完成的 `.loci` 或已拆样本的 paired RAD reads + 新 WGS reads | 独立 R1/R2 arm、恢复状态、验证后的严格矩阵 | 不推断 R1/R2 中间插入区；不直接证明 allele dropout |

## 推荐路线

```bash
# 从已完成 ipyrad .loci 建立可复用 probe
cli/geneminer2 rad-probe --ipyrad-loci assembly.loci -o rad_probe

# 为新 WGS 样本恢复独立的两条 arm
cli/geneminer2 rad --rad-probe rad_probe/rad_reference \
  -f wgs_samples.tsv -o rad_out -p 8

# 写入独立验证后的 strict 矩阵
cli/geneminer2 rad-validate --rad-probe rad_probe/rad_reference \
  --rad-recovery rad_out/rad_recovery -o rad_validate_out
```

## 输入与 probe 构建

`rad-probe` 支持以下输入路径：

- `--ipyrad-loci FILE`：已完成的 `.loci`，推荐；
- `--ipyrad-params FILE`：调用 ipyrad，默认执行 `1--7`，再读取 `project_dir/<assembly>_outfiles/<assembly>.loci`；可用 `--ipyrad-executable` 和 `--ipyrad-steps` 调整。若 ipyrad 输出被移动，可同时用 `--ipyrad-loci` 指定其路径；
- `--rad-denovo -f paired_rad_samples.tsv`：已经 demultiplex 的 paired RAD reads。

`rad-probe` 不重写 demultiplex、酶切识别、去接头或复杂样本内/样本间 clustering；这些应由 ipyrad 参数明确控制。de novo 模式是保守的候选 probe builder，不等价于 ipyrad 的严格 RAD 矩阵。

```text
paired_rad_samples.tsv
sample_id<TAB>R1.fastq.gz<TAB>R2.fastq.gz
```

原生 de novo 模式按样本计数 canonical k-mer、保留 solid k-mer，并按成对 arm 的 minimizer 建立候选 stack。它只输出满足最小 reads 深度、最短长度和跨样本支持的 paired arm。每个样本—stack 最多保留 128 对确定性抽样 reads 建 consensus，同时记录总支持 pair 数。默认 `k=31`、solid 深度 `3`、跨样本支持 `2`、最短 arm `60 bp`。`--rad-overhang` 只检查 R1；仅当 R2 从已知第二酶切端开始时才传 `--rad-overhang-r2`。

## 恢复与验证

两种 probe 都复用 MainFilter、refilter 和 `original-rust` 组装。`rad` 仅接受新的 WGS sample：若 sample 已在 RAD 参考中，R1/R2 sample 集合不一致，或出现重复归一化名称与损坏 FASTA，流程会停止而不会静默合并 header。

`--rad-min-arm-breadth` 默认 `0.80`。`rad` 生成的是 recovery matrix；应随后运行 `rad-validate`。验证将每个 candidate arm 与其 locus 的多等位 bait 比较，并与所有其他 locus 的同 arm bait 竞争。默认要求 query 与参考覆盖均至少 `0.80`、identity 至少 `0.90`、本 locus 得分相对最佳外源 locus 高至少 `5%`。只有 R1、R2 都通过的 WGS sample 才进入 strict 矩阵。可用 `--rad-validate-min-breadth`、`--rad-validate-min-identity`、`--rad-validate-min-delta` 调整门槛。

## 输出与解释

- `rad_reference/arms/`：每 locus 的 R1/R2 多等位 bait FASTA；
- `rad_recovery/`：每个 WGS sample 的招募、refilter 与组装结果；
- `rad_matrix/rad_sample_locus.tsv`：sample × locus 的 arm 与联合恢复状态；
- `rad_matrix/recovered_arms/`：保留有证据单 arm 的未比对 FASTA；
- `rad_matrix/paired_arms/`：两个 arm 都恢复的 WGS sample，尚未 validate；
- `rad_validated/rad_validation.tsv`：逐 sample × locus × arm 的验证指标和状态；
- `rad_validated/strict_arms/`：原始 RAD bait 加上两个 arm 都验证通过的 WGS sample。

R1/R2 永远是独立观察，流程不会桥接中间未测序区。所有 arm FASTA 都未比对；应在确定缺失数据策略后再进行 MSA。`rad-validate` 不重新招募或组装 reads，也不会修改原始 `rad_recovery/` 或 `rad_matrix/`。

`rad_missing_wgs_recovered` 只表示输入 RAD 矩阵缺失的 locus 已由 WGS 恢复；它本身不是 restriction-site allele dropout 的直接证据。若需作该解释，还应有酶切位点和跨 locus WGS 证据。
