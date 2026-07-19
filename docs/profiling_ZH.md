# 3. Profiling

[English version](profiling_EN.md)

Profiling 是从 WGS 或 metagenome reads 中免组装获取并定量任意扩增子 marker 的流程。它报告 marker group 的相对信号，不生成 contig，也不是生物体或细胞比例。

## 流程

```text
一次 GeneMiner2 k-mer 招募 → Themisto 伪比对 → mSWEEP group 估计
```

不会运行 `refilter`、`assemble`、`combine` 或 `tree`。`--profile-kmer-size` 同时设定招募与 Themisto 的 k-mer，必须是 15–31 的奇数。

## 必需输入

- 通过 `-r` 直接提供一个 `.fa` 或 `.fasta` marker 参考库。
- 两列 `--profile-group-map` TSV：`reference_id<TAB>group`。
- 含 WGS 或 metagenome reads 的样本表。

reference ID 是 FASTA 标题第一个空白前字段。每条参考必须映射至一个 group；重复行仅在 group 完全相同时允许。`themisto` 与 `mSWEEP` 须位于 `PATH`，或以 `--profile-themisto`、`--profile-msweep` 显式提供。可选的 `--profile-decoy` FASTA 用于表示可能的非目标序列。

## 运行与 cache

```bash
cli/geneminer2 profiling \
  -f samples.tsv -r marker_reference.fasta \
  --profile-group-map marker_groups.tsv \
  -o output -p 8 --profile-decoy non_target.fasta
```

Themisto index 使用内容寻址 cache。用 `--profile-index-dir` 在多次运行间共享；只有需要重建时才加 `--profile-force-rebuild`。

## 输出与解释

每个样本写入 `marker_profile/`：

- `marker_group_abundance.tsv`：group、证据、检出状态和相对信号。
- `marker_qc.tsv`：伪比对、target/decoy、证据、mSWEEP 和参数计数。
- `marker_reference_metadata.tsv`：reference ID、Themisto color 与 group 映射。

`evidence_queries` 与 `exclusive_queries` 按单条 FASTA/FASTQ query record 计数，并非 paired fragment。`relative_proportion` 在独占证据规则后重新归一化，仍是未校准的 marker 信号比例。

## QC 与校准

解释检出结果前，应检查伪比对 query 数、target-decoy 共享证据、decoy 独占证据和各 group 的独占支持。应以匹配 marker 与测序深度的阴性对照、混合样与 downsampling 设定 `--profile-min-evidence`。跨样本比较时，固定参考库、group map、k-mer、阈值与 decoy 策略。

参见[Filter](filter_ZH.md)、[输出说明](../manual/ZH_CN/output.md)和[命令行指南](../manual/ZH_CN/command_line.md)。
