# TE / repeatome 模式

`geneminer2 te` 是面向 genome-skimming 或 WGS 短读长的参考无关 repeatome 工作流。它输出保守的重复单元、注释证据和样本丰度；**不是**完整 TE 注释、插入位点检测或 TE 系统树流程。

## 运行

```bash
cli/geneminer2 te -f te_samples.tsv -o te_out -p 32
```

TE 使用独立样本表，不需要 `-r`：

```text
taxon_id  sample_id  read1  read2
Taxon_A   A01        /data/A01_R1.fq.gz  /data/A01_R2.fq.gz
Taxon_B   B01        /data/B01.fq.gz
```

第四列 `read2` 可省略。可选的 `--te-read-ledger FILE` 接受 `sample_id read_id` 两列，用于排除已被高置信目标序列占用的 reads。

## 工作流

```text
discover → curate → annotate → quantify
```

- **discover**：以 taxon 为单位均衡抽样，发现 atomic seed groups。
- **curate**：只运行一次 MainFilter，保留 candidate reads 和 paired-end linkage 证据，建立 exact-equivalence（EQ）library。
- **annotate**：为每个 EQ 生成长度受限、受 reads 支持的 fragment；识别 tandem/satellite 特征，并可进行保守的 library 同源注释。
- **quantify**：复用 curate 的 candidate reads，输出样本 RPM 与存在状态；不再次运行 MainFilter。

默认 `--te-stage all` 依次完成四步。也可用 `--te-stage discover|curate|annotate|quantify` 重跑单步；后续阶段会验证上游 manifest，输入或参数不一致时必须从相应上游阶段重跑。

## 可选 library 注释

```bash
cli/geneminer2 te -f te_samples.tsv -o te_out -p 32 \
  --te-library curated_te.fa
```

`--te-library` 是可选的已分类 TE FASTA；header 使用 `name#Class/Subclass`，例如 `RTE1#DNA/TcMar`。高置信分类同时要求足够的 fragment 长度和 read-pair 支持、identity、coverage，以及相对第二类别的分数差。阈值可用 `--te-annotate-min-*` 调整。

注释不会合并或删除 EQ。PE bridges 与 read transitions 始终只作为 `linked_not_merged` 证据；短、低支持或冲突 fragment 保留为 `unknown_repeat` 或 `unknown_interspersed_repeat`。内置比较是保守的 rolling-k-mer/diagonal 筛选，不替代需要基因组 assembly 的 EDTA 或 RepeatModeler。

## 输出与解释

```text
te_out/
├── 01_discover/                 atomic catalog、seeds、manifest
├── 02_curate/                   EQ library、candidate reads、linkage、manifest
├── 03_annotate/
│   ├── fragments/EQ00001.fasta  read-supported fragment（若可生成）
│   ├── annotation_evidence.tsv  fragment 与同源证据
│   ├── annotated_catalog.tsv    最终类别、置信度和决策
│   └── manifest.tsv
└── 04_quantify/
    ├── repeat_signal.tsv        每 sample × EQ 的定量与注释
    └── taxon_repeat_matrix.tsv  taxon 中位 RPM 与状态
```

`signal_rpm = 1,000,000 × specific_pairs / effective_pairs`。`PRESENT` 要求至少 100 个 effective pairs、3 个 specific pairs，且 specific 支持比例不少于 0.70。对于 UCE capture 数据，结果只能解释为 off-target repeat signal，不代表无偏的全基因组 TE 含量。
