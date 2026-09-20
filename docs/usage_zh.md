# 使用说明

## 1. 环境安装

```bash
mamba env create -f environment.yml
conda activate rnaseq-pipeline
pip install -e .
```

如果只想在已有环境中安装 Python 部分：

```bash
pip install -e .
```

在 Windows 上建议使用 WSL2（Ubuntu）运行完整流程。WSL 中 D/E 盘路径分别为 `/mnt/d`、`/mnt/e`，所有命令都在 WSL 终端中执行。

## 2. 建立项目

```bash
rnaseq-pipeline init -o morchella_project
cd morchella_project
```

编辑 `samples.tsv`：

```

如果样本表是逗号分隔的 CSV，将 `config.yaml` 中的 `samples_delimiter` 改为 `","` 即可。
sample	group	replicate	fastq_1	fastq_2
control_1	control	1	/data/control_1_R1.fq.gz	/data/control_1_R2.fq.gz
treat4_1	treat4	1	/data/treat4_1_R1.fq.gz	/data/treat4_1_R2.fq.gz
```

编辑 `config.yaml` 中的参考基因组、GTF、注释表和比较组。

## 3. 配置检查

```bash
rnaseq-pipeline check -c config.yaml
rnaseq-pipeline check -c config.yaml --strict
```

`--strict` 会同时检查 FASTQ、参考基因组、GTF、注释表是否存在。

## 4. 运行

```bash
# 完整流程
rnaseq-pipeline run -c config.yaml

# 只运行指定步骤
rnaseq-pipeline run -c config.yaml --steps fastp,index,align,count
rnaseq-pipeline run -c config.yaml --steps deseq,annotate,categories

# 只根据已有结果出图
rnaseq-pipeline figures -c config.yaml
rnaseq-pipeline figures -c config.yaml --formats png,pdf

# 只查看将要执行的命令
rnaseq-pipeline run -c config.yaml --dry-run
```

如果 FASTQ 已经完成质控，在 `config.yaml` 中设置：

```yaml
project:
  input_is_clean: true
```

流程会自动跳过 fastp，直接从 HISAT2 比对开始。已有 BAM 时可以从
`--steps count,deseq,...` 开始；已有 featureCounts 矩阵时可以从
`--steps deseq,annotate,categories,enrichment,figures` 开始。

可用步骤：

`fastp, index, align, count, stringtie, deseq, annotate, categories, enrichment, figures, multiqc`

## 5. 输出说明

| 目录/文件 | 内容 |
|---|---|
| `results/qc/` | fastp HTML/JSON 报告 |
| `results/clean/` | 质控后的 FASTQ（`keep_intermediate: false` 时自动删除） |
| `results/bam/` | 排序并建索引的 BAM |
| `results/counts/featureCounts_counts.txt` | 基因表达矩阵 |
| `results/DE/` | DESeq2 全部结果、DEG 列表、PCA、DEG 汇总 |
| `results/DE_annotated/` | 注释后的 DEG 表 |
| `results/category_analysis/` | 自定义基因集分类表和计数汇总 |
| `results/enrichment/` | GO/KEGG 富集结果 |
| `results/figures/` | 自动生成的 DEG、火山、功能基因集和富集气泡图（PNG/PDF） |

富集气泡图默认只画 `p.adjust < 0.05` 的条目；可通过
`analysis.figures.enrichment_threshold` 改为 `nominal`（原始 p < 0.05）或
`all`（全部条目）。
| `results/pipeline_summary.md` | 样本和结果摘要 |

## 6. 换成自己的物种

只需要修改 `config.yaml`：

1. `reference.genome_fasta` 和 `reference.gtf` 换成自己的参考文件；
2. 如果没有 HISAT2 索引，保持 `hisat2_index: ""` 并设置 `build_index: true`；
3. `annotation.table` 换成自己的功能注释表；
4. `analysis.contrasts` 换成自己的分组比较；
5. `analysis.categories` 换成自己的 GO/KEGG 条目和关键词。

例如把抗氧化基因集换成植物 APX/CAT/SOD：

```yaml
antioxidant:
  label: Antioxidant
  go: ["GO:0004601", "GO:0004784", "GO:0042744"]
  keywords: ["ascorbate peroxidase", "catalase", "superoxide dismutase"]
```

## 7. 注意事项

- `featureCounts` 的链特异性通过 `analysis.strandness` 配置：`0` 非链特异性，`1` forward，`2` reverse。
- 其他物种可以直接修改 `reference.annotation_format`（GTF/GFF）、`reference.feature_type`（exon/gene）、`reference.gene_attribute`（gene_id/ID）和 `analysis.strandness`（0/1/2）。
- 注释表如果不是制表符分隔，把 `annotation.delimiter` 改成 `","`；列名大小写不敏感。
- 如果 `annotation_format: GFF` 且需要运行 StringTie，请额外提供一个 GTF 注释文件（StringTie 只接受 GTF）。
- 如果没有功能注释文件，可以设置 `annotation.table: ""`、`annotation.from_gtf: true`，程序会从 GTF/GFF 属性中提取 gene_name/product；此时 GO/KEGG 富集需要额外映射文件，否则会自动跳过。
- 差异表达阈值默认 `padj < 0.05` 且 `|log2FC| >= 1`，可在 `analysis.deg` 中修改。
- 若样本缺少重复，仍可运行 DESeq2，但统计功效会下降，应在论文中说明。
- GO/KEGG 富集需要提供 `go_gene_table` 和/或 `kegg_gene_table`；没有映射文件时会自动跳过富集步骤。
