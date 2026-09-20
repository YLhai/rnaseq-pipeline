# rnaseq-pipeline

一个面向真菌/植物转录组的可复用命令行流程，覆盖：

```
raw FASTQ
  -> fastp 质控与接头去除
  -> HISAT2 比对
  -> samtools 排序
  -> featureCounts 基因定量
  -> DESeq2 差异表达分析
  -> 差异基因注释
  -> 抗氧化/能量代谢/萜类代谢等自定义基因集分类
  -> GO / KEGG 富集（clusterProfiler）
  -> MultiQC 汇总（可选）
```

## 1. 安装

推荐使用 conda/mamba 创建环境：

```bash
git clone <your-repository-url>
cd rnaseq-pipeline
mamba env create -f environment.yml
conda activate rnaseq-pipeline
pip install -e .
```

也可以基于已有环境安装 Python 包：

```bash
conda install -c conda-forge -c bioconda fastp hisat2 samtools subread stringtie multiqc r-base bioconductor-deseq2 bioconductor-clusterprofiler bioconductor-go.db bioconductor-annotationdbi bioconductor-keggrest
pip install -e .
```

### Windows 用户

HISAT2、featureCounts 等工具在 Windows 原生环境支持有限，建议在 WSL2（Ubuntu）中安装和使用：

```bash
# 在 WSL 中
cd /mnt/e/codex/rnaseq-pipeline
mamba env create -f environment.yml
conda activate rnaseq-pipeline
pip install -e .
```

WSL 中访问 D/E 盘分别为 `/mnt/d` 和 `/mnt/e`，`config.yaml` 与 `samples.tsv` 中使用 Linux 路径即可。

## 2. 初始化项目

```bash
rnaseq-pipeline init -o my_morchella_project
cd my_morchella_project
```

会生成：

```
config.yaml          # 主配置文件
samples.tsv          # 样本表
```

编辑 `samples.tsv`，每行一个样本：

```
sample	group	replicate	fastq_1	fastq_2
control_1	control	1	/path/control_1_R1.fq.gz	/path/control_1_R2.fq.gz
treat4_1	treat4	1	/path/treat4_1_R1.fq.gz	/path/treat4_1_R2.fq.gz
```

## 3. 检查配置

```bash
rnaseq-pipeline check -c config.yaml
```

会检查样本表、参考基因组、GTF、注释文件、R 包和外部软件是否就绪。

## 4. 运行

```bash
# 全流程
rnaseq-pipeline run -c config.yaml

# 分析完成后单独重新出图
rnaseq-pipeline figures -c config.yaml

# 只跑到某一阶段
rnaseq-pipeline run -c config.yaml --steps fastp,index,align,count

# 先查看将要执行什么
rnaseq-pipeline run -c config.yaml --dry-run
```

### 如果原始数据已经做过质控

如果 `samples.tsv` 里的 FASTQ 已经是去接头、过滤低质量后的 clean reads，
在 `config.yaml` 中设置：

```yaml
project:
  input_is_clean: true
```

然后正常运行即可，fastp 会自动跳过，流程直接从 HISAT2 比对开始：

```bash
rnaseq-pipeline run -c config.yaml
```

也可以显式指定：

```bash
rnaseq-pipeline run -c config.yaml --steps index,align,count,deseq,annotate,categories,enrichment,figures
```

如果 clean reads 已经完成比对、你只有 BAM 或 featureCounts 表达矩阵，
可以分别从 `align` 或 `count`/`deseq` 步骤开始：

```bash
rnaseq-pipeline run -c config.yaml --steps count,deseq,annotate,categories,enrichment,figures
```

## 5. 主要输出

```
results/
├── clean/                 # fastp 输出（可通过 keep_intermediate 控制是否保留）
├── qc/                    # fastp 报告
├── bam/                   # 排序 BAM
├── counts/                # featureCounts 定量结果
├── DE/                    # DESeq2 结果
├── DE_annotated/          # 差异基因注释
├── category_analysis/     # 自定义功能基因集分类表
├── enrichment/            # GO/KEGG 富集
├── figures/               # 自动生成的论文图（PNG/PDF）
└── pipeline_summary.md    # 流程摘要
```

## 6. 定制功能基因集

`config.yaml` 中的 `analysis.categories` 用于自定义基因集，例如抗氧化、能量代谢和萜类代谢：

```yaml
analysis:
  categories:
    antioxidant:
      label: Antioxidant
      go: ["GO:0004601", "GO:0004784", "GO:0004364"]
      keywords: ["catalase", "superoxide dismutase", "glutathione"]
    terpenoid:
      label: Terpenoid
      kegg: ["ko00900", "ko00909"]
      keywords: ["farnesyl", "squalene", "mevalonate"]
```

只要换成自己的 GO/KEGG 条目和关键词，就能复用于其他物种或其他研究主题。

### 自动出图

分析完成后运行：

```bash
rnaseq-pipeline figures -c config.yaml
# 或
rnaseq-pipeline run -c config.yaml --steps figures
```

输出目录为 `results/figures/`，默认同时生成 PNG 和 PDF，当前包含：

- `DEG_counts`：各比较组上下调 DEG 数量
- `Volcano_*`：每个比较组的火山图
- `Category_DEG_counts`：自定义功能基因集的 DEG 数量
- `Category_*_heatmap`：功能基因集在各比较组的 log2FC 热图
- `Enrichment_bubble_<比较组>_up/down_GO`：上调和下调分别出图，GO BP/MF/CC 用分面分类；`..._up/down_KEGG`：KEGG 气泡图。x=Gene ratio，点大小=Count，颜色=p.adjust 红蓝渐变；比较组和方向只在文件名中体现，图内不加标题
- `Sample_correlation`：样本间 log2 CPM 相关性热图

可用 `--formats png` 或 `--formats pdf` 只输出单一格式。

富集气泡图默认使用严格显著结果（`p.adjust < 0.05`）。没有通过校正阈值的条目时，
对应的气泡图不会生成。如需使用名义显著或全部条目：

```yaml
analysis:
  figures:
    enrichment_threshold: strict   # strict | nominal | all
```

## 7. 适配其他物种

本流程不包含任何物种专属的硬编码逻辑，换成其他物种只需要修改 `config.yaml`：

1. 替换 `reference.genome_fasta` 和 `reference.gtf`；
2. 如果注释是 GFF3，设置 `reference.annotation_format: GFF`；
3. 根据注释文件设置 `reference.feature_type` 和 `reference.gene_attribute`
   （GTF 常用 `exon`/`gene_id`，GFF3 常用 `gene`/`ID`）；
4. 链特异性文库把 `analysis.strandness` 设为 `1`（forward）或 `2`（reverse）；
5. 替换 `annotation.table`，并用 `annotation.delimiter` 指定制表符或逗号；
6. 替换 `analysis.contrasts` 的分组名称；
7. 按研究主题替换 `analysis.categories`；
8. 有 GO 映射时可直接提供 `annotation.go_gene_table`；没有时流程会从注释表的 GO 列自动生成；
9. 有 KEGG 映射时提供 `annotation.kegg_gene_table`、`ko_pathway_table` 和 `pathway_names`，没有则自动跳过 KEGG。

如果暂时没有功能注释文件，可以只跑核心分析：

```yaml
annotation:
  table: ""
  from_gtf: true     # 从 GTF/GFF 的 gene_name/product/description 自动提取
analysis:
  enrichment: false  # 没有 GO/KEGG 映射时关闭富集
```

没有注释时仍会输出表达矩阵、DESeq2 差异结果、DEG 数量图和火山图；
功能基因集分类只能依据 GTF 中的 gene_name/product 关键词，GO/KEGG 富集无法进行。

功能注释表、GO 映射表和 KEGG 映射表的完整格式见
[docs/annotation_formats.md](docs/annotation_formats.md)，`examples/` 目录下有可直接照抄的示例文件。

StringTie 步骤默认关闭，默认使用 featureCounts 的表达矩阵进行 DESeq2 分析；
如需启动，将 `analysis.stringtie` 设为 `true`。

## 8. 依赖的软件

| 软件 | 用途 |
|---|---|
| fastp | 质控、接头去除 |
| HISAT2 | 剪接感知比对 |
| samtools | BAM 排序与索引 |
| featureCounts (subread) | 基因水平定量 |
| StringTie | 可选转录本组装 |
| DESeq2 | 差异表达分析 |
| clusterProfiler / GO.db / AnnotationDbi / KEGGREST | GO/KEGG 富集 |
| MultiQC | 质控汇总（可选） |

## 9. 引用

如果使用本流程，请引用相应工具：

- fastp, HISAT2, samtools, featureCounts, StringTie, DESeq2
- clusterProfiler (Xu et al., 2024; Wu et al., 2021)
- GO.db, AnnotationDbi, KEGGREST

本项目的 `docs/references/` 目录提供了 EndNote 可导入的 RIS 文件。

## 10. 许可

MIT License。可自由修改、分发和用于科研项目。
