# 功能注释文件格式说明

流程支持两种方式获得功能注释：

1. 提供一个功能注释表（`annotation.table`）；
2. 没有注释表时，设置 `annotation.from_gtf: true`，从 GTF/GFF 属性中提取基因名/产物。

下面重点说明第 1 种方式。

## 1. 主注释表 `annotation.table`

### 基本要求

- 文本文件，默认制表符分隔；CSV 请设置 `annotation.delimiter: ","`；
- 必须有一行表头；
- 每一行是一个基因的记录；
- 基因 ID 必须和 featureCounts 使用的 ID 一致，即与 `reference.gene_attribute`
  指向的 GTF/GFF 属性一致（默认 `gene_id`）；
- 同一个基因可以出现多行，程序会自动合并；
- 列的顺序不限，列名大小写、空格、下划线不敏感；
- 额外列会被保留到 `DE_annotated/*_DEGs_annotated.csv`。

### 必需或推荐列

| 配置项 | 作用 | 示例列名 |
|---|---|---|
| `gene_id_column` | 基因 ID | GeneID、gene_id、gene id |
| `product_column` | 产物/功能描述 | Product、Description、Annotation |
| `go_column` | GO 注释文本 | GO Terms、GO、go_terms |

如果目标只是完成差异表达分析，`product_column` 和 `go_column` 都可以为空；
但功能基因集分类和 GO 富集需要它们。

### 最小示例

文件：`examples/annotation.example.tsv`

```
GeneID	Product	GO Terms
g1	catalase	GO:0004096 - catalase activity; GO:0006979 - response to oxidative stress
g2	superoxide dismutase	GO:0004784 - superoxide dismutase activity
g3	hypothetical protein	
```

### GO 列的写法

程序用正则表达式 `GO:\d+` 提取 GO ID，所以下面几种写法都可以：

```
GO:0004096
GO:0004096 - catalase activity
GO:0004096; GO:0006979
GO_function: GO:0004784 - superoxide dismutase activity [Evidence IEA]
```

提取出的 ID 会自动生成 gene-to-GO 映射，无需再单独准备 `go_gene_table`。

## 2. GO 映射表 `annotation.go_gene_table`

如果主注释表里没有 GO 列，可以单独提供。支持有表头或无表头，两列或三列：

```
GeneID	GO
g1	GO:0004096
g1	GO:0006979
g2	GO:0004784
```

或带 ontology 的格式：

```
GeneID	GO	Ontology
g1	GO:0004096	MF
g1	GO:0006979	BP
```

文件示例：`examples/gene_go.example.tsv`

## 3. KEGG 基因到 KO 表 `annotation.kegg_gene_table`

建议使用有表头的两列文件，列名可以是 `GeneID` 和 `KO`：

```
GeneID	KO
g1	K00001
g2	K00121
```

程序也兼容没有表头的两列表，第一列视为基因，第二列视为 KO。KO 可以写成
`K00001` 或 `ko:K00001`。文件示例：`examples/gene_ko.example.tsv`

## 4. KO 到通路表 `annotation.ko_pathway_table`

两列，无表头，第一列 KO，第二列通路：

```
ko:K00001	path:ko00010
ko:K00001	path:ko00620
ko:K00121	path:ko00010
```

也可以使用不带前缀的写法：

```
K00001	00010
K00001	00620
```

文件示例：`examples/ko_pathway.example.tsv`

## 5. 通路名称表 `annotation.pathway_names`

两列，无表头：

```
ko00010	Glycolysis / Gluconeogenesis
ko00020	Citrate cycle (TCA cycle)
ko00900	Terpenoid backbone biosynthesis
```

文件示例：`examples/pathway_names.example.tsv`

## 6. 没有外部注释文件时

在 `config.yaml` 中写：

```yaml
annotation:
  table: ""
  from_gtf: true
```

程序会解析 GTF/GFF 属性中的 `gene_id`、`gene_name`、`product`、`Name`、
`description`，生成：

```
results/annotation/annotation_from_gtf.tsv
```

该文件可以直接作为后续分析的注释表使用，但没有 GO/KEGG 信息时仍无法做富集。

## 7. 常见错误

- 基因 ID 与 featureCounts 的 ID 不一致，导致注释全部为空；
- 表头列名和 config 不匹配，例如表里是 `gene_id`，配置写 `GeneID`，
  一般情况下可以自动匹配，但请不要写成完全无关的名字；
- 注释表是 CSV 但 `annotation.delimiter` 仍然是 `"\t"`；
- 同一个基因的多行信息用逗号分隔且未加引号，CSV 解析会错列；
- 直接把 Excel 文件保存成 `.xls/.xlsx` 而没有导出为 TSV/CSV。
