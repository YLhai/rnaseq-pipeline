#!/usr/bin/env Rscript
suppressMessages(library(DESeq2))

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (grepl("^--[^=]+=", a)) {
      key <- sub("^--([^=]+)=.*$", "\\1", a)
      val <- sub("^--[^=]+=", "", a)
      out[[key]] <- val
      i <- i + 1
    } else if (grepl("^--", a)) {
      key <- sub("^--", "", a)
      out[[key]] <- args[[i + 1]]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  out
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
need <- c("counts", "metadata", "contrasts", "outdir")
missing <- need[!need %in% names(opt)]
if (length(missing)) {
  stop("missing required arguments: ", paste(missing, collapse = ", "))
}

counts_file <- opt$counts
metadata_file <- opt$metadata
contrasts_file <- opt$contrasts
outdir <- opt$outdir
design_var <- if (!is.null(opt$`design-var`)) opt$`design-var` else "group"
padj_cut <- if (!is.null(opt$padj)) as.numeric(opt$padj) else 0.05
lfc_cut <- if (!is.null(opt$log2fc)) as.numeric(opt$log2fc) else 1.0
covariates <- character(0)
if (!is.null(opt$covariates) && nzchar(opt$covariates)) {
  covariates <- strsplit(opt$covariates, ",", fixed = TRUE)[[1]]
  covariates <- trimws(covariates[nzchar(trimws(covariates))])
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

first_line <- readLines(counts_file, n = 1, warn = FALSE)
skip_lines <- if (grepl("^\\s*#", first_line)) 1 else 0
cts <- read.delim(counts_file, row.names = 1, header = TRUE,
                  check.names = FALSE, comment.char = "",
                  skip = skip_lines)
annotation_cols <- intersect(c("Chr", "Start", "End", "Strand", "Length"),
                             colnames(cts))
if (length(annotation_cols)) {
  cts <- cts[, setdiff(colnames(cts), annotation_cols), drop = FALSE]
} else if (ncol(cts) > 5) {
  # Fall back to dropping the five standard featureCounts annotation columns.
  cts <- cts[, -(1:5), drop = FALSE]
}

meta <- read.delim(metadata_file, header = TRUE, check.names = FALSE,
                   stringsAsFactors = FALSE)
if (!"sample" %in% colnames(meta)) stop("metadata must contain a 'sample' column")
if (!design_var %in% colnames(meta)) {
  stop("design variable '", design_var, "' not found in metadata")
}

# featureCounts may use the full BAM path as the column name, for example
# "/data/bam/control_1.sorted.bam". Convert it back to the sample name.
sample_from_column <- function(x) {
  x <- basename(x)
  x <- sub("\\.bam$", "", x, ignore.case = TRUE)
  x <- sub("\\.sorted$", "", x, ignore.case = TRUE)
  x <- sub("\\.cram$", "", x, ignore.case = TRUE)
  x
}
count_samples <- sample_from_column(colnames(cts))
if (!anyDuplicated(count_samples) &&
    all(count_samples %in% meta$sample)) {
  colnames(cts) <- count_samples
}
if (!all(colnames(cts) %in% meta$sample)) {
  stop(
    "count columns do not match metadata$sample.\n",
    "count columns: ", paste(head(colnames(cts), 10), collapse = ", "), "\n",
    "metadata sample: ", paste(head(meta$sample, 10), collapse = ", ")
  )
}
meta <- meta[match(colnames(cts), meta$sample), , drop = FALSE]
if (anyNA(meta$sample)) {
  stop("metadata and count columns do not match; check sample names")
}

coldata <- data.frame(row.names = meta$sample)
for (col in setdiff(colnames(meta), "sample")) {
  values <- meta[[col]]
  numeric_values <- suppressWarnings(as.numeric(values))
  if (!anyNA(numeric_values) && length(unique(numeric_values)) > 2) {
    coldata[[col]] <- numeric_values
  } else {
    coldata[[col]] <- factor(values)
  }
}
coldata[[design_var]] <- factor(coldata[[design_var]])

design_terms <- c(design_var, covariates)
design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))

count_mat <- round(as.matrix(cts))
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = coldata,
                              design = design_formula)
dds <- DESeq(dds)

vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
pdf(file.path(outdir, "pca.pdf"), width = 7, height = 6)
print(plotPCA(vsd, intgroup = design_var))
dev.off()

contrasts <- read.delim(contrasts_file, header = TRUE,
                        stringsAsFactors = FALSE)
summary_table <- data.frame()

for (i in seq_len(nrow(contrasts))) {
  nm <- contrasts$name[i]
  num <- contrasts$numerator[i]
  den <- contrasts$denominator[i]
  res <- results(dds, contrast = c(design_var, num, den), alpha = padj_cut)
  res <- res[order(res$padj), ]
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  write.csv(df, file.path(outdir, paste0(nm, "_all_results.csv")),
            row.names = FALSE)

  deg <- subset(df, !is.na(padj) & padj < padj_cut & abs(log2FoldChange) >= lfc_cut)
  write.csv(deg[order(deg$padj), ],
            file.path(outdir, paste0(nm, "_DEGs.csv")), row.names = FALSE)
  summary_table <- rbind(summary_table, data.frame(
    comparison = nm,
    up = sum(deg$log2FoldChange > 0),
    down = sum(deg$log2FoldChange < 0),
    total = nrow(deg)
  ))
}

write.csv(summary_table, file.path(outdir, "DEG_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(outdir, "sessionInfo.txt"))
cat("DESeq2 summary:\n")
print(summary_table)
