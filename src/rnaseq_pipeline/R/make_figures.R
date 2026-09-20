#!/usr/bin/env Rscript
suppressMessages({
  library(ggplot2)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
cfg_file <- NULL
for (i in seq_along(args)) {
  if (args[[i]] == "--config" && i < length(args)) cfg_file <- args[[i + 1]]
}
if (is.null(cfg_file)) stop("usage: make_figures.R --config figures_config.json")
cfg <- fromJSON(cfg_file, simplifyVector = FALSE)

outdir <- cfg$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
formats <- if (is.null(cfg$formats) || !length(cfg$formats)) c("png", "pdf") else cfg$formats

save_plot <- function(p, base, width, height) {
  if ("pdf" %in% formats) {
    tryCatch(
      ggsave(paste0(base, ".pdf"), p, width = width, height = height,
             device = grDevices::cairo_pdf),
      error = function(e) ggsave(paste0(base, ".pdf"), p,
                                 width = width, height = height, device = "pdf")
    )
  }
  if ("png" %in% formats) {
    tryCatch(
      ggsave(paste0(base, ".png"), p, width = width, height = height,
             dpi = 300, device = "png", type = "cairo"),
      error = function(e) ggsave(paste0(base, ".png"), p,
                                 width = width, height = height, dpi = 300)
    )
  }
}

read_counts <- function(path) {
  first <- readLines(path, n = 1, warn = FALSE)
  skip <- if (grepl("^\\s*#", first)) 1 else 0
  d <- read.delim(path, skip = skip, check.names = FALSE,
                  comment.char = "", stringsAsFactors = FALSE)
  ann <- intersect(c("Chr", "Start", "End", "Strand", "Length"), colnames(d))
  if (length(ann)) d <- d[, setdiff(colnames(d), ann), drop = FALSE]
  d
}

contrast_names <- vapply(cfg$contrasts, function(x) as.character(x$name), "")

# ---------------- 1. DEG count overview ----------------
deg_file <- file.path(cfg$de_dir, "DEG_summary.csv")
if (file.exists(deg_file)) {
  deg <- read.csv(deg_file, stringsAsFactors = FALSE)
  deg_long <- rbind(
    data.frame(comparison = deg$comparison, count = deg$up, direction = "Up"),
    data.frame(comparison = deg$comparison, count = deg$down, direction = "Down")
  )
  deg_long$direction <- factor(deg_long$direction, levels = c("Up", "Down"))
  p <- ggplot(deg_long, aes(x = comparison, y = count, fill = direction)) +
    geom_col(position = "stack", width = 0.62, color = "grey20", linewidth = 0.2) +
    geom_text(aes(label = ifelse(count > 0, count, "")),
              position = position_stack(vjust = 0.5), size = 3.4) +
    scale_fill_manual(values = c(Up = "#B2182B", Down = "#2166AC")) +
    labs(title = "Differentially expressed genes",
         subtitle = paste0("padj < ", cfg$padj, " and |log2FC| >= ", cfg$log2fc),
         x = NULL, y = "Number of DEGs", fill = "Regulation") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_plot(p, file.path(outdir, "DEG_counts"), 7, 4.8)
}

# ---------------- 2. Volcano plots ----------------
for (nm in contrast_names) {
  f <- file.path(cfg$de_dir, paste0(nm, "_all_results.csv"))
  if (!file.exists(f)) next
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("log2FoldChange", "pvalue", "padj") %in% colnames(d))) next
  d <- d[!is.na(d$pvalue) & !is.na(d$log2FoldChange), ]
  d$group <- "Not significant"
  sig <- !is.na(d$padj) & d$padj < cfg$padj & abs(d$log2FoldChange) >= cfg$log2fc
  d$group[sig & d$log2FoldChange > 0] <- "Up"
  d$group[sig & d$log2FoldChange < 0] <- "Down"
  d$group <- factor(d$group, levels = c("Up", "Down", "Not significant"))
  p <- ggplot(d, aes(log2FoldChange, -log10(pvalue), color = group)) +
    geom_point(size = 0.7, alpha = 0.65) +
    scale_color_manual(values = c(Up = "#B2182B", Down = "#2166AC",
                                  "Not significant" = "grey70")) +
    labs(title = paste0("Volcano plot: ", nm),
         x = "log2 fold change", y = "-log10 p-value", color = NULL) +
    theme_bw(base_size = 12)
  save_plot(p, file.path(outdir, paste0("Volcano_", nm)), 6.4, 5.4)
}

# ---------------- 3. Category DEG counts ----------------
cat_summary <- file.path(cfg$category_dir, "category_summary.csv")
if (file.exists(cat_summary)) {
  cs <- read.csv(cat_summary, stringsAsFactors = FALSE)
  long <- rbind(
    data.frame(category = cs$category, comparison = cs$comparison,
               count = cs$up, direction = "Up"),
    data.frame(category = cs$category, comparison = cs$comparison,
               count = cs$down, direction = "Down")
  )
  long$direction <- factor(long$direction, levels = c("Up", "Down"))
  p <- ggplot(long, aes(x = comparison, y = count, fill = direction)) +
    geom_col(position = "stack", width = 0.6, color = "grey25", linewidth = 0.2) +
    geom_text(aes(label = ifelse(count > 0, count, "")),
              position = position_stack(vjust = 0.5), size = 2.8) +
    facet_wrap(~category, scales = "free_y") +
    scale_fill_manual(values = c(Up = "#B2182B", Down = "#2166AC")) +
    labs(title = "DEGs in custom functional gene sets",
         x = NULL, y = "Number of DEGs", fill = "Regulation") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_plot(p, file.path(outdir, "Category_DEG_counts"), 8, 4.8)
}

# ---------------- 4. Category log2FC heatmaps ----------------
cat_files <- list.files(cfg$category_dir, pattern = "^combined_.*_matrix\\.csv$",
                        full.names = TRUE)
for (f in cat_files) {
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  log2_cols <- grep("_log2FC$", colnames(d), value = TRUE)
  deg_cols <- grep("_DEG$", colnames(d), value = TRUE)
  if (!length(log2_cols) || !length(deg_cols)) next
  keep <- rowSums(d[deg_cols] == "yes", na.rm = TRUE) > 0
  d <- d[keep, , drop = FALSE]
  if (!nrow(d)) next
  if (nrow(d) > 70) {
    score <- apply(abs(d[log2_cols]), 1, max, na.rm = TRUE)
    d <- d[order(score, decreasing = TRUE)[1:70], , drop = FALSE]
  }
  labels <- d$gene
  if ("product" %in% colnames(d)) {
    prod <- ifelse(is.na(d$product) | d$product == "" |
                     grepl("hypothetical", d$product, ignore.case = TRUE),
                   "", paste0(" | ", substr(d$product, 1, 40)))
    labels <- paste0(d$gene, prod)
  }
  long <- data.frame()
  for (i in seq_len(nrow(d))) {
    for (col in log2_cols) {
      value <- suppressWarnings(as.numeric(d[[col]][i]))
      if (is.na(value)) next
      deg_col <- sub("_log2FC$", "_DEG", col)
      star <- if (deg_col %in% deg_cols && !is.na(d[[deg_col]][i]) &&
                  d[[deg_col]][i] == "yes") "*" else ""
      long <- rbind(long, data.frame(
        gene = labels[i],
        comparison = sub("_log2FC$", "", col),
        value = value, star = star, stringsAsFactors = FALSE
      ))
    }
  }
  if (!nrow(long)) next
  lim <- max(1, ceiling(max(abs(long$value), na.rm = TRUE)))
  p <- ggplot(long, aes(x = comparison, y = gene, fill = value)) +
    geom_tile(color = "grey88", linewidth = 0.3) +
    geom_text(aes(label = star), size = 3, na.rm = TRUE) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0, limits = c(-lim, lim), name = "log2FC") +
    labs(title = sub("\\.csv$", "", basename(f)),
         subtitle = "* padj < threshold and |log2FC| >= threshold",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 7.5),
          axis.text.x = element_text(angle = 30, hjust = 1))
  h <- max(4, 1.6 + 0.28 * length(unique(long$gene)))
  base_name <- sub("^combined_", "", sub("\\.csv$", "", basename(f)))
  base_name <- sub("_matrix$", "", base_name)
  save_plot(p, file.path(outdir, paste0("Category_", base_name, "_heatmap")), 9, h)
}

# ---------------- 5. Enrichment bubble plots ----------------
old_enrichment <- list.files(
  outdir, pattern = "^(Enrichment_top_terms|Enrichment_bubble_).*\\.(png|pdf)$",
  full.names = TRUE
)
if (length(old_enrichment)) file.remove(old_enrichment)

parse_enrichment_name <- function(path) {
  nm <- sub("\\.csv$", "", basename(path))
  db <- NA_character_
  for (candidate in c("GO_BP", "GO_MF", "GO_CC", "KEGG")) {
    if (grepl(candidate, nm, fixed = TRUE)) {
      db <- candidate
      break
    }
  }
  if (is.na(db)) return(NULL)
  direction <- "All"
  if (grepl("(^|_)up(_|$)", nm)) direction <- "Up"
  if (grepl("(^|_)down(_|$)", nm)) direction <- "Down"
  comparison <- nm
  comparison <- gsub("_(GO_BP|GO_MF|GO_CC|KEGG)", "", comparison)
  comparison <- gsub("_(up|down)", "", comparison)
  comparison <- gsub("_(nominal|all)$", "", comparison)
  comparison <- gsub("_+", "_", comparison)
  comparison <- gsub("^_|_$", "", comparison)
  list(comparison = comparison, direction = direction, db = db)
}

parse_gene_ratio <- function(x) {
  vapply(as.character(x), function(s) {
    parts <- strsplit(s, "/", fixed = TRUE)[[1]]
    if (length(parts) != 2) return(NA_real_)
    num <- suppressWarnings(as.numeric(parts[1]))
    den <- suppressWarnings(as.numeric(parts[2]))
    if (is.na(num) || is.na(den) || den <= 0) return(NA_real_)
    num / den
  }, numeric(1))
}

wrap_terms <- function(x, width = 48) {
  vapply(as.character(x), function(s) {
    paste(strwrap(s, width = width), collapse = "\n")
  }, character(1))
}

en_files <- list.files(cfg$enrichment_dir, pattern = "\\.csv$", full.names = TRUE)
en_files <- en_files[!grepl("enrichment_summary", basename(en_files))]
enrichment_threshold <- if (is.null(cfg$enrichment_threshold)) {
  "strict"
} else {
  as.character(cfg$enrichment_threshold)
}
if (length(en_files)) {
  base_names <- basename(en_files)
  is_nominal <- grepl("_nominal\\.csv$", base_names)
  is_all <- grepl("_all\\.csv$", base_names)
  is_strict <- !is_nominal & !is_all
  if (enrichment_threshold == "strict") {
    en_files <- en_files[is_strict]
    if (!length(en_files)) {
      cat("No enrichment results with p.adjust < 0.05; bubble plots skipped.\n")
    }
  } else if (enrichment_threshold == "all") {
    rank <- ifelse(is_all, 1, ifelse(is_strict, 2, 3))
    ord <- order(rank, na.last = TRUE)
    en_files <- en_files[ord]
  } else {
    rank <- ifelse(is_nominal, 1, ifelse(is_strict, 2, 3))
    ord <- order(rank, na.last = TRUE)
    en_files <- en_files[ord]
  }
  parsed <- lapply(en_files, parse_enrichment_name)
  keys <- vapply(parsed, function(x) {
    if (is.null(x)) NA_character_
    else paste(x$comparison, x$direction, x$db, sep = "|")
  }, character(1))
  keep <- !is.na(keys) & !duplicated(keys)
  en_files <- en_files[keep]
}
if (length(en_files)) {
  enrichment_rows <- list()
  for (f in en_files) {
    info <- parse_enrichment_name(f)
    if (is.null(info)) next
    d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    if (!nrow(d) || !"pvalue" %in% colnames(d)) next
    d <- d[order(d$pvalue), , drop = FALSE]
    d <- head(d, 60)
    term <- if ("Description" %in% colnames(d)) d$Description else d$ID
    ratio <- if ("GeneRatio" %in% colnames(d)) parse_gene_ratio(d$GeneRatio) else NA_real_
    count <- if ("Count" %in% colnames(d)) d$Count else NA_real_
    enrichment_rows[[length(enrichment_rows) + 1]] <- data.frame(
      comparison = info$comparison,
      direction = info$direction,
      db = info$db,
      ontology = if (info$db %in% c("GO_BP", "GO_MF", "GO_CC")) {
        sub("^GO_", "", info$db)
      } else {
        ""
      },
      term = wrap_terms(term),
      ratio = ratio,
      pvalue = d$pvalue,
      padj = if ("p.adjust" %in% colnames(d)) d$p.adjust else NA_real_,
      count = count,
      stringsAsFactors = FALSE
    )
  }
  if (length(enrichment_rows)) {
    all_enrichment <- do.call(rbind, enrichment_rows)
    all_enrichment$group_db <- ifelse(
      all_enrichment$db %in% c("GO_BP", "GO_MF", "GO_CC"), "GO", "KEGG"
    )
    all_enrichment$group_key <- paste(all_enrichment$comparison,
                                      all_enrichment$group_db, sep = "|")
    has_direction <- tapply(all_enrichment$direction, all_enrichment$group_key,
                            function(x) any(x %in% c("Up", "Down")))
    drop_all <- all_enrichment$direction == "All" &
      has_direction[all_enrichment$group_key]
    all_enrichment <- all_enrichment[!drop_all, , drop = FALSE]
    group_keys <- unique(paste(all_enrichment$comparison,
                               all_enrichment$group_db,
                               all_enrichment$direction, sep = "\r"))
    for (key in group_keys) {
      parts <- strsplit(key, "\r", fixed = TRUE)[[1]]
      sub_data <- all_enrichment[
        all_enrichment$comparison == parts[[1]] &
          all_enrichment$group_db == parts[[2]] &
          all_enrichment$direction == parts[[3]], ,
        drop = FALSE
      ]
      if (parts[[2]] == "GO") {
        strata <- sub_data$ontology
        top_n <- 8
      } else {
        strata <- rep("KEGG", nrow(sub_data))
        top_n <- 12
      }
      sub_data <- do.call(rbind, lapply(split(sub_data, strata),
                                        function(d) head(d[order(d$pvalue), ], top_n)))
      sub_data <- sub_data[order(sub_data$pvalue), , drop = FALSE]
      sub_data <- head(sub_data, if (parts[[2]] == "GO") 36 else 12)
      if (all(!is.na(sub_data$ratio))) {
        sub_data$x <- sub_data$ratio
        x_label <- "Gene ratio"
      } else {
        sub_data$x <- -log10(sub_data$pvalue)
        x_label <- "-log10 P"
      }
      sub_data$size_value <- if (any(!is.na(sub_data$count))) {
        ifelse(is.na(sub_data$count), min(sub_data$count, na.rm = TRUE),
               sub_data$count)
      } else {
        1
      }
      sub_data$term <- factor(sub_data$term,
                              levels = rev(unique(sub_data$term)))
      if (any(!is.na(sub_data$padj) & sub_data$padj > 0)) {
        sub_data$color_value <- ifelse(
          !is.na(sub_data$padj) & sub_data$padj > 0,
          sub_data$padj,
          sub_data$pvalue
        )
        color_label <- "p.adjust"
      } else {
        sub_data$color_value <- sub_data$pvalue
        color_label <- "P value"
      }
      p <- ggplot(sub_data,
                  aes(x = x, y = term, size = size_value,
                      color = color_value)) +
        geom_point(alpha = 0.9) +
        scale_color_gradient(low = "#D73027", high = "#4575B4",
                             name = color_label) +
        scale_size_continuous(range = c(2.2, 7), name = "Count") +
        labs(x = x_label, y = NULL) +
        theme_bw(base_size = 11) +
        theme(
          axis.text.y = element_text(size = 7.5),
          panel.grid.minor = element_blank()
        )
      if (parts[[2]] == "GO") {
        p <- p + facet_grid(ontology ~ ., scales = "free_y", space = "free_y")
      }
      direction_tag <- tolower(parts[[3]])
      base_name <- paste0("Enrichment_bubble_", parts[[1]], "_",
                          direction_tag, "_", parts[[2]])
      save_plot(p, file.path(outdir, base_name),
                8.5, min(30, max(4.5, 0.34 * nrow(sub_data) + 2)))
    }
  }
}

# ---------------- 6. Sample correlation heatmap ----------------
if (!is.null(cfg$counts_file) && file.exists(cfg$counts_file)) {
  counts <- read_counts(cfg$counts_file)
  if (ncol(counts) >= 2) {
    mat <- as.matrix(counts[, -1, drop = FALSE])
    suppressWarnings(storage.mode(mat) <- "double")
    mat[is.na(mat)] <- 0
    cpm <- log2(t(t(mat) / pmax(colSums(mat), 1)) * 1e6 + 1)
    cor_mat <- cor(cpm)
    cor_long <- data.frame(
      sample1 = rep(rownames(cor_mat), each = ncol(cor_mat)),
      sample2 = rep(colnames(cor_mat), times = nrow(cor_mat)),
      value = as.vector(cor_mat)
    )
    p <- ggplot(cor_long, aes(sample1, sample2, fill = value)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = sprintf("%.2f", value)), size = 2.6) +
      scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                           midpoint = mean(cor_mat), name = "Pearson r") +
      labs(title = "Sample correlation (log2 CPM)", x = NULL, y = NULL) +
      theme_minimal(base_size = 11) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1))
    save_plot(p, file.path(outdir, "Sample_correlation"),
              7, 6)
  }
}

cat("figures written to", outdir, "\n")
