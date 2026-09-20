#!/usr/bin/env Rscript
suppressMessages({
  library(clusterProfiler)
  library(GO.db)
  library(AnnotationDbi)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
cfg_file <- NULL
for (i in seq_along(args)) {
  if (args[[i]] == "--config" && i < length(args)) cfg_file <- args[[i + 1]]
}
if (is.null(cfg_file)) stop("usage: run_enrichment.R --config enrichment_config.json")
cfg <- fromJSON(cfg_file, simplifyVector = FALSE)

outdir <- cfg$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

background <- readLines(cfg$background, warn = FALSE)
background <- unique(background[nzchar(background)])

deg_sets <- list()
for (nm in names(cfg$deg_files)) {
  f <- cfg$deg_files[[nm]]
  if (!file.exists(f)) next
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"gene" %in% colnames(d)) next
  deg_sets[[nm]] <- d
}

go_map <- NULL
if (nzchar(cfg$go_map) && file.exists(cfg$go_map)) {
  raw <- read.delim(cfg$go_map, header = FALSE, stringsAsFactors = FALSE,
                    check.names = FALSE)
  if (nrow(raw) && any(grepl("gene|GO", raw[1, ], ignore.case = TRUE))) {
    raw <- raw[-1, , drop = FALSE]
  }
  if (ncol(raw) >= 2) {
    go_map <- data.frame(gene = raw[[1]], go = raw[[2]],
                         stringsAsFactors = FALSE)
    if (ncol(raw) >= 3) {
      go_map$ontology <- raw[[3]]
    } else {
      info <- suppressMessages(AnnotationDbi::select(
        GO.db, keys = unique(go_map$go),
        columns = c("TERM", "ONTOLOGY"), keytype = "GOID"))
      go_map <- merge(go_map, data.frame(go = info$GOID,
                                         term = info$TERM,
                                         ontology = info$ONTOLOGY,
                                         stringsAsFactors = FALSE),
                      by = "go", all.x = TRUE)
    }
  }
}

kegg_map <- NULL
if (nzchar(cfg$kegg_gene) && file.exists(cfg$kegg_gene) &&
    nzchar(cfg$ko_pathway) && file.exists(cfg$ko_pathway)) {
  first_line <- readLines(cfg$kegg_gene, n = 1, warn = FALSE)
  has_header <- grepl("gene|ko", first_line, ignore.case = TRUE)
  gene_ko <- read.delim(cfg$kegg_gene, header = has_header,
                        stringsAsFactors = FALSE, check.names = FALSE)
  if (!has_header && ncol(gene_ko) >= 2) {
    colnames(gene_ko)[1:2] <- c("GeneID", "KO")
  }
  gene_col <- if ("GeneID" %in% colnames(gene_ko)) "GeneID" else colnames(gene_ko)[1]
  ko_col <- if ("KO" %in% colnames(gene_ko)) "KO" else colnames(gene_ko)[2]
  ko_path <- read.delim(cfg$ko_pathway, header = FALSE, stringsAsFactors = FALSE)
  ko_path[[1]] <- sub("^ko:", "", ko_path[[1]])
  ko_path[[2]] <- sub("^path:ko", "", ko_path[[2]])
  gene_ko <- data.frame(gene = gene_ko[[gene_col]],
                        ko = sub("^ko:", "", gene_ko[[ko_col]]),
                        stringsAsFactors = FALSE)
  merged <- merge(gene_ko, ko_path, by.x = "ko", by.y = 1)
  kegg_map <- data.frame(term = merged[[3]], gene = merged$gene,
                         stringsAsFactors = FALSE)
  if (nzchar(cfg$pathway_names) && file.exists(cfg$pathway_names)) {
    pn <- read.delim(cfg$pathway_names, header = FALSE, stringsAsFactors = FALSE)
    pn[[1]] <- sub("^ko", "", pn[[1]])
    names_map <- setNames(pn[[2]], pn[[1]])
  } else {
    names_map <- setNames(unique(kegg_map$term), unique(kegg_map$term))
  }
}

summary_rows <- data.frame()

run_one <- function(genes, tag) {
  if (length(genes) < 3) return(invisible(NULL))
  if (!is.null(go_map)) {
    for (ont in c("BP", "MF", "CC")) {
      sub_map <- go_map[!is.na(go_map$ontology) & go_map$ontology == ont,
                        c("go", "gene")]
      sub_name <- unique(go_map[go_map$ontology == ont, c("go", "term")])
      if (!nrow(sub_map)) next
      colnames(sub_map) <- c("term", "gene")
      colnames(sub_name) <- c("term", "name")
      bg <- intersect(background, unique(sub_map$gene))
      hit <- intersect(genes, unique(sub_map$gene))
      if (length(hit) < 3) next
      ego <- enricher(hit, universe = bg, TERM2GENE = sub_map,
                      TERM2NAME = sub_name,
                      pAdjustMethod = cfg$padj_method,
                      pvalueCutoff = cfg$pvalue_cutoff,
                      qvalueCutoff = cfg$qvalue_cutoff,
                      minGSSize = cfg$min_gene_set_size,
                      maxGSSize = cfg$max_gene_set_size)
      if (!is.null(ego) && nrow(as.data.frame(ego))) {
        res <- as.data.frame(ego)
        res <- res[order(res$pvalue), ]
        write.csv(res, file.path(outdir, paste0(tag, "_GO_", ont, "_all.csv")),
                  row.names = FALSE)
        nom <- res[!is.na(res$pvalue) & res$pvalue < 0.05, ]
        strict <- nom[!is.na(nom$p.adjust) & nom$p.adjust < 0.05, ]
        if (nrow(nom)) {
          write.csv(nom, file.path(outdir, paste0(tag, "_GO_", ont, "_nominal.csv")),
                    row.names = FALSE)
          summary_rows <<- rbind(summary_rows, data.frame(
            set = tag, database = paste0("GO_", ont, "_nominal"), n = nrow(nom)))
        }
        if (nrow(strict)) {
          write.csv(strict, file.path(outdir, paste0(tag, "_GO_", ont, ".csv")),
                    row.names = FALSE)
          summary_rows <<- rbind(summary_rows, data.frame(
            set = tag, database = paste0("GO_", ont), n = nrow(strict)))
        }
      }
    }
  }
  if (!is.null(kegg_map)) {
    sub_map <- unique(kegg_map)
    colnames(sub_map) <- c("term", "gene")
    bg <- intersect(background, unique(sub_map$gene))
    hit <- intersect(genes, unique(sub_map$gene))
    if (length(hit) >= 3) {
      eko <- enricher(hit, universe = bg, TERM2GENE = sub_map,
                      TERM2NAME = data.frame(term = names(names_map),
                                             name = as.character(names_map)),
                      pAdjustMethod = cfg$padj_method,
                      pvalueCutoff = cfg$pvalue_cutoff,
                      qvalueCutoff = cfg$qvalue_cutoff,
                      minGSSize = cfg$min_gene_set_size,
                      maxGSSize = cfg$max_gene_set_size)
      if (!is.null(eko) && nrow(as.data.frame(eko))) {
        res <- as.data.frame(eko)
        res <- res[order(res$pvalue), ]
        write.csv(res, file.path(outdir, paste0(tag, "_KEGG_all.csv")),
                  row.names = FALSE)
        nom <- res[!is.na(res$pvalue) & res$pvalue < 0.05, ]
        strict <- nom[!is.na(nom$p.adjust) & nom$p.adjust < 0.05, ]
        if (nrow(nom)) {
          write.csv(nom, file.path(outdir, paste0(tag, "_KEGG_nominal.csv")),
                    row.names = FALSE)
          summary_rows <<- rbind(summary_rows, data.frame(
            set = tag, database = "KEGG_nominal", n = nrow(nom)))
        }
        if (nrow(strict)) {
          write.csv(strict, file.path(outdir, paste0(tag, "_KEGG.csv")),
                    row.names = FALSE)
          summary_rows <<- rbind(summary_rows, data.frame(
            set = tag, database = "KEGG", n = nrow(strict)))
        }
      }
    }
  }
}

for (nm in names(deg_sets)) {
  d <- deg_sets[[nm]]
  genes <- unique(d$gene)
  run_one(genes, paste0(nm, "_all"))
  if ("log2FoldChange" %in% colnames(d)) {
    up <- unique(d$gene[!is.na(d$log2FoldChange) & d$log2FoldChange > 0])
    down <- unique(d$gene[!is.na(d$log2FoldChange) & d$log2FoldChange < 0])
    run_one(up, paste0(nm, "_up"))
    run_one(down, paste0(nm, "_down"))
  }
}

write.csv(summary_rows, file.path(outdir, "enrichment_summary.csv"),
          row.names = FALSE)
cat("enrichment done\n")
