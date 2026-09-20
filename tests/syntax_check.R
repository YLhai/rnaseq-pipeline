files <- c(
  "src/rnaseq_pipeline/R/run_deseq2.R",
  "src/rnaseq_pipeline/R/run_enrichment.R",
  "src/rnaseq_pipeline/R/make_figures.R"
)
for (f in files) {
  parse(file = f)
  cat("SYNTAX OK:", f, "\n")
}
