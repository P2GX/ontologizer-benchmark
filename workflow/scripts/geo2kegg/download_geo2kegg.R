library(GSEABenchmarkeR)
library(EnrichmentBrowser)
library(org.Hs.eg.db)

alpha_threshold <- 0.05
log_fc_threshold <- 1.0
output_tmp <- file.path(getwd(), "results/geo2kegg")


cat("Loading geo2kegg compendium...\n")
geo2kegg <- loadEData("geo2kegg")
names(geo2kegg)

geo2kegg <- maPreproc(geo2kegg)
geo2kegg <- runDE(geo2kegg, de.method = "limma", padj.method = "flexible")

# Define the dataset name explicitly (since loop is commented out)
for (ds_name in names(geo2kegg)) {
  output_dir <- file.path(output_tmp, ds_name)
  dir.create(output_dir, showWarnings = FALSE)
  cat(sprintf("Analyzing dataset: %s\n", ds_name))
  se <- geo2kegg[[ds_name]]

  all_entrez_ids <- rownames(se)

  all_symbols <- mapIds(
    org.Hs.eg.db,
    keys = all_entrez_ids,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )

  de_results <- as.data.frame(rowData(se))
  de_results$EntrezID <- all_entrez_ids
  de_results$Symbol <- mapIds(
    org.Hs.eg.db,
    keys = all_entrez_ids,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )

  study_mask <- de_results$ADJ.PVAL < alpha_threshold &
    abs(de_results$FC) > log_fc_threshold
  de_significant <- de_results[which(study_mask), ]

  study_file <- file.path(output_dir, sprintf("%s_study.txt", ds_name))
  pop_file <- file.path(output_dir, sprintf("%s_population.txt", ds_name))

  write(na.omit(de_significant$Symbol), file = study_file)
  write(na.omit(de_results$Symbol), file = pop_file)
}
