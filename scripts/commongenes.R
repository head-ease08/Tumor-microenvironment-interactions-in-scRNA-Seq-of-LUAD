library(recount3)
projects <- available_projects()

tcga <- subset(
  projects,
  project_home == "data_sources/tcga"
)

luad_info <- subset(tcga, project == "LUAD")

rse <- create_rse(luad_info)

stype <- colData(rse)$tcga.gdc_cases.samples.sample_type

rse_tumor <- rse[, stype == "Primary Tumor"]

expr <- assays(rse_tumor)$raw_counts
expr <- expr[rowSums(expr) > 0, ]
expr <- log2(expr + 1)

rd <- as.data.frame(rowData(rse_tumor))

expr_ens <- sub("\\..*", "", rownames(expr))

rd_ens <- rd$ensembl_clean

common <- intersect(expr_ens, rd_ens)

map <- rd[rd$ensembl_clean %in% common, c("ensembl_clean", "gene_name")]

expr_filtered <- expr[expr_ens %in% common, ]
rownames(expr_filtered) <- map$gene_name[match(expr_ens[expr_ens %in% common], map$ensembl_clean)]

genes = readLines("../datafor_gepia2.txt")
genes <- trimws(genes)
common_genes <- intersect(rownames(expr_filtered), genes)
writeLines(common_genes, "../data/genes_filtered_with_recount.txt")
