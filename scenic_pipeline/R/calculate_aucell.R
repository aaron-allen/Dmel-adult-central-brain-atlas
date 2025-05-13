# Calculate AUCell scores for high confidence regulons

library(optparse)
library(tidyverse)
library(Seurat)
library(AUCell)

option_list <- list(
    make_option(c("--sample_name"), default = "merged-all",
                help = "name of sample"),
    make_option(c("--data_type"), default = "raw",
                help = "was pySCENIC run on raw or normalised counts?"),
    make_option(c("--seed"), default = 123,
                help = "random seed for AUCell_buildRankings function"),
    make_option(c("--filtered_expression"), default = "filtered_expression.dir/raw.dir/merged-all.dir/filtered-expression.csv",
                help = "filtered raw counts matrix generated in filter_csv.py"),
    make_option(c("--auc_threshold"), default = 0.05,
                help = "threshold to calculate the AUC")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (interactive()) {
    opt <- list()
    opt$sample_name <- "merged-all"
    opt$data_type <- "raw"
    opt$seed <- 123
    opt$filtered_expression <- "filtered_expression.dir/raw.dir/merged-all.dir/filtered-expression.csv"
    opt$auc_threshold <- 0.05
}

message("sample_name: ", opt$sample_name,
        "\ndata_type: ", opt$data_type,
        "\nfiltered_expression: ", opt$filtered_expression,
        "\nauc_threshold:", opt$auc_threshold)

base_dir <- paste0("results.dir/aggregated.dir/", opt$data_type, ".dir/",
                   opt$sample_name, ".dir/")

# Read in data ----

high_confidence_regulons <- read_csv(
    paste0(base_dir, "high_confidence_regulons.csv"))

counts <- read.csv(opt$filtered_expression, row.names = 1)

# AUCell ----

# Rank genes by their expression across all cells
set.seed(opt$seed)
cells_rankings <- AUCell_buildRankings(as.matrix(counts), plotStats = FALSE)
saveRDS(cells_rankings, paste0(base_dir, "aucell_rankings.rds"))

if (!dir.exists(paste0(base_dir, "plots.dir"))) {
    dir.create(paste0(base_dir, "plots.dir"))
}

pdf(paste0(base_dir, "plots.dir/AUCell_number_genes_detected_per_cell.pdf"))
plotGeneCount(as.matrix(counts))
dev.off()

# Calculate the AUCell scores for the regulons
gene_sets <- split(high_confidence_regulons,
                   high_confidence_regulons$regulon_name)
gene_sets <- lapply(gene_sets, function(x) x$gene_name)

cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings,
                            aucMaxRank = round(nrow(cells_rankings) * as.numeric(opt$auc_threshold)))

# Save AUCell scores
aucell_scores <- cells_AUC@assays@data$AUC
aucell_scores <- as.data.frame(aucell_scores) %>%
    rownames_to_column(var = "Regulon")
write_csv(aucell_scores, paste0(base_dir, "aucell.csv"))

# Explore AUCell thresholds
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist = FALSE,
                                             nCores = 1, assign = TRUE)
saveRDS(cells_assignment, paste0(base_dir, "aucell_cells_assignment.rds"))

# Extract thresholds for the regulons
thresholds <- c()
for (gene_set in names(gene_sets)) {
    threshold <- cells_assignment[[gene_set]]$aucThr$selected
    thresholds[gene_set] <- threshold
}

# Plotting the AUCell histogram of the gene sets
for (gene_set in names(thresholds)) {
    pdf(paste0(base_dir, "plots.dir/", gene_set, "_AUCell_thresholds.pdf"))
    AUCell_plotHist(cells_AUC[gene_set, ], aucThr = thresholds[gene_set])
    abline(v = thresholds[gene_set])
    dev.off()
 }
