# Generate a file containing high confidence regulons (regulons that occur in > n% of runs, and genes that are found in those regulons on > 80% of runs)

# Author: Lucy Garner
#  modifications : Aaron M. Allen

library(optparse)
library(tidyverse)

option_list <- list(
    make_option(c("--file_path"),
                help = "path to results folder where different pySCENIC run outputs are"),
    make_option(c("--sample_name"), default = "merged-all",
                help = "name of sample"),
    make_option(c("--data_type"), default = "raw",
                help = "was pySCENIC run on raw or normalised counts?"),
    make_option(c("--number_runs"), default = 100,
                help = "how many times pySCENIC workflow has been run"),
    make_option(c("--number_genes"), default = 1,
                help = "number of genes required for high confidence regulons"),
    make_option(c("--percent_occurance"), default = 50,
                help = "percent occurance required for high confidence regulons")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (interactive()) {
    opt <- list()
    opt$file_path <- "results/server/pyscenic_multi"
    opt$sample_name <- "merged-all"
    opt$data_type <- "raw"
    opt$number_runs <- 100
    opt$number_genes <- 1
    opt$percent_occurance <- 50
}

message("file_path: ", opt$file_path,
        "\nsample_name: ", opt$sample_name,
        "\ndata_type: ", opt$data_type,
        "\nnumber_runs: ", opt$number_runs,
        "\nnumber_genes: ", opt$number_genes,
        "\npercent_occurance: ", opt$percent_occurance)

# Define functions ---

# Convert regulons.csv table into long format
convert_to_long <- function(regulon_table) {
    # Max number of genes in regulon
    regulon_sizes <- regulon_table$genes %>%
        str_count(pattern = ", ")
    regulon_sizes <- regulon_sizes + 1
    names(regulon_sizes) <- regulon_table$regulon_name
    max_n_genes_in_regulon <- max(regulon_sizes)

    regulons_long_gene <- regulon_table %>%
        separate(genes, into = paste0("gene", seq(1, max_n_genes_in_regulon)),
                 sep = ", ") %>%
        pivot_longer(seq(3, ncol(.) - 3), names_to = "gene", values_to = "gene_name",
                     values_drop_na = TRUE) %>%
        select(regulon_name, gene_name)

    regulon_table$weights <- gsub("\\]", "", regulon_table$weights)
    regulons_long_weight <- regulon_table %>%
        select(-genes, -score, -context) %>%
        separate(weights, into = paste0("weight", seq(1, max_n_genes_in_regulon)),
                 sep = ", ") %>%
        pivot_longer(seq(3, ncol(.)), names_to = "weightnumber", values_to = "weight",
                     values_drop_na = TRUE) %>%
        select(regulon_name, weight)

    regulons_long <- cbind(regulons_long_gene, regulons_long_weight$weight)
    colnames(regulons_long)[3] <- "weight"

    return(regulons_long)
}

# Read data ----

regulon_list <- list()
for (i in 1:opt$number_runs) {
    file_path <- paste0(opt$file_path, "/", "run_", i, ".dir/", opt$data_type,
                        ".dir/", opt$sample_name, ".dir/")
    regulons <- read_csv(paste0(file_path, "regulons.csv"))
    regulon_list[[i]] <- regulons
}
saveRDS(regulon_list,
        paste0(opt$file_path, "/aggregated.dir/", opt$data_type,
               ".dir/", opt$sample_name, ".dir/all_runs_regulon_list.rds"))

# Convert regulons tables into long format
regulons_long <- lapply(regulon_list, convert_to_long)
regulons_long <- lapply(1:opt$number_run, function(x) {
    regulons_long[[x]] <- regulons_long[[x]] %>%
        mutate(run_number = x)
})
regulons_long_df <- do.call("rbind", regulons_long)
write_csv(regulons_long_df,
          paste0(opt$file_path, "/aggregated.dir/", opt$data_type,
                 ".dir/", opt$sample_name, ".dir/all_runs_regulons_long_df.csv"))

# High confidence regulons ---

# Count how many times each regulon is identified
regulons <- unique(unlist(sapply(regulon_list, function(x) x$regulon_name)))    # vector of all regulons

regulon_occurrence <- c()
for (i in 1:length(regulons)) {
    regulon_occurrence[regulons[i]] <- sum(sapply(
        1:opt$number_run,
        function(x) regulons[i] %in% regulon_list[[x]]$regulon_name))
}

regulon_occurrence <- sort(regulon_occurrence, decreasing = TRUE)
regulon_occurrence <- regulon_occurrence/opt$number_runs * 100
saveRDS(regulon_occurrence,
        paste0(opt$file_path, "/aggregated.dir/", opt$data_type,
               ".dir/", opt$sample_name, ".dir/regulon_occurrence.rds"))

if (!dir.exists(paste0(opt$file_path, "/aggregated.dir/", opt$data_type,
                       ".dir/", opt$sample_name, ".dir/plots.dir"))) {
    dir.create(paste0(opt$file_path, "/aggregated.dir/", opt$data_type,
                      ".dir/", opt$sample_name, ".dir/plots.dir"))
}
pdf(paste0(opt$file_path, "/aggregated.dir/", opt$data_type,
           ".dir/", opt$sample_name, ".dir/plots.dir/regulon_occurrence.pdf"))
hist(regulon_occurrence, breaks = 100)
dev.off()

# Define high confidence regulons as those identified on > n% of runs
# Will filter further later to retain only genes identified within a regulon on > n% of runs and only regulons with >= opt$number_genes genes
high_conf_regulons <- names(regulon_occurrence)[regulon_occurrence >= opt$percent_occurance]

# Filter to keep only high confidence regulons
# Determine the percentage of runs in which each gene is found within a given regulon
regulon_gene_occurrence <- regulons_long_df %>%
    filter(regulon_name %in% high_conf_regulons) %>%
    group_by(regulon_name, gene_name) %>%
    tally() %>%
    mutate(percent_occurrence = n/opt$number_runs * 100)

# Keep genes within a regulon if they are found in that regulon in > 20% of runs
high_conf_regulons_genes <- regulon_gene_occurrence %>%
    # filter(percent_occurrence > 80)
    filter(percent_occurrence > 20)

# Count the number of genes in each high confidence regulon
number_genes_regulon <- high_conf_regulons_genes %>%
    group_by(regulon_name) %>%
    tally()

# Filter high confidence regulons to retain only those that contain >= opt$number_genes genes
high_conf_regulons <- number_genes_regulon$regulon_name[number_genes_regulon$n >= opt$number_genes]
message("Number of high confidence regulons is ", length(high_conf_regulons))

high_conf_regulons_genes <- high_conf_regulons_genes %>%
    filter(regulon_name %in% high_conf_regulons)
write_csv(high_conf_regulons_genes %>%
              select(regulon_name, gene_name),
          paste0(opt$file_path, "/aggregated.dir/", opt$data_type,
                 ".dir/", opt$sample_name, ".dir/high_confidence_regulons.csv"))

regulon_gene_occurrence <- regulon_gene_occurrence %>%
    filter(regulon_name %in% high_conf_regulons)
write_csv(regulon_gene_occurrence,
          paste0(opt$file_path, "/aggregated.dir/", opt$data_type,
                 ".dir/", opt$sample_name,
                 ".dir/regulon_gene_occurrence_high_confidence_regulons.csv"))
