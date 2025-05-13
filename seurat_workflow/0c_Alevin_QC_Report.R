#!/usr/bin/env Rscript


# Author: Aaron M Allen
# Date: 2023.01.24
#
#

# Description:
#
# Create an html report of the alevin alignments using the `alevinQC` package.
#
#
#
#
#



# Packages
library(alevinQC)

# Setup
job_id <- commandArgs(trailingOnly = TRUE)[2]
params_file <- commandArgs(trailingOnly = TRUE)[3]
if (is.na(params_file)) {
    params_file <- "pipeline_R_params.R"
}
params_path <- paste0("src/param_files/", params_file)
source(params_path)

local_path_mod <- ""
curr_machine <- Sys.info()["nodename"]
if (curr_machine == "mentok") {
    local_path_mod <- "/mnt/data/backups/cluster_backups/cbrg"
}
raw_path <- paste0(local_path_mod, raw_path)
my_sample <- commandArgs(trailingOnly = TRUE)[4]
output_dir <- paste0(raw_path, my_sample)


# Generate Report
alevinQCReport(baseDir = paste0(raw_path, my_sample),
               sampleId = my_sample,
               outputFile = paste0(my_sample, "-job_id_", job_id, "-alevin_report.html"),
               outputDir = output_dir,
               outputFormat = "html_document",
               showCode = TRUE,
               forceOverwrite = TRUE)
