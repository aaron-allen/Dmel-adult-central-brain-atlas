


curr_dir <- getwd()
print(paste0("starting in ", curr_dir))
setwd(paste0(curr_dir, "/src/"))
curr_dir <- getwd()
print(paste0("now in ", curr_dir))

params_file <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(params_file)) {
    params_file <- "pipeline_R_params.R"
}
params_path <- paste0("param_files/", params_file)
print(paste0("The params file is :   ", params_path))
source(params_path)

dir.create("../analyses/logs", recursive = TRUE, showWarnings = FALSE)
dir.create("../analyses/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("../analyses/inputs", recursive = TRUE, showWarnings = FALSE)
dir.create("../analyses/markers", recursive = TRUE, showWarnings = FALSE)
dir.create("../analyses/rds_files", recursive = TRUE, showWarnings = FALSE)
dir.create("../docs", recursive = TRUE, showWarnings = FALSE)


if (remove_doublets) {
    dir.create(doubletfinder_output_path, recursive = TRUE, showWarnings = FALSE)
    dir.create(scrublet_solo_prep_output_path, recursive = TRUE, showWarnings = FALSE)
    dir.create(scrublet_calls_path, recursive = TRUE, showWarnings = FALSE)
    dir.create(solo_calls_path, recursive = TRUE, showWarnings = FALSE)
# for (i in seq_along(sample_names)) {
# 	dir.create(paste0(scrublet_solo_prep_output_path, sample_names[[i]]), recursive = TRUE, showWarnings = FALSE)
# }
}

if (remove_ambient) {
    if (ambient_method == "soupx") {
        dir.create(soupx_output_path, recursive = TRUE, showWarnings = FALSE)
    }
    if (ambient_method == "decontx") {
        dir.create(decontx_output_path, recursive = TRUE, showWarnings = FALSE)
    }
}
