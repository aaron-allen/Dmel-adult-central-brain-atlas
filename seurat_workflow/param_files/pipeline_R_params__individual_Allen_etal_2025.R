#!/usr/bin/env Rscript


# Author: Aaron M Allen
# Date: ?
#     updated: 2022.12.15
#
#

# Description:

## This is a set of parameters for the scRNA-seq pipeline that will be sourced
## at the start of all the Rmarkdown files.
##
##
## Broadly speaking, the parameters script is divided into 7 sections:
##    1. Setup - file paths, sample names, etc.
##    2. Initial filtering - trimming by number of cells, number of features, etc.
##    3. Doublet removal - using doubletfinder, scrublet, and solo
##    4. Ambient RNA removal - using SoupX or decontX
##    5. Normalization and integration - lognorm or SCT, and Seurat FindAnchors etc or Harmony
##    6. Clustering - ... ya know, clustering ...
##    7. Annotation - using cell type markers, reference sc data, and bulk sequencing data
##


## Caveats:
## I have not rigorously tested all processing options of this "pipeline". I primarily
## use alevin for alignment (not CellRanger), log-normalization (and not SCT),
## SoupX (and not decontX), and Harmony (and not Seurat's FindAnchors etc). So there
## might still be a few bugs in the other options. I also have not tested this with
## very small datasets, so there might be some issues with that as well.








###########################################################
###########################################################
###########################################################
###   Setup   #############################################
###########################################################
###########################################################
###########################################################


# data loading

# rename to something like "dataset" or "experiment"
dataset <- "cns"

# need to rethink this parameter doesn't quite work in new workflow
input_type <- "alevin"     # options - "cellranger","alevin"

# paths etc
raw_path <- file.path("/project/cncb/shared/proj068/analyses/alevin/fru_iso/Alevin_r_6.30_allcb/")
raw_names <- list("FemaleBrainRep1_rna", "FemaleBrainRep2_rna","MaleBrainRep1_rna", "MaleBrainRep2_rna")
sample_names <- list("brain_fem_rep1", "brain_fem_rep2","brain_male_rep1", "brain_male_rep2")

objects_path <- file.path("../analyses/rds_files/")


# Run number parameters.
# Some of the scripts are run multiple times depending on whether doublets
# are removed, whether ambient is removed, and whether subclustering is done.
# Not entirely sure if I still need this, but it's here for now.
run_number_1a <- 1
run_number_1b <- 1
run_number_4a <- 1
run_number_4b <- 1
run_number_4c <- 1
run_number_4d <- 1


# Round the data to integers.
# Alevin, Soupx, and decontX all generate non-integer values.
# If this is set to then the data are stochastically/probabilistically
# rounded to integers (as in SoupX's "roundToInt" option)
round_the_data <- TRUE


# Parallization parameters
parallelize <- TRUE
futures_plan <- "multisession"      # still not sure what the best option is for the cluster
futures_n_cores <- 4
futures_ram <- 20
doubletfinder_n_cores <- 10
tsne_n_cores <- 28


# Cell types and Markers:
gmt_file_path <- "/project/cncb/aallen/gene_lists/gmt/D_melanogaster_celltype_markers.gmt"
# Cell types for which the first marker in the gmt file will be plotted on umaps/tsnes to visualize embedding and clustering.
qc_plot_cell_type <- c("Cholinergic", "Glutamatergic", "GABAergic", "Monoaminergic", "Kenyon_cell", "Glia")


# Adding meta data:

# sexing cells
insilico_sex <- TRUE          # use lncRNA::roX1 and lncRNA::roX2 to sex the cells
sample_sex <- FALSE           # grep the sample names for female/male

# add tissue/dataset meta data if multiple things are combined
# will extract the characters to first underscore in sample name
multi_tissue <- FALSE








###########################################################
###########################################################
###########################################################
###   Initial filetering   ################################
###########################################################
###########################################################
###########################################################


# Filtering params
min_tot_cells <- 200   # used to remover replicates with too few cells
min_cells <- 0         # used while creating seurat object
min_features <- 0      # used while creating seurat object
min_umi  <- 500        # used to filter seurat object
min_genes <- 300       # used to filter seurat object
max_umi <- Inf         # used to filter seurat object
max_genes <- Inf       # used to filter seurat object
max_mito <- 15         # used to filter seurat object
max_rrna <- 100        # used to filter seurat object
max_rprot <- 100       # used to filter seurat object
max_hsp <- 15          # used to filter seurat object

# which params get plotted
qc_params_plot <- c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_hsp", "percent_rProt", "percent_rRNA")
qc_params_plot_thesh <- c(4000, 10000, 15, 5, 40, 2)

# 10x Chromium v3 seems to get a lot of PolIII features, if you want to remove them set this to TRUE.
remove_poliii_genes <- TRUE

# if alignments were run with a custom gene model with different isoforms split off, and if
# you don't want them included in the genes used for PCA etc., then list them here.
# Also list any other (trans)genes you may want to be ignored for the clustering etc.
iso_genes_remove <- c("fru-plus", "fru-dP1", "fru-dP1-plus", "fruA", "fruB", "fruC",
                      "fruD1", "fruD2", "fruABint", "fruBCint", "fruCint", "fruPr1Ex1",
                      "fruPr1Ex2", "fruCOM", "dsxM", "dsxF",
                      "EGFP", "GFP", "mCherry", "Gal4", "GAL4",
                      "lncRNA:roX1", "lncRNA:roX2"
                      )





###########################################################
###########################################################
###########################################################
###   Doublet Removal    ##################################
###########################################################
###########################################################
###########################################################


remove_doublets <- FALSE  # TRUE, FALSE

###### doublet comparisons
doublet_path <- file.path("../analyses/doublets/")

##### doubletfinder
doubletfinder_output_path <- file.path("../analyses/doublets/doubletfinder/")
doubletfinder_calls <- file.path(paste0(doubletfinder_output_path,
                                        dataset, "_",
                                        input_type, "_",
                                        "doubletfinder_calls.csv")
)

##### scrublet solo prep
scrublet_solo_prep_output_path <- file.path("../analyses/doublets/filtered_inputs/")
scrublet_calls_path <- file.path("../analyses/doublets/scrublet/")
solo_calls_path <- file.path("../analyses/doublets/solo/")


# Cells to remove:
## How many methods need to call a cell a doublet for you to want to remove them?
## 0 = don't remove any cells (but then why are you running this ..?..)
## 1 = union of all methods
## 2 = at least 2 methods need to think it is a doublet
## 3 = strict intersection of all 3 methods
cells_remove <- 2


# pre-process individual samples
use_umi <- TRUE

# Note: we aren't using vars.to.regress for the initial processing as what was regressed isn't written to the object,
#   and DoubletFinder only looks at the object for when it reruns the Seurat analysis internally





###########################################################
###########################################################
###########################################################
###   Ambient Removal    ##################################
###########################################################
###########################################################
###########################################################

##### soupx and decontX

remove_ambient <- FALSE     #  TRUE or FALSE
ambient_method <- "soupx"  # "soupx" or "decontX"

##### soupx
soupx_output_path <- file.path("../analyses/ambient/soupx/")
##### decontX
decontx_output_path <- file.path("../analyses/ambient/decontx/")

# min and max UMI for soupx to use as empty droptlets
soup_min <- 0
soup_max <- 100

# total number of UMI of the sum of all markers with a group to consider "real" expression when plotting Euler diagrams
# and Upset plots
coexpression_thresh <- 2

# Used in ambient removal scripts. These cell types need to be truely mutually exclusive, and the marker genes for
# these cell types in the gmt file must only be in these cell types.
mutually_exclusive_markers_celltypes <- c("Cholinergic", "Glutamatergic", "GABAergic", "Glia")





###########################################################
###########################################################
###########################################################
###   Normalize and Integrate      ########################
###########################################################
###########################################################
###########################################################


# Upon interations of running the pipeline, you may find that some cells that you want to remove.
# Write out a list of the cells that you want to keep, and it can be read in here.
user_cells_remove <- FALSE
user_cells_remove_file <- "../analyses/cells_to_keep.csv"


neighbours <- 30             # used for RunUMAP and FindNeighbors
n_var_features <- 4000       # number of variable features to use for PCA

# Normalization and variable feature selection
normalization <- "LogNormalize"           # options "SCT", "LogNormalize"
integration <- "harmony"                  # options "harmony", "cca"
var_gene_method <- "mean.var.plot"        # options "mean.var.plot", "vst" ,"dispersion"
min_mean_cutoff <- 0.001                  # cutoffs for var_gene_method <- "mean.var.plot"
max_mean_cutoff <- Inf                    # cutoffs for var_gene_method <- "mean.var.plot"
min_dispersion_cutoff <- 0.1              # cutoffs for var_gene_method <- "mean.var.plot"
max_dispersion_cutoff <- Inf              # cutoffs for var_gene_method <- "mean.var.plot"

# number of PCs to use for clustering
curr_pc <- 60                # "ideal" pc
min_pc <- 10                 # min pc to test in range of pcs
max_pc <- 90                 # max pc to test in range of pcs
step_pc <- 10                # step size from min to max pc to test in range of pcs

# FindIntegrationAnchors params
nn_method <- "rann"                       # options "annoy", "rann"
nn_eps <- 0                               # used if nn_method <- "rann"
nn_trees <- 200                           # used if nn_method <- "annoy"

# RunHarmony params
harmony_group <- c("orig.ident", "sex")
harmony_theta <- rep(2, length(harmony_group))     # Diversity clustering penalty parameter
harmony_lambda <- rep(1, length(harmony_group))    # Ridge regression penalty parameter
harmony_nclust <- 200                              # Number of clusters in model
harmony_tau <- 30                                  # expected number of cells per cluster

# ScaleData regression
vars_to_regress <- c("orig.ident", "nCount_RNA", "percent_mt", "percent_hsp", "sex")
sct_vars_to_regress <- setdiff(vars_to_regress, c("orig.ident", "nCount_RNA", "sex"))





###########################################################
###########################################################
###########################################################
#########   Clustering    #################################
###########################################################
###########################################################
###########################################################


# Which clustering algorithm to use?
# 1 = original Louvain algorithm;
# 2 = Louvain algorithm with multilevel refinement;
# 3 = SLM algorithm;
# 4 = Leiden algorithm
clustering_algorithm <- 4

# series of cluster resolutions to test
resolution_test <- c(seq(1, 5, 0.5), seq(6, 20, 2), seq(25, 50, 5))
curr_res <- 10               # "ideal" cluster resolution

# FindClusters params
leiden_graph_method <- "igraph"
group_singletons <- FALSE

# cluster markers params
test_use <- "wilcox"
avg_logFC_thres <- 2
p_val_adj_thres <- 0.05
lim_n_sig <- 5
min_pct <- 0.25
min_cluster_size <- 50
only_pos <- TRUE





###########################################################
###########################################################
###########################################################
##############   Annotation    ############################
###########################################################
###########################################################
###########################################################


# Annotation parameters for modulescore and k-means
n_crt_features <- 100
annotate_by_K <- TRUE
my_k <- 2
my_k_itr <- 100

# Core cell types that will be annotated. The way we do this is in a "competitive" fashion so these groups should be
# mostly exclusive, but more importantly every cell should belong to one of the groups. This approach doesn't work as
# well when there are "gaps" in the annotations. It also should be noted that these "cell types" should robustly cluster
# separately from each other. If they co-cluster, the cluster level annotations will suffer.
annotation_cell_type <- c("Cholinergic", "Glutamatergic", "GABAergic", "Monoaminergic", "Kenyon_cell", "Glia")

# To annotate the cells, we can filter module scores > 0, > k-means derived threshold, or manually supplied threshold.
# These values should be on the l2-normalized scale.
manual_thresh <- FALSE
manual_thresholds <- c(0.003, 0.006, 0.005, 0.015, 0.007, 0.01)

# Cell types that shouldn't really be in the data, but might be present due to contamination during dissection
contamination_cell_types <- c("Salivary_gland", "Sperm_etc", "Fat_body", "Trachea", "Oenocyte", "Hemocyte", "Epithelia", "Muscle")

# "annotation_broad_celltype", "annotation_broad_cell_type_avg", "annotation_broad_cell_type_ext"
subcluster_by <- "annotation_broad_cell_type_ext_fix"

# sublcustering the following - "Cholinergic","Glutamatergic","GABAergic","Monoaminergic","Kenyon_cell","Motor_neurons","Glia","NB_INP"
subcluster_pcs <- c(50, 40, 40, 30, 30, 30, 30, 30)
subcluster_res <- c(4, 3, 3, 2, 2, 2, 2, 2)

# Annotate from other scRNA-seq data:
# Path to rds objects of reference rds objects
reference_object_path <- file.path("/project/cncb/aallen/reference_datasets/single_cell/")

# Annotate from bulk:
# Path to FBgn to gene symbol csv.
fbgn_gene_path <- "/project/cncb/aallen/genomes/dmel/FBgn_to_symbol.csv"

# Path to bulk seq data sets.
turner_evans_path <- "/project/cncb/aallen/reference_datasets/bulk/GSE155329_Turner-Evans_2021_CX/GSE155329_Turner-Evans_BulkRNASeq_updated.xlsx"
konstantinides_path <- "/project/cncb/aallen/reference_datasets/bulk/GSE103772_Konstantinides_2018_optic_lobe/FACS_sorted_cell_types_RNASeq.txt"
davis_path <- "/project/cncb/aallen/reference_datasets/bulk/GSE116969_Davis_2020_optic_lobe_KC/GSE116969_dataTable7a.genes_x_cells_TPM.modeled_genes.txt.gz"
