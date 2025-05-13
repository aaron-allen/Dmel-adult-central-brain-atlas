library(BiocParallel)
library(DropletUtils)
library(cowplot)
library(ggplot2)
library(tximport)
library(fishpond)
library(optparse)



theme_set(theme_bw())
multicoreParam <- MulticoreParam(workers = 8)
multicoreParam

options(stringsAsFactors = FALSE)

option_list <- list(
    make_option(c("--input_file"), default = "../../Alevin_r_6.30_allcb/mem_G1_rep2_rna/alevin/quants_mat.gz",
                help ="quants.mat.gz  filepath from alevin for a sample"),
    make_option(c("--output_dir"), default = "../../Alevin_r_6.30_allcb/mem_G1_rep2_rna/alevin",
                help ="sample alevin output directory"),
    make_option(c("--sample_label"), default = "mem_G1_rep2",
                help ="sample label"),
    make_option(c("--filter_fdr"), default = 0.001,
                type = 'numeric',
                help = "FDR filter for removal of barcodes with no cells"),
    make_option(
        c( "--filter_empty"),
        action = "store",
        default = 1,
        help = "Should barcodes estimated to have no cells be removed from the output object?"),
    
    make_option(
        c( "--filter_genes"),
        action = "store",
        default = 0,
        help = "should rRNA, rProt, mito and hsp genes be filtered out from the counts matrix before running empty drops"),
    make_option(
        c( "--plot_empty_pval"),
        action = "store",
        default = 1,
        help = "plot distirbution of p-values for assumed empty droplets"),
    
    make_option(
        c( "--lower"),
        action = "store",
        default = 100,
        type = 'numeric',
        help = "A numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets."
    ),
    make_option(
        c( "--niters"),
        action = "store",
        default = 10000,
        type = 'numeric',
        help = "An integer scalar specifying the number of iterations to use for the Monte Carlo p-value calculations."
    ),
    make_option(
        c( "--test_ambient"),
        action = "store",
        default = 0,
        help = "A logical scalar indicating whether results should be returned for barcodes with totals less than or equal to lower."
    ),
    make_option(
        c( "--ignore"),
        action = "store",
        default = 0,
        type = 'numeric',
        help = "A numeric scalar specifying the lower bound on the total UMI count, at or below which barcodes will be ignored."
    ),
    make_option(
        c( "--retain"),
        action = "store",
        default = 0,
        type = 'numeric',
        help = "A numeric scalar specifying the threshold for the total UMI count above which all barcodes are assumed to contain cells."
    ),
    make_option(c("--figures_dir"), default = "../../figures.dir/alevin_droplet_utils",
                help ="directory path for plots"),
    
    make_option(c("--barcode_plot_function"), default = "./Alevin_dropletBarcodePlot.R",
                help ="Rscript for barcode plot function")
)

opt <- parse_args(OptionParser(option_list = option_list, add_help_option = FALSE))

message("Running with options:")
print(opt)
sample <- opt$sample_label
cat("Running sample:", sample , fill = T)


output_dir_raw <- paste0(opt$output_dir,"/", "raw_cellranger_output")
output_dir_filt <- paste0(opt$output_dir, "/", "filtered_cellranger_output")

if (!dir.exists(output_dir_raw)) {
    dir.create(output_dir_raw, recursive = T)
}

if (!dir.exists(output_dir_filt)) {
    dir.create(output_dir_filt, recursive = T)
}

if (!dir.exists(opt$figures_dir)) {
    dir.create(opt$figures_dir, recursive = T)
}

opt_table <- data.frame(value=unlist(opt), stringsAsFactors = FALSE)

if(opt$test_ambient == 0){
    test_ambient <- FALSE
} else {
    test_ambient <- TRUE
}

if(opt$ignore == 0){
    ignore <- NULL
} else {
    ignore <- opt$ignore
}

if(opt$retain == 0){
    retain <- NULL
} else {
    retain <- opt$retain
}

cat( "modified option param test_ambient values are:", test_ambient, fill =T)
cat( "modified option param ignore values are:", ignore, fill =T)
cat( "modified option param retain values are:", retain, fill =T)


# Read in allCB alevin quants matrix as as sparse matrix using tximport
file.exists(opt$input_file)
txi <- tximport(opt$input_file, type = "alevin", sparse = T)

# write raw counts into 10x CellRanger version 3 format using write10xcounts

write10xCounts(x=txi$counts, path= output_dir_raw,type = "sparse", version = "3", overwrite = T)


# Run empty drops to get empty bead calls
sce <- read10xCounts(output_dir_raw, sample.names = sample, col.names = T, version= "3")
sce

if(opt$filter_genes==1){
    counts <- counts(sce)
    mito_genes <- grep(x = row.names(counts), pattern ="^mt:" , value = T)
#    rRNA <- grep(x= row.names(counts), pattern=  "(S|-)rRNA", value =T)
    rProt <- grep(x=row.names(counts), pattern="Rp(L|S)", value=T)
#    Hsp <- grep(x=row.names(counts), pattern="^Hsp", value=T)
    
    # qc is a list of character vectors with mitochondrial and ribosomal gene names
    keep <- !row.names(counts) %in% c(mito_genes, rProt)
    print(" mito and  rProt genes have been filtered from the counts matrix provided to empty drops")
    # emptyDrops performs Monte Carlo simulations to compute p-values,
    # so we need to set the seed to obtain reproducible results.
    set.seed(57676)
    #set.seed(100)
    #e.out <- DropletUtils::emptyDrops(counts[keep, ])
    e.out <- emptyDrops( m = counts[keep,], lower = opt$lower, niters = opt$niters, test.ambient = test_ambient, ignore = ignore, retain = retain)
    
    # See ?emptyDrops for an explanation of why there are NA values.
    print(summary(e.out$FDR <= opt$filter_fdr))
    rm(counts)
} else {
    print("no mito and  rProt  genes have been fiiltered out from the counts matrix provided to empty drops")
    counts <- counts(sce)
    
    
    # emptyDrops performs Monte Carlo simulations to compute p-values,
    # so we need to set the seed to obtain reproducible results.
    set.seed(57676)
    #set.seed(100)
    #e.out <- DropletUtils::emptyDrops(counts[keep, ])
    e.out <- emptyDrops( m = counts, lower = opt$lower, niters = opt$niters, test.ambient = test_ambient, ignore = ignore, retain = retain)
    
    # See ?emptyDrops for an explanation of why there are NA values.
    print(summary(e.out$FDR <= opt$filter_fdr)) 
    rm(counts)
}
#Removing mitochondrial and ribosomal genes so that they don't get used when determining if droplets are different from ambient droplets:




#The number of Monte Carlo iterations determines the lower bound for the p-values (Phipson and Smyth 2010). 
#The Limited field in the output indicates whether or not the computed  p-value for a particular barcode is bounded by the number of iterations.
#If any non-significant barcodes are TRUE for Limited, we may need to increase the number of iterations. 
#A larger number of iterations will result in a lower  p-value for these barcodes, which may allow them to be detected after correcting for multiple testing.


print(table(Sig=e.out$FDR <= opt$filter_fdr, Limited=e.out$Limited))








#As mentioned above, emptyDrops() assumes that barcodes with low total UMI counts are empty droplets. 
#Thus, the null hypothesis should be true for all of these barcodes. We can check whether the hypothesis testing procedure holds its size by examining the distribution of  
#p-values for low-total barcodes with test.ambient=TRUE

set.seed(106890)
if(opt$plot_empty_pval == 1){
    all.out <- emptyDrops(counts(sce), lower= opt$lower, test.ambient= TRUE) 
    pdf(file = paste0( opt$figures_dir,"/",opt$sample_label,"_","distribution_pvalues_assumed_emptydroplets.pdf"), height = 5, width = 7)
    hist(all.out$PValue[all.out$Total <= opt$lower & all.out$Total > 0],
         xlab="P-value", main="", col="grey80") 
    dev.off()
} else {
    print("plot empty pval option set to NO, so no plots produced")
}

# other diagnostic plots
# Droplets detected as cells should show up with large negative log-probabilities or 
#very large total counts (based on the knee point reported by barcodeRanks



is.cell <- e.out$FDR <= opt$filter_fdr & ! is.na(e.out$FDR)
cat(c(
    paste0(
        'At an FDR of ', opt$filter_fdr,', estimate that ',
        sum(is.cell, na.rm = TRUE),
        ' barcodes have cells.'
    ),
    ifelse( opt$filter_empty == 1, paste(" Workflow Will filter to",  sum(is.cell, na.rm = TRUE), 'barcodes.'), ''),
    '\nParameter values:',
    capture.output(print(opt_table))
), sep = '\n')

pdf(file = paste0( opt$figures_dir,"/",opt$sample_label,"_","total_counts_vs_neg_log_probability.pdf"), height = 5, width = 7)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
     xlab="Total UMI count", ylab="-Log Probability")
dev.off()



# plot standard barcode rank and barcod edensity plot for raw data
source(opt$barcode_plot_function)

bc_sce <- data.frame(V1 = colData(sce)$Barcode, V2=colSums(assays(sce)[[1]]))
bc_sce <- bc_sce[order(bc_sce$V2, decreasing = TRUE), ]

# Get the roryk cutoff
roryk_count_cutoff <- pick_roryk_cutoff(bc_sce$V2)

# Run dropletUtils' barcodeRanks to get knee etc
br.out <- barcodeRanks(counts(sce), lower = opt$lower)

# UMI count threshols
dropletutils_knee <- metadata(br.out)$knee
dropletutils_inflection <- metadata(br.out)$inflection
emptydrops_retain <- metadata(e.out)$retain
#e.out_is.cell <- e.out[which(is.cell),]
#emptydrops_cells <- max(e.out_is.cell$Total)


# plot the standard barcode rank and tthe density plot side by side
plots <-list(
    dropletutils = barcode_rank_plot(br.out, roryk_count_cutoff, dropletutils_knee, dropletutils_inflection,emptydrops_retain, name = sample),
    roryk = barcode_density_plot(bc_sce$V2, roryk_count_cutoff, dropletutils_knee, dropletutils_inflection,emptydrops_retain,name = sample)
)
plot_noleg <- plot_grid( plots$dropletutils,plots$roryk + theme(legend.position = "none"), nrow = 1)

# extract the legend from one of the plots
legend <- get_legend(
    # create some space to the left of the legend
    plots$roryk + theme(legend.box.margin = margin(0, 0, 0,15))
)


pdf(width = 9, height = 5, file= paste0( opt$figures_dir,"/",opt$sample_label,"_","barcode_rank_density_plots.pdf"))
plot_grid( plot_noleg, legend, nrow = 1, rel_widths = c(1.6,0.6))
dev.off()








# Filter empty cells if specified

# Subsetting the matrix to the cell-containing droplets.
# (using 'which()' to handle NAs smoothly).
if(opt$filter_empty == 1) {
   # colnames(e.out) <- paste0('empty', colnames(e.out))
    #colData(sce)[,colnames(e.out)] <- e.out
    sce_cells <- sce[,which(e.out$FDR <= opt$filter_fdr), drop=FALSE]
    
    count_cells <- counts(sce_cells)
    # Write 10x format files for detected cells for samples as filtered cellranger output
    
    write10xCounts(x=count_cells, path= output_dir_filt,type = "sparse", version = "3", overwrite = T)
    
} else {
    print("No filtered cellranger output produced")
}



# save empty drops object as rds within alevin folder for each sample
saveRDS(e.out, file = paste0(opt$output_dir,"/" ,sample, "_emptydrops_res",".rds"))



