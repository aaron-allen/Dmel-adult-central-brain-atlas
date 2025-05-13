setwd("/project/cncb/shared/proj136/analyses/scenic/scenic_multi_dsx/")
library(Seurat)
library(tidyverse)
library(data.table)

dsx_data <- readRDS(file="data.dir/...")
write_rds(x = dsx_data, file = "data.dir/dsx_scenic.rds")

# dsx_data_normalised <- GetAssayData(object= dsx_data,slot="data")
# dsx_data_normalised <- as.matrix(dsx_data_normalised)
# write.table(dsx_data_normalised, file="data.dir/dsx_normalised-expression.csv", sep=",",quote=F)
# # fwrite(dsx_data_normalised, "data.dir/dsx_normalised-expression.csv", sep=",",quote=F)


dsx_data_raw <- GetAssayData(object=dsx_data, slot="counts", assay="RNA")
dsx_data_raw <- as.matrix(dsx_data_raw)
write.table(dsx_data_raw, file="data.dir/dsx_raw-expression.csv", sep=",",quote=F)
# fwrite(dsx_data_raw, file="data.dir/dsx_raw-expression.csv", sep=",",quote=F)

metadata <- dsx_data@meta.data

meta_ann <- data.frame(rownames(metadata), metadata$cell_type)
colnames(meta_ann)[1] <- "cellid"
colnames(meta_ann)[2] <- "cell_type"
write.table(meta_ann, file="data.dir/dsx_cell_type-annotation.csv", sep=",",quote=F, row.names = F)
# fwrite(meta_ann, file="data.dir/dsx_cell_type-annotation.csv", sep=",",quote=F, row.names = F)


meta_sex <- data.frame(rownames(metadata), metadata$sex)
colnames(meta_sex)[1] <- "cellid"
colnames(meta_sex)[2] <- "sexed"
write.table(meta_sex, file="data.dir/dsx_sexed-annotation.csv", sep=",",quote=F, row.names = F)
# fwrite(meta_sex, file="data.dir/dsx_sexed-annotation.csv", sep=",",quote=F, row.names = F)


metadata$meta_sex <- paste(metadata$cell_type,"_", metadata$sex,sep = "")
meta_sex <- data.frame(rownames(metadata), metadata$meta_sex)
colnames(meta_sex)[1] <- "cellid"
colnames(meta_sex)[2] <- "cell_type_sexed"
write.table(meta_sex, file="data.dir/dsx_cell_type_sexed-annotation.csv", sep=",",quote=F, row.names = F)
# fwrite(meta_sex, file="data.dir/dsx_cell_type_sexed-annotation.csv", sep=",",quote=F, row.names = F)


# rm(dsx_data, dsx_data_normalised, )
