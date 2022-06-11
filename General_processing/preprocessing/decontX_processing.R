# loading the libaries 

library(celda)
library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)

args = commandArgs(trailingOnly=TRUE)
h5ad_path = args[1]
h5seurat_path = args[2]
tumor = args[3]

Convert(h5ad_path, dest = "h5seurat", overwrite = TRUE)
scBC  <- LoadH5Seurat(h5seurat_path)
scBC

scBC_to_scX <- as.SingleCellExperiment(scBC)
decontX_output <- decontX(scBC_to_scX)

decontX_output_counts <- assay(decontX_output,'decontXcounts')
write.csv(decontX_output_counts, paste0("decontx_counts_output/decontX_" ,tumor, '.csv'))