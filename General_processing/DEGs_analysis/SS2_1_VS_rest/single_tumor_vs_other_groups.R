library(Seurat)
library(SeuratDisk)

args = commandArgs(trailingOnly=TRUE)
h5ad_path = args[1]
h5seurat_path = args[2]
tumor = args[3]
met_group = args[4]


print("converting scanpy to seurat")
Convert(h5ad_path, dest = "h5seurat", overwrite = TRUE)
scBC  <- LoadH5Seurat(h5sseurat_path)
scBC

print("running MAST")
Idents(scBC) <- "metastatic_potential_group"
scBC$tumor_name_numeric <- as.factor(scBC$Tumor_ID)
scBC$tumor_name_numeric <- as.numeric(scBC$tumor_name_numeric)
DEG_df <- FindMarkers(scBC, test.use = "MAST",ident.1 = met_group, ident.2 = NULL,latent.vars = c("tumor_name_numeric"),
                      logfc.threshold = 0,min.pct=0)

print("saving csv")
write.csv(DEG_df, paste0("1_VS_rest_" ,tumor, '.csv'))
print("Done")
