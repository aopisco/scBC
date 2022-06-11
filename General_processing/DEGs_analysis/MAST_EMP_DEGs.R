
library(Seurat)
library(SeuratDisk)

Convert("SS2_processed.h5ad", dest = "h5seurat", overwrite = TRUE)
scBC_SS2 <- LoadH5Seurat("SS2_processed.h5seurat")

Convert("MULTI_processed.h5ad", dest = "h5seurat", overwrite = TRUE)
scBC_MULTI <- LoadH5Seurat("MULTI_processed.h5seurat")

Idents(scBC_SS2) <- "EMP_stage"
Idents(scBC_MULTI) <- "EMP_stage"

EMT_list <- c('EMP Intermediate', 'Epithelial-like', 'Mesenchymal-like')

scBC_SS2$tumor_name_numeric <- as.factor(scBC_SS2$Tumor_ID)
scBC_SS2$tumor_name_numeric <- as.numeric(scBC_SS2$tumor_name_numeric)

scBC_MULTI$tumor_name_numeric <- as.factor(scBC_MULTI$Tumor_ID)
scBC_MULTI$tumor_name_numeric <- as.numeric(scBC_MULTI$tumor_name_numeric)

for (EMT in EMT_list){
  DEG_df <- FindMarkers(scBC_SS2, test.use = "MAST",ident.1 = EMT, ident.2 = NULL,latent.vars = c("tumor_name_numeric"), logfc.threshold = 0, min.pct = 0)
  write.csv(DEG_df, paste0("SS2_global_" ,EMT, '.csv'))
}


for (EMT in EMT_list){
  DEG_df <- FindMarkers(scBC_MULTI, test.use = "MAST",ident.1 = EMT, ident.2 = NULL,latent.vars = c("tumor_name_numeric"), logfc.threshold = 0, min.pct = 0)
  write.csv(DEG_df, paste0("MAST_global_" ,EMT, '.csv'))
}


