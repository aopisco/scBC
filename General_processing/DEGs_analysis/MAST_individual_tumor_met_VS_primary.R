
library(Seurat)
library(SeuratDisk)

Convert("SS2_processed.h5ad", dest = "h5seurat", overwrite = TRUE)
scBC_SS2 <- LoadH5Seurat("SS2_processed.h5seurat")

Idents(scBC_SS2) <- "Tumor_ID"

DEG_df <- FindMarkers(scBC_SS2, test.use = "MAST",ident.1 = "Metastatic", ident.2 = "Tumor",logfc.threshold = 0,min.pct=0)

write.csv(DEG_df, 'MAST_global_met_vs_primary.csv')

tumor_list <- c('H3204',
                'H4272',
                'H5097',
                'H5471',
                'HCI001',
                'HCI005',
                'HCI009',
                'HCI010',
                'HCI011',
                'J2036',
                'J53353',
                'J55454')

for (tumor in tumor_list){
  tumor_sub = subset(scBC_SS2, idents = c(tumor))
  Idents(tumor_sub) <- "sort"
  DEG_df <- FindMarkers(tumor_sub, test.use = "MAST",ident.1 = "Metastatic", ident.2 = "Tumor",logfc.threshold = 0,min.pct=0)
  write.csv(DEG_df, paste0("met_vs_primary" ,tumor, '.csv'))
  
}