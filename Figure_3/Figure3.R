#This script includes analysis and visualization of scRNAseq data published in Winkler et al. 2022

library(Seurat)
library(ggplot2)

col_tumor <- c("HCI010"="#1598cc","J53353"="#80b1d3","J2036"="#8dd3c7",
               "HCI001"="#fdb462","HCI009"="#fb8072","H5471"="#b3de69",
               "H4272"="#ccebc5","H3204"="#ffed6f","J55454"="#e8ee88",
               "HCI002"="#d9d9d9","HCI011"="#bebada","H5097"="#bc80bd",
               "HCI005"="#fccde5") 

#Figure 3D----------------

#complex heatmap of sleceted genes with average expression
library(ComplexHeatmap)
#see here: https://jokergoo.github.io/ComplexHeatmap-reference/book/
library(RColorBrewer)
library(circlize)

#complex heatmap for SS2 and 10x combined

#take only primary tumor of ss2
Idents(seu_ss2) <- seu_ss2@meta.data$sort
table(Idents(seu_ss2))
seu_ss2_PT <- seu_ss2[,Idents(seu_ss2)=="Tumor"]
table(seu_ss2_PT$Tumor_ID, seu_ss2_PT$sort)
table(seu_ss2_PT@meta.data$sort)

Idents(seu_10Xv2) <- seu_10Xv2@meta.data$sort_updated
table(Idents(seu_10Xv2))
seu_10Xv2_PT <- seu_10Xv2[,Idents(seu_10Xv2)=="tumor"]
table(seu_10Xv2_PT$Tumor_ID_updated, seu_10Xv2_PT$sort_updated)
table(seu_10Xv2_PT@meta.data$sort_updated)

#calculate average expression
select_genes <-c("HLA-A", "HLA-B", "HLA-C","HLA-E","B2M", "TAP1",
                 "NFKBIA", "PSMB3", "SQSTM1", "LAMP2",
                 "IFI6", "IFI35",
                 "TPM2", "BGN", "GJA1", "EMP3", "CTHRC1", "PCOLCE", "VIM",
                 "S100A4", "MUC1")

seu_av_tumorID_ss2 <- AverageExpression(seu_ss2_PT, group.by = "Tumor_ID", 
                                        features = select_genes, 
                                        return.seurat = TRUE)

seu_av_tumorID_10x <- AverageExpression(seu_10Xv2_PT, group.by = "Tumor_ID_updated", 
                                        features = select_genes, 
                                        return.seurat = TRUE)

seu_av_tumorID_10x@meta.data$Tumor_ID <- colnames(seu_av_tumorID_10x)
seu_av_tumorID_10x@meta.data$method <- "MULTI"
seu_av_tumorID_ss2@meta.data$Tumor_ID <- colnames(seu_av_tumorID_ss2)
seu_av_tumorID_ss2@meta.data$method <- "SS2"

#rename cells into unique cell names
seu_av_tumorID_10x <- RenameCells(seu_av_tumorID_10x, new.names = paste(colnames(seu_av_tumorID_10x),"10x", sep="_"))
seu_av_tumorID_ss2 <- RenameCells(seu_av_tumorID_ss2, new.names = paste(colnames(seu_av_tumorID_ss2),"ss2", sep="_"))

#create combined data set of ss2 and 10x for complex heatmap
k<-cbind(seu_av_tumorID_10x@assays$RNA@scale.data, seu_av_tumorID_ss2@assays$RNA@scale.data)

meta<-rbind(seu_av_tumorID_10x@meta.data, seu_av_tumorID_ss2@meta.data)

#add metastatic potential group
meta$metastatic_potential_group <- "Low"
meta$metastatic_potential_group[meta$Tumor_ID %in% c("HCI011", "HCI009","HCI001")] <- "Moderate"
meta$metastatic_potential_group[meta$Tumor_ID %in% c("H5097", "J53353", "J2036", "HCI010")] <- "High"

#order colnames in k accordingly 
head(rownames(meta))
head(colnames(k))

#set order for plotting correctly
myorder<-c("HCI002_10x", "J55454_10x","J55454_ss2", "H5471_ss2", "HCI005_10x", "HCI005_ss2", "H3204_ss2",
           "H4272_10x","H4272_ss2","HCI009_ss2", "HCI011_10x", "HCI011_ss2", "HCI001_10x", "HCI001_ss2",
           "H5097_10x", "H5097_ss2", "J2036_10x", "J2036_ss2", "J53353_10x", "J53353_ss2", "HCI010_10x","HCI010_ss2")
meta<-meta[myorder,]
k<-k[,myorder]
head(rownames(meta))
head(colnames(k))

#draw heatmap

ha = HeatmapAnnotation(tumor_ID = meta$Tumor_ID,
                       method=meta$method,
                       col = list(tumor_ID=col_tumor, method=col_method),
                       border = FALSE,gap = unit(2, "points"),
                       simple_anno_size = unit(0.2, "cm"),
                       show_annotation_name = FALSE, show_legend = FALSE,
                       annotation_name_side = "left")

ht_list <- Heatmap(k, name = "expression",  cluster_rows = FALSE, cluster_columns = FALSE,
                   show_column_names = FALSE,show_row_dend = FALSE, show_column_dend =FALSE, 
                   show_heatmap_legend = FALSE,
                   col = colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7","#B2182B")),
                   top_annotation = ha,
                   row_names_gp = gpar(fontsize = 7))
draw(ht_list)
dev.off()

#Supplementary Figure 3G--------

#add proliferation score
#https://www.nature.com/articles/s41523-021-00216-w
proliferation <- list(c("ANLN", "CCNE1", "CDC20", "CDC6", "CDCA1", "CENPF", "CEP55", "EXO1", "KIF2C", "KNTC2", "MELK",
                        "MKI67", "ORC6L", "PTTG1", "RRM2", "TYMS", "UCE2C","UBE2T"))

seu_10Xv2 <- AddModuleScore(seu_10Xv2, features = proliferation, name="proliferation")

#use primary tumor only!

p<-RidgePlot(seu_10xv2_PT, features = "proliferation1", group.by = "mets")+ 
  NoLegend()+NoAxes()+
  scale_fill_manual(name="mets", values = c("yes"="#666666", "no"="#BCBBBC"))+
  theme(plot.title = element_blank())

#barplot mets and cell cycle phase
data.plot <- as.data.frame(table(seu_10xv2_PT@meta.data$phase, seu_10xv2_PT@meta.data$mets))
colnames(data.plot) <- c("phase", "mets", "cellcount")

p<-ggplot(data = data.plot, aes(x=mets, y=cellcount, fill=phase))+
  geom_bar(stat="identity", position = "fill")+
  scale_fill_brewer(palette = "Set2")+
  theme_minimal() +
  theme(panel.grid=element_blank(),
        panel.border = element_blank(),
        title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.2, linetype = "solid"),
        legend.position="none",
        axis.line  = element_line(colour = "black", size = 0.2, linetype = "solid"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

#Supplementary Figure 3F------------

#add signature score
Myc_Lee2022 <- read.table(file = "Goga_Mouse_Sig3_biomaRtMapped_2020_08_03.txt", 
                       header = TRUE, sep = "",  dec = ".")
Myc_Lee2022  <-list(Myc_Lee2022$HGNC.symbol)
antigen_own <- list(c("HLA-A", "HLA-B", "HLA-C","HLA-E","B2M", "TAP1",
                      "NFKBIA", "PSMB3", "SQSTM1", "LAMP2",
                      "IFI6", "IFI35"))

seu_10Xv2_PT<- AddModuleScore(object=seu_10Xv2_PT, features = Myc_Lee2022, ctrl = 10, name = "myc_Lee")
seu_10Xv2_PT<- AddModuleScore(object=seu_10Xv2_PT, features = antigen_own, ctrl = 10, name = "antigen_own")

#visualize in scatter plot
p<-FeatureScatter(seu_10Xv2_PT, feature1 = "antigen_own1", feature2 = "myc_joyce1", 
                  group.by = "Tumor_ID_updated", cols = col_tumor, pt.size = 0.1)+ NoLegend()+
  theme_minimal() +
  theme(panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.ticks =element_line(size = 0.25),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position="none",
        axis.line  = element_line(size = 0.25))


