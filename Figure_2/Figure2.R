##This script includes analysis and visualization of scRNAseq data published in Winkler et al. 2022

library(Seurat)
library(ggplot2)


#Pathway enrichment analysis for Figure 2B----------------
#pathway enrichment analysis using fgsea for DEGs between primary tumor and metastatic cells
#smart seq2

#negative avg_log2FC means upregulated in primary tumor
#positive avg_log2FC means upregulated in metastasis
#I am ordering avg_log2FC increasing (starting with negative avg_log2FC)
#that means
#positive NES means enriched in metastasis
#negative NES means enriched in primary tumor

# see this reference: https://bioconductor.org/packages/release/bioc/manuals/fgsea/man/fgsea.pdf
# This is also a more thorough vignette:
# https://davetang.org/muse/2018/01/10/using-fast-preranked-gene-set-enrichment-analysis-fgsea-package/


library(fgsea)
library(org.Hs.eg.db)
library(reactome.db)
library(ggplot2)
library(msigdbr)

#gene set group from msigdbr
#### Option 1 - Run fgsea against the genesets from the msigdb 
# This page gives the category and subcategory descriptions: https://www.gsea-msigdb.org/gsea/msigdb/annotate.jsp
# Choose the sets you want to test, use the command "?msigdbr" for list of category and subcategory choices


m_df = msigdbr(species = "Homo sapiens", category = "H" )#HALLMARK
#m_df = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP" )#GO biological process
#m_df = msigdbr(species = "Homo sapiens", category = "C5" )#GO 

m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

#GLOBAL-----------------------
DGE_global <- read.csv("MAST_global_sort.csv", row.names = 1)
dim(DGE_global)#12307

#enrichment analysis
DGE_global <- DGE_global[order(DGE_global[,"avg_log2FC"]),]
genelist <- DGE_global$avg_log2FC
names(genelist) <-rownames(DGE_global)
#tail(genelist)
#length(genelist)
#head(genelist)

fgsea_res <- fgsea(pathways = m_list, 
                   stats = genelist,
                   minSize=15,
                   maxSize=500,
                   nperm=10000)

fgsea_H_DGE_global_nosubset <-fgsea_res
fgsea_H_DGE_global <- subset(fgsea_res[fgsea_res$pval <=0.05,])
fgsea_H_DGE_global <- fgsea_H_DGE_global[order(fgsea_H_DGE_global[,"NES"], decreasing = TRUE),]

#fgsea_GO_DGE_global_nosubset <-fgsea_res
#fgsea_GO_DGE_global <- subset(fgsea_res[fgsea_res$pval <0.05,])
#fgsea_GO_DGE_global <- fgsea_GO_DGE_global[order(fgsea_GO_DGE_global[,"ES"], decreasing = TRUE),]

##Figure 2B visualize pathway enrichment of global DEGs primary tumor vs metastasis -------
#plot top enriched pathways

data.plot <- fgsea_H_DGE_global[fgsea_H_DGE_global$pathway %in% topPathways,]
#data.plot <- fgsea_GO_DGE_global[fgsea_GO_DGE_global$pathway %in% topPathways,] #Supplementary Figure 2B

#remove HALLMARK
data.plot$pathway <- gsub("HALLMARK_", "", data.plot$pathway)

#top 10 significant enriched pathways bar plot for figure 2B
ggplot(data=data.plot, 
          aes(x=NES, y=reorder(pathway, NES, sum), fill=cut(NES, c(-Inf,0, Inf))))+
  geom_bar(stat="identity")+
  theme_minimal()+
  scale_fill_manual(name="NES", values=c("(-Inf,0]"="#ff7f00", "(0, Inf]"="#377eb8"))+
  theme(panel.grid=element_blank(),
        axis.line  = element_line(colour = "black", size = 0.2, linetype = "solid"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=5, color="black"),
        axis.text.x = element_text(size=5, color="black"),
        axis.text.y = element_text(size=5, color="black"),
        axis.ticks = element_line(colour = "black", size = 0.2, linetype = "solid"),
        legend.position="none") 

#Supplementary Figure 2C----------
#complex heat map for genes involved in top enriched pathways 

fgsea_H_DGE_global[order(fgsea_H_DGE_global[,"NES"]), "pathway"][[1]]
genes.tumor <-fgsea_H_DGE_global[order(fgsea_H_DGE_global[,"NES"]), "leadingEdge"][[1]]
genes.tumor.hypoxia <- genes.tumor[[1]]
genes.met <-fgsea_H_DGE_global[order(-fgsea_H_DGE_global[,"NES"]), "leadingEdge"][[1]]
genes.met.MYC <- genes.met[[1]]
fgsea_H_DGE_global[order(-fgsea_H_DGE_global[,"NES"]), "pathway"][[1]]

library(ComplexHeatmap)
library(circlize)          

j<-as.matrix(seu_ss2@assays$RNA@scale.data)
met.MYC <- j[rownames(j) %in% genes.met.MYC,]
tumor.hypoxia <- j[rownames(j) %in% genes.tumor.hypoxia,]

meta<-seu_ss2@meta.data

col_tumor <- c("HCI010"="#1598cc","J53353"="#80b1d3","J2036"="#8dd3c7",
               "HCI001"="#fdb462","HCI009"="#fb8072","H5471"="#b3de69",
               "H4272"="#ccebc5","H3204"="#ffed6f","J55454"="#e8ee88",
               "HCI002"="#d9d9d9","HCI011"="#bebada","H5097"="#bc80bd",
               "HCI005"="#fccde5")
col_tissue <- c("Tumor" = "#ff7f00", "Metastatic" = "#377eb8")
meta$sort <- factor(meta$sort, levels = c("Tumor", "Metastatic"))

ha = HeatmapAnnotation(tumor_ID = meta$Tumor_ID, type = meta$sort,
                       col = list(type = col_tissue, tumor_ID=col_tumor),
                       annotation_name_side = "left",
                       border = TRUE,gap = unit(3, "points"),
                       simple_anno_size = unit(0.5, "cm"))

#heatmap for MYC genes up in met
heatmap_plot <- Heatmap(met.MYC, name = "expression",  cluster_rows = TRUE, cluster_columns = TRUE,
                   show_column_names = FALSE,
                   show_row_dend = FALSE, show_column_dend =FALSE,
                   column_split=factor(meta$sort, levels = c("Tumor", "Metastatic")),
                   col = colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7","#B2182B")),top_annotation = ha)
draw(heatmap_plot)

#heatmap for hypoxia genes up in primary tumor
heatmap_plot <- Heatmap(tumor.hypoxia, name = "expression",  cluster_rows = TRUE, cluster_columns = TRUE,
                        show_column_names = FALSE,
                        show_row_dend = FALSE, show_column_dend =FALSE,
                        column_split=factor(meta$sort, levels = c("Tumor", "Metastatic")),
                        col = colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7","#B2182B")),top_annotation = ha)

draw(heatmap_plot)
dev.off()

#Top enriched pathways enrichment plots
#change the index below to plot any of the results
#top tumor enriched (Hypoxia) - negative NES 
p<- plotEnrichment(m_list[[fgsea_res[order(NES), ][1,]$pathway]],genelist) + 
  labs(title=fgsea_res[order(NES), ][1,]$pathway)+
  theme(title = element_text(size=12),
        axis.text = element_text(size=12))
#top metastasis enriched (MYC) - positive NES
p<- plotEnrichment(m_list[[fgsea_res[order(-NES), ][1,]$pathway]],genelist) + 
  labs(title=fgsea_res[order(-NES), ][1,]$pathway)+
  theme(title = element_text(size=12),
        axis.text = element_text(size=12))

#Pathway enrichment analysis for Figure 2E----------------
#perform fgsea for all tumors individually  
#J55454---------------------
# too few cells

#HCI005---------------------------------
DGE_HCI005 <- read.csv("MAST_HCI005.csv", row.names = 1)
dim(DGE_HCI005)#12362
head(DGE_HCI005)

DGE_HCI005 <- DGE_HCI005[order(DGE_HCI005[,"avg_log2FC"]),]

genelist <- DGE_HCI005$avg_log2FC
names(genelist) <-rownames(DGE_HCI005)

fgsea_res <- fgsea(pathways = m_list, 
                   stats = genelist,
                   minSize=15,
                   maxSize=500,
                   nperm=10000)

#fgsea_H_DGE_HCI005_nosubset <-fgsea_res
#fgsea_H_DGE_HCI005 <- subset(fgsea_res[fgsea_res$pval <0.05,])
#fgsea_H_DGE_HCI005 <- fgsea_H_DGE_HCI005[order(fgsea_H_DGE_HCI005[,"ES"], decreasing = TRUE),]

fgsea_GO_DGE_HCI005_nosubset <-fgsea_res
fgsea_GO_DGE_HCI005 <- subset(fgsea_res[fgsea_res$pval <0.05,])
fgsea_GO_DGE_HCI005 <- fgsea_GO_DGE_HCI005[order(fgsea_GO_DGE_HCI005[,"ES"], decreasing = TRUE),]


#H4272-----------------------
DGE_H4272 <- read.csv("MAST_H4272.csv", row.names = 1)
dim(DGE_H4272)#12470

DGE_H4272 <- DGE_H4272[order(DGE_H4272[,"avg_log2FC"]),]
genelist <- DGE_H4272$avg_log2FC
names(genelist) <-rownames(DGE_H4272)

fgsea_res <- fgsea(pathways = m_list, 
                   stats = genelist,
                   minSize=15,
                   maxSize=500,
                   nperm=10000)

#fgsea_H_DGE_H4272_nosubset <-fgsea_res
#fgsea_H_DGE_H4272 <- subset(fgsea_res[fgsea_res$pval <0.05,])
#fgsea_H_DGE_H4272 <- fgsea_H_DGE_H4272[order(fgsea_H_DGE_H4272[,"ES"], decreasing = TRUE),]

fgsea_GO_DGE_H4272_nosubset <-fgsea_res
fgsea_GO_DGE_H4272 <- subset(fgsea_res[fgsea_res$pval <0.05,])
fgsea_GO_DGE_H4272 <- fgsea_GO_DGE_H4272[order(fgsea_GO_DGE_H4272[,"ES"], decreasing = TRUE),]


#H3204-----------------------
DGE_H3204 <- read.csv("MAST_H3204.csv", row.names = 1)
dim(DGE_H3204)#12372

DGE_H3204 <- DGE_H3204[order(DGE_H3204[,"avg_log2FC"]),]
genelist <- DGE_H3204$avg_log2FC
names(genelist) <-rownames(DGE_H3204)

fgsea_res <- fgsea(pathways = m_list, 
                   stats = genelist,
                   minSize=15,
                   maxSize=500,
                   nperm=10000)

#fgsea_H_DGE_H3204_nosubset <-fgsea_res
#fgsea_H_DGE_H3204 <- subset(fgsea_res[fgsea_res$pval <0.05,])
#fgsea_H_DGE_H3204 <- fgsea_H_DGE_H3204[order(fgsea_H_DGE_H3204[,"ES"], decreasing = TRUE),]

fgsea_GO_DGE_H3204_nosubset <-fgsea_res
fgsea_GO_DGE_H3204 <- subset(fgsea_res[fgsea_res$pval <0.05,])
fgsea_GO_DGE_H3204 <- fgsea_GO_DGE_H3204[order(fgsea_GO_DGE_H3204[,"ES"], decreasing = TRUE),]


#H5471------------------
#too few cells

#HCI009-----------------------
DGE_HCI009 <- read.csv("MAST_HCI009.csv", row.names = 1)
dim(DGE_HCI009)#10892

DGE_HCI009 <- DGE_HCI009[order(DGE_HCI009[,"avg_log2FC"]),]
genelist <- DGE_HCI009$avg_log2FC
names(genelist) <-rownames(DGE_HCI009)

fgsea_res <- fgsea(pathways = m_list, 
                   stats = genelist,
                   minSize=15,
                   maxSize=500,
                   nperm=10000)

#fgsea_H_DGE_HCI009_nosubset <-fgsea_res
#fgsea_H_DGE_HCI009 <- subset(fgsea_res[fgsea_res$pval <0.05,])
#fgsea_H_DGE_HCI009 <- fgsea_H_DGE_HCI009[order(fgsea_H_DGE_HCI009[,"ES"], decreasing = TRUE),]

fgsea_GO_DGE_HCI009_nosubset <-fgsea_res
fgsea_GO_DGE_HCI009 <- subset(fgsea_res[fgsea_res$pval <0.05,])
fgsea_GO_DGE_HCI009 <- fgsea_GO_DGE_HCI009[order(fgsea_GO_DGE_HCI009[,"ES"], decreasing = TRUE),]

#HCI001-----------------------
DGE_HCI001 <- read.csv("MAST_HCI001.csv", row.names = 1)
dim(DGE_HCI001)#10048

DGE_HCI001 <- DGE_HCI001[order(DGE_HCI001[,"avg_log2FC"]),]
genelist <- DGE_HCI001$avg_log2FC
names(genelist) <-rownames(DGE_HCI001)

fgsea_res <- fgsea(pathways = m_list, 
                   stats = genelist,
                   minSize=15,
                   maxSize=500,
                   nperm=10000)

#fgsea_H_DGE_HCI001_nosubset <-fgsea_res
#fgsea_H_DGE_HCI001 <- subset(fgsea_res[fgsea_res$pval <0.05,])
#fgsea_H_DGE_HCI001 <- fgsea_H_DGE_HCI001[order(fgsea_H_DGE_HCI001[,"ES"], decreasing = TRUE),]

fgsea_GO_DGE_HCI001_nosubset <-fgsea_res
fgsea_GO_DGE_HCI001 <- subset(fgsea_res[fgsea_res$pval <0.05,])
fgsea_GO_DGE_HCI001 <- fgsea_GO_DGE_HCI001[order(fgsea_GO_DGE_HCI001[,"ES"], decreasing = TRUE),]

#HCI011-----------------------
DGE_HCI011 <- read.csv("MAST_HCI011.csv", row.names = 1)
dim(DGE_HCI011)#10158

DGE_HCI011 <- DGE_HCI011[order(DGE_HCI011[,"avg_log2FC"]),]
genelist <- DGE_HCI011$avg_log2FC
names(genelist) <-rownames(DGE_HCI011)

fgsea_res <- fgsea(pathways = m_list, 
                   stats = genelist,
                   minSize=15,
                   maxSize=500,
                   nperm=10000)

#fgsea_H_DGE_HCI011_nosubset <-fgsea_res
#fgsea_H_DGE_HCI011 <- subset(fgsea_res[fgsea_res$pval <0.05,])
#fgsea_H_DGE_HCI011 <- fgsea_H_DGE_HCI011[order(fgsea_H_DGE_HCI011[,"ES"], decreasing = TRUE),]

fgsea_GO_DGE_HCI011_nosubset <-fgsea_res
fgsea_GO_DGE_HCI011 <- subset(fgsea_res[fgsea_res$pval <0.05,])
fgsea_GO_DGE_HCI011 <- fgsea_GO_DGE_HCI011[order(fgsea_GO_DGE_HCI011[,"ES"], decreasing = TRUE),]


#H5097-----------------------
DGE_H5097 <- read.csv("MAST_H5097.csv", row.names = 1)
dim(DGE_H5097)#9861

DGE_H5097 <- DGE_H5097[order(DGE_H5097[,"avg_log2FC"]),]
genelist <- DGE_H5097$avg_log2FC
names(genelist) <-rownames(DGE_H5097)
#tail(genelist)
#length(genelist)
#head(genelist)

fgsea_res <- fgsea(pathways = m_list, 
                   stats = genelist,
                   minSize=15,
                   maxSize=500,
                   nperm=10000)

#fgsea_H_DGE_H5097_nosubset <-fgsea_res
#fgsea_H_DGE_H5097 <- subset(fgsea_res[fgsea_res$pval <0.05,])
#fgsea_H_DGE_H5097 <- fgsea_H_DGE_H5097[order(fgsea_H_DGE_H5097[,"ES"], decreasing = TRUE),]

fgsea_GO_DGE_H5097_nosubset <-fgsea_res
fgsea_GO_DGE_H5097 <- subset(fgsea_res[fgsea_res$pval <0.05,])
fgsea_GO_DGE_H5097 <- fgsea_GO_DGE_H5097[order(fgsea_GO_DGE_H5097[,"ES"], decreasing = TRUE),]

#J2036-----------------------
DGE_J2036 <- read.csv("MAST_J2036.csv", row.names = 1)
dim(DGE_J2036)#11484

DGE_J2036 <- DGE_J2036[order(DGE_J2036[,"avg_log2FC"]),]
genelist <- DGE_J2036$avg_log2FC
names(genelist) <-rownames(DGE_J2036)
#tail(genelist)
#length(genelist)
#head(genelist)

fgsea_res <- fgsea(pathways = m_list, 
                   stats = genelist,
                   minSize=15,
                   maxSize=500,
                   nperm=10000)

#fgsea_H_DGE_J2036_nosubset <-fgsea_res
#fgsea_H_DGE_J2036 <- subset(fgsea_res[fgsea_res$pval <0.05,])
#fgsea_H_DGE_J2036 <- fgsea_H_DGE_J2036[order(fgsea_H_DGE_J2036[,"ES"], decreasing = TRUE),]

fgsea_GO_DGE_J2036_nosubset <-fgsea_res
fgsea_GO_DGE_J2036 <- subset(fgsea_res[fgsea_res$pval <0.05,])
fgsea_GO_DGE_J2036 <- fgsea_GO_DGE_J2036[order(fgsea_GO_DGE_J2036[,"ES"], decreasing = TRUE),]

#J53353-----------------------
DGE_J53353 <- read.csv("MAST_J53353.csv", row.names = 1)
dim(DGE_J53353)#12950
head(DGE_J53353)

DGE_J53353 <- DGE_J53353[order(DGE_J53353[,"avg_log2FC"]),]
genelist <- DGE_J53353$avg_log2FC
names(genelist) <-rownames(DGE_J53353)
#tail(genelist)
#length(genelist)
#head(genelist)

fgsea_res <- fgsea(pathways = m_list, 
                   stats = genelist,
                   minSize=15,
                   maxSize=500,
                   nperm=10000)

#fgsea_H_DGE_J53353_nosubset <-fgsea_res
#fgsea_H_DGE_J53353 <- subset(fgsea_res[fgsea_res$pval <0.05,])
#fgsea_H_DGE_J53353 <- fgsea_H_DGE_J53353[order(fgsea_H_DGE_J53353[,"ES"], decreasing = TRUE),]

fgsea_GO_DGE_J53353_nosubset <-fgsea_res
fgsea_GO_DGE_J53353 <- subset(fgsea_res[fgsea_res$pval <0.05,])
fgsea_GO_DGE_J53353 <- fgsea_GO_DGE_J53353[order(fgsea_GO_DGE_J53353[,"ES"], decreasing = TRUE),]


#HCI010-----------------------
DGE_HCI010 <- read.csv("MAST_HCI010.csv", row.names = 1)
dim(DGE_HCI010)#14449

DGE_HCI010 <- DGE_HCI010[order(DGE_HCI010[,"avg_log2FC"]),]
genelist <- DGE_HCI010$avg_log2FC
names(genelist) <-rownames(DGE_HCI010)
#tail(genelist)
#length(genelist)
#head(genelist)

fgsea_res <- fgsea(pathways = m_list, 
                   stats = genelist,
                   minSize=15,
                   maxSize=500,
                   nperm=10000)

#fgsea_H_DGE_HCI010_nosubset <-fgsea_res
#fgsea_H_DGE_HCI010 <- subset(fgsea_res[fgsea_res$pval <0.05,])
#fgsea_H_DGE_HCI010 <- fgsea_H_DGE_HCI010[order(fgsea_H_DGE_HCI010[,"ES"], decreasing = TRUE),]

fgsea_GO_DGE_HCI010_nosubset <-fgsea_res
fgsea_GO_DGE_HCI010 <- subset(fgsea_res[fgsea_res$pval <0.05,])
fgsea_GO_DGE_HCI010 <- fgsea_GO_DGE_HCI010[order(fgsea_GO_DGE_HCI010[,"ES"], decreasing = TRUE),]




#combine enrichment analysis of individual tumors  ---------
#take unfiltered enrichment analysis
#HALLMARK
fgsea_H_DGE_global_nosubset$Tumor <-"global"
fgsea_H_DGE_H3204_nosubset$Tumor <-"H3204"
fgsea_H_DGE_H4272_nosubset$Tumor <-"H4272"
fgsea_H_DGE_H5097_nosubset$Tumor <-"H5097"
fgsea_H_DGE_HCI001_nosubset$Tumor <-"HCI001"
fgsea_H_DGE_HCI010_nosubset$Tumor <-"HCI010"
fgsea_H_DGE_HCI005_nosubset$Tumor <-"HCI005"
fgsea_H_DGE_HCI009_nosubset$Tumor <-"HCI009"
fgsea_H_DGE_HCI011_nosubset$Tumor <-"HCI011"
fgsea_H_DGE_J2036_nosubset$Tumor <-"J2036"
fgsea_H_DGE_J53353_nosubset$Tumor <-"J53353"

fgsea_H_DGE_all_nosubset <- rbind(fgsea_H_DGE_global_nosubset, fgsea_H_DGE_H3204_nosubset, fgsea_H_DGE_H4272_nosubset, fgsea_H_DGE_H5097_nosubset,
                                  fgsea_H_DGE_HCI001_nosubset, fgsea_H_DGE_HCI005_nosubset, fgsea_H_DGE_HCI009_nosubset, fgsea_H_DGE_HCI010_nosubset,
                                  fgsea_H_DGE_HCI011_nosubset, fgsea_H_DGE_J2036_nosubset, fgsea_H_DGE_J53353_nosubset)


fgsea_H_DGE_all_nosubset$leadingEdge <- vapply(fgsea_H_DGE_all_nosubset$leadingEdge, paste, collapse = ", ", character(1L))


#subset on shared pathways between tumor models-----------
#subset on shared between more than 4 tumor models based on significant p-value
Hsubset_met <- c("HALLMARK_MYC_TARGETS_V1",
                 "HALLMARK_E2F_TARGETS",
                 "HALLMARK_G2M_CHECKPOINT",
                 "HALLMARK_APICAL_JUNCTION")

#subset on shared between more than 4 tumor models based on significant p-value
Hsubset_tumor <- c("HALLMARK_GLYCOLYSIS",
                   "HALLMARK_HYPOXIA",
                   "HALLMARK_MTORC1_SIGNALING",
                   "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                   "HALLMARK_ESTROGEN_RESPONSE_LATE",
                   "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                   "HALLMARK_UV_RESPONSE_UP")
shared_H<- unique(c(Hsubset_met, Hsubset_tumor))


#subset on shared 
test<-fgsea_H_DGE_all_nosubset[which(fgsea_H_DGE_all_nosubset$pathway %in% shared_H),]
#remove "HALLMARK" from pathway
test$pathway <- gsub("HALLMARK_", "", test$pathway)
test$pathway <- factor(test$pathway, levels=c(gsub("HALLMARK_", "", Hsubset_tumor), 
                                              gsub("HALLMARK_", "", Hsubset_met[!Hsubset_met %in% Hsubset_tumor])))


##Figure 2E Visualization bubble plot of shared pathways ----
#using test as data input
test$Tumor <- factor(test$Tumor, levels=c("global",
                                          "HCI005", 
                                          "H3204", 
                                          "H4272",
                                          "HCI009",  
                                          "HCI011", 
                                          "HCI001", 
                                          "H5097", 
                                          "J2036", 
                                          "J53353",
                                          "HCI010"))

p <-ggplot(data=test, aes(x=Tumor,y=pathway, size = -log10(pval), color=NES)) +
  geom_point() +
  geom_point(data=test[which(test$pval <= 0.05),],
             aes(x=Tumor,y=pathway), shape=21, color="black", stroke=0.5) +
  scale_color_gradient2(low="#ff7f00", high="#377eb8", midpoint=0, 
                        name="NES", limits = c(-2.6,2.6)) +
  scale_size(range = c(0.01, 5), name="p-value",labels=c("0.1", "0.01","0.001"), breaks= c(1,2,3))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.1),
        # panel.spacing.x = unit(0.1, "cm"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        # axis.text.x = element_text(size=10, color="black", angle=90, hjust=1),
        axis.text.x =element_blank(), 
        axis.text.y = element_text(size=11, color="black"),
        legend.key = element_blank(), 
        legend.position = "none",
        legend.title = element_text(size=6, color="black"),
        legend.text = element_text(size=6, color="black"),
        axis.ticks = element_line(colour = "black", size = 0.1, linetype = "solid"),
        axis.line = element_line(colour = "black", size = 0.1, linetype = "solid"))
