#This script includes analysis and visualization of scRNAseq data published in Winkler et al. 2022

library(Seurat)
library(ggplot2)

col_tumor <- c("HCI010"="#1598cc","J53353"="#80b1d3","J2036"="#8dd3c7",
               "HCI001"="#fdb462","HCI009"="#fb8072","H5471"="#b3de69",
               "H4272"="#ccebc5","H3204"="#ffed6f","J55454"="#e8ee88",
               "HCI002"="#d9d9d9","HCI011"="#bebada","H5097"="#bc80bd",
               "HCI005"="#fccde5") 

#Figure 1-----------------------

#Supplementary Figure 1G-------------------
#extract normalized expression of ESR1, PGR, ERBB2
data_receptor<- cbind(t(as.data.frame(
  seu_ss2@assays$RNA@data[rownames(seu_ss2@assays$RNA@data) %in% c("ESR1", "PGR", "ERBB2"),])),
  seu_ss2@meta.data)

#make average expression of ESR1, PGR, ERBB2 per tumor model and tissue 
meanESR1<-aggregate(data_receptor$ESR1, list(data_receptor$sort, data_receptor$Tumor_ID), FUN=mean) 
meanPGR<-aggregate(data_receptor$PGR, list(data_receptor$sort, data_receptor$Tumor_ID), FUN=mean) 
meanERBB2<-aggregate(data_receptor$ERBB2, list(data_receptor$sort, data_receptor$Tumor_ID), FUN=mean) 

#calculate Pearson correlation
#ESR1
cortest <- cor.test(meanESR1[grep("Tumor", meanESR1$Group.1), "x"], 
                    meanESR1[grep("Metastatic", meanESR1$Group.1), "x"]) #Pearson R=0.9937307  ,  p-value =  7.547e-11
cortest$estimate^2#0.9875008 
#PGR
cortest <- cor.test(meanPGR[grep("Tumor", meanPGR$Group.1), "x"], 
                    meanPGR[grep("Metastatic", meanPGR$Group.1), "x"]) #Pearson R=0.9803243  ,  p-value =  2.247e-08
cortest$estimate^2#0.9610358
#ERBB2
cortest <- cor.test(meanERBB2[grep("Tumor", meanERBB2$Group.1), "x"], 
                    meanERBB2[grep("Metastatic", meanERBB2$Group.1), "x"]) #Pearson R=0.7580199 ,  p-value =  0.00428
cortest$estimate^2#0.5745942

#visualize correlation in scatter plot
#ESR1
data.plot <- cbind(meanESR1[grep("Tumor", meanESR1$Group.1),],
                   meanESR1[grep("Metastatic", meanESR1$Group.1),])
colnames(data.plot) <- c("Tumor", "Tumor_model", "Tumor_mean", "Metastatis", "Tumor_model2","Metastasis_mean")

ggplot(data.plot, aes(x=Tumor_mean, y= Metastasis_mean, color=Tumor_model)) + 
  geom_smooth(method=lm , color="black", fill="lightgrey", se=TRUE,size = 0.25)+
  geom_point(size=1) +
  theme_minimal() +
  labs(title="ESR1")+
  scale_color_manual(values=col_tumor)+
  theme(panel.grid=element_blank(),
        axis.ticks =element_line(size = 0.25),
        axis.line  = element_line(size = 0.25))
#PGR
data.plot <- cbind(meanPGR[grep("Tumor", meanPGR$Group.1),],
                   meanPGR[grep("Metastatic", meanPGR$Group.1),])
colnames(data.plot) <- c("Tumor", "Tumor_model", "Tumor_mean", "Metastatis", "Tumor_model2","Metastasis_mean")

ggplot(data.plot, aes(x=Tumor_mean, y= Metastasis_mean, color=Tumor_model)) + 
  geom_smooth(method=lm , color="black", fill="lightgrey", se=TRUE,size = 0.25)+
  geom_point(size=1) +
  theme_minimal() +
  labs(title="PGR")+
  scale_color_manual(values=col_tumor)+
  theme(panel.grid=element_blank(),
        axis.ticks =element_line(size = 0.25),
        axis.line  = element_line(size = 0.25))

#ERBB2
data.plot <- cbind(meanERBB2[grep("Tumor", meanERBB2$Group.1),],
                   meanERBB2[grep("Metastatic", meanERBB2$Group.1),])
colnames(data.plot) <- c("Tumor", "Tumor_model", "Tumor_mean", "Metastatis", "Tumor_model2","Metastasis_mean")

ggplot(data.plot, aes(x=Tumor_mean, y= Metastasis_mean, color=Tumor_model)) + 
  geom_smooth(method=lm , color="black", fill="lightgrey", se=TRUE,size = 0.25)+
  geom_point(size=1) +
  theme_minimal() +
  labs(title="ERBB2")+
  scale_color_manual(values=col_tumor)+
  theme(panel.grid=element_blank(),
        axis.ticks =element_line(size = 0.25),
        axis.line  = element_line(size = 0.25))

#Supplementary Figure 1I-------------------
DimPlot(seu_ss2, group.by = "plate_ID", reduction = "umap", pt.size=0.1)

#Supplementary Figure 1J-------------------
DimPlot(seu_ss2, group.by = "Animal_ID", reduction = "umap", pt.size=0.1)
 



