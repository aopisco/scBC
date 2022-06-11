#This script includes analysis and visualization of scRNAseq data published in Winkler et al. 2022

library(Seurat)
library(ggplot2)
library(viridis)

#custom colors
col_tumor <- c("HCI010"="#1598cc","J53353"="#80b1d3","J2036"="#8dd3c7",
               "HCI001"="#fdb462","HCI009"="#fb8072","H5471"="#b3de69",
               "H4272"="#ccebc5","H3204"="#ffed6f","J55454"="#e8ee88",
               "HCI002"="#d9d9d9","HCI011"="#bebada","H5097"="#bc80bd",
               "HCI005"="#fccde5") 
col_EMT <- c("EMT Intermediate"='#BC80B7','Epithelial-like' = "#81B1D3", "Mesenchymal-like"="#F47F72")

#Figure 4E----------------

#EMT stage line plot

p <-ggplot(seu_ss2@meta.data, aes(x=reorder(rownames(seu_ss2@meta.data), EMT_score), y=EMT_score, color=EMT_stage)) + 
  geom_point(size=0.1)+
  geom_rug(sides="l", alpha=0.03)+
  theme_minimal() +
  scale_color_manual(values = col_EMT)+
  geom_hline(yintercept=c(-0.2,0.2), size=0.25, linetype="dashed")+
  theme(panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.25),
        axis.line  = element_line(size = 0.25))

#Figure 4F--------------
#stacked barplot
data.plot<- as.data.frame(table(seu_ss2$EMT_stage, seu_ss2$Tumor_ID))
colnames(data.plot) <- c("EMT", "Tumor_model", "Cellcount")
#order by mesenchymal proportions
data.plot$Tumor <- factor(data.plot$Tumor, levels=c("HCI005","HCI011","HCI009","H5097", "H4272" ,"H5471","H3204",
                                            "J2036","HCI001","HCI010", "J55454","J53353"))
data.plot$EMT <- factor(data.plot$EMT, levels=c( "Epithelial-like", "EMT Intermediate", "Mesenchymal-like"))

p <-ggplot(
  data=data.plot, 
  aes(x=Tumor_model, y = Cellcount, fill=EMT)) +
  geom_bar(stat="identity", position = position_fill((vjust = 0.5))) +
  labs(y="Cell fraction") +
  scale_fill_manual(values = col_EMT)+
  theme_minimal()+
  theme(plot.title = element_blank(),panel.grid=element_blank(),
        axis.line = element_line(colour = "black", size = 0.25, linetype = "solid"),
        axis.ticks.y = element_line(colour = "black", size = 0.25, linetype = "solid"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

#Figure 4G--------------

#do this for Smart-Seq2 data and MULTI-Seq data (Supplementary Figure 4I)

#Smart-Seq2 data:
data_TF <- cbind(seu_ss2@meta.data, 
                 t(as.data.frame(seu_ss2@assays$RNA@data[rownames(seu_ss2@assays$RNA@data) %in% c("SNAI1","SNAI2","TWIST1","ZEB1","ZEB2"),])))

data_TF$EMT_stage <- factor(data_TF$EMT_stage, levels=c( "Epithelial-like", "EMT Intermediate", "Mesenchymal-like"))

#Visualize TF expressing positive cells per EMP cell stage
#Violin plot
#do this for c("SNAI1","SNAI2","TWIST1","ZEB1","ZEB2")
p <-ggplot(data_TF[data_TF$TWIST1 >0,], aes(x=EMT_stage, y=TWIST1, fill=EMT_stage)) + 
  geom_violin(color = NA) +
  scale_fill_manual(values =col_EMT)+
  geom_jitter(shape=16, position=position_jitter(0.1), size=0.2)+
  theme_minimal() +
  theme(panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_line(size = 0.25),
        legend.position="none",
        axis.line  = element_line(size = 0.25))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

#bar diagram of fraction positive expressing cells per EMP cell stage
#do this for c("SNAI1","SNAI2","TWIST1","ZEB1","ZEB2")
TWIST1_pos_E <- c(sum(data_TF$TWIST1[data_TF$EMT_stage == "Epithelial-like"] >0),
                  sum(data_TF$TWIST1[data_TF$EMT_stage == "Epithelial-like"] == 0))
TWIST1_pos_M <- c(sum(data_TF$TWIST1[data_TF$EMT_stage == "Mesenchymal-like"] >0),
                  sum(data_TF$TWIST1[data_TF$EMT_stage == "Mesenchymal-like"] ==0))
TWIST1_pos_I <- c(sum(data_TF$TWIST1[data_TF$EMT_stage == "EMT Intermediate"] >0),
                  sum(data_TF$TWIST1[data_TF$EMT_stage == "EMT Intermediate"] ==0))
TWIST_pos <- c(TWIST1_pos_E, TWIST1_pos_M, TWIST1_pos_I)
pos_neg <- c("pos", "neg","pos", "neg","pos", "neg")
EMT_stage<- c("Epithelial-like", "Epithelial-like",
              "Mesenchymal-like","Mesenchymal-like",
              "EMT Intermediate","EMT Intermediate")
counts_pos <- data.frame(TWIST_pos, pos_neg, EMT_stage)
counts_pos$EMT_stage <- factor(counts_pos$EMT_stage, levels=c( "Epithelial-like", "EMT Intermediate", "Mesenchymal-like"))

p <-ggplot(counts_pos, aes(x=EMT_stage, y=TWIST_pos, fill=pos_neg, color=EMT_stage)) + 
  geom_bar(stat = "identity", position = "fill", size=0.5) +
  scale_fill_manual(values=c("pos"="Grey", "neg"="White"))+
  scale_color_manual(values=col_EMT)+
  theme_minimal() +
  theme(panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_line(size = 0.25),
        legend.title = element_blank(),
        legend.position="none",
        axis.line  = element_line(size = 0.25))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

#Supplementary Figure 4A--------------------
#expression of EMT markers Smart-Seq2
data.plot <- cbind(seu_ss2@meta.data, 
              t(as.data.frame(seu_ss2@assays$RNA@data[rownames(seu_ss2@assays$RNA@data) %in% c("EPCAM", "CDH1","VIM", "CDH2", "FN1"),])))

#VIM
p <-ggplot(data.plot, aes(x=E_score, y=M_score, color=VIM)) + 
  geom_point(size=0.1) +
  scale_color_continuous(type = "viridis")+
  labs(title = "VIM")+
  theme_minimal() +
  theme(panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_line(size = 0.25),
        axis.line  = element_line(size = 0.25))

#do this for c("EPCAM", "CDH1","VIM", "CDH2", "FN1")
#do the same using MULTI-Seq data

#Supplementary Figure 4B-----------------------------
#Calculate mean EMP score per tumor model and tissue 
mean_EMP_score<-aggregate(seu_ss2@meta.data$EMT_score, list(seu_ss2@meta.data$sort, seu_ss2@meta.data$Tumor_ID), FUN=mean) 
colnames(mean_EMP_score) <-c("Tissue", "Tumor_model", "EMP_score")

#Calculate Pearson correlation
cortest <- cor.test(mean_EMP_score[grep("Tumor", mean_EMP_score$Tissue), "EMP_score"], 
                    mean_EMP_score[grep("Metastatic", mean_EMP_score$Tissue), "EMP_score"]) #Pearson R=0.8777764,  p-value = 0.0001744
cortest$estimate^2 #0.7704914 

#visualize correlation in scatter plot
data.plot <- cbind(mean_EMP_score[grep("Tumor", mean_EMP_score$Tissue),],
                   mean_EMP_score[grep("Metastatic", mean_EMP_score$Tissue),])
colnames(data.plot) <- c("Tumor", "Tumor_model", "Primary_tumor_mean", "Metastatis", "Tumor_model2","Metastasis_mean")

p <-ggplot(data.plot, aes(x=Primary_tumor_mean, y=Metastasis_mean, color=Tumor_model)) + 
  geom_smooth(method=lm , color="black", fill="lightgrey", se=TRUE,size = 0.25)+
  geom_point(size=1) +
  theme_minimal() +
  labs(title="EMP signature")+
  scale_color_manual(values=col_tumor)+
  theme(panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.ticks =element_line(size = 0.25),
        axis.line  = element_line(size = 0.25))




#Supplementary Figure 4G----------------

#EMT stage line plot
p <-ggplot(seu_10Xv2@meta.data, aes(x=reorder(rownames(seu_10Xv2@meta.data), EMT_score), y=EMT_score, color=EMT_stage)) + 
  geom_point(size=0.1)+
  geom_rug(sides="l", alpha=0.03)+
  theme_minimal() +
  scale_color_manual(values = col_EMT)+
  geom_hline(yintercept=c(-0.2,0.2), size=0.25, linetype="dashed")+
  theme(panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.25),
        axis.line  = element_line(size = 0.25))

#Supplementary Figure 4I---------------
#do the same like in Figure 4G for MULTI-Seq data
#Figure 4G--------------

#do this for Smart-Seq2 data and MULTI-Seq data (Supplementary Figure 4I)

#Smart-Seq2 data:
data_TF <- cbind(seu_ss2@meta.data, 
                   t(as.data.frame(seu_ss2@assays$RNA@data[rownames(seu_ss2@assays$RNA@data) %in% c("SNAI1","SNAI2","TWIST1","ZEB1","ZEB2"),])))

data_TF$EMT_stage <- factor(data_TF$EMT_stage, levels=c( "Epithelial-like", "EMT Intermediate", "Mesenchymal-like"))

#Visualize TF expressing positive cells per EMP cell stage
#Violin plot
#do this for c("SNAI1","SNAI2","TWIST1","ZEB1","ZEB2")
p <-ggplot(data_TF[data_TF$TWIST1 >0,], aes(x=EMT_stage, y=TWIST1, fill=EMT_stage)) + 
  geom_violin(color = NA) +
  scale_fill_manual(values =col_EMT)+
  geom_jitter(shape=16, position=position_jitter(0.1), size=0.2)+
  theme_minimal() +
  theme(panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_line(size = 0.25),
        legend.position="none",
        axis.line  = element_line(size = 0.25))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

#bar diagram of fraction positive expressing cells per EMP cell stage
#do this for c("SNAI1","SNAI2","TWIST1","ZEB1","ZEB2")
TWIST1_pos_E <- c(sum(data_TF$TWIST1[data_TF$EMT_stage == "Epithelial-like"] >0),
                  sum(data_TF$TWIST1[data_TF$EMT_stage == "Epithelial-like"] == 0))
TWIST1_pos_M <- c(sum(data_TF$TWIST1[data_TF$EMT_stage == "Mesenchymal-like"] >0),
                  sum(data_TF$TWIST1[data_TF$EMT_stage == "Mesenchymal-like"] ==0))
TWIST1_pos_I <- c(sum(data_TF$TWIST1[data_TF$EMT_stage == "EMT Intermediate"] >0),
                  sum(data_TF$TWIST1[data_TF$EMT_stage == "EMT Intermediate"] ==0))
TWIST_pos <- c(TWIST1_pos_E, TWIST1_pos_M, TWIST1_pos_I)
pos_neg <- c("pos", "neg","pos", "neg","pos", "neg")
EMT_stage<- c("Epithelial-like", "Epithelial-like",
              "Mesenchymal-like","Mesenchymal-like",
              "EMT Intermediate","EMT Intermediate")
counts_pos <- data.frame(TWIST_pos, pos_neg, EMT_stage)
counts_pos$EMT_stage <- factor(counts_pos$EMT_stage, levels=c( "Epithelial-like", "EMT Intermediate", "Mesenchymal-like"))

p <-ggplot(counts_pos, aes(x=EMT_stage, y=TWIST_pos, fill=pos_neg, color=EMT_stage)) + 
  geom_bar(stat = "identity", position = "fill", size=0.5) +
  scale_fill_manual(values=c("pos"="Grey", "neg"="White"))+
  scale_color_manual(values=col_EMT)+
  theme_minimal() +
  theme(panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_line(size = 0.25),
        legend.title = element_blank(),
        legend.position="none",
        axis.line  = element_line(size = 0.25))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
