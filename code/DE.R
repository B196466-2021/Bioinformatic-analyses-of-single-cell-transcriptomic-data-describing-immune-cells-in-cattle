#! /usr/bin/Rscript
DefaultAssay(compare) <- "RNA"
compare$label<-as.character(compare$integrated_snn_res.0.8)

# Define dendritic cells group 1 
compare$label<-ifelse(compare$label %in% c(1,3,4,6,13),"DC_g1",compare$label)
compare$celltype<-paste(compare$label,compare$type,sep="_")
Idents(compare)<-"celltype"
# Find differentially expressed genes
DC_g1.response<- FindMarkers(compare, ident.1 = "DC_g1_Post_BCG", ident.2 = "DC_g1_Pre_BCG", verbose = FALSE)
head(DC_g1.response, n = 15)
write.csv(DC_g1.response,file="table1.csv")

# Volcano plot
png("volcano.png")
EnhancedVolcano(DC_g1.response,lab = rownames(DC_g1.response),x = "avg_log2FC",y = "p_val_adj",pCutoff = 0.05,FCcutoff = 0.5)
dev.off

# Define dendritic cells group 1 subset
Idents(compare)<-compare$label
DC1<- subset(compare, idents=c("DC_g1"))
#FeaturePlot
FeaturePlot<-FeaturePlot(DC1,features = c("RNF19B"), split.by = "type", max.cutoff ="q90", cols = c("grey", "red"))
FeaturePlot
ggsave("FeaturePlot.png",FeaturePlot)

plots <- VlnPlot(DC1, features = c("RNF19B"), group.by = "type",pt.size = 0.01, combine = FALSE,assay = "RNA",log=TRUE)
wrap_plots(plots = plots, ncol = 1)

plots <- VlnPlot(DC, features = c("RNF19B"), split.by = "type",pt.size = 0.01, combine = FALSE,assay = "RNA",log=TRUE)
wrap_plots(plots = plots, ncol = 1)

# Define dendritic cells group 2 
compare$label<-ifelse(compare$label %in% c(8,16),"DC_g2",compare$label)
compare$celltype<-paste(compare$label,compare$type,sep="_")
Idents(compare)<-"celltype"
# Find differentially expressed genes
DC_g2.response<- FindMarkers(compare, ident.1 = "DC_g2_Post_BCG", ident.2 = "DC_g2_Pre_BCG", verbose = FALSE)
head(DC_g2.response, n = 15)
write.csv(DC_g2.response,file="table2.csv")

# Volcano plot
png("volcano_g2.png")
EnhancedVolcano(DC_g2.response,lab = rownames(DC_g2.response),x = "avg_log2FC",y = "p_val_adj",pCutoff = 0.05,FCcutoff = 0.5)
dev.off

# Define dendritic cells group 3 
compare$label<-ifelse(compare$label %in% c(2),"DC_g3",compare$label)
compare$celltype<-paste(compare$label,compare$type,sep="_")
Idents(compare)<-"celltype"
# Find differentially expressed genes
DC_g3.response<- FindMarkers(compare, ident.1 = "DC_g3_Post_BCG", ident.2 = "DC_g3_Pre_BCG", verbose = FALSE)
head(DC_g3.response, n = 15)
write.csv(DC_g3.response,file="table3.csv")

# Volcano plot
png("volcano_g3.png")
EnhancedVolcano(DC_g3.response,lab = rownames(DC_g3.response),x = "avg_log2FC",y = "p_val_adj",pCutoff = 0.05,FCcutoff = 0.5)
dev.off

# Define dendritic cells group 3 subset
Idents(compare)<-compare$label
DC3<- subset(compare, idents=c("DC_g3"))
# FeaturePlot
FeaturePlot<-FeaturePlot(DC3,features = c("RBBP6"), split.by = "type", max.cutoff ="q90", cols = c("grey", "red"))
FeaturePlot
ggsave("FeaturePlot.png",FeaturePlot)

plots <- VlnPlot(DC3, features = c("RBBP6"), group.by = "type",pt.size = 0.01, combine = FALSE,assay = "RNA",log=TRUE)
wrap_plots(plots = plots, ncol = 1)

plots <- VlnPlot(DC, features = c("RBBP6"), split.by = "type",pt.size = 0.01, combine = FALSE,assay = "RNA",log=TRUE)
wrap_plots(plots = plots, ncol = 1)

# Define dendritic cells group 4 
compare$label<-ifelse(compare$label %in% c(0,10,11),"DC_g4",compare$label)
compare$celltype<-paste(compare$label,compare$type,sep="_")
Idents(compare)<-"celltype"
# Find differentially expressed genes
DC_g4.response<- FindMarkers(compare, ident.1 = "DC_g4_Post_BCG", ident.2 = "DC_g4_Pre_BCG", verbose = FALSE)
head(DC_g4.response, n = 15)
write.csv(DC_g4.response,file="table4.csv")

# Volcano Plot
png("volcano_g4.png")
EnhancedVolcano(DC_g4.response,lab = rownames(DC_g4.response),x = "avg_log2FC",y = "p_val_adj",pCutoff = 0.05,FCcutoff = 0.5)
dev.off

# Define dendritic cells group 4 subset
Idents(compare)<-compare$label
DC4<- subset(compare, idents=c("DC_g4"))
#FeaturePlot
FeaturePlot<-FeaturePlot(DC4,features = c("SRSF5"), split.by = "type", max.cutoff ="q90", cols = c("grey", "red"))
FeaturePlot
ggsave("FeaturePlot.png",FeaturePlot)

plots <- VlnPlot(DC4, features = c("SRSF5"), group.by = "type",pt.size = 0.01, combine = FALSE,assay = "RNA",log=TRUE)
wrap_plots(plots = plots, ncol = 1)

plots <- VlnPlot(DC, features = c("SRSF5"), split.by = "type",pt.size = 0.01, combine = FALSE,assay = "RNA",log=TRUE)
wrap_plots(plots = plots, ncol = 1)

# Define dendritic cells group 5 
compare$label<-ifelse(compare$label %in% c(18),"DC_g5",compare$label)
compare$celltype<-paste(compare$label,compare$type,sep="_")
Idents(compare)<-"celltype"
# Find differentially expressed genes
DC_g5.response<- FindMarkers(compare, ident.1 = "DC_g5_Post_BCG", ident.2 = "DC_g5_Pre_BCG", verbose = FALSE)
head(DC_g5.response, n = 15)
write.csv(DC_g5.response,file="table5.csv")

# Volcano Plot
png("volcano_g5.png")
EnhancedVolcano(DC_g5.response,lab = rownames(DC_g5.response),x = "avg_log2FC",y = "p_val_adj",pCutoff = 0.05,FCcutoff = 0.5)
dev.off

# Define macrophage group
compare$label<-ifelse(compare$label %in% c(17),"Macrophage",compare$label)
compare$celltype<-paste(compare$label,compare$type,sep="_")
Idents(compare)<-"celltype"
# Find differentially expressed genes
Macrophage.compare.response<- FindMarkers(compare, ident.1 = "Macrophage_Post_BCG", ident.2 = "Macrophage_Pre_BCG", verbose = FALSE)
head(Macrophage.compare.response, n = 15)
write.csv(Macrophage.compare.response,file="Macrophage.csv")

# Volcano Plot
png("Macrophage.png")
EnhancedVolcano(Macrophage.compare.response,lab = rownames(Macrophage.compare.response),x = "avg_log2FC",y = "p_val_adj",pCutoff = 0.05,FCcutoff = 0.5)
dev.off

# Define macrophage group subset
Idents(compare)<-compare$label
macrophage<- subset(compare, idents=c("Macrophage"))
#FeaturePlot
FeaturePlot<-FeaturePlot(macrophage, features = c("ARHGDIB"), split.by = "type", max.cutoff = "q90",cols = c("grey", "red"))
FeaturePlot
ggsave("FeaturePlot.png",FeaturePlot)

plots <- VlnPlot(macrophage, features = c("ARHGDIB"), group.by = "type",pt.size = 0.01,combine = FALSE,assay = "RNA",log=TRUE)
wrap_plots(plots = plots, ncol = 1)


# Define the whole dendritic cells subset
Idents(compare)<-compare$label
DC<- subset(compare, idents=c("DC_g1","DC_g2","DC_g3","DC_g4","DC_g5"))
# Find all markers 
all.markers <- FindAllMarkers(object = DC)

# find the top 20 genes
all.markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top20
# make heatmap
DefaultAssay(DC) <- "integrated"
Heatmap<-DoHeatmap(DC, features = top20$gene) 
ggsave("Heatmap.png",Heatmap,height=10)
# create table
top20
write.csv(top20,file="top20.csv")
