#! /usr/bin/Rscript

# Cell type annotation using SingleR
monaco.ref <- celldex::MonacoImmuneData()
compare.sce <- as.SingleCellExperiment(compare)
monaco.main <- SingleR(test = compare.sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
compare@meta.data$monaco.main <- monaco.main$pruned.labels
compare.annotation <- SetIdent(compare, value = "monaco.main")
DimPlot(compare.annotation, label = T , repel = T, label.size = 3)



# Show dendritic cells groups location
FeaturePlot<-FeaturePlot(compare, features = c("nFeature_RNA","nCount_RNA"),ncol = 4,pt.size = 0.1,cols = c("grey", "red")) 
ggsave("FeaturePlot1.png",FeaturePlot,width=20)
FeaturePlot<-FeaturePlot(compare, features = c("percent.mt","percent.rb"),ncol = 4,pt.size = 0.1,cols = c("grey", "red")) 
ggsave("FeaturePlot2.png",FeaturePlot,width=20)
cell_cycle<-FeaturePlot(compare, features = c("S.Score", "G2M.Score"),ncol = 4, pt.size = 0.1,cols = c("grey", "red"))
ggsave("cell_cycle.png",cell_cycle,width=20)
LYZ<-FeaturePlot(compare, features = c("LYZ"),cols = c("grey", "red"))
ggsave("LYZ.png",LYZ)
CD1B<-FeaturePlot(compare, features = c("CD1B"),cols = c("grey", "red"))
ggsave("CD1B.png",CD1B)
RBP4<-FeaturePlot(compare, features = c("RBP4"),cols = c("grey", "red"))
ggsave("RBP4.png",RBP4)
DNASE1L3<-FeaturePlot(compare, features = c("DNASE1L3"),cols = c("grey", "red"))
ggsave("DNASE1L3.png",DNASE1L3)
SIRPA<-FeaturePlot(compare, features = c("SIRPA"),cols = c("grey", "red"))
ggsave("SIRPA.png",SIRPA)
EPCAM<-FeaturePlot(compare, features = c("EPCAM"),cols = c("grey", "red"))
ggsave("EPCAM.png",EPCAM)
RUNX2<-FeaturePlot(compare, features = c("RUNX2"),cols = c("grey", "red"))
ggsave("RUNX2.png",RUNX2)
CSF1R<-FeaturePlot(compare, features = c("CSF1R"),cols = c("grey", "red"))
ggsave("CSF1R.png",CSF1R)
S100B<-FeaturePlot(compare, features = c("S100B"),cols = c("grey", "red"))
ggsave("S100B.png",S100B)
CD14<-FeaturePlot(compare, features = c("CD14"),cols = c("grey", "red"))
ggsave("CD14.png",CD14)
