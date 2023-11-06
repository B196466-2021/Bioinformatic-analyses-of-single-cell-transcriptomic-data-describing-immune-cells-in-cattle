#! /usr/bin/Rscript
#########################
# Load libraries
#########################
library(Seurat)
library(Matrix)
library(SeuratData)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scater)
library(SeuratDisk)
library(SingleCellExperiment)
library(scDblFinder)
library(BiocManager)
library(celldex)
library(SingleR)
library(cowplot)
library(DoubletFinder)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
library(stringr)
library(enrichplot)
library(topGO)
library(GSVA)
library(pathview)
library(DOSE)
library(GSEABase)
library(EnhancedVolcano)
library(enrichR)
library(garnett)
library(rlang)

# Define the data and project information
data_paths <- c(
  "~/Desktop/11858HJPool01-N__11858HJ0001_cellranger_filtered_feature_bc_matrix",
  "~/Desktop/11858HJPool01-N__11858HJ0002_cellranger_filtered_feature_bc_matrix",
  "~/Desktop/11858HJPool01-N__11858HJ0003_cellranger_filtered_feature_bc_matrix",
  "~/Desktop/11858HJPool01-N__11858HJ0004_cellranger_filtered_feature_bc_matrix",
  "~/Desktop/11858HJPool02-N__11858HJ0005_cellranger_filtered_feature_bc_matrix",
  "~/Desktop/11858HJPool02-N__11858HJ0006_cellranger_filtered_feature_bc_matrix",
  "~/Desktop/11858HJPool02-N__11858HJ0007_cellranger_filtered_feature_bc_matrix",
  "~/Desktop/11858HJPool02-N__11858HJ0008_cellranger_filtered_feature_bc_matrix"
)

project_names <- c(
  "602791_Pre_BCG_DC",
  "602791_Post_BCG_DC",
  "303215_Pre_BCG_DC",
  "303215_Post_BCG_DC",
  "402740_Pre_BCG_DC",
  "402740_Post_BCG_DC",
  "Pre_BCG_Lymph_node",
  "Post_BCG_Lymph_node"
)

# Load data and create Seurat objects using a loop
seurat_objects <- list()

for (i in 1:length(data_paths)) {
  data <- Read10X(data.dir = data_paths[i])
  result <- CreateSeuratObject(counts = data, project = project_names[i])
  result$type <- project_names[i]  # Set type
  seurat_objects[[i]] <- result
}


# Merge dataset into one single seurat object
pre<-merge(result1, c(result3,result5,result7), add.cell.ids = c("602791_Pre_BCG_DC", "303215_Pre_BCG_DC", "402740_Pre_BCG_DC","Pre_BCG_Lymph_node"))
post<-merge(result2, c(result4,result6,result8), add.cell.ids = c("602791_Post_BCG_DC","303215_Post_BCG_DC","402740_Post_BCG_DC","Post_BCG_Lymph_node"))
pre$type <- "Pre_BCG"
post$type <- "Post_BCG"
compare<-merge(pre, c(post), add.cell.ids = c("Pre_BCG", "Post_BCG"))
compare

# Clear unused variables
rm(data1,result1,data2,result2,data3,result3,data4,result4,data5,result5,data6,result6,data7,result7,data8,result8,pre,post)
gc()


##############################
## doubletFinder 
##############################
compare = ScaleData(compare,verbose =F)
compare = RunPCA(compare, verbose = F, npcs = 20)
compare = RunUMAP(compare, dims = 1:20, reduction = "pca",verbose = F)
compare = FindNeighbors(compare, dims = 1:20,verbose = FALSE)
compare = FindClusters(compare, resolution = 0.8, verbose = FALSE)

## Find the optimal pK value
sweep.res.list <- paramSweep_v3(compare, PCs = 1:20, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## Exclude homologous doublets that cannot be detected and optimize the expected number of doublets
DoubletRate = ncol(compare)*8*1e-6                     
homotypic.prop <- modelHomotypic(compare$seurat_clusters) 
nExp_poi <- round(DoubletRate*ncol(compare)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Identify doublets using defined parameters
compare<- doubletFinder_v3(compare, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F, sct = T)

# Predict doublet
DF.name = colnames(compare@meta.data)[grepl("DF.classification", colnames(compare@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(compare, group.by = "orig.ident") + NoAxes(),DimPlot(compare, group.by = DF.name) + NoAxes())

# Remove doublet
compare = compare[, compare@meta.data[, DF.name] == "Singlet"]
dim(compare)
# Remove
rm(sweep.res.list,sweep.stats,bcmvn)
gc()

#####################################
##  Integrate datasets ##############
#####################################

compare <- SplitObject(compare, split.by = "orig.ident")
compare <- lapply(X = compare, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = compare)
compare.anchors <- FindIntegrationAnchors(object.list = compare, anchor.features = features)
compare <- IntegrateData(anchorset = compare.anchors)
DefaultAssay(compare) <- "integrated"



#####################################
## Quality control ##################
#####################################
# Quality control
compare[["percent.mt"]] <- PercentageFeatureSet(compare, pattern = "^MT")
compare[["percent.rb"]] <- PercentageFeatureSet(compare, pattern = "^RP[SL]")
# Visualize QC metrics as a violin plot
VlnPlot(compare, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),group.by = "orig.ident",ncol = 4,pt.size = 0.1) 

# Scatter Plot
Scatter1<-FeatureScatter(compare, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident", pt.size = 0.5)
Scatter1
ggsave("Scatter1.png",Scatter1)
Scatter2<-FeatureScatter(compare, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident", pt.size = 0.5)
Scatter2
ggsave("Scatter2.png",Scatter2)
Scatter3<-FeatureScatter(compare, feature1 = "nCount_RNA", feature2 = "percent.rb",group.by = "orig.ident", pt.size = 0.5)
Scatter3
ggsave("Scatter3.png",Scatter3)
Scatter4<-FeatureScatter(compare, feature1 = "percent.rb", feature2 = "percent.mt",group.by = "orig.ident", pt.size = 0.5)
Scatter4
ggsave("Scatter4.png",Scatter4)

# Detection-based filtering
# consider cells with at least 200 detected genes and genes need to be expressed in at least 3 cells
selected_c <- WhichCells(compare, expression =  nFeature_RNA >200 & nFeature_RNA < 4000)
selected_f <- rownames(compare)[Matrix::rowSums(compare) > 3]
compare <- subset(compare, features = selected_f, cells = selected_c)


# Mito/Ribo filtering
selected_mito <- WhichCells(compare, expression = percent.mt < 0.5 & percent.mt >0.05)
selected_ribo <- WhichCells(compare, expression = percent.rb < 60 & percent.rb > 5)
selected_RNA <- WhichCells(compare, expression = nCount_RNA < 39000)
# and subset the object to only keep those cells
compare <- subset(compare, cells = selected_mito)
compare <- subset(compare, cells = selected_ribo)
compare <- subset(compare, cells = selected_RNA)
dim(compare)
table(compare$orig.ident)

# Plot filtered QC
VlnPlot(compare, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),group.by = "orig.ident",ncol = 4,pt.size = 0.1) 


# Calculate cell-cycle scores
compare <- NormalizeData(compare)
compare <- CellCycleScoring(object = compare, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)

# Cell cycle plot
cell_cycle<-VlnPlot(compare, features = c("S.Score", "G2M.Score"), group.by = "orig.ident",ncol = 4, pt.size = 0.1)
cell_cycle
ggsave("cell_cycle.png",cell_cycle,width=20, height=10)




#######################################
# normalize data and find feature genes
#######################################
compare <- NormalizeData(compare)
compare <- FindVariableFeatures(compare, verbose = F,selection.method = "vst", nfeatures = 2000)
# Run the standard workflow for visualization and clustering
# All genes are treated with a mean of 0 and a standard deviation of 1. Highly expressed genes and lowly expressed genes are treated equally. Standardization
# scaladata cell cycle If there are differences, just pick what you want
compare <- ScaleData(compare,verbose = F, vars.to.regress = c("nFeature_RNA", "percent.mt"))
#Keep the variation and process 2000 to 20 features
compare <- RunPCA(compare,npcs =20, verbose = FALSE)
# UMAP TSNE 20 dimensions and then 2 dimensions UMAP has many cells, fast, global accuracy PC1 PC2
compare <- RunUMAP(compare, reduction = "pca", dims = 1:20,verbose = FALSE)

# PCA plot
dimplot<- DimPlot(compare, reduction = "pca",group.by = "orig.ident",dims = 1:2)
# Elbow Plot
elbow<- ElbowPlot(compare, ndims = 30)


compare <- FindNeighbors(compare, dims = 1:20,verbose = FALSE)
compare <- FindClusters(compare, resolution = 0.8, verbose = FALSE)
#table(compare@meta.data$RNA_snn_res.0.8)
compare <- RunUMAP(compare, dims = 1:20, do.fast = TRUE)
DimPlot(compare,reduction = "umap",label=T)
saveRDS(compare, file="compare_seurat_1.Rds")

