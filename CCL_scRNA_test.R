library(dplyr)
library(Seurat)
library(metap)

folders=list.files('./')
folders

sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(paste0(folder,"/filtered_feature_bc_matrix")),  project = folder)
})

for(i in 1:length(folders)){
  sceList[[i]][["percent.mt"]] <- PercentageFeatureSet(sceList[[i]], pattern = "^MT-")
  sceList[[i]] <- subset(sceList[[i]],subset=nFeature_RNA>20 & nFeature_RNA <2500 & percent.mt<10)
}

names(sceList)  = folders
sum(sapply(sceList, function(x)ncol(x@assays$RNA@counts)))

# Set up control object
nc <- merge(x = sceList[[8]], y = sceList[[9]], add.cell.ids = names(sceList)[8:9] , project = "NC")
nc$label <- "NC"
nc <- NormalizeData(nc, verbose = FALSE)
nc <- FindVariableFeatures(nc, selection.method = "vst", nfeatures = 2000)

# Set up case object
ccl <- merge(x = sceList[[1]], y = sceList[2:6], add.cell.ids = names(sceList)[1:6] , project = "CCL")
ccl$label <- "CCL"
ccl <- NormalizeData(ccl, verbose = FALSE)
ccl <- FindVariableFeatures(ccl, selection.method = "vst", nfeatures = 2000)

###normalize, integration and filter
ccl.anchors <- FindIntegrationAnchors(object.list = list(nc, ccl), dims = 1:20)
ccl.combined <- IntegrateData(anchorset = ccl.anchors, dims = 1:20)
DefaultAssay(ccl.combined) <- "integrated"

#######UMAP

ccl.combined <- ScaleData(ccl.combined, verbose = FALSE)
ccl.combined <- RunPCA(ccl.combined, npcs = 30, verbose = FALSE)
ccl.combined <- RunUMAP(ccl.combined, reduction = "pca", dims = 1:30)
ccl.combined <- FindNeighbors(ccl.combined, reduction = "pca", dims = 1:30)
ccl.combined <- FindClusters(ccl.combined, resolution = 0.5)

p1 <- DimPlot(ccl.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(ccl.combined, reduction = "umap", label = TRUE, repel = TRUE)

p1 + p2 

DimPlot(ccl.combined, reduction = "umap", split.by = "label")

#### Identify conserved cell type markers

DefaultAssay(ccl.combined) <- "RNA"
con.markers <- FindConservedMarkers(ccl.combined, ident.1 = 0, grouping.var = "label", verbose = FALSE)
head(con.markers)

FeaturePlot(ccl.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", 
    "CCL2", "PPBP"), min.cutoff = "q9")


ccl.markers <- FindAllMarkers(ccl.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ccl.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

IL7R, CCR7

VlnPlot(ccl.combined, features = c("IL7R", "CCR7"))