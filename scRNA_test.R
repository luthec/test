library(dplyr)
library(Seurat)

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

sceList <- lapply(X = sceList, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  })

###normalize, integration and filter
normal.features <- SelectIntegrationFeatures(object.list = sceList[8:9])
normal.anchors <- FindIntegrationAnchors(object.list = sceList[8:9], anchor.features = normal.features)
normal.combined <- IntegrateData(anchorset = normal.anchors)
DefaultAssay(normal.combined) <- "normal_integrated"

ccl.features <- SelectIntegrationFeatures(object.list = sceList[1:6])
ccl.anchors <- FindIntegrationAnchors(object.list = sceList[1:6], anchor.features = ccl.features)
ccl.combined <- IntegrateData(anchorset = ccl.anchors)
DefaultAssay(ccl.combined) <- "ccl_integrated"


#######UMAP

normal.combined <- ScaleData(normal.combined, verbose = FALSE)
normal.combined <- RunPCA(normal.combined, npcs = 30, verbose = FALSE)
normal.combined <- RunUMAP(normal.combined, reduction = "pca", dims = 1:30)
normal.combined <- FindNeighbors(normal.combined, reduction = "pca", dims = 1:30)
normal.combined <- FindClusters(normal.combined, resolution = 0.5)


ccl.combined <- ScaleData(ccl.combined, verbose = FALSE)
ccl.combined <- RunPCA(ccl.combined, npcs = 30, verbose = FALSE)
ccl.combined <- RunUMAP(ccl.combined, reduction = "pca", dims = 1:30)
ccl.combined <- FindNeighbors(ccl.combined, reduction = "pca", dims = 1:30)
ccl.combined <- FindClusters(ccl.combined, resolution = 0.5)


p1 <- DimPlot(normal.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(normal.combined, reduction = "umap", label = TRUE, repel = TRUE)


p3 <- DimPlot(ccl.combined, reduction = "umap", group.by = "orig.ident")
p4 <- DimPlot(ccl.combined, reduction = "umap", label = TRUE, repel = TRUE)

p1 + p2 + p3 + p4







table(sce.m$orig.ident)

table(sce.m_f$orig.ident)

save(sce.m_f,file = 'sce.ccl.Rdata')




##############directly integrate

sce.m <- merge(sceList[[1]], 
                 y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],
                       sceList[[6]],sceList[[7]],sceList[[8]],sceList[[9]]), 
                 add.cell.ids = folders, 
                 project = "ccl_test")
sce.m
