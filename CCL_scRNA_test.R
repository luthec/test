library(dplyr)
library(Seurat)
library(SeuratData)
# library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2)


SCT_intergration <- function(assay_list,n_dim=20) {  
features <- SelectIntegrationFeatures(object.list = assay_list, nfeatures = 3000)
datasets <- PrepSCTIntegration(object.list = assay_list, anchor.features = features, verbose = TRUE)
anchors <- FindIntegrationAnchors(object.list = datasets, normalization.method = "SCT",anchor.features = features, verbose = TRUE, reference=1,reduction = "cca")
IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = TRUE, k.weight = 46)
}

Clustering_Cells <- function(x,n_dim=20) {
# x <- NormalizeData(x)
# x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
DefaultAssay(x) <- "integrated"
x <- ScaleData(x, verbose = FALSE)
x <- RunPCA(x, npcs = 30, verbose = FALSE)
x <- RunUMAP(x, reduction = "pca", dims = 1:n_dim)
x <- FindNeighbors(x, reduction = "pca", dims = 1:n_dim)
FindClusters(x, resolution = 0.5, algorithm = 2)
}

##########start


folders=list.files('./')
folders

sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder),  project = folder)
})

for(i in 1:length(folders)){
  print(i)
  sceList[[i]] <- PercentageFeatureSet(sceList[[i]], pattern = "^MT-", col.name = "percent.mt") 
  sceList[[i]] <- subset(sceList[[i]],subset = nCount_RNA >200 & percent.mt < 20) %>% 
  SCTransform(vars.to.regress = "percent.mt", verbose = FALSE,vst.flavor = "v2")
}

names(sceList)  = folders

##sum(sapply(sceList, function(x)ncol(x@assays$RNA@counts)))

for(i in 1:10) { sceList[[i]] $label <-"CLL"}
for(i in 11:17) { sceList[[i]] $label <-"NC"}


objs <- merge(x = sceList[[1]], y = sceList[2:17], add.cell.ids = names(sceList), project = "CLL")

DefaultAssay(objs) <- 'RNA'

objs <- objs %>% NormalizeData() %>%
        FindVariableFeatures() %>%
        ScaleData() %>% 
        RunPCA(npcs = 30, verbose = F) 
        
DefaultAssay(objs) <- 'RNA'      

objs_sct <- objs  %>% SCTransform() 
objs_sct <- IntegrateLayers(
            object = objs_sct,
            method = RPCAIntegration,
            normalization.method = "SCT",
            verbose = F)


objs <- IntegrateLayers(
  object = objs, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)



 %>%    FindNeighbors(dims = 1:30, reduction = "integrated.dr") %>%
        FindClusters(resolution = 2)








objs <-SCT_intergration(sceList[1:17])

DefaultAssay(objs) <- 'SCT'

objs.tcells <- subset(x = objs, subset =CD3D > 0 | CD3E > 0 | CD3G > 0)
tcells.cluster = Clustering_Cells(objs.tcells)

objs.bcells <- subset(x = objs, subset =MS4A1 > 0 | CD79A > 0 | CD19 > 0)
bcells.cluster = Clustering_Cells(objs.bcells)


#######figure
p1 <- DimPlot(tcells.cluster, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(tcells.cluster, reduction = "umap", label = TRUE, repel = TRUE)

p3 <- DimPlot(tcells.cluster, reduction = "umap", split.by = "label", label = TRUE)
p4 <- FeaturePlot(tcells.cluster, features = c("CD8A", "GZMK", "CCL5", "CCR7","ITGA4","BCL2"), min.cutoff = "q9")


p5 <- DimPlot(bcells.cluster, reduction = "umap", group.by = "orig.ident")
p6 <- DimPlot(bcells.cluster, reduction = "umap", label = TRUE, repel = TRUE)
p7 <- DimPlot(bcells.cluster, reduction = "umap", split.by = "label", label = TRUE)

outpdf=paste("UMAP","_all.pdf",sep='')
pdf(outpdf, width = 16, height = 10)

print(p1+p2)
print(p3)
print(p4)

VlnPlot(tcells.cluster, features = c("CD8A"))
VlnPlot(tcells.cluster, features = c("ITGA4"))
VlnPlot(tcells.cluster, features = c("BCL2"))

print(p5+p6)
print(p7)

dev.off()

tcells.cluster <- RunAzimuth(tcells.cluster, reference = "pbmcref")