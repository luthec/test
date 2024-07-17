library(dplyr)
library(Seurat)
library(SeuratData)
# library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(future)
library(patchwork)
plan("multicore", workers = 4)



SCT_intergration <- function(assay_list,n_dim=20) {  
features <- SelectIntegrationFeatures(object.list = assay_list, nfeatures = 3000)
datasets <- PrepSCTIntegration(object.list = assay_list, anchor.features = features, verbose = TRUE)
anchors <- FindIntegrationAnchors(object.list = datasets, normalization.method = "SCT",anchor.features = features, verbose = TRUE, reference=1,reduction = "cca")
IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = TRUE, k.weight = 46)
}

Clustering_Cells <- function(x,re,n_dim=20) {
# x <- NormalizeData(x)
# x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
# DefaultAssay(x) <- "integrated"
# x <- ScaleData(x, verbose = FALSE)
x <- RunPCA(x, npcs = 30, verbose = FALSE)
x <- FindNeighbors(x, reduction = re, dims = 1:n_dim)
x <- FindClusters(x,graph.name = "RNA_snn", reduction = re,resolution = 0.5, algorithm = 2)
RunUMAP(x, reduction = re, dims = 1:n_dim)
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
  sceList[[i]] <- subset(sceList[[i]],subset = nFeature_RNA > 1000 & nCount_RNA >200 & percent.mt < 20) 
  # SCTransform(vars.to.regress = "percent.mt", verbose = FALSE,vst.flavor = "v2")
}

names(sceList)  = folders

##sum(sapply(sceList, function(x)ncol(x@assays$RNA@counts)))

for(i in 1:10) { sceList[[i]] $label <-"CLL"}
for(i in 11:17) { sceList[[i]] $label <-"NC"}


objs <- merge(x = sceList[[1]], y = sceList[2:17], add.cell.ids = names(sceList), project = "CLL")

####cell annotation
objs <- subset(objs, nFeature_RNA > 1000)
# objs <- RunAzimuth(objs, reference = "pbmcref")

#pre-normalization
DefaultAssay(objs) <- 'RNA'
objs <- objs %>% NormalizeData() %>%
        FindVariableFeatures() %>%
        ScaleData() %>% 
        RunPCA(npcs = 30, verbose = F) 
        
####integration 1
DefaultAssay(objs) <- 'RNA'      
options(future.globals.maxSize = 50000 * 1024^2)
objs <- objs  %>% SCTransform() 
objs <- IntegrateLayers(
            object = objs,
            method = RPCAIntegration,
            normalization.method = "SCT",
            verbose = F)

objs <- FindNeighbors(objs, dims = 1:30, reduction = "integrated.dr")
objs <- FindClusters(objs,graph.name = "RNA_snn",resolution = 2, cluster.name = "integrated_clusters")
objs <- RunUMAP(objs, dims = 1:30, reduction = "integrated.dr", reduction.name = "umap.integrated")

####integration 2
DefaultAssay(objs) <- 'RNA'  
objs <- IntegrateLayers(
  object = objs, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)
objs <- FindNeighbors(objs, reduction = "integrated.cca", dims = 1:30)
objs <- FindClusters(objs, resolution = 2, cluster.name = "cca_clusters")
objs <- RunUMAP(objs, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

####integration 3
DefaultAssay(objs) <- 'RNA'  
objs <- IntegrateLayers(
  object = objs, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
objs <- FindNeighbors(objs, reduction = "harmony", dims = 1:30)
objs <- FindClusters(objs, resolution = 2, cluster.name = "harmony_clusters")
objs <- RunUMAP(objs, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

###
objs <- JoinLayers(objs)


p01 <- DimPlot(
  objs,
  reduction = "umap.cca",
  group.by = c("orig.ident", "cca_clusters"),
  combine = FALSE, label.size = 2
)

p02 <- DimPlot(
  objs,
  reduction = "umap.integrated",
  group.by = c("orig.ident",  "integrated_clusters"),
  combine = FALSE, label.size = 2
)

p03 <- DimPlot(
  objs,
  reduction = "umap.harmony",
  group.by = c("orig.ident",  "harmony_clusters"),
  combine = FALSE, label.size = 2
)

wrap_plots(c(p01, p02, p03), ncol = 2)



# names(alldata@graphs) 

#常见的细胞群的maker基因分别是：上皮细胞（EPCAM、KRT19、CLDN4）、基质（PECAM1、CLO1A2、VWF）、增殖性（MKI67、STMN1、PCNA）、T（CD3D、CD3E、CD2）、B（CD79A，IGHG1，MS4A1），NK（KLRD1、GNLY、KLRF1）和髓系（CSF1R、CSF3R、CD68）细胞。
DefaultAssay(objs) <- 'SCT'

objs.tcells <- subset(x = objs, subset =CD3D > 0 | CD3E > 0 | CD3G > 0 | CD2 > 0)
tcells.cluster = Clustering_Cells(objs.tcells,"integrated.cca")

objs.bcells <- subset(x = objs, subset =MS4A1 > 0 | CD79A > 0 | CD19 > 0| IGHG1 > 0)
bcells.cluster = Clustering_Cells(objs.bcells,"integrated.cca")


#######figure
p1 <- DimPlot(tcells.cluster, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(tcells.cluster, reduction = "umap", label = TRUE, repel = TRUE)

p3 <- DimPlot(tcells.cluster, reduction = "umap", split.by = "label", label = TRUE)
p4 <- FeaturePlot(tcells.cluster, reduction = "umap",features = c('CD8A','CD4','CD3D', 'CD3E', 'CD3G','CD2'), min.cutoff = "q9")
p42 <- FeaturePlot(tcells.cluster, reduction = "umap",features = c("ITGA4","BCL2"), min.cutoff = "q9")


p5 <- DimPlot(bcells.cluster, reduction = "umap", group.by = "orig.ident")
p6 <- DimPlot(bcells.cluster, reduction = "umap", label = TRUE, repel = TRUE)
p7 <- DimPlot(bcells.cluster, reduction = "umap", split.by = "label", label = TRUE)
p8 <- FeaturePlot(bcells.cluster, reduction = "umap",features = c("MS4A1","CD79A","CD19","IGHG1"), min.cutoff = "q9")

outpdf=paste("UMAP","_all.pdf",sep='')
pdf(outpdf, width = 16, height = 10)

wrap_plots(c(p01, p02, p03), ncol = 2)

print(p1+p2)
print(p3)
print(p4)
print(p42)

VlnPlot(tcells.cluster, features = c("CD8A"))
VlnPlot(tcells.cluster, features = c("CD4"))
VlnPlot(tcells.cluster, features = c("ITGA4"))
VlnPlot(tcells.cluster, features = c("BCL2"))

print(p5+p6)
print(p7)
print(p8)

dev.off()

tcells.cluster <- RunAzimuth(tcells.cluster, reference = "pbmcref")