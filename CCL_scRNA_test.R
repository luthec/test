library(dplyr)
library(Seurat)
library(ggrepel)
library(ggpubr)


SCT_intergration <- function(assay_list,n_dim=20) {  
features <- SelectIntegrationFeatures(object.list = assay_list, nfeatures = 3000)
datasets <- PrepSCTIntegration(object.list = assay_list, anchor.features = features, verbose = TRUE)
anchors <- FindIntegrationAnchors(object.list = datasets, normalization.method = "SCT",anchor.features = features, verbose = TRUE, reference=1,reduction = "cca")
IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = TRUE, k.weight = 46)
}

Clustering_Cells <- function(x,n_dim=20) {
# x <- NormalizeData(x)
# x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
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
sum(sapply(sceList, function(x)ncol(x@assays$RNA@counts)))

for(i in 1:10) { sceList[[i]] $label <-"CLL"}
for(i in 11:17) { sceList[[i]] $label <-"NC"}

objs <-SCT_intergration(sceList[1:17])
objs.tcells <- subset(x = objs, subset =CD3D > 0 | CD3E > 0 | CD3G > 0)

tcells.cluster = Clustering_Cells(objs.tcells)

# Set up control object

nc <-SCT_intergration(sceList[11:17])
nc$label <- "NC"

# DefaultAssay(nc) <- 'RNA'
# nc <- JoinLayers(nc)

cll <-SCT_intergration(sceList[1:10])
cll$label <- "CLL"


DefaultAssay(cll) <- 'SCT'
DefaultAssay(nc) <- 'SCT'

nc.bcells <- subset(x = nc, subset =MS4A1 > 0 | CD79A > 0 | CD19 > 0)
cll.bcells <- subset(x = cll, subset =MS4A1 > 0 | CD79A > 0 | CD19 > 0)

#######figure
p1 <- DimPlot(tcells.cluster, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(tcells.cluster, reduction = "umap", label = TRUE, repel = TRUE)

p3 <- DimPlot(tcells.cluster, reduction = "umap", split.by = "label", label = TRUE)
p4 <- FeaturePlot(tcells.cluster, features = c("CD8A", "GZMK", "CCL5", "CCR7","ITGA4","BCL2"), min.cutoff = "q9")

outpdf=paste("UMAP","_all.pdf",sep='')
pdf(outpdf, width = 16, height = 10)

#print(ggarrange(p1,p2, ncol = 2, nrow = 1,widths=c(1, 1)))
print(p1+p2)
print(p3)
print(p4)

VlnPlot(tcells.cluster, features = c("CD8A"))
VlnPlot(tcells.cluster, features = c("ITGA4"))
VlnPlot(tcells.cluster, features = c("BCL2"))

dev.off()


