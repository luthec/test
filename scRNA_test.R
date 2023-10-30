library(dplyr)
library(Seurat)

folders=list.files('./')
folders

sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(paste0(folder,"/filtered_feature_bc_matrix")),  project = folder)
})

for(i in 1:length(folders)){
  sceList[[i]][["percent.mt"]] <- PercentageFeatureSet(sceList[[i]], pattern = "^MT-")
  sceList[[i]] <- subset(sceList[[i]],subset=nFeature_RNA>200 & percent.mt<20)
}

names(sceList)  = folders
sum(sapply(sceList, function(x)ncol(x@assays$RNA@counts)))

sceList <- lapply(X = sceList, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })


ccl.anchors <- FindIntegrationAnchors(object.list = sceList[1:6], dims = 1:20)
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

p3 <- FeaturePlot(ccl.combined, features = c("CD8A", "IL7R", "CCR7", "S100A4", "GNLY", "NKG7", "FCGR3A", "MS4A7", "CD14"), min.cutoff = "q9")



markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5", 
    "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1", 
    "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ")

DotPlot(ccl.combined, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, ) + RotatedAxis()




0	IL7R, CCR7	Naive CD4+ T
1	CD14, LYZ	CD14+ Mono
2	IL7R, S100A4	Memory CD4+
3	MS4A1	B
4	CD8A	CD8+ T
5	FCGR3A, MS4A7	FCGR3A+ Mono
6	GNLY, NKG7	NK
7	FCER1A, CST3	DC
8	PPBP	Platelet









CD3D
CD3E
CD3G
TRAC
CD4
TCF7
CD27
IL7R
CD8A
CD8B
GNLY
NKG7
CST7
GZMK
CCL5
FCER1G
CD19
MS4A1
CD79A
CD79B
MZB1
IGHD
IGHM
NCAM1
KLRB1
KLRD1
KLRF1
KLRC1
KLRC2
KLRC3
KLRC4
CD14
FCGR3A
FCGR3B
ITGAL
ITGAM
VCAN
FCN1
S100A8
S100A9
CSF3R
CSF1R
CX3CR1
TYROBP
LYZ
S100A12
CDKN1C
MS4A7
HLA−DPB1
HLA−DPA1
HLA−DQA1
ITGAX
CD1C
CD1E
FCER1A
CLEC10A
FCGR2B
IL3RA
GZMB
JCHAIN
IRF7
TCF4
LILRA4
CLEC4C
CD38
XBP1
SLAMF7
IGHA1
IGHA2
IGHG1
IGHG2
IGHG3
IGHG4
PF4
PPBP
GP5
ITGA2B
NRGN
TUBB1
SPARC
RGS18
MYL9
GNG11
DPP4
CTSL





#####################
normal.combined <- ScaleData(normal.combined, verbose = FALSE)
normal.combined <- RunPCA(normal.combined, npcs = 30, verbose = FALSE)
normal.combined <- RunUMAP(normal.combined, reduction = "pca", dims = 1:30)
normal.combined <- FindNeighbors(normal.combined, reduction = "pca", dims = 1:30)
normal.combined <- FindClusters(normal.combined, resolution = 0.5)


normal.features <- SelectIntegrationFeatures(object.list = sceList[8:9])
normal.anchors <- FindIntegrationAnchors(object.list = sceList[8:9], anchor.features = normal.features)
normal.combined <- IntegrateData(anchorset = normal.anchors)
DefaultAssay(normal.combined) <- "integrated"



p1 <- DimPlot(normal.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(normal.combined, reduction = "umap", label = TRUE, repel = TRUE)


library(celldex)
library(SingleR)





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
