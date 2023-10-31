library(dplyr)
library(Seurat)
library(metap)
library(ggrepel)
library(ggpubr)

rm(list = ls())
gc()

folders=list.files('./')
folders

sceList = lapply(folders,function(folder){ 
  print(folder)
  obj = CreateSeuratObject(counts = Read10X(paste0(folder,"/filtered_feature_bc_matrix")),  project = folder)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj,subset = nCount_RNA >200 & percent.mt < 20)
})
names(sceList)  = folders

####pre-QC
ccl.test <- merge(x = sceList[[1]], y = sceList[2:9], add.cell.ids = names(sceList), project = "CCL")

p0s<- VlnPlot(ccl.test , features = c("CD4", "CD8A"), slot = "counts",combine = FALSE)
CombinePlots(plots = p0s, ncol = 1)    

nc <- merge(x = sceList[[8]], y = sceList[[9]], add.cell.ids = names(sceList)[8:9] , project = "NC") %>% SCTransform(verbose = FALSE)
nc$label <- "NC"
ccl <- merge(x = sceList[[2]], y = sceList[3:6], add.cell.ids = names(sceList)[2:6] , project = "CCL") %>% SCTransform(verbose = FALSE)
ccl$label <- "CCL"
jsh <- sceList[[7]] %>% SCTransform(verbose = FALSE)
jsh$label <- "JSH"

sample_list=c(nc,ccl,jsh)
features <- SelectIntegrationFeatures(object.list = sample_list, nfeatures = 2000)
datasets <- PrepSCTIntegration(object.list = sample_list, anchor.features = features, verbose = TRUE)
datasets <- lapply(X = datasets, FUN = RunPCA, verbose = FALSE, features = features)
anchors <- FindIntegrationAnchors(object.list = datasets, normalization.method = "SCT",anchor.features = features, verbose = TRUE, reference=1,reduction = "cca")
objs <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = TRUE)

objs <- RunPCA(objs, verbose = FALSE, approx = FALSE, npcs = 30)
# objs <- RunUMAP(objs, reduction = "pca", dims = 1:30, umap.method = "umap-learn", metric = "correlation")
# objs <- RunTSNE(objs, reduction = "pca",dims = 1:30)
objs <- RunUMAP(objs, reduction = "pca", dims = 1:30)
objs <- FindNeighbors(objs, reduction = "pca",dims = 1:30)
objs <- FindClusters(objs, resolution = 0.5, algorithm = 2)

DefaultAssay(objs) <- 'SCT'

p1 <- DimPlot(objs, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(objs, reduction = "umap", label = TRUE, repel = TRUE)

p3 <- DimPlot(objs, reduction = "umap", split.by = "label")
p4 <- FeaturePlot(objs, features = c("CD8A", "IL7R", "CCR7", "S100A4", "GNLY", "NKG7", "FCGR3A", "MS4A7", "CD14"), min.cutoff = "q9")



VlnPlot(objs, features = c("CD8A"))

VlnPlot(objs, features = c("CD4"))


plots <- VlnPlot(objs, features = c("CD4", "CD8A"), split.by = "label", slot = "counts",
    pt.size = 0, combine = FALSE)

markers.to.plot <- c("CD4","CD8A","CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5", 
     "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1", 
    "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ")
p5 = DotPlot(objs, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, ) + RotatedAxis()

outpdf=paste("UMAP","_all.pdf",sep='')
pdf(outpdf, width = 16, height = 10)

#print(ggarrange(p1,p2, ncol = 2, nrow = 1,widths=c(1, 1)))
print(p1+p2)
print(p3)
print(p4)
CombinePlots(plots = plots, ncol = 1)
print(p5)
dev.off()


#### Identify   cell type  



objs<- PrepSCTFindMarkers(objs,assay = "SCT", verbose = TRUE)
DEG <- FindAllMarkers(objs,

                      logfc.threshold=0.25,

                      min.diff.pct = 0.25,

                      max.cells.per.ident = 10000,

                      only.pos=T)

mark_gene <- DEG %>% mutate(avg_logFC=avg_log2FC) %>% filter(p_val_adj<0.05)

signature <- readxl::read_excel('Cell_marker_Human.xlsx')

sig_gene  <- signature %>% as.data.frame() %>% filter(`Tissue type`=="Peripheral blood" & `Cell type`== "Cancer cell") %>% mutate(V1=`Cell name`,V2=Symbol) %>% select(V1,V2)  %>% unique(.) 


