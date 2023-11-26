library(dplyr)
library(Seurat)
library(metap)
library(ggrepel)
library(ggpubr)

set.seed(0)
rm(list = ls())
gc()

SCT_intergration <- function(assay_list,n_dim=20) {  
features <- SelectIntegrationFeatures(object.list = assay_list, nfeatures = 3000)
datasets <- PrepSCTIntegration(object.list = assay_list, anchor.features = features, verbose = TRUE)
#datasets <- lapply(X = datasets, FUN = RunPCA, verbose = FALSE, features = features)
anchors <- FindIntegrationAnchors(object.list = datasets, normalization.method = "SCT",anchor.features = features, verbose = TRUE, reference=1,reduction = "cca")
objs <-IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = TRUE, k.weight = 46)

objs <- RunPCA(objs, verbose = FALSE, approx = FALSE, npcs = 30)
# objs <- RunUMAP(objs, reduction = "pca", dims = 1:30, umap.method = "umap-learn", metric = "correlation")
# objs <- RunTSNE(objs, reduction = "pca",dims = 1:30)
objs <- RunUMAP(objs, reduction = "pca", dims = 1:n_dim)
objs <- FindNeighbors(objs, reduction = "pca",dims = 1:n_dim)
FindClusters(objs, resolution = 0.5, algorithm = 2)
}


all_markers.to.plot <- c("CD4","CD8A","CD8B","CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5", 
     "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1", 
    "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ")

tcell_markers.to.plot <- c("CD4","CD8A", "CD8B",'CD40LG',"CCR7", "IL7R", "SELL", "DTHD1", "LYST", "LEF1", "CCR4", "GZMH")

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
cll.test <- merge(x = sceList[[1]], y = sceList[2:9], add.cell.ids = names(sceList), project = "CLL")

p0s<- VlnPlot(cll.test , features = c("CD4", "CD8A"), slot = "counts",combine = FALSE)
CombinePlots(plots = p0s, ncol = 1)    

nc <- merge(x = sceList[[8]], y = sceList[[9]], add.cell.ids = names(sceList)[8:9] , project = "NC") %>% SCTransform(verbose = FALSE)
nc$label <- "NC"
cll <- merge(x = sceList[[2]], y = sceList[3:6], add.cell.ids = names(sceList)[2:6] , project = "CLL") %>% SCTransform(verbose = FALSE)
cll$label <- "CLL"
jsh <- sceList[[7]] %>% SCTransform(verbose = FALSE)
jsh$label <- "JSH"

sample_list=c(nc,cll,jsh)
#features <- SelectIntegrationFeatures(object.list = sample_list, nfeatures = 2000)
#datasets <- PrepSCTIntegration(object.list = sample_list, anchor.features = features, verbose = TRUE)
#datasets <- lapply(X = datasets, FUN = RunPCA, verbose = FALSE, features = features)
#anchors <- FindIntegrationAnchors(object.list = datasets, normalization.method = "SCT",anchor.features = features, verbose = TRUE, reference=1,reduction = "cca")
#objs <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = TRUE)

#objs <- RunPCA(objs, verbose = FALSE, approx = FALSE, npcs = 30)
# objs <- RunUMAP(objs, reduction = "pca", dims = 1:30, umap.method = "umap-learn", metric = "correlation")
# objs <- RunTSNE(objs, reduction = "pca",dims = 1:30)
#objs <- RunUMAP(objs, reduction = "pca", dims = 1:30)
#objs <- FindNeighbors(objs, reduction = "pca",dims = 1:30)
objs <- SCT_intergration(sample_list,30)

DefaultAssay(objs) <- 'SCT'

p1 <- DimPlot(objs, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(objs, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(objs, reduction = "umap", split.by = "label")
# p4 <- FeaturePlot(objs, features = c("CD8A", "IL7R", "CCR7", "S100A4", "GNLY", "NKG7", "FCGR3A", "MS4A7", "CD14"), min.cutoff = "q9")

#### Identify   cell type  

library(celldex)
library(SingleR)
hpc <- celldex::HumanPrimaryCellAtlasData() 
# hpc <- get(load("../single_ref/ref_Human_all.RData"))
library(BiocParallel)
pred.scRNA <- SingleR(test = objs@assays$integrated@data, 
                      ref = hpc,
                      labels = hpc$label.main, 
                      clusters = objs@active.ident)

p33 = plotScoreHeatmap(pred.scRNA, clusters=pred.scRNA@rownames, fontsize.row = 9,show_colnames = T)

pred.scRNA %>% filter(labels=="T cell")

# T-Cell Markers
signature <- readxl::read_excel('C:\\Users\\tliu05\\Desktop\\2023projects\\01_DS\\jiangsu\\Cell_marker_Human.xlsx')
sig_gene  <- signature %>% as.data.frame() %>% filter(`Tissue type`=="Peripheral blood" & `Cell type`== "Cancer cell") %>% mutate(V1=`Cell name`,V2=Symbol) %>% select(V1,V2)  %>% unique(.) 
sig_gene %>% filter(V1=="CD8 T cell" | V1=="CD4 T cell" | V1=="T cell")
sig_gene %>% filter(grepl("T cell",V1))

p4 <- FeaturePlot(objs, features = c("CD8A","CD8B","CD4","CD3D","CD3E","CD3G"),  label = TRUE)
plots0 <- VlnPlot(objs, features = c("CD4", "CD8A","CD8B"), split.by = "label", slot = "counts", assay = "RNA",pt.size = 0, combine = FALSE)
p5 = DotPlot(objs, features = rev(all_markers.to.plot), cols = c("blue", "red"), dot.scale = 8, ) + RotatedAxis()

###target genes check
#BCL2 family
# avg.cd8t.cells.10$gene[grep(genes.to.label[2],avg.cd8t.cells.10$gene)]
# c("BCL2L15","BCL2L11","BCL2L14","BCL2L2","BCL2L2-PABPN1","BCL2L10" ,"BCL2A1","BCL2","BCL2L1","BCL2L12","BCL2L13")   



cd8t.cells.10 <- subset(objs, idents = "10",subset =CD8A > 0 | CD8B > 0)
cd8t.cells.10 <- RenameIdents(cd8t.cells.10,  `10` = "CD8 T")

markers.to.plot <- c("CD3D","CD3E","CD3G", "CD8A","CD8B", "ITGA4","BCL2","BCL2A1")
DotPlot(cd8t.cells.10, features = rev(markers.to.plot),  dot.scale = 8, cols = c("blue", "red","green"), split.by = "label") + RotatedAxis()

#####average expression
aggregate_cd8t.cells.10 <- AggregateExpression(cd8t.cells.10, group.by = "label", return.seurat = TRUE)
genes.to.label = c("ITGA4","BCL2","BCL2A1")

p001 <- CellScatter(aggregate_cd8t.cells.10, "CLL", "NC", highlight = genes.to.label) + ggtitle("CD8+ T Cells")
p001 <- LabelPoints(plot = p001, points = genes.to.label, repel = TRUE)

plots1 = VlnPlot(cd8t.cells.10, features = genes.to.label, split.by = "label", slot = "counts", assay = "RNA",pt.size = 0, combine = FALSE)
CombinePlots(plots = plots1, ncol = 1)


# Create Seurat Object of T-Cell Clusters
DefaultAssay(objs) <- "RNA"
Tcells <- subset(x = objs, subset =CD3D > 0 | CD3E > 0 | CD3G > 0)
#Tcells = subset(objs, cells = WhichCells(objs, idents= c("0", "5", "2", "11", "13", "6")) )
tcell_list2 <- SplitObject(Tcells, split.by = "label")
tcell_list2 <- lapply(tcell_list2, SCTransform, vars.to.regress = "percent.mt")

tcell_combined2 <-SCT_intergration(tcell_list2)

pt1 <- DimPlot(tcell_combined2, group.by = "orig.ident")
pt2 <- DimPlot(tcell_combined2, label = TRUE, repel = TRUE)
pt3 <- DimPlot(tcell_combined2, label=TRUE, reduction = "umap", split.by = "label")
pt4 <- DoHeatmap(tcell_combined2, features = tcell_markers.to.plot , assay = "RNA", slot = "data", angle = 90) + scale_fill_gradientn(colors = c("white", "red"))
plotst <- VlnPlot(tcell_combined2, features = c("CD4", "CD8A","CD8B"), split.by = "label", slot = "counts", assay = "RNA",pt.size = 0, combine = FALSE)
pt5 = DotPlot(tcell_combined2, features = rev(all_markers.to.plot), cols = c("blue", "red"), dot.scale = 8, ) + RotatedAxis()


outpdf=paste("UMAP","_all.pdf",sep='')
pdf(outpdf, width = 16, height = 10)

#print(ggarrange(p1,p2, ncol = 2, nrow = 1,widths=c(1, 1)))
print(p1+p2)
grid::grid.newpage()
print(p33)
print(p4)
CombinePlots(plots = plots0, ncol = 1)
print(p3)
print(p5)

print(pt1+pt2)
print(pt4)
CombinePlots(plots = plotst, ncol = 1)
print(pt3)
print(pt5)
dev.off()


#### Identify  T cell type  

Tcells = subset(objs, cells = WhichCells(objs, idents= c("0", "5", "2", "11", "13", "6")) )


objs<- PrepSCTFindMarkers(objs,assay = "SCT", verbose = TRUE)
DEG <- FindAllMarkers(objs,

                      logfc.threshold=0.25,

                      min.diff.pct = 0.25,

                      max.cells.per.ident = 10000,

                      only.pos=T)

mark_gene <- DEG %>% mutate(avg_logFC=avg_log2FC) %>% filter(p_val_adj<0.05)




