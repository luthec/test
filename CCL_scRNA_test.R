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
library(MAST)
library(irGSEA)
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

Specific_OBJ <- function(spe.obj,tar.gene,tar.cells='T_Cell') {
#poscells <- WhichCells(cd8t.cells.10, expression = ITGA4 > 0)
expr <- FetchData(object = spe.obj, vars = tar.gene ,layer = "counts", assay = "RNA")
poscells <- spe.obj[, which(x = expr >= 1)] %>% colnames()

spe.obj$TAR_exp<- paste0(spe.obj$label,"_",tar.cells,"_",tar.gene,ifelse(colnames(spe.obj) %in% poscells, "+", "-"))
print(table(spe.obj$TAR_exp))
spe.obj<- PrepSCTFindMarkers(spe.obj)
Idents(spe.obj) <- "TAR_exp"
subset(spe.obj,features=setdiff(rownames(spe.obj),tar.gene ))
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



# names(objs.tcells@graphs) 
# names(objs.tcells@reductions)[5]="umap"

#常见的细胞群的maker基因分别是：上皮细胞（EPCAM、KRT19、CLDN4）、基质（PECAM1、CLO1A2、VWF）、增殖性（MKI67、STMN1、PCNA）、T（CD3D、CD3E、CD2）、B（CD79A，IGHG1，MS4A1），NK（KLRD1、GNLY、KLRF1）和髓系（CSF1R、CSF3R、CD68）细胞。
DefaultAssay(objs) <- 'SCT'

objs.tcells <- subset(x = objs, subset =CD3D > 0 | CD3E > 0 | CD3G > 0 | CD2 > 0)
tcells.cluster = Clustering_Cells(objs.tcells,"integrated.cca")

objs.bcells <- subset(x = objs, subset =MS4A1 > 0 | CD79A > 0 | CD19 > 0| IGHG1 > 0)
bcells.cluster = Clustering_Cells(objs.bcells,"integrated.cca")

##cell annotation
library(HGNChelper)
source("https://raw.githubusercontent.com/kris-nader/sc-type/master/R/sctype_wrapper.R"); 

tcells.cluster <- run_sctype(tcells.cluster, known_tissue_type="Immune system",custom_marker_file="C:/Users/tliu05/Desktop/2023projects/01_DS/jiangsu/cell_annotation/ScTypeDB_full.xlsx",name="sctype_classification",plot=TRUE)

## grep(pattern = "^IL", x = rownames(bcells.cluster), value = TRUE)

##cell chat

library(scriabin)

ITGA4.exp.Tcells.spe.names = subset(tcells.cluster[,tcells.cluster$label == "CLL"], idents = "5",subset =( CD8A > 0 | CD8B > 0) & ITGA4 > 0 ) %>% colnames()

ITGA4.nonexp.Tcells.spe.names = subset(tcells.cluster[,tcells.cluster$label == "CLL"], idents = "5",subset =( CD8A > 0 | CD8B > 0) & ITGA4 <= 0 ) %>% colnames()

Immune_supress.Bcells.spe.names = subset(bcells.cluster[,bcells.cluster$label == "CLL"], idents = c("5","6","8","9"), invert = TRUE)  %>% 
                                  subset(,subset = TGFB1 > 0) %>%
                                  subset(, downsample = 100) %>%
                                  colnames()

objs$celltype<- case_when(colnames(objs) %in% ITGA4.exp.Tcells.spe.names ~ paste0(objs$label,"_CD8+ITGA4+_Tcells"),
                          colnames(objs) %in% ITGA4.nonexp.Tcells.spe.names ~ paste0(objs$label,"_CD8+ITGA4-_Tcells"),     
                          colnames(objs) %in% Immune_supress.Bcells.spe.names ~ paste0(objs$label,"_TGFB1+_Bcells"))


# print(table(objs2$celltype))
CCIM <- function(obj,senders.names,receivers.names){
        obj_ccim <- GenerateCCIM(obj, 
                             senders = senders.names,
                             receivers = receivers.names)

                    NormalizeData(obj_ccim) %>% 
                    ScaleData() %>%
                    FindVariableFeatures() %>%
                    RunPCA() %>% 
                    RunUMAP(dims = 1:10) %>%
                    subset(subset = nCount_CCIM > 1) %>% 
                    FindNeighbors(dims = 1:10) %>%
                    FindClusters(resolution = 0.2)
}

objs_ccim = CCIM(subset(objs, subset = nCount_RNA > 100),ITGA4.exp.Tcells.spe.names,Immune_supress.Bcells.spe.names)




DimPlot(objs_ccim, group.by = "receiver_celltype")
DimPlot(objs_ccim, label = T, repel = T) + NoLegend()


CCIMFeaturePlot(objs_ccim, seu = objs, features = c("ITGA4"), type_plot = "sender")


####cellchat
library(CellChat)


CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# objs_cc = objs[ ,!is.na(objs$celltype)]
# cellchat <- createCellChat(object = objs_cc, meta = objs_cc@meta.data, group.by = "celltype", assay = "RNA")
cellchat <- createCellChat(object = objs[ ,!is.na(objs$celltype)], group.by = "celltype", assay = "RNA")


##DB setup
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB$interaction$annotation %>% unique()
# [1] "Secreted Signaling" "ECM-Receptor"       "Cell-Cell Contact" 

CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact", key = "annotation") # use Secreted Signaling

# Only uses the Secreted Signaling from CellChatDB v1
#  CellChatDB.use <- subsetDB(CellChatDB, search = list(c("Secreted Signaling"), c("CellChatDB v1")), key = c("annotation", "version"))

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)


# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat,features = CellChatDB.use$geneInfo$Symbol)


# DB_ligands = dplyr::glimpse(CellChatDB$interaction)$ligand %>% unique()
# DB_receptors = dplyr::glimpse(CellChatDB$interaction)$receptor %>% unique()
# cellchat <- subsetData(cellchat,features = c("ITGA4",DB_ligands,DB_receptors))

cellchat <- updateCellChat(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
#识别过表达配体受体对
cellchat <- identifyOverExpressedInteractions(cellchat)

#project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- projectData(cellchat, PPI.human)
cellchat@data.project[1:4,1:4]

cellchat@LR

cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

#all the inferred cell-cell communications at the level of ligands/receptors
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "cell-cell_communications.all.csv")

#access the the inferred communications at the level of signaling pathways
df.net1 <- subsetCommunication(cellchat,slot.name = "netP")

#gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
levels(cellchat@idents)
df.net2 <- subsetCommunication(cellchat, sources.use = c("Epi"), targets.use = c("Fibroblast" ,"T")) 

#gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
df.net3 <- subsetCommunication(cellchat, signaling = c("CCL", "TGFb"))

#计算每个信号通路相关的所有配体-受体相互作用的通信结果
cellchat <- computeCommunProbPathway(cellchat)
#计算整合的细胞类型之间通信结果
cellchat <- aggregateNet(cellchat)

#######figure
p1 <- DimPlot(tcells.cluster, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(tcells.cluster, reduction = "umap", label = TRUE, repel = TRUE)

p3 <- DimPlot(tcells.cluster, reduction = "umap", split.by = "label", label = TRUE)
p4 <- FeaturePlot(tcells.cluster, reduction = "umap",features = c('CD8A','CD4','CD3D', 'CD3E', 'CD3G','CD2'), min.cutoff = "q9")
p42 <- FeaturePlot(tcells.cluster, reduction = "umap",features = c("ITGA4","BCL2"), min.cutoff = "q9")
p43 <- DimPlot(tcells.cluster, reduction = "umap", group.by = "sctype_classification")



p5 <- DimPlot(bcells.cluster, reduction = "umap", group.by = "orig.ident")
p6 <- DimPlot(bcells.cluster, reduction = "umap", label = TRUE, repel = TRUE)
p7 <- DimPlot(bcells.cluster, reduction = "umap", split.by = "label", label = TRUE)

p8 <- FeaturePlot(bcells.cluster, reduction = "umap",features = c("MS4A1","CD79A","CD19","IGHG1"), min.cutoff = "q9")
p9 <- FeaturePlot(bcells.cluster, reduction = "umap",features = c("TGFB1","IL10"), min.cutoff = "q9")


outpdf=paste("UMAP","_all.pdf",sep='')
pdf(outpdf, width = 16, height = 10)

wrap_plots(c(p01, p02, p03), ncol = 2)

print(p1+p2)
print(p3)
print(p4)
print(p42)
print(p43)

VlnPlot(tcells.cluster, features = c("CD8A"))
VlnPlot(tcells.cluster, features = c("CD4"))
VlnPlot(tcells.cluster, features = c("ITGA4"))
VlnPlot(tcells.cluster, features = c("BCL2"))

print(p5+p6)
print(p7)
print(p8)
print(p9)

dev.off()

