library(dplyr)
library(Seurat)
library(Azimuth)

set.seed(0)
rm(list = ls())
gc()

##Read in
folders=list.files('./')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder),  project = folder)
})

for(i in 1:length(folders)){
  print(i)
  ####cell filtering
  sceList[[i]] <- PercentageFeatureSet(sceList[[i]], pattern = "^MT-", col.name = "percent.mt") 
  sceList[[i]] <- subset(sceList[[i]],subset = nFeature_RNA > 1000 & nCount_RNA >200 & percent.mt < 20) 
}
names(sceList)  = folders

### label
for(i in 1:2) { sceList[[i]] $label <-"HC"}
for(i in 3:4) { sceList[[i]] $label <-"LTBI"}
for(i in 5:7) { sceList[[i]] $label <-"TB"}


objs <- merge(x = sceList[[1]], y = sceList[2:7], add.cell.ids = names(sceList), project = "TB")

#pre-normalization
DefaultAssay(objs) <- 'RNA'
objs <- objs %>% NormalizeData() %>%
        FindVariableFeatures() %>%
        ScaleData() %>% 
        RunPCA(npcs = 30, verbose = F) 

####integration 
DefaultAssay(objs) <- 'RNA'  
objs <- IntegrateLayers(
  object = objs, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
objs <- FindNeighbors(objs, reduction = "harmony", dims = 1:30)
objs <- FindClusters(objs, resolution = 2, cluster.name = "harmony_clusters")
objs <- RunUMAP(objs, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
objs <- JoinLayers(objs)

###cell type annotation
options(timeout=1000)
objs <- RunAzimuth(objs, reference = "pbmcref")

###DE

cell_type = objs$predicted.celltype.l1 %>% table() %>% names()
objs$celltype.condition <- paste(objs$predicted.celltype.l1, objs$label, sep = "_")
Idents(objs) <- "celltype.condition"

##MAST
for (i in 1:length(cell_type)){
  print(cell_type[i])
  FindMarkers(objs, ident.1 = paste0(cell_type[i],"_TB"), ident.2 = paste0(cell_type[i],"_HC"), verbose = FALSE,test.use = "MAST") %>% 
  rownames_to_column(var = "gene_names_col") %>%
  filter(gene_names_col %in% c("TNFAIP8","TNFAIP8L1","TNFAIP8L2","TNFAIP8L3"))
}

##
# pseudobulk the counts based on donor-condition-celltype
pseudo_objs <- AggregateExpression(objs, assays = "RNA", return.seurat = T, group.by = c("label", "orig.ident", "predicted.celltype.l1"))
head(Cells(pseudo_objs))

pseudo_objs$celltype.condition <- paste(pseudo_objs$predicted.celltype.l1, pseudo_objs$label, sep = "_")
Idents(pseudo_objs) <- "celltype.condition"

print("TB vs LTBI")
for (i in 1:length(cell_type)){
  print(cell_type[i])
  FindMarkers(pseudo_objs, ident.1 = paste0(cell_type[i],"_TB"), ident.2 = paste0(cell_type[i],"_LTBI"), min.cells.group = 2, verbose = FALSE, test.use = "DESeq2") %>% 
  rownames_to_column(var = "gene_names_col") %>%
  filter(gene_names_col %in% c("TNFAIP8","TNFAIP8L1","TNFAIP8L2","TNFAIP8L3")) %>% print()
}

print("TB vs HC")
for (i in 1:length(cell_type)){
  print(cell_type[i])
  FindMarkers(pseudo_objs, ident.1 = paste0(cell_type[i],"_TB"), ident.2 = paste0(cell_type[i],"_HC"), min.cells.group = 2, verbose = FALSE, test.use = "DESeq2") %>% 
  rownames_to_column(var = "gene_names_col") %>%
  filter(gene_names_col %in% c("TNFAIP8","TNFAIP8L1","TNFAIP8L2","TNFAIP8L3")) %>% print()
}

write.table(pseudo_objs@active.ident, file='Convert_UMI_Label.tsv', quote=FALSE, sep='\t', col.names = TRUE)
write.table(pseudo_objs@assays[["RNA"]]$counts, file='Gene_Count_per_Cell.tsv', quote=FALSE, sep='\t', col.names = TRUE)
###figure


p1 <- DimPlot(objs,group.by = "predicted.celltype.l1", split.by = "label", label = TRUE)
p2 <- FeaturePlot(objs,split.by = "label",features = c('TNFAIP8','TNFAIP8L1'), min.cutoff = "q9")
p3 <- FeaturePlot(objs,split.by = "label",features = c('TNFAIP8L2','TNFAIP8L3'), min.cutoff = "q9")

p4 <-VlnPlot(objs,group.by = "predicted.celltype.l1",split.by = "label", features = c("TNFAIP8"))
p5 <-VlnPlot(objs,group.by = "predicted.celltype.l1",split.by = "label", features = c("TNFAIP8L1"))
p6 <-VlnPlot(objs,group.by = "predicted.celltype.l1",split.by = "label", features = c("TNFAIP8L2"))
p7 <-VlnPlot(objs,group.by = "predicted.celltype.l1",split.by = "label", features = c("TNFAIP8L3"))


outpdf=paste("TB","_TNFAIP8.pdf",sep='')
pdf(outpdf, width = 16, height = 10)

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)

dev.off()