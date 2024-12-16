library(dplyr)
library(Seurat)
library(metap)
library(ggrepel)
library(ggpubr)
library(MAST)
library(irGSEA)

set.seed(0)
rm(list = ls())


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