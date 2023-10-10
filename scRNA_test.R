library(dplyr)
library(Seurat)

folders=folders=list.files('./')
folders

sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(paste0(folder,"/filtered_feature_bc_matrix")),  project = folder)
})

sce.m <- merge(sceList[[1]], 
                 y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],
                       sceList[[6]],sceList[[7]],sceList[[8]],sceList[[9]]), 
                 add.cell.ids = folders, 
                 project = "ccl_test")
sce.m


sce.m[["percent.mt"]]<-PercentageFeatureSet(sce.m,pattern = "^MT-")

sce.m_f<-subset(sce.m,subset=nFeature_RNA>20 & nFeature_RNA <2500 & percent.mt<5)


table(sce.m$orig.ident)

table(sce.m_f$orig.ident)

save(sce.m_f,file = 'sce.ccl.Rdata')




##############

for(i in 1:length(dirs)){
  scelist[[i]][["percent.mt"]] <- PercentageFeatureSet(scelist[[i]], pattern = "^MT-")
  scelist[[i]] <- subset(scelist[[i]], subset = percent.mt < 10)
}
names(scelist)  = paste0("R",1:4)
sum(sapply(scelist, function(x)ncol(x@assays$RNA@counts)))


ccl0.data <- Read10X(data.dir="C:/Users/tliu05/Desktop/raw.data/CLL/filtered_feature_bc_matrix")

ccl0 <- CreateSeuratObject(counts = ccl0.data, project = "ccl_test" , min.cells = 3,min.features = 200)

ccl0[["percent.mt"]]<-PercentageFeatureSet(ccl0,pattern = "^MT-")

ccl0_f<-subset(ccl0,subset=nFeature_RNA>20 & nFeature_RNA <2500 & percent.mt<5)

head(ccl0_f@meta.data, 5)
