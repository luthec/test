library(dplyr)
library(Seurat)
library(scater)
library(BiocParallel)
library(EnsDb.Hsapiens.v86)

rm(list = ls())
folders=list.files('./')
folders

sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(paste0(folder,"/filtered_feature_bc_matrix")),  project = folder)
})

names(sceList)  = folders
sum(sapply(sceList, function(x)ncol(x@assays$RNA@counts)))

ccl <- merge(x = sceList[[1]], y = sceList[2:9], add.cell.ids = names(sceList), project = "CCL")
ccl <- as.SingleCellExperiment(ccl)

#############filter
rowData(ccl)$Chr <- mapIds(EnsDb.Hsapiens.v86, keys=rownames(ccl),column="SEQNAME", keytype="SYMBOL")

bpp <- SnowParam(8)
ccl <- unfiltered <- addPerCellQC(ccl, BPPARAM=bpp,subsets=list(Mito=which(rowData(ccl)$Chr=="MT")))

qc <- quickPerCellQC(colData(ccl), batch=ccl$ident,sub.fields="subsets_Mito_percent")
ccl <- ccl[,!qc$discard]
unfiltered$discard <- qc$discard

gridExtra::grid.arrange(
    plotColData(unfiltered, x="ident", y="sum", colour_by="discard") +
        scale_y_log10() + ggtitle("Total count"),
    plotColData(unfiltered, x="ident", y="detected", colour_by="discard") +
        scale_y_log10() + ggtitle("Detected features"),
    plotColData(unfiltered, x="ident", y="subsets_Mito_percent",
        colour_by="discard") + ggtitle("Mito percent"),
    ncol=2
)

plotColData(unfiltered, x="sum", y="subsets_Mito_percent", 
    colour_by="discard") + scale_x_log10()

##normalization
ccl <- logNormCounts(ccl, size_factors = ccl$sum)
summary(sizeFactors(ccl))

#####variance model
library(scran)
set.seed(1010010101)
dec.ccl <- modelGeneVarByPoisson(ccl, 
    block=ccl$ident, BPPARAM=bpp)
top.ccl <- getTopHVGs(dec.ccl, n=5000)

par(mfrow=c(5,2))
blocked.stats <- dec.ccl$per.block
for (i in colnames(blocked.stats)) {
    current <- blocked.stats[[i]]
    plot(current$mean, current$total, main=i, pch=16, cex=0.5,
        xlab="Mean of log-expression", ylab="Variance of log-expression")
    curfit <- metadata(current)
    curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)
}

###integration
library(batchelor)
library(BiocNeighbors)

set.seed(1010001)
merged.ccl <- fastMNN(ccl, batch = ccl$ident, subset.row = top.ccl,
     BSPARAM=BiocSingular::RandomParam(deferred = TRUE), 
     BNPARAM=AnnoyParam(),
     BPPARAM=bpp)

reducedDim(ccl, 'MNN') <- reducedDim(merged.ccl, 'corrected')
metadata(merged.ccl)$merge.info$lost.var


set.seed(01010100)
ccl <- runUMAP(ccl, dimred="MNN",
    external_neighbors=TRUE, 
    BNPARAM=AnnoyParam(),
    BPPARAM=bpp,
    n_threads=bpnworkers(bpp))


library(bluster)

set.seed(1000)
colLabels(ccl) <- clusterRows(reducedDim(ccl, "MNN"),
    TwoStepParam(KmeansParam(centers=1000), NNGraphParam(k=5)))

table(colLabels(ccl))

tab <- table(Cluster=colLabels(ccl), Donor=ccl$ident)
library(pheatmap)
pheatmap(log10(tab+10), color=viridis::viridis(100))


scrambled <- sample(ncol(ccl))
gridExtra::grid.arrange(
    plotUMAP(ccl, colour_by="label", text_by="label"),
    plotUMAP(ccl[,scrambled], colour_by="ident"), nrow=2
)

########DE & celltype
markers.ccl <- findMarkers(ccl, block = ccl$ident, 
    direction = 'up', lfc = 1, BPPARAM=bpp)

top.markers <- markers.ccl[["4"]]
best <- top.markers[top.markers$Top <= 10,]
lfcs <- getMarkerEffects(best)

library(pheatmap)
pheatmap(lfcs, breaks=seq(-5, 5, length.out=101))


se.aggregated <- sumCountsAcrossCells(ccl, id=colLabels(ccl), BPPARAM=bpp)

library(celldex)
hpc <- HumanPrimaryCellAtlasData()

library(SingleR)
anno.single <- SingleR(se.aggregated, ref = hpc, labels = hpc$label.main,
    assay.type.test="sum")
anno.single