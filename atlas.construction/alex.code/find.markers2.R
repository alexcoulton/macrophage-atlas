source('~/work/ucl/scripts/misc/functions.R')
library(Seurat)
library(ggsci)
library(rhdf5)
library(Seurat)
library(patchwork)
library(dplyr)
library(Matrix)

i = commandArgs(trailingOnly=TRUE)
i = as.numeric(i[[1]])

print(i)
print(paste0('finding markers for cluster ', i))

mac.int = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/all.dat.w.norm.rpca.integrated.rds')
mac.clus = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/all.dat.w.norm.w.cluster.meta.data.rds')

mac.int@meta.data = mac.clus
Idents(mac.int) = 'integrated_snn_res.0.525'

#DefaultAssay(mac.int) = 'SCT'
#mac.int = PrepSCTFindMarkers(mac.int)

DefaultAssay(mac.int) = 'RNA'
mac.int <- NormalizeData(mac.int, normalization.method = "LogNormalize", scale.factor = 10000)

clusmarkers = FindMarkers(mac.int, ident.1 = i, min.pct = 0.25)

#cluster 23 failed because it only has a few cells in it.
#clusmarkers = FindMarkers(mac.int, ident.1 = 23, min.pct = 0.25)
#which(Idents(mac.int) == '23')

saveRDS(
    clusmarkers,
    paste0('~/work/ucl/data/ucl.projects/macrophage/markers/lognorm.assay.all.dat.w.norm.clus', i, 'markers.rds')
)

#for(i in 10:26){
    #clusmarkers = FindMarkers(mac.int, ident.1 = i, min.pct = 0.25)
    #saveRDS(
        #clusmarkers,
        #paste0('~/work/ucl/data/ucl.projects/macrophage/markers/clus', i, 'markers.rds')
    #)
#}

