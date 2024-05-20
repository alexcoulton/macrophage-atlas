source('~/work/ucl/scripts/misc/functions.R')
library(rhdf5)
library(Seurat)
library(patchwork)
library(dplyr)
library(Matrix)

#mac = readRDS('~/work/ucl/data/ucl.projects/macrophage/all.mac.comb.w.meta.w.norm.rds')

#mac.split = SplitObject(mac, split.by = 'study')

#mac.split = lapply(mac.split, function(x){
    #x = SCTransform(x, variable.features.n = 2000)
#})

#features <- SelectIntegrationFeatures(object.list = mac.split, nfeatures = 2000)
#mac.split = PrepSCTIntegration(object.list = mac.split, anchor.features = features)

#mac.split = lapply(mac.split, function(x){
    #x = RunPCA(x, features = features, verbose = F)
    #x
#})

#anchors <- FindIntegrationAnchors(
    #object.list = mac.split,
    #anchor.features = features,
    #normalization.method = 'SCT',
    #reduction = 'rpca',
    #dims = 1:50
#)

#saveRDS(anchors, '~/work/ucl/bigdata/ucl.projects/macrophage/anchors.new.w.norm.rpca.rds')

source('~/work/ucl/scripts/ucl.projects/macrophage/jun.scripts/integrate_functions_mod.R')
anchors = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/anchors.new.w.norm.rpca.rds')

mac.integrated <- IntegrateDataW(anchorset = anchors, dims = 1:50, normalization.method = 'SCT')
saveRDS(mac.integrated, '~/work/ucl/bigdata/ucl.projects/macrophage/all.dat.w.norm.rpca.integrated.rds')
