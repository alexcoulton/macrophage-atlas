source('~/work/ucl/scripts/misc/functions.R')
library(Seurat)
library(ggsci)
library(harmony)
library(rhdf5)
library(Seurat)
library(patchwork)
library(dplyr)
library(Matrix)

mac.int = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/all.dat.w.norm.rpca.integrated.rds')

DefaultAssay(mac.int)
mac.int <- RunPCA(mac.int, npcs = 50, verbose = FALSE)

pl(ElbowPlot(mac.int, ndims = 50))

print('run UMAP')
mac.int <- RunUMAP(mac.int, reduction = "pca", dims = 1:50)

plng(DimPlot(mac.int), 1000, 1000)

mac.int

x = Embeddings(mac.int, reduction = 'umap')

plot1 = ggplot(as.data.frame(x), aes(x = UMAP_1, y = UMAP_2)) +
    geom_point()

plng(plot1, 1000, 1000)

print('find neighbours')
mac.int <- FindNeighbors(mac.int, reduction = "pca", dims = 1:30)
mac.int <- FindClusters(mac.int, resolution = 0.525)

Reductions(mac.int)

head(mac.int@meta.data)
nunique(mac.int@meta.data$integrated_snn_res.0.525)




