source('~/work/ucl/scripts/misc/functions.R')
library(rhdf5)
library(Seurat)
library(gridExtra)
library(dplyr)
library(Matrix)

#kim.dat = read.tsv('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/kim/GSE131907_Lung_Cancer_raw_UMI_matrix.txt')

#genes1 = kim.dat[, 1]
#cells1 = colnames(kim.dat)

#kim.dat = kim.dat[, -1]

#rownames(kim.dat) = genes1

#kim.dat[1:10, 1:10]

#kim.dat = kim.dat[, gsub('.*?_(.*)_.*', '\\1', colnames(kim.dat)) == 'LUNG']

#kim.dat = as.matrix(kim.dat)

#write.tsv(kim.dat, '~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/kim/only.lung.tsv')


#kim.dat2 = Matrix(kim.dat, sparse = T)

#saveRDS(kim.dat2, '~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/kim/only.lung.sparse.rds')
kim.dat2 = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/kim/only.lung.sparse.rds')

kim.seurat = CreateSeuratObject(
    kim.dat2,
    project = 'test',
    min.cells = 3,
    min.features = 200
)

kim.seurat

