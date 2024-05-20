source('~/work/ucl/scripts/misc/functions.R')
library(rhdf5)
library(Seurat)
library(dplyr)
library(Matrix)

macro.data.dir = '~/work/ucl/bigdata/ucl.projects/macrophage/'

khaliq.dat = read.csv('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/khaliq/GSE200997/suppl/GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv')

rownames(khaliq.dat) = khaliq.dat$X
khaliq.dat = khaliq.dat[, -1]
khaliq.dat = as.matrix(khaliq.dat)
khaliq.dat2 = Matrix(khaliq.dat, sparse = T)

'CD68' %in% rownames(khaliq.dat2)

khaliq.anno = read.csv('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/khaliq/CRC-Single-Cell/Data/Myeloid/Myeloid_Metadata.csv')

head(khaliq.anno)

table(khaliq.anno$seurat_clusters)

head(khaliq.anno)

table(khaliq.anno$seurat_clusters)

khaliq.macro = khaliq.anno %>%
    filter(grepl('MACRO', seurat_clusters))

head(khaliq.macro)

dim(khaliq.macro)
sum(colnames(khaliq.dat) %in% khaliq.macro$X)

khaliq.macro.dat = khaliq.dat2[, colnames(khaliq.dat2) %in% khaliq.macro$X]


dim(khaliq.macro.dat)
saveRDS(khaliq.macro.dat, '~/work/ucl/bigdata/ucl.projects/macrophage/khaliq/khaliq.macrophages.RDS')

mac = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/khaliq/khaliq.macrophages.RDS')

############################
#REMOVE ADJACENT NORMAL SAMPLES 
############################

mac = mac[, !grepl('^B_cac[0-9]*', colnames(mac))]
saveRDS(mac, '~/work/ucl/bigdata/ucl.projects/macrophage/khaliq/khaliq.macrophages2.RDS')

