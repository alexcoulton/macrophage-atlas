#has macrophage annotations already
source('~/work/ucl/scripts/misc/functions.R')
library(rhdf5)
library(Seurat)
library(dplyr)
library(Matrix)

macro.data.dir = '~/work/ucl/bigdata/ucl.projects/macrophage/'

wu.files = list.files('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/wu/GSE176078/suppl/', pattern = 'CID', full.names = T)

wu.dat = lapply(wu.files, function(x){
    f2 = list.files(x, full.names = T)
    mtx.path = f2[grepl('sparse.mtx', f2)]
    barcodes.path = f2[grepl('barcodes.tsv', f2)]
    genes.path = f2[grepl('genes.tsv', f2)]
    g = ReadMtx(mtx = mtx.path, cells = barcodes.path, features = genes.path, feature.column = 1)
    g
})

wu.dat.comb = do.call(cbind, wu.dat)
head(wu.dat.comb)

wu.meta = read.csv('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/wu/Whole_miniatlas_meta.csv')
head(wu.meta)
unique(wu.meta$celltype_minor)

wu.macro = wu.meta %>%
    filter(celltype_minor == 'Macrophage') %>%
    .[, 1]

dim(wu.dat.comb)
length(wu.macro)

wu.dat.comb2 = wu.dat.comb[, colnames(wu.dat.comb) %in% wu.macro]
saveRDS(wu.dat.comb2, '~/work/ucl/bigdata/ucl.projects/macrophage/wu/wu.macrophages.RDS')


mac = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/wu/wu.macrophages.RDS')

unique(gsub('(.*)_.*', '\\1', colnames(mac)))

length(wu.macro)
ncol(wu.dat.comb) - length(wu.macro)
