source('~/work/ucl/scripts/misc/functions.R')
library(rhdf5)
library(Seurat)
library(dplyr)
library(Matrix)

biermann.files = list.files('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/biermann/GSE200218/suppl/counts/', pattern = '.csv', full.names = T)

biermann.dat = lapply(biermann.files, read.csv)
biermann.dat[[1]][, 1]

colnames(biermann.dat[[1]]) %in%
    colnames(biermann.dat[[2]])

nunique(do.call(c, lapply(biermann.dat, colnames)))
sum(uulapply(biermann.dat, ncol))
uulapply(biermann.dat, nrow)

biermann.meta = read.csv('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/biermann/GSE200218/suppl/GSE200218_sc_sn_metadata.csv')

head(biermann.meta)
unique(biermann.meta$cell_type_fine)
unique(biermann.meta$cell_type_main)

sort(table(biermann.meta$cell_type_fine))

biermann.meta.macro = biermann.meta %>%
    filter(grepl('MDM', cell_type_fine))

h2(biermann.meta.macro)

biermann.dat[[2]][, 1]
biermann.dat.comb = do.call(cbind, biermann.dat)

head(biermann.dat.comb)
biermann.dat.comb[1:10, 1:10]
ncol(biermann.dat.comb) - nrow(biermann.meta.macro)



head(biermann.dat.comb[, 1])

biermann.meta.macro$X = gsub('-', '.', biermann.meta.macro$X)

b.macro.dat = biermann.dat.comb[, colnames(biermann.dat.comb) %in% biermann.meta.macro$X]

rownames(b.macro.dat)

rownames(biermann.dat[[1]])

b.macro.dat[, 1]
dim(b.macro.dat)

rownames(b.macro.dat) = biermann.dat.comb[, 1]

head(b.macro.dat)
b.macro.dat = as.matrix(b.macro.dat)
b.macro.dat2 = Matrix(b.macro.dat, sparse = T)

head(b.macro.dat2)

saveRDS(b.macro.dat2, '~/work/ucl/bigdata/ucl.projects/macrophage/biermann/biermann.macrophages.RDS')

############################
#METADATA 
############################

biermann.meta = read.csv('~/work/ucl/bigdata/ucl.projects/macrophage/biermann/GSE200218/suppl/GSE200218_sc_sn_metadata.csv')
mac = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/biermann/biermann.macrophages.RDS')

mac.names = unique(gsub('[A-Z]+\\.(.*)', '\\1', colnames(mac)))

head(biermann.meta)
head(biermann.meta$X)
colnames(biermann.meta)[1:10]

biermann.meta$X[1:10]

bm2 = biermann.meta %>%
    mutate(id1 = gsub('[A-Z]+-(.*)', '\\1', X)) %>%
    tibble %>%
    select(id1, orig.ident) %>%
    distinct %>%
    filter(id1 %in% mac.names) 


mac.names %in% bm2$id1
bm2$id1 %in% mac.names

table(gsub('[A-Z]+\\.(.*)', '\\1', colnames(mac)))

write.tsv(bm2, '~/work/ucl/data/ucl.projects/macrophage/new.study.metadata/biermann.meta.tsv')


