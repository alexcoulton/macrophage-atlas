source('~/work/ucl/scripts/misc/functions.R')
library(rhdf5)
library(Seurat)
library(ggsci)
library(dplyr)
library(Matrix)

macro.data.dir = '~/work/ucl/bigdata/ucl.projects/macrophage/'

antunes.files.1 = list.files('/camp/home/coultoa/work/ucl/bigdata/ucl.projects/macrophage/antunes/newly.diagnosed.gbm.tam/filtered_feature_bc_matrix_HumanNewlyDiagnGBM/filtered_feature_bc_matrix', full.names = T)
antunes.anno = read.csv('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/antunes/newly.diagnosed.gbm.tam/annot_Human_ND_GBM_TAM.csv')

#derive number of macs and number of non-macs from annotation file
length(which(antunes.anno$cluster == 'Unknown'))
length(which(antunes.anno$cluster != 'Unknown'))

head(antunes.anno)
antunes.anno$cluster

antunes.dat1 = ReadMtx(mtx = antunes.files.1[[3]], cells = antunes.files.1[[1]], features = antunes.files.1[[2]])
antunes.macro.dat1 = antunes.dat1[, colnames(antunes.dat1) %in% antunes.anno$cell]

antunes.files.2 = list.files('/camp/home/coultoa/work/ucl/bigdata/ucl.projects/macrophage/antunes/recurrent.gbm.tam/filtered_feature_bc_matrix_HumanRecurrentGBM/filtered_gene_bc_matrices/GRCh38/', full.names = T)

antunes.dat2 = ReadMtx(mtx = antunes.files.2[[3]], cells = antunes.files.2[[1]], features = antunes.files.2[[2]])
antunes.anno2 = read.csv('~/work/ucl/bigdata/ucl.projects/macrophage/antunes/recurrent.gbm.tam/annot_Human_R_GBM_TAM.csv')


antunes.anno2$key = paste0(antunes.anno2$cell, '_', antunes.anno2$sample)
antunes.anno$key = paste0(antunes.anno$cell, '_', antunes.anno$sample)

antunes.macro.dat2 = antunes.dat2[, colnames(antunes.dat2) %in% antunes.anno2$cell]

colnames(antunes.macro.dat1) = 
    paste0(colnames(antunes.macro.dat1), '_ND')

colnames(antunes.macro.dat2) = 
    paste0(colnames(antunes.macro.dat2), '_R')

all(rownames(antunes.macro.dat1) == rownames(antunes.macro.dat2))

antunes.macro.comb = cbind(antunes.macro.dat1, antunes.macro.dat2)

saveRDS(antunes.macro.comb, '~/work/ucl/bigdata/ucl.projects/macrophage/antunes/antunes.macrophage.RDS')



#antunes.mac = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/antunes/antunes.macrophage.RDS')
#head(antunes.mac)
#unique(gsub('.*-', '', colnames(antunes.mac)))


############################
#CROSS-CHECKING ANTUNES ANNO WITH OUR CLUSTERING 
############################

source('~/work/ucl/scripts/ucl.projects/macrophage/mac.util.R')
antunes.anno = read.csv('~/work/ucl/bigdata/ucl.projects/macrophage/antunes/newly.diagnosed.gbm.tam/annot_Human_ND_GBM_TAM.csv')
antunes.anno2 = read.csv('~/work/ucl/bigdata/ucl.projects/macrophage/antunes/recurrent.gbm.tam/annot_Human_R_GBM_TAM.csv')
load.mac.clus()

head(mac.clus)
head(antunes.anno)
head(antunes.anno2)

antunes.anno$cell = paste0(antunes.anno$cell, '_ND')
antunes.anno2$cell = paste0(antunes.anno2$cell, '_R')

antunes.anno$ident = as.character(antunes.anno$ident)
antunes.anno2$ident = as.character(antunes.anno2$ident)

antunes.anno.all = bind_rows(antunes.anno, antunes.anno2)
antunes.anno.all$cell = paste0('antunes-', antunes.anno.all$cell)

mac.clus.antunes = mac.clus[mac.clus$cellid %in% antunes.anno.all$cell, ]

antunes.anno.all$cellid = antunes.anno.all$cell


match(mac.clus.antunes$cellid)

mac.clus$antunes.cluster = antunes.anno.all$cluster[match(mac.clus$cellid, antunes.anno.all$cellid)]

m1 = mac.clus %>%
    filter(study == 'antunes')

head(m1)

antunes.plot = ggplot(m1, aes(x = UMAP_1, y = UMAP_2, color = antunes.cluster)) +
    geom_point() +
    scale_color_igv() +
    facet_wrap(.~antunes.cluster)

plng(antunes.plot, 1000, 1000)

