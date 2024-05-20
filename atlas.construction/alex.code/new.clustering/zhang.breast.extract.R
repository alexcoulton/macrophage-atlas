source('~/work/ucl/scripts/misc/functions.R')
library(rhdf5)
library(Seurat)
library(dplyr)
library(Matrix)

macro.data.dir = '~/work/ucl/bigdata/ucl.projects/macrophage/'

zhang.files = sort(list.files(paste0(macro.data.dir, '/zhang.breast/GSE169246/suppl'), full.names = T))
zhang.files
zhang.dat2 = ReadMtx(mtx = zhang.files[[5]], cells = zhang.files[[4]], features = zhang.files[[6]], feature.column = 1)

zhang.meta = read.csv('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/zhang.breast/zhang.breast.cell.annotation.csv')

num.mac = sum(grepl('macro', zhang.meta$Cluster))
nrow(zhang.meta) - num.mac

zhang.meta = zhang.meta %>%
    filter(Tissue != 'blood')

zhang.macro = zhang.meta[grepl('macro', zhang.meta$Cluster), ]

colnames(zhang.macro)[1] = 'barcode'
zhang.macro.dat = zhang.dat2[, colnames(zhang.dat2) %in% zhang.macro$barcode]
zhang.non.macro.dat = zhang.dat2[, !colnames(zhang.dat2) %in% zhang.macro$barcode]

saveRDS(zhang.macro.dat, '~/work/ucl/bigdata/ucl.projects/macrophage/zhang.breast/zhang.breast.macrophages.RDS')

zhang.breast = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/zhang.breast/zhang.breast.macrophages.RDS')

#num patients
unique(gsub('.*_(P[0-9]*)_.*', '\\1', sort(unique(gsub('[A-Z]*\\.', '', colnames(zhang.breast))))))

sort(unique(gsub('[A-Z]*\\.', '', colnames(zhang.breast))))

############################
#BENCHMARK 
############################

macro.sig = 
    c('MARCO', 'CXCL5', 'SCG5', 'SULT1C2', 
    'MSR1', 'CTSK', 'PTGDS', 'COLEC12', 
    'GPC4', 'PCOLCE2', 'CHIT1', 'KAL1', 
    'CLEC5A', 'ME1', 'DNASE2B', 'CCL7',
    'CD163', 'FN1', 'GM2A', 'SCARB2', 'BCAT1',
    'RAI14', 'COL8A2', 'APOE', 'CHI3L1', 
    'ATG7', 'CD84', 'FDX1', 'MS4A4A', 'SGMS1',
    'EMP1', 'CYBB', 'CD68')

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

macro.sig2 = macro.sig[macro.sig %in% rownames(zhang.macro.dat)]

x = apply(zhang.macro.dat[macro.sig2, ], 2, gm_mean)
x2 = apply(zhang.non.macro.dat[macro.sig2, ], 2, gm_mean)

mean(x)
mean(x2)

#pl(hist(x, breaks = 100, xlim = c(0, 8)))
#pl(hist(x2, breaks = 100, xlim = c(0, 8)))



