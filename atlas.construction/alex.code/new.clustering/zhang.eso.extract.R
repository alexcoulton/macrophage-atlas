library(parallel)

macro.sig = 
    c('MARCO', 'CXCL5', 'SCG5', 'SULT1C2', 
    'MSR1', 'CTSK', 'PTGDS', 'COLEC12', 
    'GPC4', 'PCOLCE2', 'CHIT1', 'KAL1', 
    'CLEC5A', 'ME1', 'DNASE2B', 'CCL7',
    'CD163', 'FN1', 'GM2A', 'SCARB2', 'BCAT1',
    'RAI14', 'COL8A2', 'APOE', 'CHI3L1', 
    'ATG7', 'CD84', 'FDX1', 'MS4A4A', 'SGMS1',
    'EMP1', 'CYBB', 'CD68')

zhangeso.dat = read.delim('~/work/ucl/bigdata/ucl.projects/macrophage/zhang.eso/GSE160269/suppl/GSE160269_CD45pos_UMIs.txt', sep = ' ')
zhangeso.dat = as.matrix(zhangeso.dat)
zhangeso.dat2 = Matrix(zhangeso.dat, sparse = T)

zhangeso.meta = read.csv('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/zhang.eso/zhang.eso.cell.annotation.csv')

head(zhangeso.meta)
unique(zhangeso.meta$celltype)

dim(zhangeso.meta)
sum(grepl('TAM', zhangeso.meta$celltype))
nrow(zhangeso.meta) - sum(grepl('TAM', zhangeso.meta$celltype))

zhang.macro = zhangeso.meta %>%
    filter(celltype %in% c('TAM01', 'TAM02', 'TAM03', 'TAM04'))

mac.bars = zhang.macro[, 1]
mac.bars[!mac.bars %in% colnames(zhangeso.dat2)]

#Here the authors had discordant cell names between their annotations
#and their expression data, the annotations not indicating tumour or normal (e.g. P130T, P130N)
#I checked there were no duplicates and then matched them up

mac.bars.split = strsplit(mac.bars, '\\.') %>%
    lapply(., function(x){
        data.frame(
            p1 = x[[1]],
            p2 = x[[2]],
            p3 = x[[3]]
        )
    }) %>%
    bind_rows

mac.bars.split$row = 1:nrow(mac.bars.split)
mac2 = split(mac.bars.split, mac.bars.split$row)
mac2[[1]]

count = 1
mac.match = mclapply(mac2, function(x){
    print(count)
    mac.orig = paste0(x[[1]], '.', x[[2]], '.', x[[3]])
    col.tag = colnames(zhangeso.dat)[
        grepl(paste0(x[[1]], '[A-Z]*\\.'), colnames(zhangeso.dat)) &
        grepl(x[[2]], colnames(zhangeso.dat)) &
        grepl(x[[3]], colnames(zhangeso.dat))
    ]
    #tryCatch(data.frame(
        #mac.tag = mac.orig,
        #col.tag = col.tag
    #), error = function(e) browser())
    count <<- count + 1
    data.frame(
            mac.tag = mac.orig,
            col.tag = col.tag
        )
}, mc.cores = 8) 

mac.match = bind_rows(mac.match)

zhangeso.macro.dat = zhangeso.dat2[, colnames(zhangeso.dat2) %in% mac.match$col.tag]
zhangeso.non.macro.dat = zhangeso.dat2[, !colnames(zhangeso.dat2) %in% mac.match$col.tag]

#saveRDS(zhangeso.macro.dat, '~/work/ucl/bigdata/ucl.projects/macrophage/zhang.eso/zhang.eso.macrophage.RDS')

mac1 = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/zhang.eso/zhang.eso.macrophage.RDS')
#head(colnames(mac1))
#unique(gsub('(.*?\\..*\\.).*', '\\1', colnames(mac1)))

############################
#BENCHMARKING MACRO SIG 
############################

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

macro.sig2 = macro.sig[macro.sig %in% rownames(zhangeso.macro.dat)]

x = apply(zhangeso.macro.dat[macro.sig2, ], 2, gm_mean)
x2 = apply(zhangeso.non.macro.dat[macro.sig2, ], 2, gm_mean)

mean(x)
mean(x2)

pl(hist(x, breaks = 100, xlim = c(0, 8)))
pl(hist(x2, breaks = 100, xlim = c(0, 8)))

