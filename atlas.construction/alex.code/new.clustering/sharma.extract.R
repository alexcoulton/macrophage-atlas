sharma = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/sharma/GSE156625/suppl/GSE156625_HCCF_integrated.RData')

head(sharma@meta.data)
sort(table(sharma@meta.data$CT_mye))
dim(sharma@meta.data)
sum(grepl('TAM', sharma@meta.data$CT_mye))
nrow(sharma@meta.data) - sum(grepl('TAM', sharma@meta.data$CT_mye))

dim(GetAssayData(sharma))

macro.sig %in% rownames(GetAssayData(sharma))

unique(sharma[[]]$patientno)

sharma[[]] %>%
    group_by(PatientID, condition) %>%
    summarise(
        n = n()
    ) %>%
    as.data.frame

table(sharma[[]]$CT_mye)
table(sharma[[]]$condition)

sharma.macro.barcodes = sharma[[]] %>%
    filter(grepl('TAM', CT_mye) & condition == 'Tumor') %>%
    rownames

sharma.macro.barcodes

sharma.macro = GetAssayData(sharma)[, colnames(GetAssayData(sharma)) %in%
    sharma.macro.barcodes] 

g1 = '~/work/ucl/bigdata/ucl.projects/macrophage/sharma/GSE156625/suppl/GSE156625_HCCmatrix.mtx.gz'
g2 = '~/work/ucl/bigdata/ucl.projects/macrophage/sharma/GSE156625/suppl/GSE156625_HCCgenes.tsv.gz'
g3 = '~/work/ucl/bigdata/ucl.projects/macrophage/sharma/GSE156625/suppl/GSE156625_HCCbarcodes.tsv.gz'

sharma.raw = ReadMtx(mtx = g1, cells = g3, features = g2)

dim(sharma.raw)

length(colnames(sharma.macro))
macro.raw.barcodes = gsub('^HCC_', '', colnames(sharma.macro))
colnames(sharma.raw) %in%
    gsub('^HCC_', '', colnames(sharma.macro)) %>% sum # check if there are any duplicates for the macrophage columns

dim(sharma.raw)

sharma.raw.non.macro = sharma.raw[, !colnames(sharma.raw) %in% macro.raw.barcodes]
sharma.raw.macro = sharma.raw[, colnames(sharma.raw) %in% macro.raw.barcodes]

saveRDS(sharma.raw.macro, '~/work/ucl/bigdata/ucl.projects/macrophage/sharma/sharma.macro.RDS')

mac = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/sharma/sharma.macro.RDS')

samps.multi = sharma[[]] %>%
    filter(condition == 'Tumor') %>%
    select(PatientID) %>%
    mutate(cellnum = gsub('.*-(.*)', '\\1', rownames(.))) %>%
    distinct %>%
    mutate(new.id = paste0(PatientID, '-s', cellnum))

for(i in 1:nrow(samps.multi)){
    cellnum = samps.multi$cellnum[[i]]
    new.id = samps.multi$new.id[[i]]
    colnames(samps.multi)
    coord = grepl(paste0('-', samps.multi$cellnum[[i]], '$'), colnames(mac))
    colnames(mac)[coord] = gsub(cellnum, new.id, colnames(mac)[coord])
}

saveRDS(mac, '~/work/ucl/bigdata/ucl.projects/macrophage/sharma/sharma.macro2.RDS')


############################
#BENCHMARKING 
############################

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

macro.sig = macro.sig[macro.sig %in% rownames(sharma.raw.macro)]
x = apply(sharma.raw.macro[macro.sig, ], 2, gm_mean)
x2 = apply(sharma.raw.non.macro[macro.sig, ], 2, gm_mean)

mean(x)
x
mean(x2)

h2(sharma.raw.macro)
#macro.sig

#dim(sharma.raw.macro)
#macro.sig
#pl(hist(sharma.raw.macro['CD68', ]))
#pl(hist(sharma.raw.non.macro['CD68', ]))

g4 = '~/work/ucl/bigdata/ucl.projects/macrophage/sharma/GSE156625/suppl/GSE156625_HCCFmatrix.mtx.gz'
g5 = '~/work/ucl/bigdata/ucl.projects/macrophage/sharma/GSE156625/suppl/GSE156625_HCCFgenes.tsv.gz'
g6 = '~/work/ucl/bigdata/ucl.projects/macrophage/sharma/GSE156625/suppl/GSE156625_HCCFbarcodes.tsv.gz'

sharma.raw2 = ReadMtx(mtx = g4, cells = g6, features = g5)

sharma.raw2.macro = sharma.raw2[, colnames(sharma.raw2) %in% macro.raw.barcodes]

dim(sharma.raw2.macro)
x3 = apply(sharma.raw.macro[macro.sig, ], 2, gm_mean)

sharma.raw2.macro



