############################
#NOTES 
############################

# jun.dat.comb.raw.rds has now gone through the gene name conversion process
# new.macrophage.seurat.comb.rds has also gone through the gene name conversion process
# all.mac.comb.w.meta.rds has also gone through the gene name conversion process

############################
#SCRIPT 
############################

source('~/work/ucl/scripts/misc/functions.R')
library(rhdf5)
library(Seurat)
library(dplyr)
library(org.Hs.eg.db)
library(Matrix)
library(HGNChelper)
#library(MatrixExtra)

macro.data.dir = '~/work/ucl/bigdata/ucl.projects/macrophage/'

path1 = '/camp/home/coultoa/work/ucl/bigdata/ucl.projects/macrophage/'

authors = c(
    'durante',
    'antunes',
    'zhang.breast',
    'zheng',
    'khaliq',
    'zhang.eso',
    'pu',
    'zhang.ovar',
    'xu',
    'sharma',
    'lu',
    'becker',
    'che',
    'biermann',
    'wu'
)

macro.paths = uulapply(authors, function(x){
    list.files(paste0(path1, x), pattern = 'macro.*RDS', full.names = T)
})

#writeLines(macro.paths, '~/work/ucl/data/ucl.projects/macrophage/mac.paths.txt')
macro.paths = readLines('~/work/ucl/data/ucl.projects/macrophage/mac.paths.txt')

macro.dat = lapply(macro.paths, function(x){
    print(1)
    readRDS(x)
})

lapply(macro.dat, dim)

############################
#GENE NAME CORRECTION 
############################

current.map = getCurrentHumanMap()

count = 1
macro.genes = lapply(macro.dat, function(x){
    study = authors[[count]]
    count <<- count + 1
    data.frame(
        study = study,
        gene = rownames(x)
    )
}) %>% bind_rows

table(sort(table(macro.genes$gene)))

g = sort(table(macro.genes$gene))

g[g == 1]

head(macro.genes)

gene.syms = checkGeneSymbols(macro.genes$gene, map = current.map)
head(gene.syms)

to.correct = gene.syms %>%
    filter(Approved == F & !is.na(Suggested.Symbol))

to.correct$present = to.correct$Suggested.Symbol %in% macro.genes$gene
to.correct = distinct(to.correct)

dim(to.correct)

macro.dat = lapply(macro.dat, function(x){
    head(to.correct)
    update = rownames(x)[rownames(x) %in% to.correct$x]
    update2 = to.correct %>% filter(x %in% update) %>% reset.rownames
    print(paste0('corrected: ', nrow(update2)))

    rownames(x)[match(update2$x, rownames(x))] = update2$Suggested.Symbol
    x
})

head(macro.dat[[1]])

############################
#CORRECT ENSEMBL IDS 
############################

head(macro.genes)
ens.ids = unique(macro.genes$gene[grepl('ENSG', macro.genes$gene)])
ens.mapped = mapIds(org.Hs.eg.db, keys = ens.ids, keytype = 'ENSEMBL', column = 'SYMBOL', multiVals = 'first')
ens.mapped = ens.mapped[!is.na(ens.mapped)]
ens.mapped = as.data.frame(ens.mapped) %>%
    mutate(ens.id = rownames(.)) %>%
    reset.rownames


macro.dat = lapply(macro.dat, function(x){
    head(ens.mapped)
    update = rownames(x)[rownames(x) %in% ens.mapped$ens.id]
    if(length(update) == 0){
        print('corrected: 0')
        return(x)
    }
    update2 = ens.mapped %>% filter(ens.id %in% update) %>% reset.rownames
    print(paste0('corrected: ', nrow(update2)))

    rownames(x)[match(update2$ens.id, rownames(x))] = update2$ens.mapped
    x
})

#all.gene.names = unique(uulapply(macro.dat, rownames))

#writeLines(sort(all.gene.names), '~/work/ucl/tmp/all.gene.names.txt')

############################
#MERGING 
############################

count = 1
macro.dat = Map(function(x, y){
    print(count)
    colnames(x) = paste0(y, '-', colnames(x))
    x
}, macro.dat, authors)

macro.seurat = lapply(macro.dat, function(x){
    x.seurat = CreateSeuratObject(
        x,
        project = 'test',
    )
})

all.seurat.merge = merge(macro.seurat[[1]], macro.seurat[2:15])
all.seurat.merge[[]]$cellname 
head(all.seurat.merge[[]])

all.seurat.merge@meta.data$cellname = rownames(all.seurat.merge@meta.data)
head(all.seurat.merge@meta.data)
all.seurat.merge@meta.data$study = gsub('([a-z]+)-.*', '\\1', all.seurat.merge@meta.data$cellname)

rownames(all.seurat.merge[[]])

x = all.seurat.merge
y = x@meta.data 

############################
#INTEGRATE DURANTE META 
############################

durante.meta = read.csv('~/work/ucl/data/ucl.projects/macrophage/new.study.metadata/durante.meta.csv', header = F)

y$sample.id = ""
head(y)

y[y$study == 'durante', "sample.id"] = 
    gsub(
        '[a-z]*-(.*)-.*-.*',
        '\\1',
        y[y$study == 'durante', "cellname"]
    )

head(y)
head(durante.meta)
durante.meta

macro.paths[grepl('2', macro.paths)]

head(y, n = 2)
head(durante.meta, n = 2)

durante.meta = durante.meta[c('V1', 'V3')]
colnames(durante.meta) = c('sample.id', 'prim.met.other')
durante.meta$study = 'durante'
durante.meta$pre.or.post.treatment = NA
durante.meta$cancer.type = 'melanoma'
durante.meta$cancer.subtype = 'uveal melanoma'
durante.meta$patient = paste0('DU', formatC(1:nrow(durante.meta), digits = 1, flag = "0"))

head(y)
head(durante.meta)

durante.meta$sample.id %in% y$sample.id
y = left_join(y, durante.meta, by = c('study', 'sample.id'), keep = F)


y %>%
    filter(study == 'durante') %>%
    head(n = 20)

############################
#INTEGRATE ANTUNES META 
############################

y$sample.id[y$study == 'antunes'] = 
    y %>%
        filter(study == 'antunes') %>%
        pull(cellname) %>%
        gsub('antunes-[A-Z]*-(.*)', '\\1', .) 

antunes.meta = y %>%
    filter(study == 'antunes') %>%
    pull(cellname) %>%
    gsub('antunes-[A-Z]*-(.*)', '\\1', .) %>%
    unique

antunes.meta = data.frame(
    sample.id = antunes.meta
)

antunes.meta$patient =  paste0('AN', formatC(1:nrow(antunes.meta), digits = 1, flag = "0"))

y %>%
    filter(study == 'antunes') %>%
    head(n = 10)

y$cancer.type[y$study == 'antunes'] = 'GBM'
y$cancer.subtype[y$study == 'antunes'] = 'GBM'
y$prim.met.other[y$study == 'antunes'] = 'primary'

y$patient[y$study == 'antunes'] = 
    antunes.meta$patient[
        match(
            y$sample.id[y$study == 'antunes'],
            antunes.meta$sample.id
        )
        ]

############################
#INTEGRATE ZHANG BREAST META 
############################
gsub('zhang.breast-[A-Z]*\\.(.*)', '\\1', y$cellname[y$study == 'zhang.breast'])

y$sample.id[y$study == 'zhang.breast'] = 
    gsub('zhang.breast-[A-Z]*\\.(.*)', '\\1', y$cellname[y$study == 'zhang.breast']) 

y$patient[y$study == 'zhang.breast'] = 
    gsub('zhang.breast-[A-Z]*\\.(.*)', '\\1', y$cellname[y$study == 'zhang.breast']) %>%
        gsub('[A-Za-z]*_(.*)_t', '\\1', .)

y$pre.or.post.treatment[y$study == 'zhang.breast'] = 
    gsub('zhang.breast-[A-Z]*\\.(.*)', '\\1', y$cellname[y$study == 'zhang.breast']) %>%
        gsub('([A-Za-z]*)_.*', '\\1', .)

y$cancer.type[y$study == 'zhang.breast'] = 'breast'
y$cancer.subtype[y$study == 'zhang.breast'] = 'TNBC'

zhang.meta = read.csv('~/work/ucl/bigdata/ucl.projects/macrophage/zhang.breast/zhang.breast.cell.annotation.csv')

zhang.meta = zhang.meta %>%
    filter(Tissue != 'blood')
colnames(zhang.meta)[1] = 'cellid'
zhang.meta$cellid = paste0('zhang.breast-', zhang.meta$cellid)

primary.cellids = zhang.meta %>%
    filter(Tissue == 'breast') %>%
    pull(cellid)

met.cellids = zhang.meta %>%
    filter(Tissue != 'breast') %>%
    pull(cellid)

y$prim.met.other[y$study == 'zhang.breast' & y$cellname %in% primary.cellids] = 'primary'
y$prim.met.other[y$study == 'zhang.breast' & y$cellname %in% met.cellids] = 'met'

#any(y$cellname %in% met.cellids)
#met.cellids

#y %>%
    #filter(study == 'zhang.breast' & is.na(prim.met.other))

#y %>%
    #filter(study == 'zhang.breast') %>%
    #pull(prim.met.other) 

############################
#INTEGRATE ZHENG META 
############################
y$sample.id[y$study == 'zheng'] = gsub('zheng-(.*)-.*-.*', '\\1', y$cellname[y$study == 'zheng'])

head(y[y$study == 'zheng', ])

zheng.df = data.frame(
    sample.id = unique(gsub('zheng-(.*)-.*-.*', '\\1', y$cellname[y$study == 'zheng'])),
    patient = paste0('ZHENG', formatC(1:4, digits = 1, flag = "0"))
)

y$patient[y$study == 'zheng'] = 
    zheng.df$patient[match(y$sample.id[y$study == 'zheng'], zheng.df$sample.id)]


y$cancer.type[y$study == 'zheng'] = 'CRC'
y$cancer.subtype[y$study == 'zheng'] = 'CRC'
y$prim.met.other[y$study == 'zheng'] = 'primary'

y[y$study == 'zheng', ]

############################
#INTEGRATE KHALIQ METADATA 
############################

head(y[y$study == 'khaliq', ])
y$sample.id[y$study == 'khaliq'] = gsub('khaliq-(.*_.*)_.*', '\\1', y$cellname[y$study == 'khaliq'])
y$cancer.type[y$study == 'khaliq'] = 'CRC'
y$cancer.subtype[y$study == 'khaliq'] = 'CRC'
y$cancer.subtype[y$study == 'khaliq'] = 'CRC'
y$prim.met.other[y$study == 'khaliq'] = 'primary'
y$patient[y$study == 'khaliq'] = paste0('KHA', gsub('khaliq-(.*_.*)_.*', '\\1', y$cellname[y$study == 'khaliq']))

############################
#ZHANG.OVARIAN 
############################

head(y[y$study == 'zhang.ovar', ])
y$patient[y$study == 'zhang.ovar'] = gsub('zhang.ovar-[A-Z]*\\.(.*)_.*', '\\1', y$cellname[y$study == 'zhang.ovar'])

zo.meta = read.tsv('~/work/ucl/bigdata/ucl.projects/macrophage/zhang.ovar/GSE165897/suppl/GSE165897_cellInfo_HGSOC.tsv')

head(zo.meta)
unique(zo.meta$treatment_phase)
zo.meta$prim.met.other = ""
zo.meta$prim.met.other[zo.meta$treatment_phase == 'post-NACT'] = 'metastatic'
zo.meta$prim.met.other[zo.meta$treatment_phase == 'treatment-naive'] = 'primary/laparoscopy'

zo.meta = zo.meta[match(gsub('\\.', '-', gsub('zhang.ovar-', '', y$cellname[y$study == 'zhang.ovar'])), zo.meta$cell), ]

y$sample.id[y$study == 'zhang.ovar'] = zo.meta[match(gsub('\\.', '-', gsub('zhang.ovar-', '', y$cellname[y$study == 'zhang.ovar'])), zo.meta$cell), ] %>%
    pull(sample)

y[y$study == 'zhang.ovar', ]$pre.or.post.treatment = zo.meta$treatment_phase
y$cancer.type[y$study == 'zhang.ovar'] = 'HGSOC'
y$cancer.subtype[y$study == 'zhang.ovar'] = 'HGSOC'
y$prim.met.other[y$study == 'zhang.ovar'] = zo.meta$prim.met.other

############################
#PU 
############################

head(y[y$study == 'pu', ])

y$sample.id[y$study == 'pu'] = 
    gsub('(pu-s[0-9]*)-[A-Z]*-.*', '\\1', y$cellname[y$study == 'pu'])

y$patient[y$study == 'pu'] = 
    y$sample.id[y$study == 'pu']

y$cancer.type[y$study == 'pu'] = 'thyroid'
y$cancer.subtype[y$study == 'pu'] = 'thyroid'
y$prim.met.other[y$study == 'pu'] = 'primary'

y$prim.met.other[y$study == 'pu' & y$sample.id == 'pu-s4'] = 'metastasis'
y$prim.met.other[y$study == 'pu' & y$sample.id == 'pu-s9'] = 'metastasis'

############################
#ZHANG ESO 
############################

head(y[y$study == 'zhang.eso', ])
y$sample.id[y$study == 'zhang.eso'] = gsub('zhang.eso-(.*?\\..*)\\..*', '\\1', y$cellname[y$study == 'zhang.eso'])
y$patient[y$study == 'zhang.eso'] = y$sample[y$study == 'zhang.eso']
y$prim.met.other[y$study == 'zhang.eso'] = 'primary'
y$prim.met.other[y$study == 'zhang.eso' & grepl('N', y$sample)] = 'normal'
y$cancer.type[y$study == 'zhang.eso'] = 'ESO'
y$cancer.subtype[y$study == 'zhang.eso'] = 'ESO'

y[y$study == 'zhang.eso', ]$cancer.type


############################
#XU 
############################

head(y[y$study == 'xu', ])
y$cancer.type[y$study == 'xu'] = 'breast'
y$cancer.subtype[y$study == 'xu'] = 'breast'
y$prim.met.other[y$study == 'xu'] = 'primary'
y$sample.id[y$study == 'xu'] = gsub('(xu-s[0-9])-.*', '\\1', y$cellname[y$study == 'xu'])
y$patient[y$study == 'xu'] = y$sample.id[y$study == 'xu']

############################
#SHARMA 
############################

head(y[y$study == 'sharma', ])
y$sample.id[y$study == 'sharma'] = gsub('sharma-[A-Z]*-(.*-.*)', '\\1', y$cellname[y$study == 'sharma'])
y$patient[y$study == 'sharma'] = 
    paste0('SHARMA_', gsub('sharma-[A-Z]*-(.*)-.*', '\\1', y$cellname[y$study == 'sharma']))

y$cancer.type[y$study == 'sharma'] = 'liver'
y$cancer.subtype[y$study == 'sharma'] = 'liver'
y$prim.met.other[y$study == 'sharma'] = 'primary'

############################
#LU 
############################

head(y[y$study == 'lu', ])
y$sample.id[y$study == 'lu'] = gsub('(lu-.*)_.*', '\\1', y$cellname[y$study == 'lu'])

y$patient[y$study == 'lu'] = 
    y$sample.id[y$study == 'lu']

y$prim.met.other[y$study == 'lu'] = 'primary'
y$cancer.type[y$study == 'lu'] = 'liver'
y$cancer.subtype[y$study == 'lu'] = 'liver'

############################
#BECKER 
############################

head(y[y$study == 'becker', ])
y$sample.id[y$study == 'becker'] = gsub('(becker-s[0-9]*)-.*-.*', '\\1', y$cellname[y$study == 'becker'])
y$prim.met.other[y$study == 'becker'] = 'primary'
y$cancer.type[y$study == 'becker'] = 'CRC'
y$cancer.subtype[y$study == 'becker'] = 'CRC'

y$patient[y$study == 'becker'] = y$sample.id[y$study == 'becker']

############################
#CHE 
############################

head(y[y$study == 'che', ])
y$sample.id[y$study == 'che'] = gsub('che-[A-Z]*_(.*_.*)', '\\1', y$cellname[y$study == 'che'])

y$prim.met.other[y$study == 'che' & grepl('CRC', y$sample.id)] = 'primary'
y$prim.met.other[y$study == 'che' & grepl('LM', y$sample.id)] = 'metastasis'

y$prim.met.other[y$study == 'che']

y$patient[y$study == 'che'] = gsub('che-[A-Z]*_(.*_.*)', '\\1', y$cellname[y$study == 'che']) %>%
    gsub('(.*)_(.*)', '\\1', .) 

y$cancer.type[y$study == 'che'] = 'CRC'
y$cancer.subtype[y$study == 'che'] = 'CRC'

############################
#BIERMANN 
############################

head(y[y$study == 'biermann', ])

biermann.meta = read.tsv('~/work/ucl/data/ucl.projects/macrophage/new.study.metadata/biermann.meta.tsv')
head(biermann.meta)
biermann.meta$sample.id = gsub('(.*_.*)_.*', '\\1', biermann.meta$orig.id)

y$sample.id[y$study == 'biermann'] = biermann.meta$sample.id[
    match(gsub('biermann-.*\\.(.*_.*)', '\\1', y$cellname[y$study == 'biermann']), biermann.meta$id1)
    ]

y$prim.met.other[y$study == 'biermann'] = 'metastasis'
y$cancer.type[y$study == 'biermann'] = 'Melanoma'
y$cancer.subtype[y$study == 'biermann'] = 'Melanoma'

y$patient[y$study == 'biermann'] = gsub('_sn|_sc', '', y$sample.id[y$study == 'biermann'])

############################
#WU 
############################

head(y[y$study == 'wu', ])
y$sample.id[y$study == 'wu'] = gsub('wu-(.*)_.*', '\\1', y$cellname[y$study == 'wu'])
y$patient[y$study == 'wu'] = y$sample.id[y$study == 'wu']
y$cancer.type[y$study == 'wu'] = 'breast'
y$cancer.subtype[y$study == 'wu'] = 'breast'
y$prim.met.other[y$study == 'wu'] = 'primary'

y = y[, !colnames(y) == 'sample']


x@meta.data = y

############################
#JUN DATA 
############################

x
#data/ucl.projects/macrophage/saveRDS(x, '~/work/ucl/data/ucl.projects/macrophage/new.macrophage.seurat.comb.rds')

jun.dat = readRDS('~/work/ucl/data/ucl.projects/macrophage/integrated_s17_macro_220606.rds')
jun.meta = jun.dat@meta.data

jun.meta

jun.meta %>%
    filter(is.na(primary_met2)) %>%
    pull(study) %>%
    unique

length(which(is.na(jun.meta$primary_met2)))
length(which(is.na(jun.meta$primary_met)))
jun.meta$primary_met
nrow(jun.meta)
head(jun.meta)
jun.meta$primary_met2
colnames(jun.meta)
jun.meta[!is.na(jun.meta$MetastasisStatus), ]
unique(jun.meta$MetastasisStatus)


jun.dat.raw = CreateSeuratObject(
    GetAssayData(jun.dat, slot = 'count'),
    meta.data = jun.dat@meta.data,
    project = 'jun'
)

head(jun.dat.raw@meta.data)

############################
#MORE GENE NAME CORRECTION 
############################

jun.dat.raw2 = GetAssayData(jun.dat.raw)

length(which(rownames(jun.dat.raw2) %in% rownames(x))) 
length(which(!rownames(jun.dat.raw2) %in% rownames(x))) 
length(which(!rownames(x) %in% rownames(jun.dat.raw2))) 

ens.ids = unique(rownames(jun.dat.raw2)[grepl('ENSG', rownames(jun.dat.raw2))])
ens.mapped = mapIds(org.Hs.eg.db, keys = ens.ids, keytype = 'ENSEMBL', column = 'SYMBOL', multiVals = 'first')
ens.mapped = ens.mapped[!is.na(ens.mapped)]
ens.mapped = as.data.frame(ens.mapped) %>%
    mutate(ens.id = rownames(.)) %>%
    reset.rownames

rownames(jun.dat.raw2)[match(ens.mapped$ens.id, rownames(jun.dat.raw2))] = 
    ens.mapped$ens.mapped

gene.syms.jun = checkGeneSymbols(rownames(jun.dat.raw2), map = current.map)
head(gene.syms.jun)

to.correct.jun = gene.syms.jun %>%
    filter(x != Suggested.Symbol & !is.na(Suggested.Symbol)) 

rownames(jun.dat.raw2)[match(to.correct.jun$x, rownames(jun.dat.raw2))] = 
    to.correct.jun$Suggested.Symbol

jun.dat.raw = CreateSeuratObject(
    jun.dat.raw2,
    meta.data = jun.dat@meta.data,
    project = 'jun'
)

############################
#CONTINUE 
############################

#saveRDS(jun.dat.raw, '~/work/ucl/data/ucl.projects/macrophage/jun.dat.comb.raw.rds')
jun.dat.raw = readRDS('~/work/ucl/data/ucl.projects/macrophage/jun.dat.comb.raw.rds')

############################
#FIX JUN MET INFORMATION 
############################

jun.meta = jun.dat.raw@meta.data

colnames(jun.meta)
g = jun.meta[c(colnames(jun.meta)[grep('met|Met', colnames(jun.meta))], 'study')]

g2 = g %>%
    group_by(study) %>%
    summarise(
        num.all.na = sum(is.na(MetastasisStatus) &
            is.na(primary_met) &
            is.na(primary_met2)),
        total = n()
    )


############################
#QIAN 
############################

jun.meta$sample[jun.meta$study == 'QIAN'] = 
    paste0(
        jun.meta$cancer[jun.meta$study == 'QIAN'],
        '_',
        jun.meta$patient[jun.meta$study == 'QIAN']
    )

jun.meta$primary_met2[
    jun.meta$study == 'QIAN' &
        jun.meta$tissue == 'Tumor' &
        !grepl('^OV_', jun.meta$sample)
] = 'primary'

#not possible to distinguish which of the OV samples for QIAN are mets,
#will have to leave as NA


############################
#AZIZI 
############################

jun.meta %>%
    filter(study == 'AZIZI') %>%
    pull(sample) %>%
    unique


jun.meta$primary_met2[jun.meta$study == 'AZIZI']
azizi.info = gsub('.*?_(.*)', '\\1', jun.meta$sample[jun.meta$study == 'AZIZI'])
azizi.info[azizi.info == 'TUMOR'] = 'primary'
azizi.info = tolower(azizi.info)
jun.meta$primary_met2[jun.meta$study == 'AZIZI'] = azizi.info

############################
#BI 
############################

jun.meta %>%
    filter(study == 'BI') %>%
    pull(sample) %>%
    unique

bi.info = gsub('(.*?)_.*', '\\1', jun.meta$sample[jun.meta$study == 'BI'])

bi.samps = read.csv('~/work/ucl/data/ucl.projects/macrophage/jun.study.meta/bi.samps.csv', header = F)
colnames(bi.samps) = c('samp', 'site', 'met.status')

bi.info2 = data.frame(samp = bi.info)
bi.info3 = left_join(bi.info2, bi.samps)

jun.meta$primary_met2[jun.meta$study == 'BI'] = bi.info3$met.status

############################
#BRAUN 
############################

jun.meta %>%
    filter(study == 'BRAUN') %>%
    pull(sample) %>%
    unique

braun.missing = jun.meta %>%
    filter(study == 'BRAUN') %>%
    dplyr::select(study, sample, primary_met, primary_met2, MetastasisStatus) %>%
    filter(
        is.na(MetastasisStatus) &
            is.na(primary_met) &
            is.na(primary_met2)
    ) %>%
    pull(sample) 


cat(braun.missing)

df1 = data.frame(
    sample = c(
        "S10_T2",
        "S11_T",
        "S12_T",
        "S15_T",
        "S16_T",
        "S5_T",
        "S7_T"
    ),
    primary_met2 = c(
        "primary",
        "primary",
        "primary",
        "primary",
        "primary",
        "metastasis",
        "primary"
    )
)

df2 = left_join(
    data.frame(sample = jun.meta$sample[jun.meta$sample %in% df1$sample]),
    df1
)

jun.meta$primary_met2[jun.meta$sample %in% df1$sample] = df2$primary_met2

############################
#CHAN 
############################

#some of the sample IDs for chan are not present in their metadata
#for samples. I'm just going to have to leave these as NA

jun.meta %>%
    filter(study == 'CHAN') %>%
    pull(sample) %>%
    unique

jun.meta$sample[jun.meta$study == 'CHAN'] = jun.meta$patient[jun.meta$study == 'CHAN']

library(readxl)
chan.meta = read_excel('~/work/ucl/data/ucl.projects/macrophage/jun.study.meta/chan.meta.xlsx', skip = 1)
colnames(chan.meta)[1] = 'sample'

chan.info = as.data.frame(chan.meta[c('sample', 'Tissue Type')])

jun.meta$sample[jun.meta$study == 'CHAN'][
    jun.meta$sample[jun.meta$study == 'CHAN'] %in% chan.info$sample
    ]

chan.info2 = chan.info$`Tissue Type`[match(
    jun.meta$sample[jun.meta$study == 'CHAN'],
    chan.info$sample
)]

chan.info2[chan.info2 == 'Recurrence'] = 'Primary'

jun.meta$primary_met2[jun.meta$study == 'CHAN'] = chan.info2

############################
#CHENG 
############################

#pan-cancer atlas without much meta data, might be difficult

############################
#JERBY 
############################

jun.meta %>%
    filter(study == 'JERBY') %>%
    pull(sample) %>%
    unique

jerby.missing = jun.meta %>%
    filter(study == 'JERBY') %>%
    dplyr::select(study, sample, primary_met, primary_met2, MetastasisStatus) %>%
    filter(
        is.na(MetastasisStatus) &
            is.na(primary_met) &
            is.na(primary_met2)
    ) %>%
    pull(sample) %>%
    unique

jerby.missing

jun.meta$primary_met2[jun.meta$study == 'JERBY' & jun.meta$sample == 'Mel129pa'] = 'primary'
jun.meta$primary_met2[jun.meta$study == 'JERBY' & jun.meta$sample == 'Mel129pb'] = 'primary'

############################
#KIM 
############################

jun.meta %>%
    filter(study == 'KIM') %>%
    pull(sample) %>%
    unique

kim.missing = jun.meta %>%
    filter(study == 'KIM') %>%
    dplyr::select(study, sample, primary_met, primary_met2, MetastasisStatus) %>%
    filter(
        is.na(MetastasisStatus) &
            is.na(primary_met) &
            is.na(primary_met2)
    ) %>%
    pull(sample) %>%
    unique

sort(kim.missing)

jun.meta$primary_met2[jun.meta$study == 'KIM'][
    grepl('^NS', jun.meta$sample[jun.meta$study == 'KIM'])
] = 'met'

jun.meta$primary_met2[jun.meta$study == 'KIM'][
    grepl('^EFFUSION', jun.meta$sample[jun.meta$study == 'KIM'])
] = 'effusion'

jun.meta$primary_met2[jun.meta$study == 'KIM'][
    grepl('^LUNG_T', jun.meta$sample[jun.meta$study == 'KIM'])
] = 'primary'

jun.meta$primary_met2[jun.meta$study == 'KIM'][
    grepl('BRONCHO_58', jun.meta$sample[jun.meta$study == 'KIM'])
] = 'primary'

jun.meta$primary_met2[jun.meta$study == 'KIM'][
    grepl('BRONCHO_11', jun.meta$sample[jun.meta$study == 'KIM'])
] = 'met'

jun.meta$primary_met2[jun.meta$study == 'KIM'][
    grepl('BRONCHO_11', jun.meta$sample[jun.meta$study == 'KIM'])
] 

jun.meta$primary_met2[jun.meta$study == 'KIM'][
    grepl('EBUS_10|EBUS_12|EBUS_13|EBUS_15|EBUS_19|EBUS_51', jun.meta$sample[jun.meta$study == 'KIM'])
] = 'met' 

jun.meta$primary_met2[jun.meta$study == 'KIM'][
    grepl('EBUS_06|EBUS_28|EBUS_49', jun.meta$sample[jun.meta$study == 'KIM'])
] = 'primary'

############################
#KRISHNA 
############################

jun.meta %>%
    filter(study == 'KRISHNA') %>%
    pull(sample) %>%
    unique

krishna.missing = jun.meta %>%
    filter(study == 'KRISHNA') %>%
    dplyr::select(study, sample, primary_met, primary_met2, MetastasisStatus) %>%
    filter(
        is.na(MetastasisStatus) &
            is.na(primary_met) &
            is.na(primary_met2)
    ) %>%
    pull(sample) %>%
    unique

cat(krishna.missing)

krishna.meta = read.tsv('~/work/ucl/data/ucl.projects/macrophage/jun.study.meta/krishna.meta.tsv')
colnames(krishna.meta)
krishna.meta = distinct(krishna.meta[c('type', 'region', 'Sample', 'Sample2')])
head(krishna.meta)

krishna.meta$primary_met2 = krishna.meta$type

krishna.meta$primary_met2[krishna.meta$type == 'LymphNode'] = 'met'
krishna.meta$primary_met2[krishna.meta$type == 'Tumor'] = 'primary'
krishna.meta$primary_met2[krishna.meta$type == 'Normal'] = 'normal'
krishna.meta$primary_met2[krishna.meta$type == 'PBMC'] = 'blood'

jun.meta$primary_met2[jun.meta$study == 'KRISHNA'] = krishna.meta$primary_met2[
    match(jun.meta$sample[jun.meta$study == 'KRISHNA'], krishna.meta$Sample2)
    ]

############################
#LEADER 
############################

jun.meta %>%
    filter(study == 'LEADER') %>%
    pull(patient) %>%
    unique

leader.missing = jun.meta %>%
    filter(study == 'LEADER') %>%
    dplyr::select(study, sample, primary_met, primary_met2, MetastasisStatus) %>%
    filter(
        is.na(MetastasisStatus) &
            is.na(primary_met) &
            is.na(primary_met2)
    ) %>%
    pull(sample) %>%
    unique

jun.meta$sample[jun.meta$study == 'LEADER'] =
    paste0(
        jun.meta$patient[jun.meta$study == 'LEADER'], '_',
        jun.meta$tissue[jun.meta$study == 'LEADER']
    )

jun.meta$primary_met2[jun.meta$study == 'LEADER' & jun.meta$tissue != 'Normal'] =
    'primary'

############################
#LI 
############################

jun.meta %>%
    filter(study == 'LI') %>%
    pull(patient) %>%
    unique

li.missing = jun.meta %>%
    filter(study == 'LI') %>%
    dplyr::select(study, sample, primary_met, primary_met2, MetastasisStatus) %>%
    filter(
        is.na(MetastasisStatus) &
            is.na(primary_met) &
            is.na(primary_met2)
    ) %>%
    pull(sample) %>%
    unique


df1 = data.frame(
    missing.samp = sort(li.missing)
)

head(df1)
df1$patient = gsub('(.*?)-.*', '\\1', df1$missing.samp)
head(df1)

li.meta = read_excel('~/work/ucl/data/ucl.projects/macrophage/jun.study.meta/li.meta.xlsx', skip = 2)
li.patient = gsub('(.*)-.*', '\\1', li.meta$Patient_ID)
li.meta$sample = li.meta$Patient_ID
li.meta$Patient_ID = li.patient

head(df1)
head(li.meta)
li.meta = li.meta[c('Patient_ID', 'sample', 'Tumor location')]
as.data.frame(li.meta)

df1$primary_met2 = ""
df1$primary_met2[grepl('^p25|^p26', df1$missing.samp)] = 'primary'
df1$missing.samp[grepl('PB', df1$missing.samp)]

df1$primary_met2[df1$primary_met2 == ""] = 'met'

jun.meta$primary_met2[jun.meta$study == 'LI'] = 
    df1$primary_met2[match(
        jun.meta$sample[jun.meta$study == 'LI'],
        df1$missing.samp
    )]

############################
#VISH 
############################

jun.meta %>%
    filter(study == 'VISH') %>%
    pull(patient) %>%
    unique

vish.missing = jun.meta %>%
    filter(study == 'VISH') %>%
    dplyr::select(study, sample, primary_met, primary_met2, MetastasisStatus) %>%
    filter(
        is.na(MetastasisStatus) &
            is.na(primary_met) &
            is.na(primary_met2)
    ) %>%
    pull(sample) %>%
    unique


jun.meta$primary_met2[jun.meta$study == 'VISH' & jun.meta$tissue == 'Tumor'] = 'primary'
jun.meta$primary_met2[jun.meta$study == 'VISH' & jun.meta$tissue == 'Blood'] = 'blood'

############################
#WU 
############################

jun.meta %>%
    filter(study == 'WU') %>%
    pull(patient) %>%
    unique

wu.missing = jun.meta %>%
    filter(study == 'WU') %>%
    dplyr::select(study, sample, primary_met, primary_met2, MetastasisStatus) %>%
    filter(
        is.na(MetastasisStatus) &
            is.na(primary_met) &
            is.na(primary_met2)
    ) %>%
    pull(sample) %>%
    unique

jun.meta$primary_met2[jun.meta$study == 'WU'] = 'primary'

############################
#YOST 
############################

jun.meta %>%
    filter(study == 'YOST_BASAL') %>%
    pull(patient) %>%
    unique

yost.missing = jun.meta %>%
    filter(study == 'YOST_BASAL') %>%
    dplyr::select(study, sample, primary_met, primary_met2, MetastasisStatus) %>%
    filter(
        is.na(MetastasisStatus) &
            is.na(primary_met) &
            is.na(primary_met2)
    ) %>%
    pull(sample) %>%
    unique

jun.meta$tissue[jun.meta$study == 'YOST_BASAL']
jun.meta$primary_met2[jun.meta$study == 'YOST_BASAL'] = 'primary'

############################
#ZILI 
############################

unique(jun.meta$tissue[jun.meta$study == 'ZILI'])

jun.meta$sample[jun.meta$study == 'ZILI'] = 
    jun.meta$patient[jun.meta$study == 'ZILI']


jun.meta %>%
    filter(study == 'ZILI') %>%
    pull(patient) %>%
    unique

zili.missing = jun.meta %>%
    filter(study == 'ZILI') %>%
    dplyr::select(study, sample, primary_met, primary_met2, MetastasisStatus) %>%
    filter(
        is.na(MetastasisStatus) &
            is.na(primary_met) &
            is.na(primary_met2)
    ) %>%
    pull(sample) %>%
    unique

unique(jun.meta$tissue[jun.meta$study == 'ZILI'])
jun.meta$primary_met2[jun.meta$study == 'ZILI'] = 'primary'


############################
#CONTINUE 
############################

jun.dat@meta.data = jun.meta

#mac1 = readRDS('~/work/ucl/data/ucl.projects/macrophage/new.macrophage.seurat.comb.rds')

mac2 = RowMergeSparseMatrices(
    GetAssayData(jun.dat.raw),
    GetAssayData(x)
)

mac2 = CreateSeuratObject(mac2, project = 'macrophages')

head(jun.dat@meta.data)

jun.meta = jun.dat@meta.data[c(
    'nCount_RNA',
    'nFeature_RNA',
    'study',
    'cellid',
    'celltype',
    'cellgroup',
    'cellgroup2',
    'major',
    'minor',
    'patient',
    'sample',
    'cancer',
    'tissue',
    'cluster_id',
    'primary_met2'
    )]

head(jun.meta)

unique(jun.meta$cancer)
unique(jun.meta$tissue)
unique(jun.meta$primary_met2)

mac.meta = x@meta.data
mac.meta = mac.meta %>%
    rename(
        cellid = cellname,
        sample = sample.id,
        primary_met2 = prim.met.other,
        cancer = cancer.type
    )

mac.meta
unique(mac.meta$study[is.na(mac.meta$primary_met2)])

meta.comb = bind_rows(jun.meta, mac.meta)

head(meta.comb)
unique(meta.comb$cancer)
nunique(match(rownames(mac2@meta.data), meta.comb$cellid))
dim(mac2@meta.data)

mac2.meta.orig = mac2@meta.data

mac2@meta.data$nCount_RNA = meta.comb$nCount_RNA[match(rownames(mac2@meta.data), meta.comb$cellid)]
mac2@meta.data$nFeature_RNA = meta.comb$nFeature_RNA[match(rownames(mac2@meta.data), meta.comb$cellid)]
mac2@meta.data$study = meta.comb$study[match(rownames(mac2@meta.data), meta.comb$cellid)]
mac2@meta.data$cellid = meta.comb$cellid[match(rownames(mac2@meta.data), meta.comb$cellid)]
mac2@meta.data$celltype = meta.comb$celltype[match(rownames(mac2@meta.data), meta.comb$cellid)]
mac2@meta.data$cellgroup = meta.comb$cellgroup[match(rownames(mac2@meta.data), meta.comb$cellid)]
mac2@meta.data$cellgroup2 = meta.comb$cellgroup2[match(rownames(mac2@meta.data), meta.comb$cellid)]
mac2@meta.data$major = meta.comb$major[match(rownames(mac2@meta.data), meta.comb$cellid)]
mac2@meta.data$minor = meta.comb$minor[match(rownames(mac2@meta.data), meta.comb$cellid)]
mac2@meta.data$patient = meta.comb$patient[match(rownames(mac2@meta.data), meta.comb$cellid)]
mac2@meta.data$sample = meta.comb$sample[match(rownames(mac2@meta.data), meta.comb$cellid)]
mac2@meta.data$cancer = meta.comb$cancer[match(rownames(mac2@meta.data), meta.comb$cellid)]
mac2@meta.data$tissue = meta.comb$tissue[match(rownames(mac2@meta.data), meta.comb$cellid)]
mac2@meta.data$cluster_id = meta.comb$cluster_id[match(rownames(mac2@meta.data), meta.comb$cellid)]
mac2@meta.data$primary_met2 = meta.comb$primary_met2[match(rownames(mac2@meta.data), meta.comb$cellid)]
mac2@meta.data$orig.ident = meta.comb$orig.ident[match(rownames(mac2@meta.data), meta.comb$cellid)]
mac2@meta.data$pre.or.post.treatment = meta.comb$pre.or.post.treatment[match(rownames(mac2@meta.data), meta.comb$cellid)]
mac2@meta.data$cancer.subtype = meta.comb$cancer.subtype[match(rownames(mac2@meta.data), meta.comb$cellid)]


length(which(is.na(mac2@meta.data$cancer)))


head(mac2@meta.data)
table(mac2@meta.data$tissue)
sort(table(mac2$primary_met2))

#cells.to.rm = mac2@meta.data %>%
    #filter(tissue %in% c('Blood', 'Effusion', 'LymphNode', 'Normal')) %>%
    #pull(cellid)

#mac2 = mac2[, !colnames(mac2) %in% cells.to.rm]

#cells.to.rm = mac2@meta.data %>%
    #filter(primary_met2 == 'normal') %>%
    #pull(cellid)

#mac2 = mac2[, !colnames(mac2) %in% cells.to.rm]

#seems the metadata got corrupted downstream somewhere? not sure where though.
#cancer field is NA for some studies
#write.tsv(mac2@meta.data, '~/work/ucl/data/ucl.projects/macrophage/all.mac.comb.meta.only.tsv')

#head(mac2@meta.data)
#saveRDS(mac2, '~/work/ucl/data/ucl.projects/macrophage/all.mac.comb.w.meta.w.norm.rds')

#write.tsv(mac2@meta.data, '~/work/ucl/data/ucl.projects/macrophage/comb.meta.w.primary.met.status.tsv')

