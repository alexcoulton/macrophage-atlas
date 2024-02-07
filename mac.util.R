############################
#LIBRARIES 
############################

source('~/work/ucl/scripts/misc/functions.R')
#source('~/work/ucl/scripts/ucl.projects/macrophage/jun.scripts/integrate_functions_mod.R')
#library(ComplexHeatmap)
#library(randomcoloR)
#library(ggplot2)
#library(gridExtra)
#library(ggradar)
#library(org.Hs.eg.db)
#library(reactome.db)
#library(fgsea)
#library(reshape2)
library(Seurat)
#library(ggsci)
#library(rhdf5)
#library(Seurat)
#library(patchwork)
library(dplyr)
#library(Matrix)
#library(statmod)
#library(speckle)
#geo.mean = function(x) exp(mean(log(x))) 
#library(UCell)
#library(viridis)

############################
#FUNCTIONS 
############################

load.mac.clus = function(){
    mac.int.embeddings = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/all.dat.w.norm.rpca.integrated.embeddings.rds')
    mac.clus.tmp = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/all.dat.w.norm.w.cluster.meta.data.rds')
    prim.met.info = read.tsv('~/work/ucl/data/ucl.projects/macrophage/comb.meta.w.primary.met.status.tsv')
    all(mac.clus.tmp$cellid == prim.met.info$cellid)
    mac.clus.tmp$primary_met2 = prim.met.info$primary_met2
    mac.clus.tmp$cancer[mac.clus.tmp$cancer == 'thyroid'] = 'THCA'
    mac.clus.tmp$cancer[mac.clus.tmp$cancer == 'melanoma'] = 'MEL'
    mac.clus.tmp$cancer[mac.clus.tmp$cancer == 'Melanoma'] = 'MEL'
    mac.clus.tmp$cancer[mac.clus.tmp$cancer == 'breast'] = 'BC'
    mac.clus.tmp$cancer[mac.clus.tmp$cancer == 'liver'] = 'LIHC'
    mac.clus.tmp$cancer[mac.clus.tmp$cancer == 'ESO'] = 'ESCA'
    mac.clus.tmp$study[mac.clus.tmp$study == 'WU'] = 'wu.lung'
    mac.clus.tmp$study = tolower(mac.clus.tmp$study)
    mac.clus.tmp$cluster = mac.clus.tmp$integrated_snn_res.0.525
    mac.clus.tmp = bind_cols(mac.clus.tmp, mac.int.embeddings$umap)
    mac.clus.tmp$cluster = as.numeric(as.character(mac.clus.tmp$cluster))
    mac.clus.tmp$sample[is.na(mac.clus.tmp$sample)] = 
        mac.clus.tmp$patient[is.na(mac.clus.tmp$sample)]
    mac.clus.tmp$study.samp = paste0(mac.clus.tmp$study, '_', mac.clus.tmp$samp)
    #these NA values are from studies I downloaded. I already filtered these for tumour only
    mac.clus.tmp$tissue[is.na(mac.clus.tmp$tissue)] = 'Tumor'
    mac.clus.tmp$ss2 = paste0(mac.clus.tmp$study.samp, '_', mac.clus.tmp$tissue)

    mac.clus.tmp$patient[mac.clus.tmp$ss2 %in% c('jerby_Mel129pb_Tumor', 'jerby_Mel129pa_Tumor')] = 'Mel129'

    clus.name.df <<- data.frame(
        cluster = 0:23,
        short.label = c(
            '0_AlvMac',
            '1_MetM2Mac',
            '2_C3Mac',
            '3_ICIMac1',
            '4_ICIMac2',
            '5_StressMac',
            '6_SPP1AREGMac',
            '7_IFNMac',
            '8_IFNGMac',
            '9_AngioMac',
            '10_InflamMac',
            '11_MetalloMac',
            '12_MBMMac',
            '13_CalciumMac',
            '14_ProliMac',
            '15_LYZMac',
            '16_ECMHomeoMac',
            '17_IFNMac3',
            '18_ECMMac',
            '19_ClassMono',
            '20_TDoub',
            '21_HemeMac',
            '22_IFNMac4',
            '23_NA'
        )
    )
    mac.clus.tmp$primary_met3 = mac.clus.tmp$primary_met2
    mac.clus.tmp$primary_met2[mac.clus.tmp$primary_met2 == 'pri'] = 
        'primary'
    mac.clus.tmp$primary_met2[mac.clus.tmp$primary_met2 == 'primary/laparoscopy'] = 
        'primary'
    mac.clus.tmp$primary_met2[mac.clus.tmp$primary_met2 == 'norm'] = 
        'normal'
    mac.clus.tmp$primary_met2 = tolower(mac.clus.tmp$primary_met2)
    mac.clus.tmp$primary_met2[mac.clus.tmp$primary_met2 == 'met'] = 
        'metastasis'
    mac.clus.tmp$primary_met2[mac.clus.tmp$primary_met2 == 'metastatic'] = 
        'metastasis'


    mac.clus.tmp = inner_join(mac.clus.tmp, clus.name.df, by = 'cluster')
    mac.clus.tmp$short.label = factor(mac.clus.tmp$short.label, levels = clus.name.df$short.label)

    ############################
    #CHAN SAMPLE INFO 
    ############################

    #library(readxl)
    #chan.info = read_excel('~/work/ucl/data/ucl.projects/macrophage/chan.clinical.sample.info.xlsx', skip = 1)
    #colnames(chan.info)[1] = 'sample'
    #chan.info = chan.info[c('sample', 'Tissue Type')]
    #colnames(chan.info)[2] = 'primary_met3'

    #chan.info
    #head(as.data.frame(chan.info))
    #mac.clus.tmp[mac.clus.tmp$study == 'chan', ]

    #head(chan.info)

    #mac.clus %>%
        #filter(study == 'qian' & is.na(primary_met2)) %>%
        #pull(tissue) %>%
        #table


    #pelka metadata
    mac.clus.tmp$primary_met2[mac.clus.tmp$study == 'pelka' & is.na(mac.clus.tmp$primary_met2)] = 'primary'
    mac.clus.tmp$primary_met3[mac.clus.tmp$study == 'pelka' & is.na(mac.clus.tmp$primary_met3)] = 'primary'
    #qian metadata
    mac.clus.tmp$primary_met2[mac.clus.tmp$study == 'qian' & is.na(mac.clus.tmp$primary_met2)] = 'primary'
    mac.clus.tmp$primary_met3[mac.clus.tmp$study == 'qian' & is.na(mac.clus.tmp$primary_met3)] = 'primary'
    #maynard metadata
    mac.clus.tmp$primary_met2[mac.clus.tmp$study == 'maynard' & is.na(mac.clus.tmp$primary_met2)] = 'primary'
    mac.clus.tmp$primary_met3[mac.clus.tmp$study == 'maynard' & is.na(mac.clus.tmp$primary_met3)] = 'primary'


    #mac.clus %>%
        #group_by(study, primary_met2, .drop = F) %>%
        #summarise(n = n()) %>%
        #as.data.frame



    mac.clus <<- mac.clus.tmp
}







load.mac.int = function(){
    load.mac.clus()
    mac.int.embeddings = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/all.dat.w.norm.rpca.integrated.embeddings.rds')
    mac.int.tmp = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/all.dat.w.norm.rpca.integrated.rds')
    mac.int.tmp@meta.data = mac.clus
    Idents(mac.int.tmp) = 'short.label'

    mac.int.tmp[['umap']] = CreateDimReducObject(
        embeddings = mac.int.embeddings$umap,
        assay = 'integrated'
    )
    mac.int <<- mac.int.tmp
}

