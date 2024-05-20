library(fgsea)
library(DESeq2)

# main figures panels labelled PANEL
# supplementary figures labelled SUPP
source('~/work/ucl/scripts/ucl.projects/macrophage/mac.util.R')
source('~/work/ucl/scripts/misc/functions.R')

load.mac.clus()

#library(ConsensusTME)

############################
#DOCETAXEL 
############################

#cpi.counts.chemo = read.csv('~/work/ucl/data/ucl.projects/macrophage/cpi3000.tpm/counts_07122023.csv')
#cpi.chemo.meta = read.csv('~/work/ucl/data/ucl.projects/macrophage/cpi3000.tpm/CPI3000_METADATA_V1.csv')

#cpi.chemo.meta$patient = gsub('^PAT-', 'PAT.', cpi.chemo.meta$patient)
#cpi.chemo.meta$patient = gsub('^EA-', 'EA.', cpi.chemo.meta$patient)

#cpi.chemo.meta = cpi.chemo.meta[cpi.chemo.meta$treatment %in% c('Docetaxel', 'chemo'), ] 

#cpi.counts.chemo = cpi.counts.chemo[, colnames(cpi.counts.chemo) %in% cpi.chemo.meta$patient]

#cpi.chemo.meta = cpi.chemo.meta[match(colnames(cpi.counts.chemo), cpi.chemo.meta$patient), ] %>%
    #reset.rownames

#we've only got response data for 75 of the chemo treated patients
#probably not worth doing additional propeller analysis on them
#sum(cpi.chemo.meta$response_CR.PR.SD6 != "")

############################
#PROCESSING RAW COUNTS / DESEQ2 & FGSEA ANALYSIS
############################

cpi.tpm = read.csv('~/work/ucl/data/ucl.projects/macrophage/cpi3000.tpm/data.for.kevin/abundance_30112023.csv')
kim.tpm = read.csv('~/work/ucl/data/ucl.projects/macrophage/cpi3000.tpm/KIM_NMED_2018.tpm_genesymbols.csv')

cpi.counts = read.csv('~/work/ucl/data/ucl.projects/macrophage/cpi3000.tpm/counts_07122023.csv')
cpi.meta2 = read.csv('~/work/ucl/data/ucl.projects/macrophage/cpi3000.tpm/CPI3000_METADATA_V1.csv')
cpi.meta2 = cpi.meta2[cpi.meta2$CPI3000_inclusion != 'NO', ]

############################
#PROCESSING OF KIM COUNTS / COHORT 
############################
#https://www.nature.com/articles/s41591-018-0101-z#Sec27

more.cpi.counts = read.csv('~/work/ucl/data/ucl.projects/macrophage/cpi3000.tpm/data.for.alex/counts.csv')

patient.ids = readLines('~/work/ucl/data/ucl.projects/macrophage/cpi3000.tpm/data.for.alex/patient.ids.txt')
kim.patient.ids = gsub('.*salmon\\/(.*?)\\/.*', '\\1', patient.ids)
kim.patient.ids = gsub('_WTS$', '', kim.patient.ids)
kim.patient.ids = kim.patient.ids[-79]

kim1 = cpi.meta2[cpi.meta2$study == 'KIM_NMED_2018', ]
kim2 = cpi.meta2[cpi.meta2$study == 'KIM_NMED_2018', ]$patient

kim.meta2 = read.tsv('~/work/ucl/data/ucl.projects/macrophage/cpi3000.tpm/filereport_read_run_PRJEB25780_tsv.txt')

kimx = gsub('.*\\/(.*?).fastq.gz', '\\1', kim.meta2$submitted_ftp)
kim.meta2$krupa.id = gsub('(.*?_.*?_.*?_.*?)_.*', '\\1', kimx)
kim.meta3 = kim.meta2[c('krupa.id', 'sample_alias')]

kim.meta3 = kim.meta3[kim.meta3$krupa.id %in% kim.patient.ids, ] %>% reset.rownames
kim.meta3 = kim.meta3[match(kim.patient.ids, kim.meta3$krupa.id), ] %>% reset.rownames

kim.meta3$samp.trunc = gsub('_RNA$', '\\1', gsub('(.*?-.*?-.*?)-.*', '\\1', kim.meta3$sample_alias))

kim.patient.ids == kim.meta3$krupa.id
kim.meta3$samp.trunc %in% kim1$patient

kim.meta3$patient = kim.meta3$samp.trunc

kim.meta4 = left_join(kim.meta3, kim1[c('patient', 'study', 'tumor_type', 'response_CR.PR.SD6')])

normal.coords = grepl('N_RNA$', kim.meta4$sample_alias)

kim.patient.ids == kim.meta4$krupa.id

kim.meta5 = kim.meta4[!normal.coords, ]

more.cpi.counts = more.cpi.counts[, !normal.coords]
kim.tpm = kim.tpm[, !normal.coords]

colnames(more.cpi.counts) = kim.meta5$patient
colnames(kim.tpm) = kim.meta5$patient

rownames(cpi.tpm) = cpi.tpm$Gene
cpi.tpm = cpi.tpm[, -(1:2)]

############################
#CONTINUE MAIN PROCESSING 
############################

cpi.counts = bind_cols(cpi.counts, more.cpi.counts)
cpi.tpm = bind_cols(cpi.tpm, kim.tpm)

cpi.meta2$patient = gsub('^PAT-', 'PAT.', cpi.meta2$patient)
cpi.meta2$patient = gsub('^EA-', 'EA.', cpi.meta2$patient)
#add X at the start of pure number identifiers following R convention
cpi.meta2$patient = gsub('^([0-9].*?)', 'X\\1', cpi.meta2$patient)

colnames(cpi.counts) = gsub('_On$|_on$|_O$', '', colnames(cpi.counts))
colnames(cpi.tpm) = gsub('_On$|_on$|_O$', '', colnames(cpi.tpm))

#the following missing patients are either chemo treated or 
#pre treatment samples. (see Exclusion_Reason column of meta data)
colnames(cpi.counts)[!colnames(cpi.counts) %in% cpi.meta2$patient]
colnames(cpi.counts)[colnames(cpi.counts) %in% cpi.meta2$patient]

#remove duplicate patient IDs that are pre treatment
cpi.meta2 = cpi.meta2[!(grepl('pre_treatment', cpi.meta2$biopsy_time) & (cpi.meta2$patient %in% cpi.meta2$patient[duplicated(cpi.meta2$patient)])), ]

cpi.meta2 = cpi.meta2[cpi.meta2$patient %in% colnames(cpi.counts), ]
cpi.counts = cpi.counts[, colnames(cpi.counts) %in% cpi.meta2$patient]
cpi.counts = cpi.counts[, colnames(cpi.counts) %in% cpi.meta2$patient]

cpi.tpm = cpi.tpm[, colnames(cpi.tpm) %in% cpi.meta2$patient]
cpi.tpm = cpi.tpm[, colnames(cpi.tpm) %in% cpi.meta2$patient]

cpi.meta2 = cpi.meta2[match(colnames(cpi.counts), cpi.meta2$patient), ]

#some basic stats
#unique(cpi.meta2$study)
#distinct(cpi.meta2[c('study', 'tumor_type')])
#sort(table(cpi.meta2$tumor_type))

coldata = cpi.meta2[c('patient', 'response_CR.PR.SD6', 'tumor_type')]
colnames(coldata)[2] = 'response'

coldata$tumor_type[coldata$tumor_type == 'LUNG CANCER'] = 'LUNG'

coldata$response = factor(coldata$response, levels = c('no_response', 'response'))
coldata$tumor_type = factor(coldata$tumor_type)

#na.vals1 = is.na(coldata$response)

#cpi.counts = cpi.counts[, -which(na.vals1)]
cpi.counts = as.matrix(cpi.counts)
cpi.counts = round(cpi.counts)
#coldata = coldata[-which(na.vals1), ]

dds <- DESeqDataSetFromMatrix(
    countData = cpi.counts,
    colData = coldata[, 2:3],
    design= ~ response + tumor_type
)

dds = DESeq(dds)
res = results(dds)

pathways.hallmark = gmtPathways('~/work/ucl/data/bio.res/hallmark/h.all.v7.4.symbols.gmt')

res = as.data.frame(res)
res = res %>%
    dplyr::arrange(stat)

ranks1 = res$stat
names(ranks1) = rownames(res)

ranks1 = ranks1[!is.na(ranks1)]

data(examplePathways)
data(exampleRanks)

exampleRanks
examplePathways

cluster.markers = read.tsv('~/work/ucl/data/ucl.projects/macrophage/all.sig.markers.per.cluster.tsv')

head(cluster.markers)

cluster.markers = cluster.markers %>%
    filter(avg_log2FC > 0) %>%
    arrange(cluster, -avg_log2FC)

cluster.split = split(cluster.markers, cluster.markers$cluster)

head(cluster.split[[1]])

cluster.pathways = lapply(cluster.split, function(x){
    x$gene[1:50]
})

fgseaRes <- fgsea(
    pathways = cluster.pathways, 
    stats    = ranks1,
    minSize  = 15,
    maxSize  = 500
)

#saveRDS(fgseaRes, '~/work/ucl/data/ucl.projects/macrophage/cpi3000.tpm/fgseaRes.no.response.correction.rds')

#here i performed the analysis but without tumor type in the 
#DESeq2 formula to check what the effect would be on the significant clusters
fgseaRes.no.tumour.type = readRDS('~/work/ucl/data/ucl.projects/macrophage/cpi3000.tpm/fgseaRes.no.response.correction.rds')

#results from differential expression / pathway analysis controlling for tumour type 
#in DESeq2 formula
fgseaRes = readRDS('~/work/ucl/data/ucl.projects/macrophage/cpi3000.tpm/fgseaRes.rds')



fgseaRes.no.tumour.type
fgseaRes


fgseaRes = fgseaRes %>%
    filter(padj < 0.1)

fgseaRes.no.tumour.type = fgseaRes.no.tumour.type %>%
    filter(padj < 0.1)

#significant pathways that are removed after correcting for tumour type
fgseaRes.no.tumour.type$pathway[!fgseaRes.no.tumour.type$pathway %in% fgseaRes$pathway]
fgseaRes$pathway %in% fgseaRes.no.tumour.type$pathway

sort(unique(mac.clus$short.label))

fgseaRes = fgseaRes %>%
    arrange(pval)

bubbleplot.df = as.data.frame(fgseaRes)[c('pathway', 'padj', 'NES')] %>%
    filter(padj < 0.1)

gold.standard = c(5, 6, 8, 11, 17, 21, 22)

bubbleplot.df
bubbleplot.df$gold = F
bubbleplot.df$gold[bubbleplot.df$pathway %in% gold.standard] = T

cluster.names = distinct(mac.clus[c('cluster', 'short.label')]) %>%
    arrange(cluster) %>%
    rename(pathway = cluster) %>%
    mutate(pathway = as.character(pathway))

bubbleplot.df = left_join(bubbleplot.df, cluster.names)

bubbleplot = ggplot(bubbleplot.df, aes(x = NES, y = short.label, size = -log(padj), color = gold)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    coord_cartesian(xlim = c(-2.8, 2.8)) +
    scale_color_manual(values = c('black', 'orange')) +
    labs(color = 'Gold standard') +
    theme_bw()

#PANEL
pdf('~/work/ucl/plots/ucl.projects/macrophage/ici.figure.components/fgsea.bubble.pdf', 4, 3)
print(bubbleplot)
dev.off()









############################
#CPI3000 IMPORT 
############################

#tpm.new = read.csv('~/work/ucl/data/ucl.projects/macrophage/cpi3000.tpm/Total_Outputs/cpi3000_abundance_october2023.csv')
#cpi.meta = read.csv('~/work/ucl/data/ucl.projects/macrophage/cpi3000.tpm/CPI3000_METADATA_V1.csv')

#x = sort(colnames(tpm.new)[grep('^Pt', colnames(tpm.new))])
#these patients are from Riaz and have both on and pre treatment samples, remove the pre-treatment samples
#to.rm = gsub('_.*$', '', sort(x[duplicated(gsub('_.*$', '', x))]))
#tpm.new = tpm.new[, !colnames(tpm.new) %in% to.rm]

#colnames(tpm.new) = gsub('\\.', '_', colnames(tpm.new))
#colnames(tpm.new) = gsub('^X', '', colnames(tpm.new))
#colnames(tpm.new) = gsub('_On$|_on$|_O$', '', colnames(tpm.new))
#colnames(tpm.new) = gsub('PAT_67d14d57330p', 'PAT_67d14d57330', colnames(tpm.new))

#cpi.meta$patient = gsub('-', '_', cpi.meta$patient)

#rownames(tpm.new) = tpm.new[, 1]
#tpm.new = tpm.new[, -1]
#tpm.new = as.matrix(tpm.new)
#tpm.new = t(tpm.new)
#tpm.new = as.data.frame(tpm.new)
#tpm.new$patient = rownames(tpm.new)

#tpm.new = left_join(tpm.new, cpi.meta[c('patient', 'response_CR.PR.SD6', 'tumor_type', 'study', 'CPI3000_inclusion')], by = 'patient')

#tpm.new = tpm.new[tpm.new$CPI3000_inclusion != 'NO', ]
#tpm.new = tpm.new %>%
    #dplyr::rename(
        #resp = response_CR.PR.SD6,
        #case = patient,
        #cancer.type = tumor_type
    #)

#tpm.new = tpm.new[!is.na(tpm.new$resp), ]
#tpm.new.df = tpm.new

#colnames(tpm.new.df)[1000:2000]

#sort(table(tpm.new$cancer.type))

#tpm.df[tpm.df$case %in% tpm.new.df$case, ]$study

############################
#CPI1000 ANALYSIS - DATA IMPORT AND PREP
############################

#old CPI cohort code
#tpm = read.csv('~/work/ucl/data/ucl.projects/cpi/Combined_TPM_matrix_all_studies_n1003.csv')
#g2 = read.tsv('~/work/ucl/data/ucl.projects/cpi/meta_analysis_input_data.txt')
#g3 = read.csv('~/work/ucl/data/ucl.projects/cpi/snv_all_patients_v2_independentRun.csv')

#mva1 = read.tsv('~/work/ucl/data/ucl.projects/cpi/mva_data_AO.txt')
#mva2 = read.tsv('~/work/ucl/data/ucl.projects/cpi/mva_data_EVA.txt')
#t1 = read.tsv('~/work/ucl/data/ucl.projects/cpi/keynote_other_tumour_type_test_set.txt')
#t2 = read.tsv('~/work/ucl/data/ucl.projects/cpi/meta_analysis_input_data.txt')
#t3 = read.tsv('~/work/ucl/data/ucl.projects/cpi/new.tmb.calc.tsv')
#t4 = read.csv('~/work/ucl/data/ucl.projects/cpi/master_sample_sheet5.csv')
#t5 = read.tsv('~/work/ucl/data/ucl.projects/cpi/TMB_metrics_with_NMD_Sept2021.txt')
#master.sample.sheet = read.csv('~/work/ucl/data/ucl.projects/cpi/master_sample_sheet5.csv')

#ctdatabase = readLines('~/work/ucl/data/ucl.projects/inhibigens/list.of.cancer.testes.antigen.genes.ctdatabase.txt')

#tpm.mat = tpm[, 3:ncol(tpm)]
#tpm.mat = as.matrix(tpm.mat)
#rownames(tpm.mat) = tpm$gene_id
#tpm.mat = t(tpm.mat)

#resp.labels = master.sample.sheet$Response_responder[match(rownames(tpm.mat), master.sample.sheet$Patient)]
#tpm.mat = tpm.mat[!is.na(resp.labels), ]
#resp.labels = resp.labels[!is.na(resp.labels)]

#colnames(t4)[colnames(t4) == 'Patient'] = 'case'
#all.tmb.files = list(t2, mva1, mva2, t1, t3, t5)

#all.tmb2 = lapply(all.tmb.files, function(x) x[c('case', 'TMB')]) %>%
    #bind_rows %>%
    #reset.rownames

#all.tmb2 = all.tmb2[!duplicated(all.tmb2$case), ]
##all.tmb2 = all.tmb2[!is.na(all.tmb2$TMB), ]

#all.tmb2 = all.tmb2[match(rownames(tpm.mat), all.tmb2$case), ]
#all.tmb2$resp = resp.labels

#tpm.df = bind_cols(tpm.mat, all.tmb2)
#tpm.df = tpm.df[!is.na(tpm.df$resp), ] %>% reset.rownames

#tpm.df$case
#tpm.df$cancer.type = master.sample.sheet$Tum_Type[match(tpm.df$case, master.sample.sheet$Patient)]
#tpm.df$cancer.type2 = master.sample.sheet$Tum_Type2[match(tpm.df$case, master.sample.sheet$Patient)]
#tpm.df$study = master.sample.sheet$Study[match(tpm.df$case, master.sample.sheet$Patient)]

#unique(tpm.df$cancer.type)


############################
#COMBINE CPI1000 AND CPI3000 
############################

#The new TPM dataset that Krupa gave me does not include all of the 
#studies that were used previously. Liu and Van Allen are now grouped 
#into the melanoma genome sequencing project, whilst voorwerk
#we apparently do not have approval to use anymore.

#tpm2 = tpm.df[, colnames(tpm.df) %in% colnames(tpm.new.df)]
#tpmnew2 = tpm.new.df[, colnames(tpm.new.df) %in% colnames(tpm.df)]


#tpm2[, tpm2$case %in% tpmnew2$case]

#these data are accounted for in the new dataset
#tpm2 = tpm2[!tpm2$study %in% c(
    #'Liu_NatureMedicine_2019',
    #'VANALLEN_SCIENCE_2015',
    #'RIAZ_CELL_2017',
    #'MARIATHASAN_NATURE_2018',
    #'MCDERMOT_NMED_2018',
    #'HUGO_CELL_2016',
    #'VOORWERK_NMED_2019'
#), ]

#sort(unique(tpm2$study))
#sort(unique(tpmnew2$study))

#all(colnames(tpm2) == colnames(tpmnew2))

#tpm.comb = bind_rows(tpm2, tpmnew2)

############################
#CONSENSUSTME 
############################

#distinct(tpm.df[c('cancer.type', 'cancer.type2')])

#tpm.df %>%
    #group_by(cancer.type, cancer.type2) %>%
    #summarise(n = n()) %>%
    #arrange(-n)

#unique(tpm.df$cancer.type)
#tcga.acro = c(MELANOMA = 'SKCM', BLADDER = 'BLCA', RENAL = 'KIRC', BREAST = 'BRCA', GASTRIC = 'STAD', LUNG = 'LUAD')

#tpm.df$tcga = unname(tcga.acro[match(tpm.df$cancer.type, names(tcga.acro))])

#tpm.split = split(tpm.df, tpm.df$tcga)
#lapply(tpm.split, function(x){
    #x2 = x[, -(which(colnames(x) == 'case'):ncol(x))]
    #browser()
#})

############################
#CPI1000 - SIGNATURE ANALYSIS 
############################

#sig.markers = read.tsv('~/work/ucl/data/ucl.projects/macrophage/all.sig.markers.per.cluster.tsv')

#sig.markers = sig.markers %>%
    #filter(avg_log2FC > 0) %>%
    #arrange(cluster, -avg_log2FC)

#head(sig.markers)

#sig.markers.split = split(sig.markers, sig.markers$cluster)
#sig.genes1 = uulapply(sig.markers.split, function(x) x$gene[1:10])

#markers.to.rm = sig.genes1[sig.genes1 %in%
    #colnames(tpm.df)[uulapply(tpm.df, function(x) length(which(is.na(x))) > 0)]]

#sig.markers = sig.markers[!sig.markers$gene %in% markers.to.rm, ]

#tpm.comb$cancer.type[tpm.comb$cancer.type == 'BLADDER'] = 'Urothelial'
#tpm.comb$cancer.type[tpm.comb$cancer.type == 'LUNG CANCER'] = 'LUNG'
#tpm.comb = tpm.comb[!is.na(tpm.comb$cancer.type), ]

############################
#DESEQ2 ATTEMPT 
############################

#library(DESeq2)

#head(tpm.comb)
#tpm.comb = tpm.comb[!duplicated(tpm.comb$case), ]
#rownames(tpm.comb) = tpm.comb$case

#head(tpm.comb)
#tpm.mat1 = as.matrix(tpm.comb[, 1:(which(colnames(tpm.comb) == 'case') - 1)])
#tpm.mat1 = t(tpm.mat1)
#coldata = tpm.comb[c('resp', 'cancer.type')]

#head(coldata)
#coldata$resp = as.numeric(factor(coldata$resp, levels = c('no_response', 'response')))
#coldata$cancer.type = as.numeric(factor(coldata$cancer.type))

#dds <- DESeqDataSetFromMatrix(countData = tpm.mat1,
                              #colData = coldata,
                              #design= ~ cancer.type + resp)



#dds <- DESeq(dds)
#resultsNames(dds) # lists the coefficients
#res <- results(dds, name="condition_trt_vs_untrt")
## or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")

#fgseaRes

############################
#LOGISTIC REGRESSION TESTING 
############################

#all.log.tests = lapply(0:22, function(z){
    #sig1 = sig.markers %>%
        #filter(cluster == z) %>%
        #head(n = 10) %>%
        #pull(gene)
    #sig1 = sig1[sig1 %in% colnames(tpm.comb)]

    #log.test.df = bind_cols(
        #tpm.comb[c('case', 'resp', 'cancer.type')],
        #geo.mean = apply(tpm.comb[sig1], 1, geo.mean)
    #)

    #log.test.df$resp = factor(log.test.df$resp)
    #log.test.df$cancer.type = factor(log.test.df$cancer.type)

    #model1 = glm(resp ~ geo.mean + cancer.type, data = log.test.df, family = 'binomial')

    #summary(model1)
#})

#all.log.tests

#count = 0
#glm.results1 = lapply(all.log.tests, function(x){
    #g = unclass(x)$coefficients
    #df1 = data.frame(
        #sig = count,
        #p = g[2, 4]
    #)
    #count <<- count + 1
    #df1
#}) %>% bind_rows %>%
    #mutate(fdr = p.adjust(p, method = 'fdr'))

#glm.results1 %>%
    #arrange(p)

#glm.results1[glm.results1$sig %in% c(5, 21, 6, 8, 17, 11, 22), ]  %>%
    #mutate(fdr2 = p.adjust(p, method = 'fdr'))


#g = read.tsv('/nemo/project/proj-tracerx-lung/public_datasets/CPI/CPI1000_ANALYSIS/RNAseq_analysis/All_studies_expectedcounts_RNAseq.txt')

#head(g)
#dim(g)

#tpm2$case %in% colnames(g)

#all.t.test = lapply(0:22, function(z){
    #print(z)
    #cluster.sig = sig.markers %>%
        #filter(cluster == z) %>%
        #head(n = 10) %>%
        #pull(gene)

    #macro.sig = 
        #c('MARCO', 'CXCL5', 'SCG5', 'SULT1C2', 
        #'MSR1', 'CTSK', 'PTGDS', 'COLEC12', 
        #'GPC4', 'PCOLCE2', 'CHIT1', 'KAL1', 
        #'CLEC5A', 'ME1', 'DNASE2B', 'CCL7',
        #'CD163', 'FN1', 'GM2A', 'SCARB2', 'BCAT1',
        #'RAI14', 'COL8A2', 'APOE', 'CHI3L1', 
        #'ATG7', 'CD84', 'FDX1', 'MS4A4A', 'SGMS1',
        #'EMP1', 'CYBB', 'CD68')


    #macro.sig %in% colnames(tpm.df)
    #cluster.sig = cluster.sig[cluster.sig %in% colnames(tpm.df)]

    #tpm.df2 = tpm.df


    #geo.mean = function(x) exp(mean(log(x))) 
    ##some of the samples in the TPM data frame have NAs for a lot of genes, need to remove these 
    ##tpm.df2 = tpm.df2[-unique(uulapply(tpm.df2, function(x) which(is.na(x)))), ]

    #sig.ratio = apply(tpm.df2[cluster.sig], 1, geo.mean) / apply(tpm.df2[macro.sig], 1, geo.mean)
    #tpm.df2$sig.ratio = sig.ratio

    #tpm.df2[cluster.sig]
    #tpm.df2$cluster.sig = apply(tpm.df2[cluster.sig], 1, geo.mean)

    #tpm.df2 = tpm.df2[!sig.ratio == Inf, ]


    #############################
    ##SIG.RATIO T-TEST 
    #############################

    #tpm.df2.resp = tpm.df2[c('sig.ratio', 'resp')]

    #tpm.df2.resp %>%
        #group_by(resp) %>%
        #summarise(
            #mean = mean(sig.ratio),
            #sd = sd(sig.ratio)
        #)

    #head(tpm.df2.resp)

    ##ecm.cpi.plot = ggplot(tpm.df2.resp, aes(x = resp, y = sig.ratio)) +
        ##geom_jitter(width = 0.3)

    ##pl(ecm.cpi.plot)

    #resp.sigs = tpm.df2.resp %>%
        #filter(resp == 'response') %>%
        #pull(sig.ratio)

    #nresp.sigs = tpm.df2.resp %>%
        #filter(resp == 'no_response') %>%
        #pull(sig.ratio)

    #resp.sigs
    #nresp.sigs
    #colnames(tpm.df2)

    #tres = t.test(resp.sigs, nresp.sigs)
    #df1 = data.frame(
        #cluster = z,
        #p.val = tres$p.value,
        #num.removed = length(which(sig.ratio == Inf))
    #)

    #############################
    ##CLUSTER SIG T-TEST 
    #############################

    #tpm.df2 = tpm.df
    #geo.mean = function(x) exp(mean(log(x))) 
    ##some of the samples in the TPM data frame have NAs for a lot of genes, need to remove these 

    #sig.ratio = apply(tpm.df2[cluster.sig], 1, geo.mean) / apply(tpm.df2[macro.sig], 1, geo.mean)
    #tpm.df2$sig.ratio = sig.ratio

    #tpm.df2[cluster.sig]
    #tpm.df2$cluster.sig = apply(tpm.df2[cluster.sig], 1, geo.mean)


    #tpm.df2.resp = tpm.df2[c('cluster.sig', 'resp')]

    #tpm.df2.resp %>%
        #group_by(resp) %>%
        #summarise(
            #mean = mean(cluster.sig),
            #sd = sd(cluster.sig)
        #)

    #head(tpm.df2.resp)

    ##ecm.cpi.plot = ggplot(tpm.df2.resp, aes(x = resp, y = cluster.sig)) +
        ##geom_jitter(width = 0.3)

    ##pl(ecm.cpi.plot)

    #resp.sigs = tpm.df2.resp %>%
        #filter(resp == 'response') %>%
        #pull(cluster.sig)

    #nresp.sigs = tpm.df2.resp %>%
        #filter(resp == 'no_response') %>%
        #pull(cluster.sig)

    #tres = t.test(resp.sigs, nresp.sigs)
    #df2 = data.frame(
        #cluster = z,
        #p.val = tres$p.value,
        #num.removed = length(which(cluster.sig == Inf))
    #)

    #df1$cat = 'sig.ratio'
    #df2$cat = 'cluster.sig'

    #bind_rows(df1, df2)
#})

#all.t.test = all.t.test %>% bind_rows

#g1 = all.t.test %>%
    #filter(cat == 'cluster.sig') %>%
    #mutate(fdr = p.adjust(p.val, method = 'fdr'))

#g2 = all.t.test %>%
    #filter(cat == 'sig.ratio') %>%
    #mutate(fdr = p.adjust(p.val, method = 'fdr'))

#g3 = bind_rows(g1, g2)

#g3$cluster = as.factor(g3$cluster)
#sig.results = g3 %>% filter(fdr < 0.1)


############################
#CPI1000 PLOTTING 
############################

cpi.p.val.plot = ggplot(g3, aes(x = cluster, y = -log10(p.val), color = cat)) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed')

#pl(cpi.p.val.plot, 10, 5)

cpi.cluster.box.plots = lapply(0:22, function(z){
    print(z)
    cluster.sig = sig.markers %>%
        filter(cluster == z) %>%
        head(n = 10) %>%
        pull(gene)

    macro.sig = 
        c('MARCO', 'CXCL5', 'SCG5', 'SULT1C2', 
        'MSR1', 'CTSK', 'PTGDS', 'COLEC12', 
        'GPC4', 'PCOLCE2', 'CHIT1', 'KAL1', 
        'CLEC5A', 'ME1', 'DNASE2B', 'CCL7',
        'CD163', 'FN1', 'GM2A', 'SCARB2', 'BCAT1',
        'RAI14', 'COL8A2', 'APOE', 'CHI3L1', 
        'ATG7', 'CD84', 'FDX1', 'MS4A4A', 'SGMS1',
        'EMP1', 'CYBB', 'CD68')


    macro.sig %in% colnames(tpm.df)
    cluster.sig = cluster.sig[cluster.sig %in% colnames(tpm.df)]

    tpm.df2 = tpm.df
    geo.mean = function(x) exp(mean(log(x))) 
    #some of the samples in the TPM data frame have NAs for a lot of genes, need to remove these 
    tpm.df2 = tpm.df2[-unique(uulapply(tpm.df2, function(x) which(is.na(x)))), ]

    sig.ratio = apply(tpm.df2[cluster.sig], 1, geo.mean) / apply(tpm.df2[macro.sig], 1, geo.mean)
    tpm.df2$sig.ratio = sig.ratio

    tpm.df2[cluster.sig]
    tpm.df2$cluster.sig = apply(tpm.df2[cluster.sig], 1, geo.mean)

    tpm.df2 = tpm.df2[!sig.ratio == Inf, ]

    ############################
    #SIG.RATIO data
    ############################

    tpm.df2.resp = tpm.df2[c('sig.ratio', 'resp')]

    tpm.df2.resp %>%
        group_by(resp) %>%
        summarise(
            mean = mean(sig.ratio),
            sd = sd(sig.ratio)
        )

    head(tpm.df2.resp)

    #ecm.cpi.plot = ggplot(tpm.df2.resp, aes(x = resp, y = sig.ratio)) +
        #geom_jitter(width = 0.3)

    #pl(ecm.cpi.plot)

    tpm.df2.resp$cat = 'sig.ratio'
    df1 = tpm.df2.resp

    ############################
    #CLUSTER data
    ############################

    tpm.df2 = tpm.df
    geo.mean = function(x) exp(mean(log(x))) 
    #some of the samples in the TPM data frame have NAs for a lot of genes, need to remove these 
    tpm.df2 = tpm.df2[-unique(uulapply(tpm.df2, function(x) which(is.na(x)))), ]

    sig.ratio = apply(tpm.df2[cluster.sig], 1, geo.mean) / apply(tpm.df2[macro.sig], 1, geo.mean)
    tpm.df2$sig.ratio = sig.ratio

    tpm.df2[cluster.sig]
    tpm.df2$cluster.sig = apply(tpm.df2[cluster.sig], 1, geo.mean)


    tpm.df2.resp = tpm.df2[c('cluster.sig', 'resp')]

    tpm.df2.resp %>%
        group_by(resp) %>%
        summarise(
            mean = mean(cluster.sig),
            sd = sd(cluster.sig)
        )

    df2 = tpm.df2.resp

    df1$cat = 'sig.ratio'
    df2$cat = 'cluster.sig'

    list(df1, df2)
})

cpi.cluster.plots.sig = lapply(cpi.cluster.box.plots, function(x) x[[1]])
cpi.cluster.plots.clus = lapply(cpi.cluster.box.plots, function(x) x[[2]])

for(i in 1:23){
    cpi.cluster.plots.sig[[i]]$cluster = i - 1
}


for(i in 1:23){
    cpi.cluster.plots.clus[[i]]$cluster = i - 1
}

cpi.sig = bind_rows(cpi.cluster.plots.sig) %>% reset.rownames
cpi.clus = bind_rows(cpi.cluster.plots.clus) %>% reset.rownames

cpi.sig = cpi.sig[!is.na(cpi.sig$sig.ratio), ]

sig.results.ratio = sig.results %>%
    filter(cat == 'sig.ratio') %>%
    mutate(resp = 'response') %>%
    mutate(cluster = as.numeric(as.character(cluster))) %>%
    inner_join(clus.name.df, by = 'cluster')

cpi.sig = inner_join(cpi.sig, clus.name.df, by = 'cluster')

cpi.sig$short.label[grepl('^[0-9]_', cpi.sig$short.label)] =
    paste0('0', cpi.sig$short.label[grepl('^[0-9]_', cpi.sig$short.label)])


cpi.plot1 = ggplot(cpi.sig, aes(x = resp, y = log(sig.ratio), fill = resp)) +
    geom_boxplot(outlier.shape = NA) +
    geom_text(
        inherit.aes = F,
        data = sig.results.ratio,
        aes(
            x = 1.5,
            y = 6,
            label = '*'
        )
    ) +
    facet_grid(cols = vars(short.label)) +
    theme_classic() +
    rot.lab() +
    theme(
        axis.text.x = element_blank(),
        strip.text = element_text(angle = 90)
    ) +
    ggtitle('Signature ratio')

#PANEL - cpi respose plot, signature ratios
#pl(cpi.plot1, 10, 6)


sig.results.clus = sig.results %>%
    filter(cat == 'cluster.sig') %>%
    mutate(resp = 'response') %>%
    mutate(cluster = as.numeric(as.character(cluster))) %>%
    inner_join(clus.name.df, by = 'cluster') 

sig.results.clus$short.label[grepl('^[0-9]_', sig.results.clus$short.label)] = 
    paste0('0', sig.results.clus$short.label[grepl('^[0-9]_', sig.results.clus$short.label)])

cpi.clus = inner_join(cpi.clus, clus.name.df, by = 'cluster')
#cpi.clus$short.label = factor(cpi.clus$short.label, levels = clus.name.df$short.label)


cpi.clus$short.label[grepl('^[0-9]_', cpi.clus$short.label)] =
    paste0('0', cpi.clus$short.label[grepl('^[0-9]_', cpi.clus$short.label)])

cpi.clus$resp[cpi.clus$resp == 'response'] = 'R'
cpi.clus$resp[cpi.clus$resp == 'no_response'] = 'NR'

head(cpi.clus)
sig.results.clus

cpi.plot2 = ggplot(cpi.clus, aes(x = resp, y = log(cluster.sig), fill = resp)) +
    geom_boxplot(outlier.shape = NA) +
    geom_text(
        inherit.aes = F,
        data = sig.results.clus,
        aes(
            x = 1.5,
            y = 8,
            label = '*'
        )
    ) +
    facet_grid(.~short.label) +
    theme_classic() +
    rot.lab()  +
    theme(
        axis.text.x = element_blank(),
        strip.text = element_text(angle = 90),
        axis.ticks.x = element_blank()
    ) +
    xlab('Response') +
    ylab('loge(cluster signature expression)') +
    labs(fill = 'Response') 


#PANEL - cpi response plot with significance markers for cluster signatures (raw, not ratios)
#pl(cpi.plot2, 8, 5)

cpi.clus.sig = cpi.clus[cpi.clus$short.label %in% sig.results.clus$short.label, ]

cpi.plot.only.sig = ggplot(cpi.clus.sig, aes(x = resp, y = log(cluster.sig), fill = resp)) +
    geom_boxplot(outlier.shape = NA) +
    geom_text(
        inherit.aes = F,
        data = sig.results.clus,
        aes(
            x = 1.5,
            y = 8,
            label = '*'
        )
    ) +
    facet_grid(.~short.label) +
    theme_classic() +
    rot.lab()  +
    theme(
        axis.text.x = element_blank(),
        strip.text = element_text(angle = 90),
        axis.ticks.x = element_blank()
    ) +
    xlab('Response') +
    ylab('loge(cluster signature expression)') +
    labs(fill = 'Response') 

#PANEL - cpi response plot only significant
#pl(cpi.plot.only.sig, 4, 4)

#pl(grid.arrange(cpi.plot1, cpi.plot2, ncol = 1), 20, 8)

clus8plot = cpi.cluster.plots.sig[[9]]

cpi.plot.clus8 = ggplot(clus8plot, aes(x = resp, y = sig.ratio, fill = resp)) +
    geom_boxplot() +
    theme_classic() +
    rot.lab() 

#pl(cpi.plot.clus8)

clus8plot2 = cpi.cluster.plots.clus[[9]]

cpi.plot.clus8.clus = ggplot(clus8plot2, aes(x = resp, y = cluster.sig, fill = resp)) +
    geom_boxplot() +
    theme_classic() +
    rot.lab() 

#pl(cpi.plot.clus8.clus)

############################
#CPI1000 CORRELATION BETWEEN ECM SIG AND T CELL SIG 
############################

sig.markers = read.tsv('~/work/ucl/data/ucl.projects/macrophage/all.sig.markers.per.cluster.tsv')

sig.markers = sig.markers %>%
    filter(avg_log2FC > 0) %>%
    arrange(cluster, -avg_log2FC)

t.cell.sig = c(
    'PRKCQ',
    'CD3D',
    'CD28',
    'LCK',
    'TRAT1',
    'BCL11B',
    'CD2',
    'TRBC1',
    'TRAC',
    'ITM2A',
    'SH2D1A',
    'CD6',
    'CD96',
    'NCALD',
    'GIMAP5',
    'TRA',
    'CD3E',
    'SKAP1'
)

cpi.tpm[ecm.sig, ]

ecm.sig = sig.markers %>%
    filter(cluster == 18) %>%
    .[1:10, ] %>%
    pull(gene)

t.cell.sig = t.cell.sig[t.cell.sig %in% rownames(cpi.tpm)]
t.cell.sig.scores = apply(cpi.tpm[t.cell.sig, ], 2, geo.mean)
ecm.sig.scores = apply(cpi.tpm[ecm.sig, ], 2, geo.mean)

tpm.t.cell.vs.ecm = data.frame(
    t.cell.sig = t.cell.sig.scores,
    ecm.sig = ecm.sig.scores,
    patient = names(t.cell.sig.scores)
)

head(tpm.t.cell.vs.ecm)
dim(tpm.t.cell.vs.ecm)

ecm.quarts = quantile(tpm.t.cell.vs.ecm$ecm.sig)
df1 = tpm.t.cell.vs.ecm[tpm.t.cell.vs.ecm$ecm.sig <= ecm.quarts[2], ]
df2 = tpm.t.cell.vs.ecm[tpm.t.cell.vs.ecm$ecm.sig >= ecm.quarts[3], ]

df1$quartile = 'Lower'
df2$quartile = 'Upper'

df3 = bind_rows(df1, df2)
df3$quartile = factor(df3$quartile, levels = c('Lower', 'Upper'))

df3
head(df1)

wilcox.test(df1$t.cell.sig, df2$t.cell.sig)


tcell.ecm.boxplot = ggplot(df3, aes(x = quartile, y = t.cell.sig)) +
    geom_boxplot() +
    ylab('T cell signature expression') +
    xlab('ECM signature quartile') +
    theme_classic()

#pl(tcell.ecm.boxplot, 4, 4)

tcell.ecm = ggplot(tpm.t.cell.vs.ecm, aes(x = t.cell.sig, y = ecm.sig)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    theme_classic() +
    xlab('T cell signature expression') +
    ylab('ECM signature expression')

#PANEL - ECM vs general T-Cell signature (measuring infiltration)
#pl(tcell.ecm)

cor.test(tpm.t.cell.vs.ecm$t.cell.sig, tpm.t.cell.vs.ecm$ecm.sig)

############################
#COSMX ANALYSIS PANELS 
############################

#source('~/work/ucl/scripts/ucl.projects/macrophage/cosmx.analysis.R')

############################
#MANA SCORE PROPELLER 
############################

#this is the original code, which loaded in mac.clus
#I've done careful mapping of sample IDs
#between this code and Danwen's sample IDs below in
#the Map() function. load.mac.clus() won't
#work here because it lacks rownames 
#and therefore violates the mapping work
#I've done below.

mac.clus = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/all.dat.w.norm.w.cluster.meta.data.rds')
mac.clus$cancer[mac.clus$cancer == 'thyroid'] = 'THCA'
mac.clus$cancer[mac.clus$cancer == 'melanoma'] = 'MEL'
mac.clus$cancer[mac.clus$cancer == 'Melanoma'] = 'MEL'
mac.clus$cancer[mac.clus$cancer == 'breast'] = 'BC'
mac.clus$cancer[mac.clus$cancer == 'liver'] = 'LIHC'
mac.clus$cancer[mac.clus$cancer == 'ESO'] = 'ESCA'

mac.clus$study[mac.clus$study == 'WU'] = 'wu.lung'
mac.clus$study = tolower(mac.clus$study)

mac.clus = mac.clus %>%
    filter(tissue == 'Tumor')

mana = readRDS('~/work/ucl/data/ucl.projects/macrophage/danwen.mana.metadata.rds')
mana$study[mana$study == 'Wu'] = 'wu.lung'
mana$study = tolower(mana$study)
mana$study[mana$study == 'vishwakarma'] = 'vish'
mana$study[mana$study == 'yost'] = 'yost_basal'
mana$study[mana$study == 'zilionis'] = 'zili'

mana = mana[mana$study %in% mac.clus$study, ]
mac.clus = mac.clus[mac.clus$study %in% mana$study, ]

unique(mana$sample_ID)
mac.clus$sample = toupper(mac.clus$sample)
mana$sample = toupper(mana$sample_ID)
mana$sample = gsub('[A-Z]*_(.*)', '\\1', mana$sample)

mana.split = split(mana, mana$study)
mac.split = split(mac.clus, mac.clus$study)

mana.mappings = Map(function(x, y){
    print(unique(x$study))
    if(unique(x$study) != unique(y$study)) print('study mismatch')
    if(unique(x$study) == 'azizi' & unique(y$study) == 'azizi'){
        x$sample = toupper(gsub('Azizi_(.*)', '\\1', x$sample_ID))
    }
    if(unique(x$study) == 'bi' & unique(y$study) == 'bi'){
        x$sample = toupper(gsub('Bi_(.*)', '\\1', x$sample_ID))
        y$sample = gsub('(.*?)_.*', '\\1', y$sample)
    }
    if(unique(x$study) == 'chan' & unique(y$study) == 'chan'){
        unique(x$sample)

        x$sample = gsub('([A-Z0-9]*)_.*', '\\1', x$sample)
        y$sample = toupper(gsub('([A-Z0-9]*)_.*', '\\1', (rownames(y))))
    }
    if(unique(x$study) == 'krishna' & unique(y$study) == 'krishna'){
        y = y[!grepl('PBMC', y$sample), ]
        y = y[!grepl('NORMAL', y$sample), ]
        y = y[!grepl('LYMPHNODE', y$sample), ]

        y$sample = gsub('([A-Z0-9]*)_.*', '\\1', y$sample)
    }
    if(unique(x$study) == 'li' & unique(y$study) == 'li'){
        unique(x$sample)
        unique(y$sample)
        y$sample = gsub('([A-Z0-9]*)-[0-9]*-.*', '\\1', y$sample)
        y$sample = gsub('-', '_', y$sample)
        x = x[x$sample %in% y$sample, ]
        y = y[y$sample %in% x$sample, ]
    }
    if(unique(x$study) == 'qian' & unique(y$study) == 'qian'){
        y = y[y$tissue == 'Tumor', ]
        y$sample = paste0(y$patient, '_TUMOR')
    }
    if(unique(x$study) == 'vish' & unique(y$study) == 'vish'){
        x$sample = gsub('_TUMOR', '', x$sample)
    }
    if(unique(x$study) == 'wu.lung' & unique(y$study) == 'wu.lung'){
        y$sample = y$patient
        y$sample
    }
    if(unique(x$study) == 'yost_basal' & unique(y$study) == 'yost_basal'){
        y$sample = gsub('\\.', '_', toupper(gsub('.*(su[0-9]*\\.[a-z]*).*', '\\1', rownames(y))))
    }
    if(unique(x$study) == 'zili' & unique(y$study) == 'zili'){
        y$sample = toupper(gsub('1$|2$', '', gsub('([a-z0-9]*)_.*', '\\1', rownames(y))))
    }

    g = x %>%
        group_by(sample) %>%
        summarise(
            mean.cxcl13.score = mean(CXCL13_score),
            mean.tcf7.score = mean(TCF7_score),
            study = unique(study),
            num.t.cells = n()
        )

    x2 = inner_join(y, g, by = c('sample', 'study'))
    x2
}, mana.split, mac.split) %>% bind_rows

distinct(mana.mappings[c('study', 'cancer')])

sum(mana.mappings$study %in% c('jerby', 'li', 'yost_basal'))
sum(!mana.mappings$study %in% c('jerby', 'li', 'yost_basal'))

mana.mappings$mana.score = mana.mappings$mean.cxcl13.score
mana.mappings$mana.score[mana.mappings$study %in% c('jerby', 'li', 'yost_basal')] = 
    mana.mappings$mean.tcf7.score[mana.mappings$study %in% c('jerby', 'li', 'yost_basal')]

#there are only 3578 macrophages in the skin datasets,
#so I'm going to just remove these from the analysis

mana.mappings = mana.mappings %>%
    filter(!study %in% c('jerby', 'li', 'yost_basal'))

mana.quarts = quantile(mana.mappings$mana.score)

mana.mappings$mana.cat[mana.mappings$mana.score < mean(mana.mappings$mana.score)] = 'low'
mana.mappings$mana.cat[mana.mappings$mana.score > mean(mana.mappings$mana.score)] = 'high'

mana.mappings = mana.mappings %>%
    filter(mana.score <= mana.quarts[[2]] | mana.score >= mana.quarts[[3]])

#mana.samp.scores = distinct(mana.mappings[c('study', 'sample', 'mana.score', 'mana.cat')])
#mana.samp.scores$study.samp = paste0(mana.samp.scores$study, '_', mana.samp.scores$samp)
#write.tsv(mana.samp.scores, '~/work/ucl/data/ucl.projects/macrophage/mana.samp.scores.tsv')

#mana.mappings$mana.quartile = cut(mana.mappings$mana.score, quantile(mana.mappings$mana.score))

colnames(mana.mappings)
mana.mappings$cluster = as.numeric(as.character(mana.mappings$integrated_snn_res.0.525))
mana.mappings$study.samp = paste0(mana.mappings$study, '_', mana.mappings$sample)

mana.score.per.cancer.type = distinct(mana.mappings[c('study.samp', 'mana.score', 'cancer')])

plot1 = ggplot(mana.score.per.cancer.type, aes(x = cancer, y = mana.score)) +
    geom_boxplot() +
    rot.lab()

#pl(plot1)

#pdf('~/work/ucl/plots/ucl.projects/macrophage/supplementary.figures/mana.score.per.cancer.type.pdf', 8, 8)
#print(plot1)
#dev.off()

unique(mana.mappings$cancer)

mana.subset = mana.mappings %>%
    filter(cancer %in% c('LUAD', 'LUSC', 'SCLC', 'NSCLC'))


#num macrophages
nrow(mana.subset)
#num t-cells
sum(distinct(mana.subset[c('study.samp', 'num.t.cells')])$num.t.cells)

#studies
unique(mana.subset$study)

prop.results.subset = propeller(
    clusters = mana.subset$cluster,
    sample = mana.subset$study.samp,
    group = mana.subset$mana.cat
)

prop.results.subset = prop.results.subset %>% filter(FDR < 0.1)
prop.results.subset

mana.mappings %>%
    group_by(cancer, mana.cat, cluster) %>%
    summarise(n = n()) %>%
    filter(cluster == 18) %>%
    arrange(mana.cat, cancer)

prop.results = propeller(
    clusters = mana.mappings$cluster,
    sample = mana.mappings$study.samp,
    group = mana.mappings$mana.cat
)

prop.results = prop.results %>% filter(FDR < 0.1)

#prop.results.quartile %>% filter(FDR < 0.05)

#LIMITED THIS ANALYSIS TO ONLY LUNG CANCER AS WE WERE GETTING STRANGE RESULTS OTHERWISE
mana.mappings.lung = mana.mappings[mana.mappings$cancer %in% c('SCLC', 'LUAD', 'LUSC', 'NSCLC'), ]

m.split = split(mana.mappings.lung, mana.mappings.lung$mana.cat)

cluster.membership = lapply(m.split, function(x){
    df1 = as.data.frame(table(x$cluster)) %>% 
        arrange(as.numeric(as.character(Var1)))
    df1$scaled.freq = df1$Freq / nrow(x)
    df1$mana.cat = unique(x$mana.cat)
    df1
}) %>% bind_rows

cluster.membership$Var1 = as.factor(as.numeric(as.character(cluster.membership$Var1)))

cluster.membership = cluster.membership %>%
    filter(Var1 %in% rownames(prop.results.subset))

load.mac.clus()

cluster.names = distinct(mac.clus[c('cluster', 'short.label')]) %>%
    arrange(cluster) %>%
    rename(Var1 = cluster) %>%
    mutate(Var1 = as.character(Var1))

cluster.membership = left_join(cluster.membership, cluster.names)

cluster.barplot = ggplot(
    cluster.membership,
    aes(
        x = short.label,
        y = scaled.freq,
        fill = mana.cat)
    ) +
    geom_col(position = 'dodge') +
    ylab('% cluster membership') +
    xlab('Cluster') +
    ggtitle('MANA score high vs low macrophage composition') +
    theme_classic() +
    rot.lab()

#PANEL - MANA barplot for w/ significantly different clusters
#I changed this, cutoffs now represent the upper quartile vs the 
#lower quartile of mana scores (for high and low)
pl(cluster.barplot, 4, 4)







############################
#ECMMAC DISTRIBUTION 
############################

#mac.clus %>%
    #group_by(short.label, cancer) %>%
    #summarise(n = n()) %>%
    #filter(short.label == '18_ECMMac') %>%
    #arrange(n) %>%
    #mutate(sum = sum(n)) %>%
    #mutate(percent = (n / sum) * 100)

#mac.clus %>%
    #group_by(short.label, cancer, ss2) %>%
    #summarise(n = n()) %>%
    #filter(short.label == '18_ECMMac') %>%
    #arrange(-n) %>%
    #as.data.frame %>%
    #head(n = 20)

#mac.clus %>%
    #filter(ss2 == 'krishna_UT1_Lower_Tumor') %>%
    #group_by(short.label) %>%
    #summarise(n = n()) %>%
    #arrange(-n) %>%
    #as.data.frame %>%
    #mutate()

#mac.clus %>%
    #filter(ss2 == 'xu_xu-s2_Tumor') %>%
    #group_by(short.label) %>%
    #summarise(n = n()) %>%
    #arrange(-n) %>%
    #as.data.frame %>%
    #head(n = 20)

