source('~/work/ucl/scripts/misc/functions.R')
library(org.Hs.eg.db)
library(HGNChelper)
library(cowplot)
library(rhdf5)
library(Seurat)
library(gridExtra)
library(dplyr)
library(Matrix)
library(hdf5r)

############################
#INITIAL PROCESSING OF LUOMA DATA
############################

#files1 = list.files('~/work/ucl/bigdata/ucl.projects/macrophage/luoma/', pattern = '.h5$', full.names = T)
#files.names = list.files('~/work/ucl/bigdata/ucl.projects/macrophage/luoma/', pattern = '.h5$')

#samp.ids = gsub('.*?_(P[0-9]+_[a-z]+-Tx).*', '\\1', files.names)

#dat.all = lapply(files1, Read10X_h5)

#for(i in 1:length(dat.all)){
    #colnames(dat.all[[i]]) = 
        #paste0(colnames(dat.all[[i]]), '_', samp.ids[[i]])
#}

#nunique(uulapply(dat.all, colnames))
#sum(uulapply(dat.all, ncol))
#colnames(dat.all[[1]])

#lapply(dat.all, dim)

#dat.comb = RowMergeSparseMatrices(
    #dat.all[[1]],
    #dat.all[2:length(dat.all)]
#)

#dim(dat.comb)

#myeloid.meta = read.tsv('~/work/ucl/bigdata/ucl.projects/macrophage/luoma/meta/GSE200996_Myeloid.tumor.single.cell.meta.data.txt')

#luoma.ids = data.frame(
    #CellType_ID = 1:9,
    #clus.label = c(
        #"CD14+ Mono",
        #"C1QB+ TAM",
        #"CXCL8+ TAM",
        #"Unassigned",
        #"CD1C+ DC",
        #"pDC",
        #"SPP1+ TAM",
        #"LAMP3+ DC",
        #"CLEC9A+ DC"
    #)
#)

#myeloid.meta = inner_join(myeloid.meta, luoma.ids, by = 'CellType_ID')

#myeloid.meta = myeloid.meta %>%
    #filter(grepl('TAM|CD14+', clus.label))


#colnames(dat.comb) = gsub('pre', 'Pre', colnames(dat.comb))
#colnames(dat.comb) = gsub('post', 'Post', colnames(dat.comb))
#colnames(dat.comb) = gsub('([A-Z]+)-[0-9]*(_P.*)', '\\1\\2', colnames(dat.comb))

#dat.mye = dat.comb[, colnames(dat.comb) %in% myeloid.meta$X]

############################
#CORRECT GENE NAMES 
############################

#current.map = getCurrentHumanMap()

#dat.genes = rownames(dat.mye)
#gene.syms = checkGeneSymbols(dat.genes, map = current.map)

#to.correct = gene.syms %>%
    #filter(Approved == F & !is.na(Suggested.Symbol))

#to.correct$present = to.correct$Suggested.Symbol %in% dat.genes
#to.correct = distinct(to.correct)

#dim(to.correct)

#dat.mye = lapply(list(dat.mye), function(x){
    #head(to.correct)
    #update = rownames(x)[rownames(x) %in% to.correct$x]
    #update2 = to.correct %>% filter(x %in% update) %>% reset.rownames
    #print(paste0('corrected: ', nrow(update2)))

    #rownames(x)[match(update2$x, rownames(x))] = update2$Suggested.Symbol
    #x
#})

#dat.mye = dat.mye[[1]]

#saveRDS(dat.mye, '~/work/ucl/bigdata/ucl.projects/macrophage/luoma/luoma.macs.extracted.rds')
luoma.macs = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/luoma/luoma.macs.extracted.rds')

luoma.macs = CreateSeuratObject(
    luoma.macs,
    project = 'projection-test',
)

############################
#PROJECTION 
############################

source('~/work/ucl/scripts/ucl.projects/macrophage/mac.util.R')
load.mac.int()

#mac.int <- RunPCA(mac.int, assay = 'integrated', features = VariableFeatures(mac.int))

#t.anchors = FindTransferAnchors(
    #reference = mac.int,
    #query = luoma.macs,
    #dims = 1:50,
    #reference.reduction = 'pca',
    #normalization.method = 'SCT'
#)

#saveRDS(t.anchors, '~/work/ucl/bigdata/ucl.projects/macrophage/luoma/luoma.transfer.anchors.rds')
t.anchors = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/luoma/luoma.transfer.anchors.rds')

############################
#PREDICTIONS BASED ON ANCHORS 
############################

predictions <- TransferData(
    anchorset = t.anchors,
    refdata = mac.int$short.label,
    dims = 1:50
)

predictions$cellid = rownames(predictions)
p2 = predictions[c('cellid', 'predicted.id')]

myeloid.meta = read.tsv('~/work/ucl/bigdata/ucl.projects/macrophage/luoma/meta/GSE200996_Myeloid.tumor.single.cell.meta.data.txt')
luoma.ids = data.frame(
    CellType_ID = 1:9,
    clus.label = c(
        "CD14+ Mono",
        "C1QB+ TAM",
        "CXCL8+ TAM",
        "Unassigned",
        "CD1C+ DC",
        "pDC",
        "SPP1+ TAM",
        "LAMP3+ DC",
        "CLEC9A+ DC"
    )
)
myeloid.meta = inner_join(myeloid.meta, luoma.ids, by = 'CellType_ID')
myeloid.meta = myeloid.meta %>%
    filter(grepl('TAM|CD14+', clus.label))


my2 = myeloid.meta[c('X', 'clus.label')]
colnames(my2)[1] = c('cellid')

p3 = inner_join(p2, my2, by = 'cellid')

luoma.heatmap.dat = p3 %>%
    group_by(clus.label) %>%
    summarise(
        freq = sort(table(predicted.id)),
        clus = names(sort(table(predicted.id)))
    ) %>%
    as.data.frame 

head(luoma.heatmap.dat)
luoma.heatmap.dat = as.data.frame(as.matrix(luoma.heatmap.dat))
luoma.heatmap.dat$freq = as.numeric(luoma.heatmap.dat$freq)

luoma.heatmap.dat = luoma.heatmap.dat %>%
    group_by(clus.label) %>%
    mutate(sum = sum(freq)) %>%
    mutate(percent = (freq / sum) * 100)

lu.heatmap = ggplot(luoma.heatmap.dat, aes(x = clus.label, y = clus, fill = percent)) +
    geom_tile() +
    scale_fill_gradient(
        low = '#ffffff',
        high = '#ff0000'
    ) +
    rot.lab() +
    xlab('Luoma et al label') +
    ylab('Atlas cluster mapping')

pl(lu.heatmap, 6, 6)

############################
#VIOLIN PLOTS 
############################

############################
#CXCL8 
############################

cxcl8.top3 = luoma.heatmap.dat %>%
    group_by(clus.label) %>%
    slice_max(freq, n = 3) %>%
    filter(clus.label == 'CXCL8+ TAM') %>%
    pull(clus) %>%
    as.character

mac.clus$cxcl8 = GetAssayData(mac.int, assay = 'SCT')['CXCL8', ]
mac.clus$cxcl8.rank = "> 3"

for(i in 1:length(cxcl8.top3)){
    mac.clus$cxcl8.rank[mac.clus$short.label %in% cxcl8.top3[[i]]] = i
}


#mac.clus$cxcl8.rank = as.numeric(mac.clus$cxcl8.rank)
mac.clus$cxcl8.rank = factor(mac.clus$cxcl8.rank, levels = c('1', '2', '3', '> 3'))

head(mac.clus)
mac.clus = mac.clus[mac.clus$short.label != '23_NA', ]

cxcl8.vio = ggplot(mac.clus, aes(x = short.label, y = cxcl8, fill = cxcl8.rank)) +
    geom_boxplot() +
    theme_classic() +
    rot.lab() +
    scale_fill_manual(
        values = c(
            '#fd0000',
            '#ff776e',
            '#ffbfbf',
            '#ffffff'
        ),
        na.value = 'transparent'
    ) +
    labs(fill = 'Mapping rank') +
    ylab('CXCL8 expression') +
    xlab('Reference cluster')

pl(cxcl8.vio, 5, 3) 

############################
#C1QB 
############################


c1qb.top3 = luoma.heatmap.dat %>%
    group_by(clus.label) %>%
    slice_max(freq, n = 3) %>%
    filter(clus.label == 'C1QB+ TAM') %>%
    pull(clus) %>%
    as.character

mac.clus$c1qb = GetAssayData(mac.int, assay = 'SCT')['C1QB', ]
mac.clus$c1qb.rank = "> 3"

for(i in 1:length(c1qb.top3)){
    mac.clus$c1qb.rank[mac.clus$short.label %in% c1qb.top3[[i]]] = i
}


#mac.clus$c1qb.rank = as.numeric(mac.clus$c1qb.rank)
mac.clus$c1qb.rank = factor(mac.clus$c1qb.rank, levels = c('1', '2', '3', '> 3'))


c1qb.vio = ggplot(mac.clus, aes(x = short.label, y = c1qb, fill = c1qb.rank)) +
    geom_boxplot() +
    theme_classic() +
    rot.lab() +
    scale_fill_manual(
        values = c(
            '#fd0000',
            '#ff776e',
            '#ffbfbf',
            '#ffffff'
        ),
        na.value = 'transparent'
    ) +
    labs(fill = 'Mapping rank') +
    ylab('C1QB expression') +
    xlab('Reference cluster')

pl(c1qb.vio, 20, 6) 

############################
#SPP1
############################

spp1.top3 = luoma.heatmap.dat %>%
    group_by(clus.label) %>%
    slice_max(freq, n = 3) %>%
    filter(clus.label == 'SPP1+ TAM') %>%
    pull(clus) %>%
    as.character

mac.clus$spp1 = GetAssayData(mac.int, assay = 'SCT')['SPP1', ]
mac.clus$spp1.rank = "> 3"

for(i in 1:length(spp1.top3)){
    mac.clus$spp1.rank[mac.clus$short.label %in% spp1.top3[[i]]] = i
}


#mac.clus$spp1.rank = as.numeric(mac.clus$spp1.rank)
mac.clus$spp1.rank = factor(mac.clus$spp1.rank, levels = c('1', '2', '3', '> 3'))

spp1.vio = ggplot(mac.clus, aes(x = short.label, y = spp1, fill = spp1.rank)) +
    geom_boxplot() +
    theme_classic() +
    rot.lab() +
    scale_fill_manual(
        values = c(
            '#fd0000',
            '#ff776e',
            '#ffbfbf',
            '#ffffff'
        ),
        na.value = 'transparent'
    ) +
    labs(fill = 'Mapping rank') +
    ylab('SPP1 expression') +
    xlab('Reference cluster' )

pl(spp1.vio, 5, 3) 


############################
#UMAP PROJECTION 
############################

#mac.int <- RunPCA(mac.int, assay = 'integrated', features = VariableFeatures(mac.int))
#mac.int <- RunUMAP(mac.int, assay = 'integrated', dims = 1:50, return.model = T)

#check if new UMAP is the same as the old one
mac.int.umap2 = Embeddings(mac.int, 'umap')
mac.int.umap2 = as.data.frame(mac.int.umap2)
mac.int.umap2$cellid = rownames(mac.int.umap2)
head(mac.int.umap2)

head(mac.int.umap2)
mm3 = inner_join(mac.int.umap2, mac.clus, by = 'cellid')
head(mm3)

mac.int.umap.plot = ggplot(
    mm3,
    aes(x = UMAP_1.x, y = UMAP_2.x, color = short.label)
    ) +
    geom_point() +
    scale_y_reverse() +
    scale_x_reverse()

plng(mac.int.umap.plot, 1000, 1000)

#load.mac.clus()

head(mac.clus[c('UMAP_1', 'UMAP_2', 'cellid')])
head(mac.int.umap2)



############################
#PERFORM UMAP PROJECTION OF NEW DATA 
############################

#luoma.umap <- MapQuery(
    #anchorset = t.anchors,
    #reference = mac.int,
    #query = luoma.macs,
    #refdata = list(short.label = "short.label"),
    #reference.reduction = "pca",
    #reduction.model = "umap"
#)

#saveRDS(luoma.umap, '~/work/ucl/bigdata/ucl.projects/macrophage/luoma/luoma.umap.rds')

luoma.umap = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/luoma/luoma.umap.rds')

luoma.umap2 = Embeddings(luoma.umap, 'ref.umap')
luoma.umap2 = as.data.frame(luoma.umap2)
luoma.umap2$cellid = rownames(luoma.umap2)
head(luoma.umap2)

luoma.umap2 = inner_join(luoma.umap2, p3, by = 'cellid')

head(luoma.umap2)
unique(luoma.umap2$predicted.id)
head(luoma.umap2)
luoma.umap2$predicted.id = factor(luoma.umap2$predicted.id, levels = clus.name.df$short.label)

luoma.umap.plot = ggplot(
    luoma.umap2,
    aes(x = refUMAP_1, y = refUMAP_2, color = predicted.id)
    ) +
    geom_point(size = 2) +
    scale_color_igv() +
    scale_y_reverse() +
    scale_x_reverse() +
    theme_classic(base_size = 25) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(
        axis.ticks = element_blank(),
        axis.text = element_blank()
    )

luoma.legend = get_legend(luoma.umap.plot)

pl(plot_grid(luoma.legend), 5, 8)


luoma.umap.plot.leg.bottom = ggplot(
    luoma.umap2,
    aes(x = refUMAP_1, y = refUMAP_2, color = predicted.id)
    ) +
    geom_point(size = 2) +
    scale_color_igv() +
    scale_y_reverse() +
    scale_x_reverse() +
    theme_classic(base_size = 25) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = 'bottom'
    )

luoma.legend2 = get_legend(luoma.umap.plot.leg.bottom)

pl(plot_grid(luoma.legend2), 8, 5)


luoma.umap.plot.no.legend = ggplot(
    luoma.umap2,
    aes(x = refUMAP_1, y = refUMAP_2, color = predicted.id)
    ) +
    geom_point(size = 2) +
    scale_color_igv() +
    scale_y_reverse() +
    scale_x_reverse() +
    theme_classic(base_size = 25) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none'
    )



plng(luoma.umap.plot.no.legend, 1000, 1000)



