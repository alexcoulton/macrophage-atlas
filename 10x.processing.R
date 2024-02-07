source('~/work/ucl/scripts/misc/functions.R')
source('~/work/ucl/scripts/ucl.projects/macrophage/mac.util.R')

library(org.Hs.eg.db)
library(jpeg)
library(HGNChelper)
library(ggsci)
library(cowplot)
library(rhdf5)
library(Seurat)
library(gridExtra)
library(dplyr)
library(Matrix)
library(hdf5r)
library(readxl)

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243275

#GSM7782696 	Breast Cancer, T2N1M0, 5' single-cell - this is from disaasociated tumour cells used for benchmarking the FPPE results
#GSM7782697 	Breast Cancer, T2N1M0, 3' single-cell - same as above but 3'
#GSM7782698 	Breast Cancer, Single Cell Gene Expression Flex - this is the FFPE sample
#GSM7782699 	Breast Cancer, Visium CytAssist Spatial Gene Expression
files1 = list.files('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/10xffpe/', pattern = '.h5$', full.names = T, recursive = T)
files.names = list.files('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/10xffpe/', pattern = '.h5$', recursive = T)

#samp.ids = gsub('.*?_(P[0-9]+_[a-z]+-Tx).*', '\\1', files.names)

dat.all = lapply(files1, Read10X_h5)

#files1 = list.files('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/10xffpe/', pattern = '.h5$', full.names = T, recursive = T)
#files.names = list.files('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/10xffpe/', pattern = '.h5$', recursive = T)

#samp.ids = gsub('.*?_(P[0-9]+_[a-z]+-Tx).*', '\\1', files.names)

#dat.all = lapply(files1, Read10X_h5)
#class(dat.all[[4]])

celltype = read_excel('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/10xffpe/GSE243275_Barcode_Cell_Type_Matrices.xlsx')

dim(dat.all[[4]])

############################
#VISIUM ANALYSIS 
############################

visium = dat.all[[4]]

sig.markers = read.tsv('~/work/ucl/data/ucl.projects/macrophage/all.sig.markers.per.cluster.tsv')

head(sig.markers)
sig2 = sig.markers %>%
    filter(avg_log2FC > 0) %>%
    arrange(cluster, -avg_log2FC)

sig.split = split(sig2, sig2$cluster)

markers = lapply(sig.split, function(x){
    x[1:10, ]$gene
})

markers[[2]]

geo.mean = function(x) exp(mean(log(x))) 

vis.seurat <- CreateSeuratObject(counts = visium)
vis.seurat <- NormalizeData(vis.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
vis.norm = GetAssayData(vis.seurat)

count1 = 1
visium.max.clus.score = apply(vis.norm, 2, function(x){
    print(count1)
    count = 0

    scores.all = lapply(markers, function(y){
        y = y[y %in% names(x)]
        df1 = data.frame(
            score = geo.mean(x[y]),
            cluster = count
        )
        count <<- count + 1
        df1
    }) %>% bind_rows

    scores.all = scores.all[!scores.all$cluster == 18, ]
    max.score.clus = scores.all[which.max(scores.all$score), ]
    count1 <<- count1 + 1
    max.score.clus
    return(max.score.clus)
})

visium2 = bind_rows(visium.max.clus.score)

visium2$barcode = colnames(vis.norm)
visium2 = visium2[visium2$score > 0, ]

sort(table(visium2$cluster))

celltype.visium = read_excel('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/10xffpe/GSE243275_Barcode_Cell_Type_Matrices.xlsx', sheet = 2)
celltype.visium = celltype.visium %>%
    rename(barcode = Barcode)

visium3 = left_join(visium2, celltype.visium)

head(visium3)

vis.sum = visium3 %>%
    pull(cluster) %>%
    table %>%
    sort %>%
    as.data.frame

colnames(vis.sum) = c('cluster', 'freq')

load.mac.clus()

x = sort(unique(mac.clus$short.label))

df1 = data.frame(
    clus.name = x,
    cluster = gsub('([0-9]*)_.*', '\\1', x)
)

vis.sum = left_join(vis.sum, df1)

vis.sum$clus.name = factor(vis.sum$clus.name, levels = vis.sum$clus.name)

vis.sum.plot = ggplot(vis.sum, aes(x = clus.name, y = freq)) +
    geom_col() +
    rot.lab()

pdf('~/work/ucl/plots/ucl.projects/macrophage/rebuttal.plots/visium.summary.plot.pdf', 8, 8)
print(vis.sum.plot)
dev.off()

visium3 %>%
    filter(cluster == 1) %>%
    pull(Annotation) %>%
    table %>%
    sort

img = readJPEG('~/work/ucl/data/ucl.projects/macrophage/10x/GSM7782699_tissue_image.jpg')

tissue.positions = read.csv('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/10xffpe/gsm7782699/spatial/tissue_positions.csv')

celltype.visium = read_excel('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/10xffpe/GSE243275_Barcode_Cell_Type_Matrices.xlsx', sheet = 2)
celltype.visium = celltype.visium %>%
    rename(barcode = Barcode)

tissue.positions = left_join(tissue.positions, celltype.visium)
tissue.positions$Cluster = factor(tissue.positions$Cluster)
tissue.positions$x_flipped = 19505 - tissue.positions$pxl_col_in_fullres

colours.10x = c(
    "b8e986",
    "1a2fcf",
    "5b2911",
    "e29d3f",
    "cf8b52",
    "a091c1",
    "6dbe5c",
    "b8e986",
    "00c5a7",
    "ced4d4",
    "00c3f1",
    "a29e1f",
    "00b5ff",
    "e39095",
    "056cda"
)

visium4 = visium3 %>%
    rename(
        mac.score = score,
        mac.cluster = cluster
    )

visium5 = visium4[c('barcode', 'mac.score', 'mac.cluster')]
visium5 = visium5 %>%
    filter(mac.cluster == 1)

tissue.positions = left_join(tissue.positions, visium5)

head(tissue.positions)

ylab1 = 'y position'
xlab1 = 'x position'

plot10x.img = ggplot(
    tissue.positions,
    aes(
        x = pxl_col_in_fullres,
        y = pxl_row_in_fullres,
        color = Cluster
        )
    ) +
    annotation_raster(
        img,
        xmin = 0,
        xmax = 19505,
        ymin = 0,
        ymax = 21571
    ) +
    coord_cartesian(xlim = c(0, 19505), ylim = c(0, 21571)) +
    theme_classic() +
    ylab(ylab1) +
    xlab(xlab1) +
    ggtitle('Tissue')

plot10x.points = ggplot(
    tissue.positions,
    aes(
        x = pxl_col_in_fullres,
        y = rev(pxl_row_in_fullres),
        color = Cluster
        )
    ) +
    geom_point() +
    coord_cartesian(xlim = c(0, 19505), ylim = c(0, 21571)) +
    scale_color_igv() +
    theme_classic() +
    ylab(ylab1) +
    xlab(xlab1) +
    ggtitle('Clusters')

plot10x.annotation = ggplot(
    tissue.positions,
    aes(
        x = pxl_col_in_fullres,
        y = rev(pxl_row_in_fullres),
        color = Annotation,
        shape = Annotation
        )
    ) +
    geom_point(size = 0.5) +
    coord_cartesian(xlim = c(0, 19505), ylim = c(0, 21571)) +
    scale_color_igv() +
    theme_classic() +
    ylab(ylab1) +
    xlab(xlab1) +
    scale_shape_manual(values = c(1:20)) +
    ggtitle('Annotations')


plot10x.points.mac = ggplot(
    tissue.positions,
    aes(
        x = pxl_col_in_fullres,
        y = rev(pxl_row_in_fullres),
        color = mac.score
        )
    ) +
    geom_point() +
    coord_cartesian(xlim = c(0, 19505), ylim = c(0, 21571)) +
    scale_color_viridis(option = 'magma', na.value = 'transparent') +
    theme_classic() +
    ylab(ylab1) +
    xlab(xlab1) +
    ggtitle('MetM2Mac Scores')

mac.score.per.anno = tissue.positions %>%
    group_by(Annotation) %>%
    summarise(
        mean.mac.score = mean(mac.score, na.rm = T),
        sd.mac.score = sd(mac.score, na.rm = T)
        ) %>%
    arrange(mean.mac.score)

mac.score.per.anno$Annotation = factor(
    mac.score.per.anno$Annotation,
    levels = unique(mac.score.per.anno$Annotation)
)

mac.anno.barplot = ggplot(mac.score.per.anno, aes(x = Annotation, y = mean.mac.score)) +
    geom_col() +
    geom_errorbar(
        aes(ymin = mean.mac.score - sd.mac.score, ymax = mean.mac.score + sd.mac.score),
        width = 0.5
    ) +
    rot.lab()

plot.10x.all = plot_grid(
    plot10x.img,
    plot10x.points,
    plot10x.points.mac,
    mac.anno.barplot,
    nrow = 2,
    ncol = 2,
    align = 'v',
    axis = 'lr'
)

#plng(plot.10x.all, 600, 600)

pdf('~/work/ucl/plots/ucl.projects/macrophage/rebuttal.plots/plot.10x.annotation.pdf', 8, 8)
print(plot10x.annotation)
dev.off()

ggsave(
    '~/work/ucl/plots/ucl.projects/macrophage/rebuttal.plots/plot.10x.all.png',
    plot.10x.all,
    width = 8,
    height = 8,
    dpi = 300
)





############################
#FIND TISSUE-SPECIFIC MACS 
############################

head(mac.clus)

clus.cancer.type.df = mac.clus %>%
    group_by(short.label, cancer) %>%
    summarise(n = n()) %>%
    arrange(short.label, n)

clus.cancer.type.df = clus.cancer.type.df %>%
    group_by(short.label) %>%
    mutate(total = sum(n)) %>%
    group_by(short.label) %>%
    mutate((n / total) * 100)

write.tsv(
    clus.cancer.type.df,
    '~/work/ucl/data/ucl.projects/macrophage/supplementary.tables/clus.cancer.type.df.tsv'
)





############################
#PROJECTION 
############################

sum(colnames(dat.all[[1]]) %in% colnames(dat.all[[2]]))
sum(colnames(dat.all[[1]]) %in% colnames(dat.all[[3]]))
sum(colnames(dat.all[[1]]) %in% colnames(dat.all[[4]]))
sum(colnames(dat.all[[2]]) %in% colnames(dat.all[[3]]))
sum(colnames(dat.all[[3]]) %in% colnames(dat.all[[4]]))

celltype$Barcode %in% colnames(dat.all[[3]])

g = celltype %>%
    filter(grepl('Macrophage', Annotation)) %>%
    pull(Barcode)

g

macs10x = dat.all[[3]][, g]

macs10x = CreateSeuratObject(
    macs10x,
    project = 'projection-test',
)

macs10x

############################
#PROJECTION 
############################

load.mac.int()

mac.int <- RunPCA(mac.int, assay = 'integrated', features = VariableFeatures(mac.int))

#t.anchors = FindTransferAnchors(
    #reference = mac.int,
    #query = macs10x,
    #dims = 1:50,
    #reference.reduction = 'pca',
    #normalization.method = 'SCT'
#)

#saveRDS(t.anchors, '~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/10xffpe/10x.transfer.anchors.rds')
t.anchors = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/10xffpe/10x.transfer.anchors.rds')

############################
#PREDICTIONS BASED ON ANCHORS 
############################

predictions <- TransferData(
    anchorset = t.anchors,
    refdata = mac.int$short.label,
    dims = 1:50
)

pred.sum = predictions %>%
    filter(prediction.score.max > 0.7) %>%
    pull(predicted.id) %>%
    table %>%
    sort %>%
    as.data.frame

colnames(pred.sum) = c('cluster', 'freq')

pred.sum

pred.sum.plot = ggplot(pred.sum, aes(x = cluster, y = freq)) +
    geom_col() +
    rot.lab()

pdf('~/work/ucl/plots/ucl.projects/macrophage/rebuttal.plots/prediction.summary.plot.pdf', 8, 8)
print(pred.sum.plot)
dev.off()


table(cut(predictions$prediction.score.max, 4))

sum(predictions$prediction.score.max > 0.7)

dim(dat.all[[4]])

head(mac.int@meta.data)

mac.int@meta.data %>%
    filter(short.label == '18_ECMMac') %>%
    group_by(cancer, short.label) %>%
    summarise(n = n())  %>%
    arrange(-n)


############################
#XENIUM 
############################

files2 = list.files(, full.names = T)
files.names = list.files('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/10xffpe/', pattern = '.h5$', recursive = T)

celltype.xen = read_excel(
    '~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/10xffpe/GSE243275_Barcode_Cell_Type_Matrices.xlsx',
    sheet = 3
)

#dim(celltype.xen)
#head(celltype.xen)
#dim(xenium[[1]])
#head(colnames(xenium[[1]]))
#celltype.xen$Barcode %in% colnames(xenium[[1]])

xenium = Read10X_h5('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/10xffpe/xenium/outs/cell_feature_matrix.h5')

length(xenium)
colnames(xenium[[1]])

rownames(xenium[[1]])
rownames(xenium[[2]])
rownames(xenium[[3]])
rownames(xenium[[4]])

xenium.genes = rownames(xenium[[1]])

m2.cluster.genes = c(
    'SELENOP',
    'SLC40A1',
    'F13A1',
    'RNASE1',
    'FOLR2',
    'PLTP',
    'STAB1',
    'PDK4',
    'LGMN',
    'DAB2',
    'MS4A6A',
    'MS4A4A',
    'MAF',
    'RBPJ',
    'CD163',
    'SLCO2B1',
    'CCL18',
    'FUCA1',
    'C1QC'
)

m2.cluster.genes %in% xenium.genes


#was confused about whether the sparse matrix Xenium files or the h5 files were the way to go
#the sparse matrix files (below) simply contain all of the h5 matrices combined.
#the first h5 matrix is the main one containing the probes from the panel.
#the others are control probes / genes 

#g1 = '~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/10xffpe/xenium/outs/cell_feature_matrix/matrix.mtx.gz'
#g2 = '~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/10xffpe/xenium/outs/cell_feature_matrix/features.tsv.gz'
#g3 = '~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/10xffpe/xenium/outs/cell_feature_matrix/barcodes.tsv.gz'
#xenium2 = ReadMtx(mtx = g1, cells = g3, features = g2)

xenium = CreateSeuratObject(
    xenium[[1]],
    project = 'projection-test',
)

head(celltype.xen)
mac.barcodes = celltype.xen[grepl('Macrophage', celltype.xen$ident), ]$Barcode
xen.macs = xenium[, mac.barcodes]

xenium.t.anchors = FindTransferAnchors(
    reference = mac.int,
    query = xen.macs,
    dims = 1:50,
    reference.reduction = 'pca',
    normalization.method = 'SCT'
)

predictions <- TransferData(
    anchorset = xenium.t.anchors,
    refdata = mac.int$short.label,
    dims = 1:50
)

head(predictions)
head(predictions)

predictions %>%
    filter(prediction.score.max > 0.9) %>%
    pull(predicted.id) %>%
    table %>%
    sort

table(cut(predictions$prediction.score.max, 4))

sum(predictions$prediction.score.max > 0.7)

dim(dat.all[[4]])



