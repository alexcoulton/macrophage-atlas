source('~/work/ucl/scripts/misc/functions.R')
library(rhdf5)
library(Seurat)
library(gridExtra)
library(dplyr)
library(Matrix)

#select only tumour or metastasis sample ids (not paratumours)
sampleids = 
    c(
        'GSM5585102',
        'GSM5585104',
        'GSM5585107',
        'GSM5585111',
        'GSM5585112',
        'GSM5585117',
        'GSM5585119',
        'GSM5585121',
        'GSM5585124'
    )


pu.dir = '/camp/home/coultoa/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/pu/GSE184362/suppl/'
g2 = list.files(pu.dir, full.names = T, recursive = T, pattern = 'GSM')

g2 = g2[uulapply(sampleids, function(x) grep(x, g2))]
pu.df = 
    data.frame(
        full.path = g2,
        file.name = basename(g2)
    )

pu.df$samp.id = multi.str.split(pu.df$file.name, '_', 1)
#filter blood and normal samples out
pu.df = pu.df %>%
    filter(!grepl('blood|normal', file.name))

pu.split = split(pu.df, pu.df$samp.id)

pu.dat = lapply(pu.split, function(x){
    x2 = x %>%
        filter(grepl('matrix', file.name))
    mtx.path = x2$full.path

    x3 = x %>%
        filter(grepl('barcodes', file.name))
    cells.path = x3$full.path 

    x4 = x %>%
        filter(grepl('features', file.name))
    genes.path = x4$full.path 

    g = ReadMtx(mtx = mtx.path, cells = cells.path, features = genes.path)
    g
})

sample.nums = paste0('s', 1:length(pu.dat))
sample.nums


pu.dat = Map(function(x, y){
    colnames(x) = paste0(y, '-', colnames(x))
    x
}, pu.dat, sample.nums)


pu.dat.comb = do.call(cbind, pu.dat)

test1 = CreateSeuratObject(
    pu.dat.comb,
    project = 'test',
    min.cells = 3,
    min.features = 200
)


test1[["percent.mt"]] <-
    PercentageFeatureSet(test1, pattern = "^MT-")

pl(
    VlnPlot(test1,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3),
20, 10)

test1 <- subset(
    test1,
    subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 15
)

test1 <- NormalizeData(test1)

test1 <- FindVariableFeatures(test1, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(test1), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(test1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#pl(grid.arrange(plot1, plot2), 20, 20)

all.genes <- rownames(test1)
test1 <- ScaleData(test1, features = all.genes)

test1 <- RunPCA(test1, features = VariableFeatures(object = test1))

#pl(ElbowPlot(test1))

test1 <- FindNeighbors(test1, dims = 1:17)
test1 <- FindClusters(test1, resolution = 1)
test1 <- RunUMAP(test1, dims = 1:20)

umapplot = DimPlot(test1, reduction = "umap") 
umapplot = LabelClusters(plot = umapplot, id = "ident")
pl(umapplot, 15, 15)

#cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)

#test1.markers <- FindAllMarkers(test1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#test1.markers %>%
    #group_by(cluster) %>%
    #slice_max(n = 10, order_by = avg_log2FC) %>%
    #as.data.frame %>%
    #select(gene, cluster)


head(GetAssayData(test1))

cd68plot = VlnPlot(test1, features = c("CD68"))
pl(cd68plot, 10, 10)

df1 = data.frame(
    expr = GetAssayData(test1)['CD68', ],
    cluster = Idents(test1)
)

head(df1)

plot.cd68 = ggplot(df1, aes(x = cluster, y = expr)) +
    geom_jitter(width = 0.2, alpha = 0.2)

pl(plot.cd68, 20, 10)

macro.sig = 
    c('MARCO', 'CXCL5', 'SCG5', 'SULT1C2', 
    'MSR1', 'CTSK', 'PTGDS', 'COLEC12', 
    'GPC4', 'PCOLCE2', 'CHIT1', 'KAL1', 
    'CLEC5A', 'ME1', 'DNASE2B', 'CCL7',
    'CD163', 'FN1', 'GM2A', 'SCARB2', 'BCAT1',
    'RAI14', 'COL8A2', 'APOE', 'CHI3L1', 
    'ATG7', 'CD84', 'FDX1', 'MS4A4A', 'SGMS1',
    'EMP1', 'CYBB', 'CD68')

macro.sig = macro.sig[macro.sig %in% rownames(GetAssayData(test1))]

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

df2 = data.frame(
    expr = apply(GetAssayData(test1)[macro.sig, ], 2, gm_mean),
    cluster = Idents(test1)
)

#pl(hist(df2$expr))
#head(df2)

plot.sig = ggplot(df2, aes(x = cluster, y = expr)) +
    geom_jitter(width = 0.2, alpha = 0.2)

pl(plot.sig, 20, 10)

threshold1 = mean(df2$expr) +
    sd(df2$expr)

cluster.means = df2 %>%
    group_by(cluster) %>%
    summarise(
        n = n(),
        mean.exp = mean(expr)
    ) %>%
    arrange(
        mean.exp
    ) %>% 
    as.data.frame() %>%
   filter(mean.exp > threshold1) %>%
   mutate(cluster = as.numeric(as.character(cluster)))

macrophages = GetAssayData(test1)[, Idents(test1) %in% cluster.means$cluster]

macrophages.raw = pu.dat.comb[, colnames(pu.dat.comb) %in%
    colnames(macrophages)]

dim(macrophages.raw)

saveRDS(macrophages.raw, '~/work/ucl/bigdata/ucl.projects/macrophage/pu/pu.macrophages.RDS')

mac = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/pu/pu.macrophages.RDS')

nrow(pu.df) / 3
colnames(mac)
pu.df

data.frame(
    samp.id = uulapply(pu.split, function(x) unique(x$samp.id)),
    s1 = paste0('s', 1:length(pu.split))
)


test1
ncol(mac)
mac[1:10, 1:10]


