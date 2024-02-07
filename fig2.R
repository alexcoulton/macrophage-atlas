# main figures panels labelled PANEL
# supplementary figures labelled SUPP
library(ape)
library(ggnewscale)
library(ggtree)
library(ggrepel)
source('~/work/ucl/scripts/ucl.projects/macrophage/mac.util.R')
load.mac.clus()
load.mac.int()

############################
#FUNCTIONS 
############################

adjust.label = function(cluster, x, y){
    label.coords$x.mean[label.coords$short.label == cluster] =
        label.coords$x.mean[label.coords$short.label == cluster] + x
    label.coords$y.mean[label.coords$short.label == cluster] = 
        label.coords$y.mean[label.coords$short.label == cluster] + y
    label.coords
}

############################
#CLUSTER UMAP 
############################

set.seed(1)
mac.clus = mac.clus[order(mac.clus$short.label), ]
#randomize order for better cluster visualization
mac.clus = mac.clus[sample(1:nrow(mac.clus), nrow(mac.clus)), ]

mac.clus$cluster = as.factor(mac.clus$cluster)

mac.clus = mac.clus %>%
    filter(cluster != 23)

############################
#GENERATE LABEL COORDS 
############################

label.coords = mac.clus %>%
    group_by(short.label) %>%
    summarize(
        x.mean = mean(UMAP_1),
        y.mean = mean(UMAP_2)
    )

label.coords$lab.text.color = "dark"
label.coords$lab.text.color[
    label.coords$short.label %in% c(
        '16_ECMHomeoMac',
        '18_ECMMac',
        '4_ICIMac2',
        '0_AlvMac',
        '7_IFNMac',
        '21_HemeMac'
    )
    ] = 'light'

#adjustments to coords
label.coords = adjust.label('0_AlvMac', -1, 0)
label.coords = adjust.label('1_MetM2Mac', 0.3, -2)
label.coords = adjust.label('4_ICIMac2', -2, 0)
label.coords = adjust.label('18_ECMMac', 0, 0.5)
label.coords = adjust.label('17_IFNMac3', 0, 1)
label.coords = adjust.label('16_ECMHomeoMac', 2, 1)
label.coords = adjust.label('3_ICIMac1', 1, 0)
label.coords = adjust.label('5_StressMac', 1, -1)
label.coords = adjust.label('12_MBMMac', -0.5, 1)
label.coords = adjust.label('13_CalciumMac', -1, 0)
label.coords = adjust.label('11_MetalloMac', 0, 1)
label.coords = adjust.label('6_SPP1AREGMac', 1.3, -0.8)
label.coords = adjust.label('2_C3Mac', 0, -0.5)
label.coords = adjust.label('10_InflamMac', 1, -1)
label.coords = adjust.label('20_TDoub', 0, 0.6)
label.coords = adjust.label('15_LYZMac', 0, 0.3)
label.coords = adjust.label('7_IFNMac', 0, 2)
label.coords = adjust.label('21_HemeMac', 3.5, -1)

clus.umap.raw = ggplot(mac.clus, aes(x = UMAP_1, y = UMAP_2, color = short.label)) +
    geom_point() +
    scale_color_igv() +
    new_scale_color() +
    geom_label(
        inherit.aes = F,
        data = label.coords,
        aes(
            label = short.label,
            x = x.mean,
            y = y.mean,
            fill = short.label,
            color = lab.text.color
            ),
        size = 8
    ) +
    scale_color_manual(values = c('#000000', '#ffffff')) +
    scale_fill_igv() +
    theme_classic() +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none',
        axis.title = element_text(size = 18)
    ) 

#PANEL - cluster UMAP
#plng(clus.umap.raw, 1000, 1000)

#label.colours = colours1$data[[2]][c('label', 'fill')]
#write.tsv(label.colours, '~/work/ucl/data/ucl.projects/macrophage/umap.label.colours.tsv')
label.colours = read.tsv('~/work/ucl/data/ucl.projects/macrophage/umap.label.colours.tsv')


clus.umap.raw.facet = ggplot(mac.clus, aes(x = UMAP_1, y = UMAP_2, color = short.label)) +
    geom_point() +
    scale_color_igv() +
    theme_classic() +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank()
    ) +
    facet_wrap(.~short.label)

#SUPP - UMAP facet 
#plng(clus.umap.raw.facet, 1000, 1000)


############################
#FGSEA 
############################

#write all reactome pathways to text file
#all_reactomeIds <- ls(reactomePATHID2EXTID)
#pathwayNames = unlist(mget(all_reactomeIds, reactomePATHID2NAME))
#human.pathways = pathwayNames[grep('Homo sapiens', pathwayNames)]
#writeLines(unname(human.pathways), '~/work/ucl/data/ucl.projects/macrophage/reactome.pathways.names.txt')


sig.markers = read.tsv('~/work/ucl/data/ucl.projects/macrophage/all.sig.markers.per.cluster.tsv')
sig.markers$entrezid = mapIds(org.Hs.eg.db, sig.markers$gene, 'ENTREZID', 'SYMBOL')
reactomePathways(sig.markers$entrezid)
sig.split = split(sig.markers, sig.markers$cluster)

#algorithm for picking pathways to test
#want to pick pathways that are in top and bottom of logFC values for each cluster
all.paths = lapply(sig.split, function(x){
    x = arrange(x, avg_log2FC)
    genes.down = head(x, n = 10)$entrezid
    genes.up = tail(x, n = 10)$entrezid
    genes.all = c(genes.down, genes.up)
    genes.all = genes.all[!is.na(genes.all)]
    genes.all
    paths.all = reactomePathways(genes.all)

    y = paths.all[order(uulapply(paths.all, length), decreasing = T)][1:4]
    names(y)
})

paths.set = sort(table(unlist(all.paths)), decreasing = T)[1:23]

#remove redundant pathways
paths.set = paths.set[!names(paths.set) %in% c(
    'Cell Cycle, Mitotic',
    'Cellular response to heat stress',
    'Gene Expression (Transcription)',
    'Immune System'
    )]

list.of.paths = list()
for(i in 0:22){
    print(paste0('running', i))
    g = sig.markers %>%
        filter(cluster == i)

    #get entrezids for the markers
    g$entrezid = mapIds(org.Hs.eg.db, g$gene, 'ENTREZID', 'SYMBOL')
    #g$rank = g$avg_log2FC * -log10(g$p_val)

    my_pathways = reactomePathways(g$entrezid)
    my_pathways = my_pathways[names(my_pathways) %in% names(paths.set)]
    #if(length(names(paths.set)[!names(paths.set) %in% names(my_pathways)]) > 0) browser()

    #paths.set[!names(paths.set) %in% names(my_pathways)]

    ranks1 = g$avg_log2FC
    names(ranks1) = g$entrezid 

    paths1 = fgsea(
        pathways = my_pathways,
        stats = ranks1,
        minSize = 15,
        maxSize = 500,
        nperm = 100000
    )


    paths1 = paths1 %>% 
        arrange(pval)

    list.of.paths = c(list.of.paths, list(paths1))

}

head(list.of.paths[[1]])
paths1 = unique(uulapply(list.of.paths, function(x) x$pathway))

############################
#PATHWAY RADAR PLOT 
############################

for(i in 1:length(list.of.paths)){
    list.of.paths[[i]]$cluster = i - 1
}

list.of.paths = bind_rows(list.of.paths)
list.of.paths = as.data.frame(list.of.paths)
list.of.paths = list.of.paths[c('cluster', 'pathway', 'pval')]
path.mat = dcast(list.of.paths, cluster ~ pathway, value.var = 'pval')
path.mat[is.na(path.mat)] = 1
path.mat = as.matrix(path.mat)
path.mat[, 2:ncol(path.mat)] = -log10(path.mat[, 2:ncol(path.mat)])

path.mat.sub = path.mat[c(6, 11, 15, 19), ]
path.mat.sub = as.data.frame(path.mat.sub)
path.mat.sub$cluster = clus.name.df$short.label[match(path.mat.sub$cluster, clus.name.df$cluster)]

#match radar plot cluster colours to UMAP
label.colours = read.tsv('~/work/ucl/data/ucl.projects/macrophage/umap.label.colours.tsv')

g1 = label.colours[label.colours$label %in% path.mat.sub$cluster, ]
g2 = g1$fill
names(g2) = g1$label

radar.plot1 = ggradar(path.mat.sub, grid.max = 5) +
    coord_cartesian(clip = 'off') +
    theme(
        plot.margin = unit(c(12, 12, 12, 12), 'lines'),
        legend.position = 'bottom'
    ) +
    scale_color_manual(values = g2)


#PANEL - radar pathway plot
pl(radar.plot1, 6, 6)



path.mat.sub = path.mat[c(7, 2, 1, 10), ]
path.mat.sub = as.data.frame(path.mat.sub)
path.mat.sub$cluster = clus.name.df$short.label[match(path.mat.sub$cluster, clus.name.df$cluster)]

g1 = label.colours[label.colours$label %in% path.mat.sub$cluster, ]
g2 = g1$fill
names(g2) = g1$label

radar.plot2 = ggradar(path.mat.sub, grid.max = 5) +
    coord_cartesian(clip = 'off') +
    theme(
        plot.margin = unit(c(12, 12, 12, 12), 'lines'),
        legend.position = 'bottom'
    ) +
    scale_color_manual(values = g2)


#PANEL - radar pathway plot
#pl(radar.plot2, 10, 10)


############################
#DOTPLOT 
############################

sig.markers = read.tsv('~/work/ucl/data/ucl.projects/macrophage/all.sig.markers.per.cluster.tsv')

top.markers = sig.markers %>%
    group_by(cluster) %>%
    summarise(
        top.markers = gene[1:2]
    ) %>%
    as.data.frame

top.markers$top.markers

DefaultAssay(mac.int) = 'SCT'

dot.plot.dat = lapply(0:23, function(x){
    exp.dat = GetAssayData(mac.int)[top.markers$top.markers, mac.clus$cluster == x]
    num.expressing = apply(exp.dat, 1, function(y) length(which(y > 0))) / ncol(exp.dat)
    mean.exp = rowMeans(exp.dat)
    data.frame(
        cluster = x,
        num.expressing = num.expressing,
        mean.exp = mean.exp,
        gene = rownames(exp.dat)
    )
}) %>% bind_rows

cluster.order = c(
    23,
    12,
    20,
    18,
    19,
    22,
    0,
    6,
    9,
    10,
    16,
    7,
    17,
    8,
    11,
    15,
    13,
    4,
    3,
    1,
    21,
    14,
    2,
    5
)

#dot.plot.dat$cluster = factor(dot.plot.dat$cluster, levels = cluster.order)
dot.plot.dat$cluster = as.numeric(as.character(dot.plot.dat$cluster))

dot.plot.dat = inner_join(dot.plot.dat, clus.name.df, by = 'cluster')
dot.plot.dat$short.label = factor(
    dot.plot.dat$short.label,
    levels = clus.name.df$short.label[cluster.order + 1]
)

dot.plot.dat = dot.plot.dat %>%
    filter(short.label != '23_NA')

head(dot.plot.dat)
unique(dot.plot.dat$short.label)

dot.plot.dat$short.label = factor(dot.plot.dat$short.label, levels = clus.name.df$short.label)


dot.plot = ggplot(dot.plot.dat, aes(x = short.label, y = gene, size = num.expressing, color = mean.exp)) +
    geom_point() +
    scale_color_gradient(
        low = '#FFFFFF',
        high = '#ff0000'
    ) +
    theme_classic() +
    rot.lab() +
    xlab('Cluster') +
    ylab('Gene') +
    labs(
        color = 'Mean expression',
        size = '% cells expressing'
    ) +
    theme(legend.position = 'bottom')


#PANEL - dot.plot
pl(dot.plot, 5, 9)

############################
#CLUSTER HCLUST 
############################

DefaultAssay(mac.int) = 'SCT'

clus.means = lapply(sort(unique(mac.clus$short.label)), function(x){
    print(x)
    rowMeans(GetAssayData(mac.int)[, mac.int@meta.data$short.label == x])
})

clus.means = bind_rows(clus.means)
clus.means = as.matrix(clus.means)
rownames(clus.means) = unique(mac.clus$short.label)
clus.dist = dist(clus.means)
#PANEL - cluster hclust
#pl(plot(hclust(clus.dist)), 10, 10)


#nn.tree = root(nj(clus.dist), outgroup = '19_ClassMono')
#nn.tree = nj(clus.dist)

#clus.dist.phylo = as.phylo(hclust(clus.dist))
#write.tree(phy = clus.dist.phylo, '~/work/ucl/data/ucl.projects/macrophage/mac.cluster.tree.newick')

clus.dist.phylo = read.tree('~/work/ucl/data/ucl.projects/macrophage/mac.cluster.tree.newick')
#clus.dist.phylo = root(clus.dist.phylo, outgroup = '19_ClassMono')


#tree.colours = c(1:nrow(label.colours))
#names(tree.colours) = label.colours$label
#clus.dist.phylo = groupOTU(clus.dist.phylo, tree.colours)

tree.colours

unclass(clus.dist.phylo)
clus.dist.phylo$tip.label = 
    gsub('(^[0-9]_)', '0\\1', clus.dist.phylo$tip.label)

clus.dist.phylo$tip.label

treeplot = ggtree(clus.dist.phylo) +
    geom_tiplab(angle = 0, hjust = 0, aes(colour = label)) +
    geom_tippoint(aes(color = label)) +
    #scale_x_reverse() +
    coord_cartesian(clip = 'off') +
    theme(
        plot.margin = unit(c(0, 6, 0, 0), 'lines'),
        legend.position = 'none'
    ) +
    scale_color_igv()
    #coord_flip(clip = 'off') 


pl(treeplot, 4, 6)







#patient.dist = hclust(dist(as.matrix(all.order.patient.mat)))
#write.tree(phy = as.phylo(patient.dist), '~/work/ucl/data/vhl/patient.hclust.newick')

