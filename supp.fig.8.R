# main figures panels labelled PANEL
# supplementary figures labelled SUPP
source('~/work/ucl/scripts/ucl.projects/macrophage/mac.util.R')
load.mac.clus()

############################
#RENDER SIZES 
############################

#sp = small panel
sp.x = 2
sp.y = 2


############################
#DATA WRANGLING 
############################

cold.tumours = mac.clus %>%
    group_by(ss2) %>%
    summarise(n = n()) %>%
    arrange(n) %>%
    filter(n < 200) %>%
    pull(ss2)

mac.clus %>%
    group_by(ss2) %>%
    summarise(n = n()) %>%
    arrange(n)

mac.clus$cold = F
mac.clus$cold[mac.clus$ss2 %in% cold.tumours] = T

############################
#HISTOGRAM OF CELL COUNTS PER SAMPLE 
############################

cell.hist = mac.clus %>%
    group_by(ss2) %>%
    summarise(n = n()) %>%
    arrange(n)

cell.hist2 = ggplot(cell.hist, aes(x = n)) +
    geom_histogram(bins = 1000) +
    theme_classic()

#PANEL - cell histogram
pl(cell.hist2, sp.x, sp.y)

############################
#COOCCURANCE CORRELATION ANALYSIS 
############################

combinations = t(combn(0:22, 2))

all.cor1 = apply(combinations, 1, function(x){
    m2 = mac.clus[mac.clus$cluster %in% x & mac.clus$cold == F & mac.clus$tissue == 'Tumor', ]

    m3 = m2 %>% 
        group_by(ss2) %>%
        summarise(
            num.clus.comparison1 = sum(cluster == x[1]),
            num.clus.comparison2 = sum(cluster == x[2])
        )

    cor1 = cor.test(m3$num.clus.comparison1, m3$num.clus.comparison2, method = 'spearman')

    data.frame(
        clus1 = x[1],
        clus2 = x[2],
        rho = cor1$estimate,
        p.val = cor1$p.value
    )

}) %>% bind_rows

all.cor1$clus1 = clus.name.df$short.label[match(all.cor1$clus1, clus.name.df$cluster)]
all.cor1$clus2 = clus.name.df$short.label[match(all.cor1$clus2, clus.name.df$cluster)]

cor.mat = matrix(nrow = 23, ncol = 23)
rownames(cor.mat) = unique(c(all.cor1$clus1, all.cor1$clus2))
colnames(cor.mat) = unique(c(all.cor1$clus1, all.cor1$clus2))

for(i in 1:23){
    for(p in 1:23){
        c1 = rownames(cor.mat)[i]
        c2 = colnames(cor.mat)[p]
        if(c1 == c2){
            rho = NA
        } else {
            rho = all.cor1 %>%
                filter(clus1 == c1 & clus2 == c2) %>%
                pull(rho)
            if(length(rho) == 0) {
                rho = all.cor1 %>%
                    filter(clus1 == c2 & clus2 == c1) %>%
                    pull(rho)
            }
        }
        cor.mat[c1, c2] = rho
    }
}

cor.mat[upper.tri(cor.mat)] = NA
diag(cor.mat) = NA
cor.mat = melt(cor.mat)


all.cor1$fdr = p.adjust(all.cor1$p.val, method = 'fdr')
top.correlations = head(arrange(all.cor1, -rho), n = 3)

head(cor.mat)
cor.mat$outline = NA
colnames(cor.mat) = c('clus1', 'clus2', 'rho', 'outline')
head(cor.mat)
cor.mat$label = ""
cor.mat$label[cor.mat$clus1 == cor.mat$clus2] = as.character(unique(cor.mat$clus1))
cor.mat$label
head(cor.mat)

for(i in 1:3){
    c1 = as.character(top.correlations$clus1[i])
    c2 = as.character(top.correlations$clus2[i])
    cor.mat$outline[cor.mat$clus1 == c2 & cor.mat$clus2 == c1] = 'Top'
}


cor.plot1 = ggplot(cor.mat, aes(x = clus1, y = clus2, fill = rho)) +
    geom_tile() +
    #geom_text(aes(label = sig)) +
    scale_fill_gradientn(
        colors = c(
            '#0006FF',
            '#ffffff',
            '#ff0000'
        ),
        limits = c(-1, 1),
        na.value = 'white'
    ) +
    #scale_color_manual(values = '#000000', na.value = 'transparent') +
    theme_classic() +
    labs(
        color = 'Rank',
        fill = "Rho"
    ) +
    geom_text(aes(label = label), color = 'black', size = 3, angle = 0, hjust = 1) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        plot.margin = unit(c(0, 0, 0, 3), 'lines'),
    ) +
    coord_cartesian(clip = 'off')
    #scale_x_discrete(position = 'top')

#PANEL - pairwise correlation matrix
pl(cor.plot1, 5, 4)



all.cor.raw = apply(combinations, 1, function(x){
    m2 = mac.clus[mac.clus$cluster %in% x & mac.clus$cold == F & mac.clus$tissue == 'Tumor', ]

    m3 = m2 %>% 
        group_by(ss2) %>%
        summarise(
            num.clus.comparison1 = sum(cluster == x[1]),
            num.clus.comparison2 = sum(cluster == x[2])
        )

    m3$comp = paste0(x[1], '_', x[2])
    m3

}) %>% bind_rows

############################
#INDIVIDUAL CORRELATION PLOTS 
############################

g.test = all.cor.raw %>%
    filter(comp == '6_10')

plot1 = ggplot(g.test, aes(x = num.clus.comparison1, y = num.clus.comparison2)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    xlab('6_SPP1AREGMac') +
    ylab('10_InflamMac') + 
    theme_classic()

#PANEL - top correlation 1
pl(plot1, sp.x, sp.y)

g.test = all.cor.raw %>%
    filter(comp == '3_16')

plot1 = ggplot(g.test, aes(x = num.clus.comparison1, y = num.clus.comparison2)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    xlab('3_ICIMac1') +
    ylab('16_ECMHomeoMac') +
    theme_classic()

#PANEL - top correlation 2
pl(plot1, sp.x, sp.y)

tail(all.cor1)

g.test = all.cor.raw %>%
    filter(comp == '0_22')

plot1 = ggplot(g.test, aes(x = num.clus.comparison1, y = num.clus.comparison2)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    xlab('0_AlvMac') +
    ylab('22_IFNMac4') +
    theme_classic()

#PANEL - negative correlation
pl(plot1, sp.x, sp.y)

############################
#RAW CORRELATION PLOT ALL SCATTERPLOTS 
############################

all.cor.raw.plot = ggplot(all.cor.raw, aes(x = num.clus.comparison1, y = num.clus.comparison2)) +
    geom_point() +
    facet_wrap(.~comp)

#SUPP - all correlations raw
#pl(all.cor.raw.plot, 40, 40)

############################
#COMPLEX HEATMAP COMPOSITION BY SAMPLE 
############################

combinations = t(combn(0:22, 2))

g2 = apply(combinations, 1, function(x){
    g = mac.clus %>%
        group_by(ss2) %>%
        summarise(co1 = all(x %in% cluster))
    length(which(g$co1))
})

mac.tab = table(mac.clus$cluster, mac.clus$ss2)
mac.tab = matrix(mac.tab, ncol = ncol(mac.tab), dimnames = dimnames(mac.tab))

#only hot
mac.tab = mac.tab[, !colnames(mac.tab) %in% unique(mac.clus$ss2[mac.clus$cold == T])]
#only tumor
#mac.tab = mac.tab[, colnames(mac.tab) %in% unique(mac.clus$ss2[mac.clus$tissue == 'Tumor'])]

#convert to proportions now only hot tumours are present
mac.tab[] = apply(mac.tab, 2, function(y) (y / sum(y)) * 100)

cancer = mac.clus$cancer[match(colnames(mac.tab), mac.clus$ss2)]
tissue = mac.clus$tissue[match(colnames(mac.tab), mac.clus$ss2)]

c1 = randomColor(nunique(cancer))
names(c1) = unique(cancer)
c2 = randomColor(nunique(tissue))
names(c2) = unique(tissue)

ggsci.ucscgb <- c(
    "chr5" = "#FF0000", "chr8" = "#FF9900", "chr9" = "#FFCC00",
    "chr12" = "#00FF00", "chr15" = "#6699FF", "chr20" = "#CC33FF",
    "chr3" = "#99991E", "chrX" = "#999999", "chr6" = "#FF00CC",
    "chr4" = "#CC0000", "chr7" = "#FFCCCC", "chr10" = "#FFFF00",
    "chr11" = "#CCFF00", "chr13" = "#358000", "chr14" = "#0000CC",
    "chr16" = "#99CCFF", "chr17" = "#00FFFF", "chr18" = "#CCFFFF",
    "chr19" = "#9900CC", "chr21" = "#CC99FF", "chr1" = "#996600",
    "chr2" = "#666600", "chr22" = "#666666", "chrY" = "#CCCCCC",
    "chrUn" = "#79CC3D", "chrM" = "#CCCC99"
)

c1 = ggsci.ucscgb[1:nunique(cancer)]
names(c1) = unique(cancer)


ggsci.npg <- c(
    "Cinnabar" = "#E64B35", "Shakespeare" = "#4DBBD5",
    "PersianGreen" = "#00A087", "Chambray" = "#3C5488",
    "Apricot" = "#F39B7F", "WildBlueYonder" = "#8491B4",
    "MonteCarlo" = "#91D1C2", "Monza" = "#DC0000",
    "RomanCoffee" = "#7E6148", "Sandrift" = "#B09C85"
)

c2 = ggsci.npg[1:nunique(tissue)]
names(c2) = unique(tissue)


rownames(mac.tab) =
    clus.name.df$short.label[match(rownames(mac.tab), clus.name.df$cluster)]
mac.tab = mac.tab[rownames(mac.tab) != '23_NA', ]

#HEATMAP PARAMETERS
hm.label.size = 15
hm.title.size = 18
dendrogram.height = 150
anno.height = 0.6

column.anno = HeatmapAnnotation(
    Cancer = cancer,
    Tissue = tissue,
    annotation_name_gp = gpar(fontsize = hm.title.size),
    col = list(
        Cancer = c1,
        Tissue = c2
    ),
    simple_anno_size = unit(anno.height, "cm"),
    annotation_legend_param = list(
        Cancer = list(
            nrow = 1,
            labels_gp = gpar(fontsize = hm.label.size),
            title_gp = gpar(fontsize = hm.title.size)
        ),
        Tissue = list(
            nrow = 1,
            labels_gp = gpar(fontsize = hm.label.size),
            title_gp = gpar(fontsize = hm.title.size)
        )
    )
)

heatmap1 = Heatmap(
    mac.tab,
    top_annotation = column.anno,
    cluster_rows = F,
    show_column_names = F,
    column_dend_height = unit(dendrogram.height, 'points'),
    row_names_gp = grid::gpar(fontsize = hm.label.size),
    heatmap_legend_param = list(
        title = '% cluster membership',
        labels_gp = gpar(fontsize = hm.label.size),
        title_gp = gpar(fontsize = hm.title.size),
        direction = 'horizontal'
    )
)

heatmap2 = draw(
    heatmap1,
    heatmap_legend_side = 'bottom',
    annotation_legend_side = 'bottom',
    merge_legend = T
)

#PANEL - mac composition heatmap
pl(heatmap2, 20, 9.5)

############################
#COOCCUR HCLUST / COMPLEX HEATMAP
############################

hclust1 = hclust(dist(t(mac.tab)))
#SUPP - mac composition hclust
#pl(plot(hclust1), 100, 100)

pco1 = cmdscale(dist(t(mac.tab)))
pco1 = as.data.frame(pco1)


pco1$ss2 = rownames(pco1)
sample.info = distinct(mac.clus[c('ss2', 'tissue', 'cancer')])
pco1 = inner_join(pco1, sample.info, by = 'ss2')

pco.plot = ggplot(pco1, aes(x = V1, y = V2, shape = tissue, color = cancer)) +
    geom_point() +
    theme_classic() +
    scale_shape_manual(values = c(15, 16, 17, 18)) +
    scale_color_ucscgb()

#PANEL - mac composition PCO
#pl(pco.plot, 10, 10)

pco.plot.facet = ggplot(pco1, aes(x = V1, y = V2, color = tissue)) +
    geom_point() +
    theme_classic() +
    scale_color_npg() +
    facet_wrap(.~cancer) +
    xlab('PCO1') +
    ylab('PCO2') +
    coord_cartesian(clip = 'off') +
    theme(
        #plot.margin = unit(c(2, 0, 2, 0), 'lines'),
        axis.text = element_blank(),
        axis.ticks = element_blank()
    ) +
    labs(
        color = 'Tissue'
    ) 

#SUPP - mac composition PCO facetted by cancer type
pl(pco.plot.facet, 6, 6)

############################
#CANCER TYPE HCLUST 
############################

mac.tab

############################
#IFN4 MAC INVESTIGATION 
############################

mac.clus %>%
    group_by(short.label, cancer) %>%
    summarise(n = n()) %>%
    arrange(short.label, n) %>%
    filter(short.label == '22_IFNMac4') %>%
    mutate(sum = sum(n)) %>%
    mutate(perc = (n / sum) * 100)


