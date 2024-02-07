# main figures panels labelled PANEL
# supplementary figures labelled SUPP

library(cowplot)
source('~/work/ucl/scripts/ucl.projects/macrophage/mac.util.R')

############################
#LOAD DATA 
############################

load.mac.clus()
mac.clus = mac.clus[mac.clus$cluster != 23, ]

############################
#COMPARISONS OF MACROPHAGE POPULATIONS 
############################
#e.g. liver mets from CRC vs primary liver tumours

crc.studies = c('becker', 'che', 'khaliq', 'lu', 'sharma', 'zheng')
m = arrange(mac.clus, study)
m = m %>% filter(study %in% crc.studies)
m$crc.comparison = m$study

m$crc.shape = NA
m$crc.shape[m$study %in% c('che')] = 'CRC Liver Met'
m$crc.shape[m$study %in% c('becker', 'khaliq', 'zheng')] = 'CRC Primary'
m$crc.shape[m$study %in% c('lu', 'sharma')] = 'LIHC Primary'
m$crc.shape = as.factor(m$crc.shape)

m$study[m$study == 'sharma'] = 'Sharma'
m$study[m$study == 'lu'] = 'Lu'
m$study[m$study == 'becker'] = 'Becker'
m$study[m$study == 'khaliq'] = 'Khaliq'
m$study[m$study == 'Zheng'] = 'Zheng'
m$study[m$study == 'Che'] = 'Che'


crc.all.umap.facet.plot.full = ggplot(m, aes(x = UMAP_1, y = UMAP_2, color = study)) +
    geom_point(size = 0.2) +
    facet_wrap(.~crc.shape) +
    theme_classic(base_size = 15) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(
        #plot.margin = unit(c(2, 0, 2, 0), 'lines'),
        axis.text = element_blank(),
        axis.ticks = element_blank()
    ) +
    labs(
        color = 'Study'
    ) 

crc.all.umap.facet.plot.18 = ggplot(m[m$cluster == 18, ], aes(x = UMAP_1, y = UMAP_2, color = study)) +
    geom_point(size = 0.2) +
    facet_wrap(.~crc.shape) +
    theme_classic(base_size = 15) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(
        #plot.margin = unit(c(2, 0, 2, 0), 'lines'),
        axis.text = element_blank(),
        axis.ticks = element_blank()
    ) +
    labs(
        color = 'Study'
    ) 

#PANEL - CRC UMAP
plng(crc.all.umap.facet.plot.full, 400, 200)

#PANEL - CRC UMAP
plng(crc.all.umap.facet.plot.18, 400, 200)

############################
#CRC PROPELLER TEST 
############################

m = mac.clus
m2 = m[c('sample', 'study', 'short.label')]
m2 = m2[m2$study %in% crc.studies, ]
colnames(m2) = c('sample', 'study', 'cluster')
m2 = reset.rownames(m2)

m2$group = ""
m2$group[m2$study %in% c('lu', 'sharma')] = 'liver_primary'
m2$group[m2$study %in% c('becker', 'khaliq', 'zheng')] = 'crc_primary'
m2$group[m2$study %in% c('che')] = 'crc_liver_met'

prop.results = propeller(clusters = m2$cluster, sample = m2$sample, group = m2$group)
prop.results %>% filter(FDR < 0.05)

############################
#CRC CLUSTER COMPOSITION PLOTTING
############################

m = mac.clus
m2 = m[c('sample', 'study', 'short.label')]

m2 = m2[m2$study %in% crc.studies, ]
colnames(m2) = c('sample', 'study', 'cluster')

m2$group = ""
m2$group[m2$study %in% c('lu', 'sharma')] = 'liver_primary'
m2$group[m2$study %in% c('becker', 'khaliq', 'zheng')] = 'crc_primary'
m2$group[m2$study %in% c('che')] = 'crc_liver_met'

m2$study[m2$study == 'sharma'] = 'Sharma'
m2$study[m2$study == 'lu'] = 'Lu'
m2$study[m2$study == 'becker'] = 'Becker'
m2$study[m2$study == 'khaliq'] = 'Khaliq'
m2$study[m2$study == 'Zheng'] = 'Zheng'
m2$study[m2$study == 'Che'] = 'Che'

m2.split = split(m2, m2$group)

cluster.membership = lapply(m2.split, function(x){
    df1 = as.data.frame(table(x$cluster)) 
    df1$group = unique(x$group)
    df1$scaled.freq = df1$Freq / nrow(x)
    df1
}) %>% bind_rows

cluster.membership = cluster.membership[cluster.membership$Var1 %in%
    rownames(prop.results %>% filter(FDR < 0.05)), ] 

head(cluster.membership)
unique(cluster.membership$group)
cluster.membership$group[cluster.membership$group == 'crc_liver_met'] = 'CRC liver met'
cluster.membership$group[cluster.membership$group == 'crc_primary'] = 'CRC primary'
cluster.membership$group[cluster.membership$group == 'liver_primary'] = 'LIHC primary'

crc.cluster.barplot = ggplot(cluster.membership, aes(x = Var1, y = scaled.freq * 100, fill = group)) +
    geom_col(position = 'dodge') +
    scale_fill_ucscgb() +
    ylab('% cells in cluster') +
    xlab('Cluster') +
    theme_classic() +
    labs(
        fill = 'Tumor type'
    ) +
    rot.lab()

crc.legend = get_legend(crc.cluster.barplot)

crc.cluster.barplot.no.legend = ggplot(cluster.membership, aes(x = Var1, y = scaled.freq * 100, fill = group)) +
    geom_col(position = 'dodge') +
    scale_fill_ucscgb() +
    ylab('% cells in cluster') +
    xlab('Cluster') +
    theme_classic() +
    labs(
        fill = 'Tumor type'
    ) +
    rot.lab() +
    theme(
        legend.position = 'none'
    )


#PANEL - CRC BARPLOT
pl(crc.cluster.barplot, 4, 4)
pl(crc.cluster.barplot.no.legend, 4, 4)
pl(plot_grid(crc.legend), 4, 4)


############################
#MELANOMA PROPELLER TEST 
############################

#without uveal
studies.mel = c('antunes', 'biermann', 'jerby', 'li')
#with uveal
#studies.mel = c('antunes', 'biermann', 'durante', 'jerby', 'li')
m = mac.clus
m = m %>% filter(study %in% studies.mel)
m2 = m[c('sample', 'study', 'short.label')]
m2 = rename(m2, cluster = short.label)

m2$group = ""
m2$group[m2$study %in% c('jerby', 'li')] = 'SKCM'
m2$group[m2$study %in% c('antunes')] = 'GBM'
m2$group[m2$study %in% c('biermann')] = 'SKCM_brain_met'
m2$group[m2$study %in% c('durante')] = 'uveal_primary'

prop.results = propeller(clusters = m2$cluster, sample = m2$sample, group = m2$group)
prop.results = prop.results %>% filter(FDR < 0.05)
prop.results = prop.results[!rownames(prop.results) == '23_NA', ]
prop.results = prop.results[1:10, ]


############################
#MELANOMA CLUSTER COMPOSITION 
############################

m.split = split(m2, m2$group)

cluster.membership = lapply(m.split, function(x){
    df1 = as.data.frame(table(x$cluster)) 
    df1$group = unique(x$group)
    df1$scaled.freq = df1$Freq / nrow(x)
    df1
}) %>% bind_rows

cluster.membership = cluster.membership[cluster.membership$Var1 %in%
    rownames(prop.results %>% filter(FDR < 0.05)), ] 

unique(cluster.membership$group)
cluster.membership$group[cluster.membership$group == 'GBM'] = 'GBM primary'
cluster.membership$group[cluster.membership$group == 'SKCM'] = 'SKCM primary'
cluster.membership$group[cluster.membership$group == 'SKCM_brain_met'] = 'SKCM brain met'

mela.cluster.barplot = ggplot(cluster.membership, aes(x = Var1, y = scaled.freq * 100, fill = group)) +
    geom_col(position = 'dodge') +
    scale_fill_ucscgb() +
    ylab('% cells in cluster') +
    xlab('Cluster') +
    theme_classic() +
    labs(
        fill = 'Tumor type'
    ) +
    rot.lab()

legend1 = get_legend(mela.cluster.barplot)

pl(plot_grid(legend1), 4, 4)

mela.cluster.barplot.no.legend = ggplot(cluster.membership, aes(x = Var1, y = scaled.freq * 100, fill = group)) +
    geom_col(position = 'dodge') +
    scale_fill_ucscgb() +
    ylab('% cells in cluster') +
    xlab('Cluster') +
    theme_classic() +
    labs(
        fill = 'Tumor type'
    ) +
    theme(legend.position = 'none') +
    rot.lab()


#PANEL - MELANOMA BARPLOT
pl(mela.cluster.barplot.no.legend, 4, 4)

############################
#MELANOMA UMAP 
############################

#without uveal
studies.mel = c('antunes', 'biermann', 'jerby', 'li')
#with uveal
#studies.mel = c('antunes', 'biermann', 'durante', 'jerby', 'li')
m = mac.clus
m = m %>% filter(study %in% studies.mel)
m2 = m[c('sample', 'study', 'short.label', 'UMAP_1', 'UMAP_2')]
m2 = rename(m2, cluster = short.label)

m2$group = ""
m2$group[m2$study %in% c('jerby', 'li')] = 'SKCM primary'
m2$group[m2$study %in% c('antunes')] = 'GBM primary'
m2$group[m2$study %in% c('biermann')] = 'SKCM brain met'
m2$group[m2$study %in% c('durante')] = 'Uveal primary'

unique(m2$study)
m2$study[m2$study == 'li'] = 'Li'
m2$study[m2$study == 'jerby'] = 'Jerby'
m2$study[m2$study == 'antunes'] = 'Antunes'
m2$study[m2$study == 'biermann'] = 'Biermann'

mela.all.umap.facet.plot.full = ggplot(m2, aes(x = UMAP_1, y = UMAP_2, color = study)) +
    geom_point(size = 0.2) +
    facet_wrap(.~group) +
    theme_classic(base_size = 12) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(
        #plot.margin = unit(c(2, 0, 2, 0), 'lines'),
        axis.text = element_blank(),
        axis.ticks = element_blank()
    ) +
    labs(
        color = 'Study'
    ) 


mela.all.umap.facet.plot.12 = ggplot(m2[m2$cluster == '12_MBMMac', ], aes(x = UMAP_1, y = UMAP_2, color = study)) +
    geom_point(size = 0.2) +
    facet_wrap(.~group) +
    theme_classic(base_size = 12) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(
        #plot.margin = unit(c(2, 0, 2, 0), 'lines'),
        axis.text = element_blank(),
        axis.ticks = element_blank()
    ) +
    labs(
        color = 'Study'
    ) 

#PANEL - MELANOMA UMAP
plng(mela.all.umap.facet.plot.full, 400, 200)
plng(mela.all.umap.facet.plot.12, 400, 200)

############################
#CLUSTER 12 COMPOSITION
############################

head(mac.clus)

mac.clus %>%
    filter(cluster == 12) %>%
    group_by(cancer) %>%
    summarise(n = n()) %>%
    arrange(n) %>%
    as.data.frame

############################
#LUNG COMPARISON 
############################

m = mac.clus[mac.clus$study %in% c('chan', 'kim', 'leader', 'maynard', 'wu.lung', 'zili'), ]

m = m[m$cancer %in% c('LUAD', 'LUSC'), ]
colnames(m)

m$sample[m$study == 'chan'] = 
    m %>%
        filter(study == 'chan') %>%
        pull(cellid) %>%
        gsub('([A-Z0-9]*_.)*_.*', '\\1', .)

m$sample[m$study == 'leader'] =
    m %>%
        filter(study == 'leader') %>%
        pull(cellid) %>%
        gsub('([0-9]*)_.*', '\\1', .)

m$sample[m$study == 'wu.lung'] = 
    m %>%
        filter(study == 'wu.lung') %>%
        pull(cellid) %>%
        gsub('(X[0-9]*)_.*', '\\1', .) 

m$sample[m$study == 'zili'] = 
    m %>%
        filter(study == 'zili') %>%
        pull(cellid) %>%
        gsub('([a-z0-9A-Z]*)_.*', '\\1', .)


m2 = m[c('study', 'cancer', 'sample', 'short.label')]

m2 = m2 %>%
    reset.rownames() %>%
    rename(cluster = short.label)

prop.results = propeller(clusters = m2$cluster, sample = m2$sample, group = m2$cancer)
prop.results %>% filter(FDR < 0.05)

m.split = split(m2, m2$cancer)

cluster.membership = lapply(m.split, function(x){
    df1 = as.data.frame(table(x$cluster)) 
    df1$scaled.freq = df1$Freq / nrow(x)
    df1$cancer = unique(x$cancer)
    df1
}) %>% bind_rows

head(cluster.membership)
cluster.membership = cluster.membership[cluster.membership$Var1 %in% rownames(prop.results[prop.results$FDR < 0.05, ]), ]

lung.cluster.barplot = ggplot(cluster.membership, aes(x = Var1, y = scaled.freq * 100, fill = cancer)) +
    geom_col(position = 'dodge') +
    scale_fill_ucscgb() +
    ylab('% cells in cluster') +
    xlab('Cluster') +
    theme_classic() +
    labs(
        fill = 'Tumor type'
    ) +
    rot.lab()

lung.legend = get_legend(lung.cluster.barplot)

lung.cluster.barplot.no.legend = ggplot(cluster.membership, aes(x = Var1, y = scaled.freq * 100, fill = cancer)) +
    geom_col(position = 'dodge') +
    scale_fill_ucscgb() +
    ylab('% cells in cluster') +
    xlab('Cluster') +
    theme_classic() +
    labs(
        fill = 'Tumor type'
    ) +
    rot.lab() +
    theme(
        legend.position = 'none'
    )


#PANEL - LUNG BARPLOT
pl(lung.cluster.barplot, 4, 4)

pl(lung.cluster.barplot.no.legend, 4, 4)
pl(plot_grid(lung.legend), 4, 4)

############################
#SAMPLE MAC COOCCURANCE 
############################

mac.split = split(mac.clus, mac.clus$study.samp)

mac.co = lapply(mac.split, function(x){
    data.frame(
        study.samp = unique(x$study.samp),
        num.cluster = nunique(x$cluster),
        cancer =  unique(x$cancer),
        num.macs = nrow(x)
    )
}) %>% bind_rows

mac.co2 = mac.co[mac.co$num.macs > 100, ]

mac.co.hist = ggplot(mac.co2, aes(x = num.cluster, fill = cancer)) +
    geom_bar() +
    scale_fill_ucscgb() +
    ylab('Sample frequency') +
    xlab('# clusters') +
    labs(
        fill = 'Cancer'
    ) +
    theme_classic()

#PANEL - HISTOGRAM SAMPLE X NUM CLUSTERS
pl(mac.co.hist, 6, 4)

mac.co %>%
    filter(cancer == 'MEL')


mac.infil = ggplot(mac.co, aes(x = num.cluster, y = num.macs)) + 
    geom_point()

#SUPP - SCATTERPLOT NUM CELLS VS NUM CLUSTERS
#pl(mac.infil, 5, 5)

co.mat = crossprod(table(mac.clus$study.samp, mac.clus$cluster))
diag(co.mat) = 0
co.melt = melt(co.mat)

co.heatmap = ggplot(co.melt, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_y_continuous(labels = 0:23, breaks = 0:23) +
    scale_x_continuous(labels = 0:23, breaks = 0:23)

pl(co.heatmap, 10, 10)

############################
#COOCCUR PACKAGE ANALYSIS 
############################

#i'm not using this package anymore because it only does a binary measure of
#co-occurance. I opted to do a correlation analysis instead.
#library(cooccur)

#all pairwise combinations of clusters occur together in some samples
#combinations = t(combn(0:22, 2))

#g2 = apply(combinations, 1, function(x){
    #g = mac.clus %>%
        #group_by(study.samp) %>%
        #summarise(co1 = all(x %in% cluster))
    #length(which(g$co1))
#})

#mac.tab = table(mac.clus$cluster, mac.clus$study.samp)
#mac.tab = matrix(mac.tab, ncol = ncol(mac.tab), dimnames = dimnames(mac.tab))

#rowSums(mac.tab)
##cooccur packages only takes binary matrix as input
#mac.tab[mac.tab > 0] = 1
#mac.tab

#head(clus.name.df)
#rownames(mac.tab) = clus.name.df$short.label
#mac.tab
#mac.tab = mac.tab[-nrow(mac.tab), ]

#co1 = cooccur(mat = mac.tab, type = 'spp_site', spp_names = T)
#co2 = co1$results

#pl(plot(co1), 20, 20)


#co2$enrich = log2(co2$obs_cooccur / co2$exp_cooccur)
#co2$enrich[co2$enrich == -Inf] = 0

#co2$sp1_name = factor(
    #co2$sp1_name,
    #levels = unique(co2$sp1_name)[order(as.numeric(gsub('([0-9]+)_.*', '\\1', unique(co2$sp1_name))))]
#)

#co2$sp2_name = factor(
    #co2$sp2_name,
    #levels = unique(co2$sp2_name)[order(as.numeric(gsub('([0-9]+)_.*', '\\1', unique(co2$sp2_name))))]
#)

#head(co2)
#max(co2$enrich)
#min(co2$enrich)

#enrich.plot = ggplot(co2, aes(x = sp1_name, y = sp2_name, fill = enrich)) +
    #geom_tile() +
    #scale_fill_gradientn(
        #colors = c(
            ##'#0006FF',
            #'#ffffff',
            #'#ff0000'
        #),
        #limits = c(0, 1),
        #na.value = 'white'
    #) +
    #rot.lab(size = 7) +
    #xlab('') +
    #ylab('') +
    #theme(
        #axis.text.y = element_text(size = 7)
    #) +
    #ggtitle('i')

#PANEL - COOCCURANCE PLOT
#pl(enrich.plot, 10, 10)


############################
#COOCCCUR SPLIT BY CANCER TYPE 
############################

#mac.clus.split = split(mac.clus, mac.clus$cancer)

#count = 1
#cooccur.by.cancer = lapply(mac.clus.split, function(z){
    #print(count)
    #mac.tab = table(z$cluster, z$study.samp)
    #mac.tab = matrix(mac.tab, ncol = ncol(mac.tab), dimnames = dimnames(mac.tab))

    ##cooccur packages only takes binary matrix as input
    #mac.tab[mac.tab > 0] = 1

    #rownames(mac.tab) = clus.name.df$short.label[sort(unique(z$cluster)) + 1]

    #co1 = cooccur(mat = mac.tab, type = 'spp_site', spp_names = T)
    #co2 = co1$results

    #co2$enrich = log2(co2$obs_cooccur / co2$exp_cooccur)
    #co2$enrich[co2$enrich == -Inf] = 0

    #co2$sp1_name = factor(
        #co2$sp1_name,
        #levels = unique(co2$sp1_name)[order(as.numeric(gsub('([0-9]+)_.*', '\\1', unique(co2$sp1_name))))]
    #)

    #co2$sp2_name = factor(
        #co2$sp2_name,
        #levels = unique(co2$sp2_name)[order(as.numeric(gsub('([0-9]+)_.*', '\\1', unique(co2$sp2_name))))]
    #)

    #count <<- count + 1

    #co2 %>%
        #filter(p_lt < 0.05)
    #co2

#})


#plots.all = lapply(cooccur.by.cancer['MEL'], function(x){
    #if(abs(min(x$enrich)) > max(x$enrich)){
        #lim1 = abs(min(x$enrich))
    #} else {
        #lim1 = abs(max(x$enrich))
    #}
    #enrich.plot = ggplot(x, aes(x = sp1_name, y = sp2_name, fill = enrich)) +
        #geom_tile() +
        #scale_fill_gradientn(
            #colors = c(
                #'#0006FF',
                #'#ffffff',
                #'#ff0000'
            #),
            #limits = c(-lim1, lim1),
            #na.value = 'white'
        #) +
        #rot.lab(size = 7) +
        #xlab('') +
        #ylab('') +
        #theme(
            #axis.text.y = element_text(size = 7)
        #)
    #enrich.plot
#})

#pl(plots.all[[1]])

#pl(do.call(grid.arrange, plots.all), 30, 30)

############################
#MAC PRESENCE BY CANCER TYPE 
############################

#cancer.cluster.presence = mac.clus %>%
    #group_by(cancer) %>%
    #summarise(
        #num.clus = nunique(short.label),
        #num.cells = n()
        #) %>%
    #arrange(num.cells) %>%
    #as.data.frame

#cc1plot = ggplot(cancer.cluster.presence, aes(x = num.clus, y = num.cells, color = cancer)) +
    #geom_point()

#SUPP
#pl(cc1plot)







