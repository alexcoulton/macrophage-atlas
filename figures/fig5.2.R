library(randomForest)
library(cowplot)
source('~/work/ucl/scripts/ucl.projects/macrophage/mac.util.R')
load.mac.clus()
load.mac.int()

DefaultAssay(mac.int) = 'SCT'

ma.genes = c(
    'APOC1',
    'APOE',
    'ACP5',
    'FABP5',
    'VEGFA',
    'SPP1',
    'LYVE1',
    'HES1',
    'FOLR2',
    'MKI67',
    'CDK1',
    'IL1B',
    'CCL3',
    'CXCL1',
    'CXCL2',
    'CXCL3',
    'CXCL5',
    'ISG15',
    'CXCL8',
    'CXCL9',
    'CXCL10',
    'ARG1',
    'MRC1',
    'CD274',
    'CX3CR1'
)

ma.exp = lapply(ma.genes, function(x){
    print(x)
    GetAssayData(mac.int)[x, ]
})

ma.exp2 = bind_cols(ma.exp)
colnames(ma.exp2) = ma.genes

mac.clus2 = bind_cols(mac.clus, ma.exp2)


############################
#LYVE1 DIVERSITY 
############################

#Of all the LYVE1+ cells, what is their distribution
#among clusters (i.e. how phenotypically diverse are they)
#and how does this distribution change with percentiles 
#of LYVE1 expression? e.g. is it more skewed towards
#one cluster at very high levels of LYVE1 expression,
#and more heterogenous at lower levels?

#head(mac.clus2)

#mean(mac.clus2$LYVE1)
#quantile(mac.clus2$LYVE1)


#lyve.quart = mac.clus2 %>%
    #filter(LYVE1 != 0) %>%
    #pull(LYVE1) %>%
    #quantile(seq(0, 1, 0.01)) %>%
    #rev


#count = 100
#lyve.clus.dist = lapply(lyve.quart, function(x){
    #df1 = mac.clus2 %>%
        #filter(LYVE1 > x) %>%
        #pull(short.label) %>%
        #table %>%
        #sort %>%
        #as.data.frame

    #df1$percentile.dat = x
    #colnames(df1)[1] = 'short.label'
    #df1$percentile = count
    #count <<- count - 1
    #df1$freq.percent = (df1$Freq / sum(df1$Freq)) * 100
    #df1
#}) %>% bind_rows

#head(lyve.clus.dist)

#lyveplot = ggplot(lyve.clus.dist, aes(x = percentile, y = freq.percent, fill = short.label)) +
    #geom_col() +
    #scale_x_reverse() +
    #scale_fill_ucscgb() +
    #ggtitle('LYVE1')

#pl(lyveplot, 10, 5)


############################
#SPP1 DIVERSITY 
############################

#head(mac.clus2)

#mean(mac.clus2$SPP1)
#quantile(mac.clus2$SPP1)

#spp1.quart = mac.clus2 %>%
    #filter(SPP1 != 0) %>%
    #pull(SPP1) %>%
    #quantile(seq(0, 1, 0.01)) %>%
    #rev

#count = 100
#spp1.clus.dist = lapply(spp1.quart, function(x){
    #df1 = mac.clus2 %>%
        #filter(SPP1 > x) %>%
        #pull(short.label) %>%
        #table %>%
        #sort %>%
        #as.data.frame

    #df1$percentile.dat = x
    #colnames(df1)[1] = 'short.label'
    #df1$percentile = count
    #count <<- count - 1
    #df1$freq.percent = (df1$Freq / sum(df1$Freq)) * 100
    #df1
#}) %>% bind_rows

#head(spp1.clus.dist)

#spp1plot = ggplot(spp1.clus.dist, aes(x = percentile, y = freq.percent, fill = short.label)) +
    #geom_col() +
    #scale_x_reverse() +
    #scale_fill_ucscgb() +
    #ggtitle('SPP1')

#pl(spp1plot, 10, 5)

############################
#MKI67 DIVERSITY 
############################

#mean(mac.clus2$MKI67)
#quantile(mac.clus2$MKI67)

#mki67.quart = mac.clus2 %>%
    #filter(MKI67 != 0) %>%
    #pull(MKI67) %>%
    #quantile(seq(0, 1, 0.01)) %>%
    #rev

#count = 100
#mki67.clus.dist = lapply(mki67.quart, function(x){
    #df1 = mac.clus2 %>%
        #filter(MKI67 > x) %>%
        #pull(short.label) %>%
        #table %>%
        #sort %>%
        #as.data.frame

    #df1$percentile.dat = x
    #colnames(df1)[1] = 'short.label'
    #df1$percentile = count
    #count <<- count - 1
    #df1$freq.percent = (df1$Freq / sum(df1$Freq)) * 100
    #df1
#}) %>% bind_rows

#head(mki67.clus.dist)

#mki67plot = ggplot(mki67.clus.dist, aes(x = percentile, y = freq.percent, fill = short.label)) +
    #geom_col() +
    #scale_x_reverse() +
    #scale_fill_ucscgb() +
    #ggtitle('MKI67')

#pl(mki67plot, 10, 5)

############################
#ALL GENES LOOP 
############################


all.clus.dist = lapply(ma.genes, function(x){
    print(x)
    mki67.quart = mac.clus2 %>%
        filter(get(x) != 0) %>%
        pull(get(x)) %>%
        quantile(seq(0, 1, 0.01)) %>%
        rev

    count = 100
    mki67.clus.dist = lapply(mki67.quart, function(y){
        df1 = mac.clus2 %>%
            filter(get(x) > y) %>%
            pull(short.label) %>%
            table %>%
            sort %>%
            as.data.frame

        df1$percentile.dat = y
        colnames(df1)[1] = 'short.label'
        df1$percentile = count
        count <<- count - 1
        df1$freq.percent = (df1$Freq / sum(df1$Freq)) * 100
        df1
    }) %>% bind_rows

    head(mki67.clus.dist)

    gene.clus.dist.plot = ggplot(mki67.clus.dist, aes(x = percentile, y = freq.percent, fill = short.label)) +
        geom_col() +
        scale_x_reverse() +
        scale_fill_ucscgb() +
        ggtitle(x)

    pl(gene.clus.dist.plot, 10, 5)
    mki67.clus.dist$gene = x
    mki67.clus.dist
}) 

all.clus2 = bind_rows(all.clus.dist)
head(all.clus2)
all.clus2$short.label = as.character(all.clus2$short.label)
all.clus2 = all.clus2[all.clus2$short.label != '23_NA', ]
all.clus2$short.label = factor(all.clus2$short.label, levels = unique(all.clus2$short.label))

all.clus.plot = ggplot(all.clus2, aes(x = percentile, y = freq.percent, fill = short.label)) +
    geom_col() +
    scale_x_reverse() +
    scale_fill_ucscgb() +
    facet_wrap(.~gene) +
    theme_classic() +
    xlab('Expression percentile') +
    ylab('% cluster proportion') +
    labs(
        fill = 'Cluster'
    ) +
    theme(legend.position = 'bottom')

clus.legend = get_legend(all.clus.plot)

pl(plot_grid(clus.legend), 5, 5)


all.clus.plot.no.legend = ggplot(all.clus2, aes(x = percentile, y = freq.percent, fill = short.label)) +
    geom_col() +
    scale_x_reverse() +
    scale_fill_ucscgb() +
    facet_wrap(.~gene) +
    theme_classic() +
    xlab('Expression percentile') +
    ylab('% cluster proportion') +
    labs(
        fill = 'Cluster'
    ) +
    theme(legend.position = 'none')


#PANEL - Ma gene dynamics plots
pl(all.clus.plot, 5, 5)

pl(all.clus.plot.no.legend, 5, 5)

############################
#UMAPS 
############################

ma.genes
mac.clus2

#PANEL - UMAP expression plots
all.gene.umaps = lapply(ma.genes, function(x){
    print(x)
    x2 = mac.clus2[c('UMAP_1', 'UMAP_2', x)]
    x2 = x2[order(x2[[x]]), ]

    plot1 = ggplot(x2,
        aes_string(x = 'UMAP_1', y = 'UMAP_2', color = x)
    ) +
    geom_point() +
    scale_color_gradient(
        low = '#ffffff',
        high = '#ff0000'
    ) +
    ggtitle(x) +
    xlab('UMAP1') +
    ylab('UMAP1') +
    theme(
        axis.ticks = element_blank(),
        axis.text = element_blank()
    )
    plng(plot1, 300, 300)
})






############################
#GRANULAR EXERSIZE 
############################

#The cluster frequency analysis above gives us a broader look
#but doesn't tell us exactly how phenotypically diverse
#macrophages expressing a particular marker are,
#as two clusters might be closely related.
#we could make some trees to assess 

#thinking more about this, some of the differences 
#would come from the influence of batch correction (or lack thereof)
#might be best to stick with the above analysis alone.




############################
#PLOTS 
############################

#lapply(ma.genes, function(x){
    ##if(x == 'LYVE1') browser()
    ##mac.clus2 = arrange(mac.clus2, by = eval(x))
    #mac.clus2 = mac.clus2[order(mac.clus2[[x]]), ]
    #plot1 = ggplot(mac.clus2, aes_string(x = 'UMAP_1', y = 'UMAP_2', color = eval(x))) + 
        #geom_point() +
        #facet_wrap(.~short.label) +
        #ggtitle(x) +
        #scale_color_gradient(
            #low = '#FFFFFF',
            #high = '#ff0000'
        #)

    #plng(plot1, 1000, 1000)
#})

############################
#OTHER EXAMPLES 
############################

#DefaultAssay(mac.int) = 'SCT'
#col1 = GetAssayData(mac.int)['CCL2', ]
#mac.clus3 = bind_cols(mac.clus2, col1)
#colnames(mac.clus3)[ncol(mac.clus3)] = 'CCL2'

#mac.clus3 = mac.clus3[order(mac.clus3$CCL2), ]

#plot2 = ggplot(mac.clus3, aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'CCL2')) + 
    #geom_point() +
    #facet_wrap(.~short.label) +
    #ggtitle('CCL2') +
    #scale_color_gradient(
        #low = '#FFFFFF',
        #high = '#ff0000'
    #)

#plng(plot2, 1000, 1000)

#sig.markers = read.tsv('~/work/ucl/data/ucl.projects/macrophage/all.sig.markers.per.cluster.tsv')
#sig.markers = sig.markers %>%
    #filter(avg_log2FC > 0)

#head(arrange(sig.markers, -avg_log2FC), n = 20)


#head(sig.markers)






############################
#HCLUST TEST / VALIDATION 
############################
#want to validate the broad distribution of some markers

#mac.clus2

#cells.to.test = mac.clus2 %>%
    #arrange(short.label, LYVE1) %>%
    #group_by(short.label) %>%
    #summarise(
        #top.lyve = tail(LYVE1, n = 1),
        #top.cell = tail(cellid, n = 1),
        #bottom.lyve = head(LYVE1, n = 1),
        #bottom.cell = head(cellid, n = 1)
    #) %>%
    #mutate(
        #top.cell2 = paste0(short.label, top.cell, '_HI'),
        #bottom.cell2 = paste0(short.label, bottom.cell, '_LO')
    #)

#cells.to.test

#DefaultAssay(mac.int) = 'integrated'
#mac.sub = GetAssayData(mac.int)[, c(cells.to.test$top.cell, cells.to.test$bottom.cell)]

#head(mac.sub)
#colnames(mac.sub) = c(cells.to.test$top.cell2, cells.to.test$bottom.cell2)

#dim(mac.sub)

#hc1 = hclust(dist(t(mac.sub)))

#pl(plot(hc1), 30, 20)






############################
#CHECK CLUSTERING 
############################

#mac.clus = mac.clus %>%
    #filter(short.label != '23_NA')

#samp.ids = mac.clus %>%
    #group_by(short.label) %>%
    #summarise(samp = sample(cellid, 100, F))

#DefaultAssay(mac.int) = 'integrated'
#mac.sub = GetAssayData(mac.int)[, samp.ids$samp]
#mac.sub = t(mac.sub)

#rownames(mac.sub) = 
    #paste0(
        #as.character(mac.clus$short.label[match(rownames(mac.sub), mac.clus$cellid)]),
        #rownames(mac.sub)
    #)

#mac.sub.dist = dist(mac.sub)
#mac.h = hclust(mac.sub.dist)

#pl(plot(mac.h), 500, 500)

#pl(mac.h, 50, 50)


#data(iris)
#head(iris)

#head(mac.sub)
#cellids = colnames(mac.sub)

#df1 = as.data.frame(t(mac.sub))
#df1$cluster = mac.clus$cluster[match(rownames(df1), mac.clus$cellid)]

#colnames(df1) = gsub('-', '_', colnames(df1))
#df1$cluster = as.numeric(as.character(df1$cluster))

#rf1 = randomForest(
    #cluster ~ .,
    #data = df1,
    #importance = T,
    #proximity = T
#)

#library(xgboost)
#library(methods)
#library(caret)

#cluster.labels = as.character(mac.clus$short.label[match(
    #rownames(mac.sub),
    #mac.clus$cellid
#)])

#match(rownames(mac.sub), )

#part1 = createDataPartition(cluster.labels, p = .8, list = F, times = 1)
#part1 = as.numeric(part1)

#train1 = mac.sub[part1, ]
#train.labels = cluster.labels[part1]

#test1 = mac.sub[-part1, ]
#test.labels = cluster.labels[-part1]
#test.labels = as.numeric(as.factor(test.labels))

#train.labels = as.numeric(factor(train.labels))

#bst = xgboost(
    #data = train1,
    #label = train.labels,
    #max.depth = 5,
    #eta = 1,
    #nthread = 2,
    #nrounds = 50,
    #objective = 'reg:squarederror'
#)


#xgboost.df = data.frame(
    #test.prediction = predict(bst, test1),
    #actual.value = test.labels
#)

#xgboost.df

#train.labels
#xgboost.df$predict.label = train.labels[round(xgboost.df$test.prediction)]

#train.labels[round(xgboost.df$test.prediction)]

#xgboost.df

############################
#RANDOM FOREST TESTING 
############################

#train1
#train1
#train2 = cbind(train1, train.labels)
#colnames(train2)

#rf.mod1 = train(
    #train.labels ~ .,
    #data = train2,
    #method = 'rf'
#)

#?train


############################
#CLUSTER UCELL SIGS 
############################

#files1 = list.files('~/work/ucl/bigdata/ucl.projects/macrophage/', pattern = 'sig.markers', full.names = T)
#files1 = files1[order(as.numeric(gsub('(^[0-9]*)\\..*', '\\1', list.files('~/work/ucl/bigdata/ucl.projects/macrophage/', pattern = 'sig.markers'))))]
#files1
#files2 = lapply(files1, readRDS)

#load.mac.clus()

#head(mac.clus)

#g1 = bind_cols(files2)
#colnames(g1) = paste0('clus', 0:23)
#mac.clus = bind_cols(mac.clus, g1)

#head(mac.clus)


#lapply(paste0('clus', 0:22), function(x){
    #print(x)
    #mac.clus = mac.clus[order(mac.clus[[x]]), ]
    #plot1 = ggplot(mac.clus,
        #aes_string(x = 'UMAP_1', y = 'UMAP_2', color = eval(x))
        #) +
        #geom_point() +
        #scale_color_viridis() +
        #facet_wrap(.~short.label) +
        #ggtitle(eval(x))
    #plng(plot1, 1000, 1000)
#})







