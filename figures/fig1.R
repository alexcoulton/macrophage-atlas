# main figures panels labelled PANEL
# supplementary figures labelled SUPP
source('~/work/ucl/scripts/ucl.projects/macrophage/mac.util.R')
load.mac.clus()
library(cowplot)
library(ggsci)

############################
#SUMMARY PLOTTING 
############################

x.lab.size = 5

head(mac.clus)
mac.clus$cancer[mac.clus$cancer == 'LUAD'] = 'LUNG'
mac.clus$cancer[mac.clus$cancer == 'LUSC'] = 'LUNG'
mac.clus$cancer[mac.clus$cancer == 'NSCLC'] = 'LUNG'
mac.clus$cancer[mac.clus$cancer == 'SCLC'] = 'LUNG'

############################
#CELLS 
############################

m.cell.sum = mac.clus %>%
    group_by(cancer, tissue) %>%
    summarise(n = n()) %>%
    arrange(-n) 

m.cell.sum2 = mac.clus %>%
    group_by(cancer) %>%
    summarise(n = n()) %>%
    arrange(-n) 

m.cell.sum$cancer = factor(m.cell.sum$cancer, levels = unique(m.cell.sum2$cancer))

cell.plot = ggplot(m.cell.sum, aes(x = cancer, y = n, fill = tissue)) + geom_col() +
    ylab('Number of cells') +
    xlab('Cancer type') + 
    theme_classic() +
    rot.lab(size = 4) +
    scale_fill_npg() +
    theme(legend.position = 'bottom')


legend1 = get_legend(cell.plot)


cell.plot.no.legend = ggplot(m.cell.sum, aes(x = cancer, y = n, fill = tissue)) + geom_col() +
    ylab('Number of cells') +
    xlab('Cancer type') + 
    theme_classic() +
    rot.lab(x.lab.size) +
    scale_fill_npg() +
    theme(legend.position = 'none')


#PANEL - number of cells per cancer type plot
#pl(cell.plot, 3, 3)

############################
#SAMPLES 
############################

m.sample.sum2 = mac.clus %>%
    group_by(cancer, tissue) %>%
    summarise(n = nunique(ss2)) %>%
    arrange(-n)

m.sample.sum3 = mac.clus %>%
    group_by(cancer) %>%
    summarise(n = nunique(ss2)) %>%
    arrange(-n)

m.sample.sum2$cancer = factor(m.sample.sum2$cancer, levels = unique(m.sample.sum3$cancer))

samples.plot = ggplot(m.sample.sum2, aes(x = cancer, y = n, fill = tissue)) + geom_col() +
    ylab('Number of samples') +
    xlab('Cancer type') + 
    theme_classic() +
    rot.lab(x.lab.size) +
    scale_fill_npg() +
    theme(legend.position = 'none')

#PANEL - number of samples per cancer type plot
#pl(samples.plot, 10, 10)

############################
#PATIENTS 
############################

m.patient.sum2 = mac.clus %>%
    group_by(cancer) %>%
    summarise(n = nunique(patient)) %>%
    arrange(-n)

m.patient.sum2$cancer = factor(m.patient.sum2$cancer, levels = unique(m.patient.sum2$cancer))

patients.plot = ggplot(m.patient.sum2, aes(x = cancer, y = n)) + geom_col() +
    ylab('Number of patients') +
    xlab('Cancer type') + 
    theme_classic() +
    rot.lab(x.lab.size) +
    scale_fill_npg()

#PANEL - number of samples per cancer type plot
#pl(patients.plot, 10, 10)

colnames(mac.clus)
unique(mac.clus$primary_met2)
unique(mac.clus[is.na(mac.clus$primary_met2), c('study')])

comb.plot1 = plot_grid(
    #legend1,
    cell.plot.no.legend,
    samples.plot,
    patients.plot,
    ncol = 3,
    axis = 'hv',
    align = 'tblr'
)

#PANEL - combined barplot horizontal
pl(comb.plot1, 8, 3)
#PANEL - horizontal legend for combined barplot
pl(plot_grid(legend1), 7, 3)

############################
#FACET 
############################

m.cell.sum$cat = 'Cell'
m.sample.sum2$cat = 'Sample'
m.patient.sum2$cat = 'Patient'

m.cc.comb = bind_rows(m.cell.sum, m.sample.sum2, m.patient.sum2)

m.cc.comb$cat = factor(m.cc.comb$cat, levels = c('Cell', 'Sample', 'Patient'))

m.cc.comb$tissue[!is.na(m.cc.comb$tissue) & m.cc.comb$tissue == 'LymphNode'] = 'Lymph node'

cc.comb.plot = ggplot(m.cc.comb, aes(x = cancer, y = n, fill = tissue)) + geom_col() +
    ylab('Count') +
    xlab('Cancer type') + 
    theme_classic() +
    rot.lab(size = 6) +
    scale_fill_npg(na.value = 'black') +
    facet_grid(cols = vars(cat), scales = 'free', space = 'free') +
    theme(
        strip.background = element_rect(fill = '#DDDDDD', color = 'transparent')
    ) +
    labs(fill = 'Tissue')

pl(cc.comb.plot, 8, 7)


#The NA values here are all from studies I downloaded. I already selected
#only the tumor samples from these studies

m.tissue.type = mac.clus %>%
    group_by(tissue) %>%
    summarise(n = n())

m.tissue.type$n2 = ""
m.tissue.type$n2[m.tissue.type$n > 50000] = 
    m.tissue.type$n[m.tissue.type$n > 50000]

m.tissue.type$tissue[m.tissue.type$tissue == 'LymphNode'] = 'Lymph node'

pieplot = ggplot(m.tissue.type, aes(x = "", y = n, fill = tissue)) +
    geom_bar(stat = 'identity', width = 1, color = 'black') +
    coord_polar('y', start = 0) +
    theme_void() +
    geom_text(
        aes(label = n2),
        position = position_stack(vjust = 0.5)
    )  +
    scale_fill_npg() +
    labs(fill = 'Tissue')


#PANEL - tissue type pie plot
pl(pieplot, 4, 4)




table(mac.clus$primary_met2)

mac.clus$primary_met3 = mac.clus$primary_met2

mac.clus$primary_met2[mac.clus$primary_met2 == 'pri'] = 
    'primary'
mac.clus$primary_met2[mac.clus$primary_met2 == 'primary/laparoscopy'] = 
    'primary'
mac.clus$primary_met2[mac.clus$primary_met2 == 'norm'] = 
    'normal'
mac.clus$primary_met2 = tolower(mac.clus$primary_met2)
mac.clus$primary_met2[mac.clus$primary_met2 == 'met'] = 
    'metastasis'
mac.clus$primary_met2[mac.clus$primary_met2 == 'metastatic'] = 
    'metastasis'

primary.met.pie = data.frame(
    Primary = sum(na.omit(mac.clus$primary_met2 == 'primary')),
    Metastasis = sum(na.omit(mac.clus$primary_met2 == 'metastasis')),
    nNA = sum(is.na(mac.clus$primary_met2))
)

pr.df2 = melt(primary.met.pie)
colnames(pr.df2) = c('Category', 'Count')

pieplot2 = ggplot(pr.df2, aes(x = "", y = Count, fill = Category)) +
    geom_bar(stat = 'identity', width = 1, color = 'black') +
    coord_polar('y', start = 0) +
    theme_void() +
    geom_text(
        aes(label = Count),
        position = position_stack(vjust = 0.5)
    )  +
    scale_fill_npg() +
    labs(fill = 'Primary / Met')


#PANEL - prim.met pie plot
pl(pieplot2, 3, 3)

############################
#CANCER TYPE DISTRIBUTION 
############################

head(mac.clus)

m.clus2 = mac.clus[c('cellid', 'cancer', 'short.label')]

m.clus3 = m.clus2 %>%
    group_by(cancer, short.label) %>%
    summarise(
        n = n()
    ) %>%
    group_by(cancer) %>%
    summarise(
        short.label = short.label,
        percent = (n / sum(n)) * 100
    )

head(m.clus3)

m.clus3 %>%
    group_by(cancer) %>%
    summarise(sum(percent))

m.clus3 = m.clus3[m.clus3$short.label != '23_NA', ]

head(m.clus3)

p1 = ggplot(m.clus3,
    aes(
        x = short.label,
        y = percent,
        color = short.label,
        fill = short.label
    )) +
    geom_col() +
    facet_grid(cols = vars(cancer)) +
    theme_bw() +
    theme(
        axis.text.x = element_blank(),
        panel.spacing.x = unit(0, 'lines')
    ) +
    scale_fill_igv() +
    scale_color_igv()

pl(p1, 10, 6)


p2 = ggplot(m.clus3,
    aes(
        x = cancer,
        y = percent,
        fill = short.label
    )) +
    geom_col() +
    scale_fill_igv() +
    theme_classic() +
    rot.lab() +
    theme(legend.position = 'bottom')

p2.no.legend = ggplot(m.clus3,
    aes(
        x = cancer,
        y = percent,
        fill = short.label
    )) +
    geom_col() +
    scale_fill_igv() +
    theme_classic() +
    rot.lab() +
    theme(legend.position = 'none') +
    xlab('Cancer type') +
    ylab('% Cluster membership')

p2.legend = get_legend(p2)
pl(plot_grid(p2.legend, 8, 4))

pl(p2.no.legend, 8, 4)

#SUPP cancer type distribution of clusters
pl(p2, 8, 8)


