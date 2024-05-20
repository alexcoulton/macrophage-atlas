### get and analyze Myeloid
### Last modified 13 Oct 2022 by Jun Murai
### 18 Feb 2022 by Jun Murai

library(Seurat)
library(SeuratDisk)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)

## build distinguishable (but not artistic) palette 
library(Polychrome)
data(palette36)
palette33 <- palette36[4:36]
names(palette33) <- NULL

setwd ("~/DC/Myeloid")
exp <- readRDS("myeloid_integrated_s17_g1700_clustered.rds")

## more detailed clustering

exp <- FindClusters(exp, resolution = 1.3)

## add full metadata to the experiment

meta.file = "~/DC/full_metadata.rds" #  "../metadata/full_metadata.rds"

meta <- readRDS(meta.file)
exp.meta <- exp@meta.data
meta <- meta%>%distinct(cellid, study, .keep_all=TRUE)  # 2068241 with 2038954, need to recheck.
exp.meta <- exp.meta %>% select(nCount_RNA, nFeature_RNA, study, integrated_snn_res.0.5)
exp.meta$cellid <- rownames(exp.meta)
exp.m2 <- exp.meta %>% left_join(meta, by=c("cellid", "study"))
identical(exp.meta$cellid, exp.m2$cellid)  # true
rownames(exp.m2) <- exp.m2$cellid
exp@meta.data <- exp.m2

### 

ct_study <- exp@meta.data %>% group_by(study, celltype) %>% summarize(num=n())
write.table(ct_study, file="Study_Celltype.txt", row.names=F, sep="\t")

saveClusterGroupPlots <- function(grouped_table, label="", num.clust=F, do.save=T, width=14) {
  tb <- summarize(grouped_table, n())
  if(label == "") label <- names(tb)[2]
  names(tb) <- c("cluster", "group", "cells")
  if (num.clust) tb$cluster <- as.numeric(as.character(tb$cluster))
  p <- ggplot(tb, aes(fill=group, y=cells, x=cluster)) + geom_bar(position="fill", stat="identity") + labs(fill=label)  # "stack"/"fill"
  if (num.clust) {
    min <- min(tb$cluster)
    max <- max(tb$cluster)
    p <- p + scale_x_continuous(breaks=min:max)
  }
  if (do.save) ggsave(file=paste0("chart_fill_",label,".png"), plot = p, width=width)
  return (p)
}

saveFPlot <- function(feature, min.cutoff=0, max.cutoff=5, ex=exp, fn="") {
  if (fn == "") fn <- feature
  fn <- paste0(fn, ".png")
  FeaturePlot(ex, features=feature, max.cutoff=max.cutoff, min.cutoff=min.cutoff, raster=F, label=T) %>% ggsave(file=fn)
}

### example to save plots
saveFPlot("ISG15") 
VlnPlot(exp, "ISG15") %>% ggsave(file="ISG15_vln.png")


# going back to exp

markers <- FindAllMarkers(exp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

### Elbow Plot
sses <- data.frame()

r <- c(0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.5, 3)
for (i in r) {
  exp  <- FindClusters(exp, resolution = i)
  sse  <- getSSEfromExp(exp)
  sses <- rbind(sses, data.frame(res = i, sse = sse))
}

clusts <- data.frame()
r <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.5, 3)
for (i in r) {
  col <- paste0("integrated_snn_res.",i)
  len <- length(unique(exp@meta.data[[col]]))
  clusts <- rbind(clusts, data.frame(res = i, cl = len))
}
sses <- sses %>% left_join(clusts)

markers55 <- FindAllMarkers(exp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# clust <- clust %>% select(cellid, starts_with("integrated_snn"))
# exp@meta.data <- exp@meta.data %>% left_join(clust)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Idents(exp) <- "integrated_snn_res.3"

cell.stats <- exp@meta.data %>% group_by(integrated_snn_res.3, cellgroup2) %>% summarize(n = n())
cell.st2 <- cell.stats %>% pivot_wider(id_cols=integrated_snn_res.3, names_from=cellgroup2, values_from=n)
cell.st2$types3 <- cell.st2$DC + cell.st2$Macro + cell.st2$Mono
cell.st2$DC.pct <- cell.st2$DC / cell.st2$types3
cell.st2$Macro.pct <- cell.st2$Macro / cell.st2$types3
cell.st2$Mono.pct  <- cell.st2$Mono / cell.st2$types3

cell.st2$judge <- ifelse(cell.st2$DC.pct>0.75, "DC", ifelse(cell.st2$Macro.pct>0.75, "Macro", ifelse(cell.st2$Mono.pct>0.75, "Mono", "")))
write.table(cell.st2, file="cell_decision.txt", sep="\t", row.names=F, col.names=TRUE)

addtypes <- read.delim("add_to_celltypes.txt")
meta <- exp@meta.data
meta$integrated_snn_res.3 <- meta$integrated_snn_res.3 %>% as.character() %>% as.integer()  
## don't as.integer() to factors. it will 
meta <- meta %>% left_join(addtypes, by=c("integrated_snn_res.3", "cellgroup2"))

meta[meta$cellgroup2 == "Macro",]$add.Mac = 1
meta[meta$cellgroup2 == "Mono",]$add.Mono = 1
meta[meta$cellgroup2 == "DC",]$add.DC = 1

rownames(meta) <- meta$cellid

saveRDS(exp@meta.data, "myeloid_metadata_10Mar.rds")
saveRDS(exp, "myeloid_integrated_update_10Mar.rds")

## resume from here 
# exp <- readRDS("myeloid_integrated_s17_g1700_clustered.rds")
# meta <- readRDS("myeloid_metadata_10Mar.rds")
# exp@meta.data <- meta


# replacing names(meta) ... integrated to myeloid
configure_subset <- function(ex) {
 ex@reductions$umap_myeloid <- ex@reductions$umap
 ex@reductions$pca_myeloid  <- ex@reductions$pca

 ### renaming metadata integrated_snn_res_* to myeloid_snn_res_*
 meta.names <- gsub("integrated_snn_res", "myeloid_snn_res", names(ex@meta.data))
 if (sum(duplicated(meta.names))) {
   stop (cat("failed: exp already have myeloid_snn_res, number of duplicates: ", sum(duplicated(meta.names))))
 }
 names(ex@meta.data) <- meta.names

 DefaultAssay(ex) <- "integrated"
 cat ("[PCA]\n");  ex <- RunPCA(ex, verbose = FALSE)
 cat ("[UMAP]\n"); ex <- RunUMAP(ex, reduction = "pca", dims = 1:50)
 cat ("[SNN]\n");  ex <- FindNeighbors(ex, reduction = "pca", dims = 1:50)
 cat ("[Clusters]\n");  ex <- FindClusters(ex, resolution = 0.5)
 
 DefaultAssay(ex) <- "SCT"
 return(ex)
}


cat ("processing DC ...\n")

exp.dc   <- subset(exp, add.DC  == 1)
exp.dc <- configure_subset(exp.dc)
saveRDS(exp.dc, "integrated_s17_dc.rds")

cat ("processing Macrophage ...\n")

exp.mac  <- subset(exp, add.Mac == 1)
exp.mac <- configure_subset(exp.mac)
saveRDS(exp.mac, "integrated_s17_macro.rds")

cat ("processing Monocyte ...\n")

exp.mono <- subset(exp, add.Mono == 1)
exp.mono <- configure_subset(exp.mono)
saveRDS(exp.mono, "integrated_s17_monocyte.rds")


