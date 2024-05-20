### Load and Plot Macrophage dataset
### last update: 29 Jul 2022 by Jun Murai
### by Jun Murai, 6 Jun, 2022

# requires R 4.1 and Seurat 4.1 
# 'conda activate R4.1' on CS cluster to use a suitable environment 

# Required memory may be around 40GB

# 1. load libraries
source ("tools.R")
source ("tools2.R")
library(tidyverse)
library(alluvial)

# 2. load data files
exp    <- readRDS("integrated_macrophage.rds")
meta   <- readRDS("integrated_metadata_s17_220729.rds")  
prj    <- readRDS("project_info.rds")
load("markers.rda")                    # read markers16, markers20, markers27
pal49 <- distantPal49()
pal80 <- distantPal80()

# 3. update experiment
if (identical(row.names(meta), row.names(exp@meta.data)) == FALSE) {
  cat ("row names of macrophage dataset and metadata is not identical.\n")
  exit(0)
}
exp@meta.data <- meta      # good to add latest metadata
Project(exp)  <- prj$title
DefaultAssay(exp) <- "SCT"

  ## exp <- PrepSCTFindMarkers(exp) : doesn't work. our integrated dataset is not compatible with FindMarker of Seurat 4.1

# 4. build base plot for blank map

dir.create("result", showWarnings=FALSE)
setwd("result")

(DimPlot(exp, cols=palette33, label=F, raster=T, raster.dpi=c(1920,1920)) + coord_fixed()) %>% ggsave(file="00_macrophage_20_2k_NL.png", width=9, height=7.2)

Idents(exp) <- exp$integrated_snn_res.0.8; 
(DimPlot(exp, cols=palette33, label=FALSE, raster=TRUE, raster.dpi=c(960,960)) + coord_fixed()) %>% ggsave(file="00_macrophage_27_NL.png", width=9, height=7.2)
(DimPlot(exp, cols=palette33, label=TRUE, raster=TRUE, raster.dpi=c(960,960)) + coord_fixed()) %>% ggsave(file="00_macrophage_27.png", width=9, height=7.2)

Idents(exp) <- exp$integrated_snn_res.1.6; 
(DimPlot(exp, cols=pal49, label=FALSE, raster=TRUE, raster.dpi=c(960,960)) + coord_fixed()) %>% ggsave(file="00_macrophage_49_NL.png", width=10, height=7.2)
(DimPlot(exp, cols=pal49, label=TRUE, raster=TRUE, raster.dpi=c(960,960)) + coord_fixed()) %>% ggsave(file="00_macrophage_49.png", width=10, height=7.2)

Idents(exp) <- exp$integrated_snn_res.3;   
(DimPlot(exp, cols=pal80, label=TRUE, raster=TRUE, raster.dpi=c(960,960)) + coord_fixed()) %>% ggsave(file="00_macrophage_80.png", width=11, height=7.2)
(DimPlot(exp, cols=pal80, label=FALSE, raster=TRUE, raster.dpi=c(960,960)) + coord_fixed()) %>% ggsave(file="00_macrophage_80_NL.png", width=11, height=7.2)

Idents(exp) <- exp$integrated_snn_res.0.5;  
(DimPlot(exp, cols=palette33, label=TRUE, raster=TRUE, raster.dpi=c(960,960)) + coord_fixed()) %>% ggsave(file="00_macrophage_20.png", width=9, height=7.2)
(DimPlot(exp, cols=palette33, label=FALSE, raster=TRUE, raster.dpi=c(960,960)) + coord_fixed()) %>% ggsave(file="00_macrophage_20_NL.png", width=9, height=7.2)

# 5. build Feature Plots

genes_nt <- c("ACP5", "APOC1", "APOE", "BAG3", "C1QA", "C1QC", "C3", "CCL18", "CCL2", "CCL20", "CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CDC20", "CENPF", "CRIP1", "CXCL1", "CXCL10", 
  "CXCL13", "CXCL9", "CYP27A1", "DNAJB1", "ENSG00000117289", "FCER1A", "FN1", "FOLR2", "HLA-DQB2", "HSPA1A", "HSPA1B", "HSPA6", "IDO1", "IFIT1", "IGKC", "IL10", "IL1B", "ISG15", "MARCO", 
  "MIF", "MKI67", "MRC1", "MT-ATP6", "MT-ND3", "MT1G", "MT1X", "MT2A", "MTRNR2L12", "NLRP3", "PLD4", "PRDM1", "S100A4", "S100A8", "SAR1A", "SPP1", "STMN1", "TGFB1", 
  "TGFB2", "TGFB3", "THBS1", "TK1", "TOP2A", "VCAN", "CCR6", "CCR7", "CD80", "CD86", "CD14")

mp_marker1 <- c("FABP4", "CES1", "CD52", "MCEMP1", "APOE", "CYP27A1", "ACP5", "APOC1", "CTSD", "SELENOP", "SLC40A1", "FOLR2", "PLTP", "C1QB", "C1QA", "SPP1", "RNASE1", "LGMN", "GPNMB", "MMP9", "A2M", "THBS1", "EREG", "AREG", "RNF144B", "IL1B", "FN1", "MARCO", "S100A4", "FCER1A", "C3", "PLD4", "HLA-DQA2", "CXCL9", "CXCL10", "GBP1", "WARS1", "CXCL11", "TAP1", "STAT1", "VCAN", "CCL2", "CCL20", "CTSL", "CTSB", "RPS20", "IGKC", "QKI", "RPL6", "RPL9")

mp_marker2 <- c("G3BP1", "MT-ND2", "MT-ND4", "MT-CO3", "STAB1", "DAB2", "ARL4C", "ISG15", "IFIT3", "IFIT1", "MT1G", "MT1X", "MT2A", "MT1H", "MKI67", "STMN1", "CENPF", "TOP2A", "CDC20", "TUBB", "CCL4L2", "CCL3L3", "CCL8", "HSPA6", "HSPA1B", "DNAJB1", "HSPA1A", "BAG3", "IER5", "JUN", "FOS", "GZMB", "CCL5", "IL32", "NKG7", "CD3D", "CD8A", "CXCL3", "FCN1", "LYZ", "RP11-1143G9.4", "S100A8", "S100A9", "AZU1", "RGS5", "IGFBP7", "SPARC", "COL1A1", "MGP", "ACTA2")

genes3 <- c("ITGAM", "ITGAX", "FABP5", "VEGFA")
genes4 <- c("PTGS2", "STAT3", "COL1A2", "COL3A1", "GBP5", "CCL3", "CXCL2", "CXCL3", "CXCL8",  "CD274", "CX3CR1", "CD4", "FOXP3")
genes_tcell <- c("CD2", "CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "FOXP3", "PDCD1")

genes <- c(genes_nt, mp_marker1, mp_marker2, genes3, genes4, genes_tcell) %>% unique()

for (i in c(20, 27, 49)) {
  info(i)
  for (g in genes) {
    saveFPlotR2(g)
  }
}

   mito <- exp@meta.data
   mito <- mito[!is.na(mito$percent_mito),]
   mito <- row.names(mito)
   fp <- FeaturePlot2(exp, features="percent_mito", raster=TRUE, raster.dpi=c(1920,1920), cols=c("#f4f4ee", "blue"), label=TRUE, order=TRUE, cells=mito   )  + coord_fixed()
   ggsave(file="percent_mito.png", fp, width=9, height=7.2)

genes3.0 <- c("LYVE1", "ACTA2", "AZU1", "CCR7", "CCR6", "CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CD80", "CDC20", "CENPF", "HES1", "ITGAL") 
genes3.0 <- c("CD69", "ARHGEF5", "PDCD1LG2", "PDCD1", "ARG1", "CD274", "CXCL5", "PIM2", "CD2")
genes3.0 <- c("FOXP3")

for (i in c(20, 27, 49)) {
  info(i)
  for (g in genes3.0) {
    saveFPlotR2(g, max.cutoff=3)
  }
}

genes2.0 <- c()

for (g in genes2.0) {
  saveFPlotR2(g, max.cutoff=2)
}

genes <- c(genes, genes3.0, genes2.0) %>% unique()

# 6. build dot plots

  # png("dotplot.png", width=1600, height=800); DotPlot(exp, features = genes) + RotatedAxis(); dev.off()

buildOneDotPlot <- function(idents, png, features) {
  currIdents <- Idents(exp)
  Idents(exp) <- idents
  png(png, width=1600, height=800); 
  print(DotPlot(exp, features = features) + RotatedAxis()) ## %>% ggsave(file=png, width=1600, height=800, units="px")
  dev.off()
  cat ("Dot Plot written to:", png, ".\n")
  Idents(exp) <- currIdents
}

buildDotPlots <- function(title, features) {
  buildOneDotPlot(idents = exp$integrated_snn_res.0.4, png = paste0(title, "16.png"), features = features)
  buildOneDotPlot(idents = exp$integrated_snn_res.0.5, png = paste0(title, "20.png"), features = features)
  buildOneDotPlot(idents = exp$integrated_snn_res.0.8, png = paste0(title, "27.png"), features = features)
  buildOneDotPlot(idents = exp$integrated_snn_res.1.6, png = paste0(title, "49.png"), features = features)
}

# non-target
buildDotPlots("markerA_", mp_marker1)
buildDotPlots("markerB_", mp_marker2)

Idents(exp) <- exp$integrated_snn_res.0.5

for (g in genes) {
  VlnPlot(exp, g) %>% ggsave(file=paste0(g, "_vln.png"))
}

# Idents(exp) <- exp$integrated_snn_res.0.5
# markers20 <- FindAllMarkers(exp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)  # 0.5

saveClusterGroupPlotsAndReturnCrossTable <- function(grouped_table, label="", num.clust=F, do.save=T, width=14, celltype="mp") {
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
  tbp <- tb %>% pivot_wider(id_cols=cluster, names_from=group, values_from=cells)
  if (do.save) ggsave(file=paste0("chart_fill_",celltype,"_",label,".png"), plot = p, width=width)
  return (tbp)
}

# 7. build cross tables

tables20 <- list()
tables20$study <- saveClusterGroupPlotsAndReturnCrossTable(exp@meta.data %>% group_by(integrated_snn_res.0.5, study), label="study", width=14, num.clust=T)

tables20$cell  <-saveClusterGroupPlotsAndReturnCrossTable(exp@meta.data %>% group_by(integrated_snn_res.0.5, cellgroup), label="cellgroup", width=14, num.clust=T)
tables20$cell2 <-saveClusterGroupPlotsAndReturnCrossTable(exp@meta.data %>% group_by(integrated_snn_res.0.5, cellgroup2), label="cellgroup2", width=14, num.clust=T)
tables20$norm_site <-saveClusterGroupPlotsAndReturnCrossTable(exp@meta.data %>% filter(tissue=="Normal") %>% group_by(integrated_snn_res.0.5, cancer), label="Normal_Tissue_Types", width=14, num.clust=T)
tables20$cancer    <-saveClusterGroupPlotsAndReturnCrossTable(exp@meta.data %>% group_by(integrated_snn_res.0.5, cancer), width=14, num.clust=T)
tables20$tissue_type <-saveClusterGroupPlotsAndReturnCrossTable(exp@meta.data %>% group_by(integrated_snn_res.0.5, tissue), width=14, num.clust=T)

write.table(tables20$study, "Mp3_CT_studies_20.txt", row.names=F, col.names=T, sep="\t", na="")
write.table(tables20$cell, "Mp3_CT_cellgroup_20.txt", row.names=F, col.names=T, sep="\t", na="")
write.table(tables20$cell2, "Mp3_CT_cellgroup2_20.txt", row.names=F, col.names=T, sep="\t", na="")
write.table(tables20$norm_site, "Mp3_CT_normal_sites_20.txt", row.names=F, col.names=T, sep="\t", na="")
write.table(tables20$cancer, "Mp3_CT_cancer_20.txt", row.names=F, col.names=T, sep="\t", na="")
write.table(tables20$tissue_type, "Mp3_CT_tissuetype_20.txt", row.names=F, col.names=T, sep="\t", na="")

markers20_top <- markers20 %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)
write.table(markers20, file="markers20_all.txt", sep="\t", row.names=F, col.names=T)
write.table(markers20_top, file="markers20_top.txt", sep="\t", row.names=F, col.names=T)

markers27_top <- markers27 %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)
write.table(markers27_top, file="markers27_top.txt", sep="\t", row.names=F, col.names=T)

# 8. Alluvial Plot between resolution

allu <- exp@meta.data %>% select(integrated_snn_res.0.4, integrated_snn_res.0.5, integrated_snn_res.0.8, cellgroup)
allu <- allu %>% group_by(integrated_snn_res.0.4, integrated_snn_res.0.5, integrated_snn_res.0.8, cellgroup) %>% summarize(n=n())
png("alluvial_plot.png")
alluvial(allu %>% select(-n), col=palette36[as.numeric(allu$integrated_snn_res.0.4)+1], freq=allu$n)
dev.off()

allu2 <- data.frame(cl16 = as.numeric(as.character(exp$integrated_snn_res.0.4)), cl20 = as.numeric(as.character(exp$integrated_snn_res.0.5)), cl27 = as.numeric(as.character(exp$integrated_snn_res.0.8)), cellgroup=exp$cellgroup)
allu2 <- allu2 %>% group_by(cl16, cl20, cl27, cellgroup) %>% summarize(n=n())
png("alluvial_plot2.png")
alluvial(allu2 %>% select(-n), col=palette36[as.numeric(allu2$cl16)+1], freq=allu2$n)
dev.off()

# 9. celltype plot for each paper

xlims <- c(min(exp@reductions$umap@cell.embeddings[,1]), max(exp@reductions$umap@cell.embeddings[,1]))
ylims <- c(min(exp@reductions$umap@cell.embeddings[,2]), max(exp@reductions$umap@cell.embeddings[,2]))

# for (i in c("BI", "BRAUN", "CHENG", "KIM", "KRISHNA", "ZILI")) {

for (i in unique(exp$study)) {
  cell.select <- exp@meta.data %>% filter(study == i) %>% rownames();   cat (i, "...");
  (DimPlot(exp, cells=cell.select, cols=palette33, group.by="celltype", shuffle=TRUE, raster=FALSE) + coord_fixed() + scale_x_continuous(limits = xlims) + scale_y_continuous(limits = ylims) + labs(title=NULL, colour=i)) %>% ggsave(file=paste0("02_",i,"_celltype.png"), width=10, height=7.2)
}
cat ("\n")

i<-"BRAUN"

  cell.select <- exp@meta.data %>% filter(study == "BRAUN") %>% rownames();   cat (i, "...");
  (DimPlot(exp, cells=cell.select, cols=palette33, group.by="minor", shuffle=TRUE, raster=FALSE) + coord_fixed() + scale_x_continuous(limits = xlims) + scale_y_continuous(limits = ylims) + labs(title=NULL, colour=i)) %>% ggsave(file=paste0("02_BRAUN_minortype_r2_NonRaster.png"), width=10, height=7.2)

# exp@reductions$umap@cell.embeddings

# 10. cell type cross table data

exp@meta.data %>% group_by(study, celltype, integrated_snn_res.0.5) %>% summarize(n()) %>% write.table("mac_study_celltype_clust20.txt", sep="\t", row.names=F)

# 11. DimPlot

saveDimPlot <- function(grp) {
  (DimPlot(exp, group.by=grp, cols=palette33, raster=TRUE, raster.dpi=c(1920, 1920), shuffle=T) + coord_fixed() + scale_x_continuous(limits = xlims) + scale_y_continuous(limits = ylims) + labs(title=NULL, colour=grp)) %>% ggsave(file=paste0("0_",grp,".png"), width=9, height=7.2)
}

saveDimPlot("study")
saveDimPlot("cellgroup2")
saveDimPlot("cellgroup")
saveDimPlot("cancer")
saveDimPlot("tissue")
saveDimPlot("cancer_site")

saveDimPlot("response_new")
saveDimPlot("platform_new")

print(gc())

# 12. save updated macrophage Seurat object file
saveRDS(exp, "integrated_macrophage_220729.rds")
