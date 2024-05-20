### Load 3 additional myeloid datasets
### Last Modified 13 Oct 2022, Jun Murai
### 18 Feb 2022, Jun Murai

# ChangeLog
#   13 Oct 2022: reformatted 
#   22 Feb 2022: removed 1 study, as it is included in another study. Now we have 3.

# What this is for: Add 3 new studies, prepare them for integration, and do integration
# 
# 1. read new 3 studies and convert to appropriate Seurat objects
# 2. load SCTransformed dataset of 14 to 'exp'
# 3. combine them, convert gene id to latest gene symbol and SCTransform new studies
# 4. save the updated datasets, prepare for integration
# 5. if all is ready, do integration 

source ("/home/murai/DC/script/tools.R")
setwd ("/home/murai/DC/scRNA_new")
exp_file <- "/home/murai/DC/Myeloid/myeloid_sctransformed.rds"
exp_out  <- "/home/murai/DC/Myeloid/myeloid_sctransformed_18.rds"
gs2new <- readRDS("/home/murai/DC/scRNA/GeneSymbol2Current.rds")

### Read 3 more studies

# 15: LAMBRECHTS
cat ("processing Lambrechts ...\n");
load("LAMBRECHTS_NATURE_2018/Myeloid.Cellview.Rds")
lamb_m <- CreateSeuratObject(log2cpm)
lamb_m <- deNormalizeSeurat(lamb_m, base2=F)
lamb_m <- convertEnsemblToGeneSymbol_KeepAll2(lamb_m)
rm(featuredata, log2cpm, tsne.data)

# 15a: build LAMBRECHTS metadata
cat ("Retrieving metadata for Lambrechts ...\n")
load("LAMBRECHTS_NATURE_2018/Allsamples.Cellview.Rds")
lamb_all <- CreateSeuratObject(log2cpm)
lamb_all <- deNormalizeSeurat(lamb_all, base2=F)
lamb_all@meta.data$nCount_RNA <- colSums(lamb_all@assays$RNA@counts)
lamb_md1 <- lamb_all@meta.data
lamb_md1$cellid <- rownames(lamb_md1)
lamb_md2 <- tsne.data
lamb_md2$cellid <- rownames(lamb_md2)

lamb_md1 <- lamb_md1 %>% dplyr::select(cellid, nCount_RNA, nFeature_RNA)
lamb_md2 <- lamb_md2 %>% dplyr::select(cellid, dbCluster)
lamb_meta <- lamb_md1 %>% left_join(lamb_md2)
rm(featuredata, log2cpm, tsne.data, lamb_all, lamb_md1, lamb_md2)

# celltype: myeloid; annotation: poor; data: count recovered; gene id: gene symbol #converted from ENSEMBL

# 16: LEADER
cat ("processing Leader ...\n");
leader <- readRDS("Leader/lung_ldm_seurat.rds")
leader_clust <- read.delim("Leader/cluster_annot_cellgroup_added.txt")
leader_clust$clust <- as.character(leader_clust$clust)
leader_clust <- leader_clust %>% dplyr::select(-norm_group) %>% dplyr::rename(major=lineage, minor=sub_lineage) # , celltype=lig_rec_group) # keep as is
lm <- leader@meta.data 
lm$id <- row.names(lm)
lm <- lm %>% left_join(leader_clust)
row.names(lm) <- lm$id
lm <- lm %>% dplyr::select(-id)
identical (row.names(leader@meta.data), row.names(lm)) # should be true
leader@meta.data <- lm
leader$isMyeloid <- leader$cellgroup2 == "Mono" | leader$cellgroup2 == "Macro" | leader$cellgroup2 == "DC"   # discard later
leader_m <- leader %>% subset(isMyeloid == TRUE)

# annotation = nCount nFeature patient clusterid +cellgroups; data = count;  gene = gene symbol

# 17: MAYNARD
cat ("processing Maynard ...\n");
maynard <- readRDS("Maynard/Maynard_immune_dieted.rds")
m2g <- read.delim("Maynard/celltype2cellgroup.txt")
m2g$isMyeloid <- m2g$cellgroup2 == "DC" | m2g$cellgroup2 == "Myeloid"
mm <- maynard@meta.data %>% dplyr::rename(celltype = immune_subtype_annotation) %>% left_join(m2g)
row.names(mm) <- mm$cell_id
# identical(row.names(mm), row.names(maynard@meta.data))
maynard@meta.data <- mm
maynard_m <- subset(maynard, isMyeloid == TRUE)

# celltype = IMMUNE / selected;  anntotation = rich (patient, cell type); gene = gene symbol; data = count

# 18: WU
cat ("processing Wu ...\n");
wu <- LoadH5Seurat("Wu_NatureComm_2021/Wu_immune.h5Seurat")
wu.annot <- read.csv("Wu_NatureComm_2021/info.csv")  
wu.annot <- wu.annot %>% dplyr::select (-study)
wu.ct <- data.frame(cancer_type = c("Lung squamous cell carcinoma(LUSC)", 
  "Lung adenocarcinoma(LUAD)", "Lung non-small cell carcinoma(NSCLC)"), 
  cancer = c("LUSC", "LUAD", "NSCLC"))
wu.annot <- wu.annot %>% left_join(wu.ct) %>% dplyr::select (-cancer_type)

wu.clust <- wu.annot %>% group_by(cluster, sub_cluster) %>% summarize()
#write.table(wu.clust, file="wu_cell_def.txt", sep="\t", quote=F)
# edit it.
wu.clust <- read.table("Wu_NatureComm_2021/wu_cell_def.txt", sep="\t")
wu.annot <- wu.annot %>% left_join(wu.clust, by=c("cluster", "sub_cluster"))
wu.annot <- wu.annot %>% dplyr::rename(major=cluster, celltype=sub_cluster)
wu.annot$isMyeloid <- wu.annot$cellgroup2 == "Mono" | wu.annot$cellgroup2 == "Macro" | wu.annot$cellgroup2 == "DC"
names(wu.annot)[1] <- "id"

wu <- addMetadataToSeurat(wu, wu.annot)
wu_m <- subset(wu, isMyeloid == TRUE)

cat ("loading other studies ...\n")
exp <- readRDS(exp_file)
### exp$LAMB    <- lamb_m   ### don't add LAMB's data. It is duplicated
exp$LEADER  <- leader_m
exp$MAYNARD <- maynard_m
exp$WU      <- wu_m

# Gene Symbol matching  for new studies

cat ("updating gene symbols ...\n")
for (i in 15:17) {
  exp[[i]] <- doConvertGeneName(exp[[i]], gs2new, shrink=T)
}

cat ("saving metadata ...\n")
## save cell metadata for new studies
annots <- list()
lamb_meta$isMyeloid <- lamb_meta$dbCluster == "Myeloid"
lamb_cnv <- data.frame(dbCluster=c("Alveolar", "B_cell", "EC", "Epi", "Fibro", "Myeloid", "T_cell", "tumor"), celltype=c("Alveolar", "B", "Endo", "Epi", "Fibro", "Myeloid", "T", "Tumor"))
annots$LAMB    <- lamb_meta
annots$LAMB <- annots$LAMB %>% left_join(lamb_cnv) %>% dplyr::select(-dbCluster)
# annots$LAMB$cellgroup <- annots$LAMB$celltype
annots$LEADER  <- leader@meta.data
annots$MAYNARD <- maynard@meta.data
annots$MAYNARD <- annots$MAYNARD %>% dplyr::select(-RNA_snn_res.0.3, -RNA_snn_res.0.5, -RNA_snn_res.0.7, -RNA_snn_res.0.1, -RNA_snn_res.0.9)
annots$MAYNARD <- annots$MAYNARD %>% dplyr::select(-RNA_snn_res.1, -nonimmune_seurat_cluster, -nonimmune_general_annotation, -immune_annotation)
annots$MAYNARD <- annots$MAYNARD %>% dplyr::select(-sort_plate_number, -pfs_over_under, -pfs_month)
annots$MAYNARD <- annots$MAYNARD %>% dplyr::select(-orig.ident, -well, -Sequence_Run1, -plate)
may_conv <- data.frame(histolgy = c("Adenocarcinoma", "Squamous"), cancer = c("LUAD", "LUSC"))   # annots$MAYNARD$histolgy  "Adenocarcinoma" "Squamous"
annots$MAYNARD <- annots$MAYNARD %>% left_join(may_conv) %>% dplyr::select (-histolgy)
annots$MAYNARD <- annots$MAYNARD %>% rename(cellid = cell_id)
annots$WU      <- wu@meta.data
annots$LEADER$cellid <- rownames(annots$LEADER)
annots$WU$cellid <- rownames(annots$WU)
saveRDS(annots, "new4_cell_metadata_v1.1.rds")

# remove large objects
rm(wu, wu.annot, wu_m, maynard_m, maynard, leader_m, leader, lamb_m, annots)

## run SCTransform for additional data

cat ("running filter cells and SCTransform ...\n")

studies <- names(exp)
for (i in 15:17) {
  study <- studies[i]
  cat ("processing ", study, "...\n")
  exp[[i]]$study <- study
  exp[[i]]$percent.mt <- PercentageFeatureSet(exp[[i]], pattern = "^MT-")
  exp[[i]] <- subset(exp[[i]], subset= nFeature_RNA > 200&nCount_RNA>300)
  exp[[i]] <- SCTransform(exp[[i]], method = "glmGamPoi", vars.to.regress = "percent.mt")
}

cat ("saving studies ...\n")
saveRDS(exp, "myeloid_sctransformed_17.rds")

# restart from here --
# exp <- readRDS("myeloid_sctransformed_17.rds")

cat ("preparing for integration (using 1700 genes) ...\n");
source ("/home/murai/DC/script/integrate_functions.R")
features <- SelectIntegrationFeatures(object.list = exp, nfeatures = 1700)     # at once.
exp      <- PrepSCTIntegration(object.list = exp, anchor.features = features)  # 41g, 31g, 1m02s
anchors  <- FindIntegrationAnchors(object.list = exp, normalization.method = "SCT", anchor.features = features) 

# check anchor stats
stat <- preCheckAnchorSet(anchors)

saveRDS(anchors, "myeloid_integrate_anchors_17_1700.rds")

PREPARED <- TRUE

if (PREPARED) {  # IF PREPARED
 # run integration, if everything is ok
 # rm(exp) - will be free some memory 
 exp <- IntegrateDataW(anchors, normalization.method = "SCT")
 saveRDS(exp, "myeloid_integrated_17studies_noclust.rds")
 
 exp <- RunPCA(exp, verbose = FALSE)
 exp <- RunUMAP(exp, reduction = "pca", dims = 1:50)
 exp <- FindNeighbors(exp, reduction = "pca", dims = 1:50)
 exp <- FindClusters(exp, resolution = 0.5)
 saveRDS(exp, "myeloid_integrated_17studies.rds")
}

