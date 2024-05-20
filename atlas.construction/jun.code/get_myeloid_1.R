### Load datasets, select myeloid cells, preprocess and integrate them 
### 13 Jan 2022, Jun Murai
### last modified 18 Feb 2022, Jun Murai

library(Seurat)
library(SeuratDisk)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)

# change log
#  18 Feb: removed codes for plotting etc., written how to run part1 and part2

# Initialization

meta <- readRDS("../metadata/myeloid_metadata.rds")

studies <- unique(meta$study)

files <- c(
  "Azizi_human.h5Seurat", # "AZIZI_CELL_2018_combined5_raw_counts.h5Seurat", 
  "BRAUN.h5Seurat", 
  "chan_counts_added.rds",   # failed to get raw count - we need it 
  "Zhanglab_Myeloid_raw_counts.h5seurat", 
  "KRISHNA_RCC.rds", 
  "KIM_counts.h5Seurat", 
  "BASAL_CELL_CARCINOMA_raw_counts.h5Seurat", 
  "Li.h5Seurat", 
  "Bi.h5seurat", 
  "Pelka_CRC.h5Seurat", 
  "JERBY-AMON_CELL_2018_raw_counts.h5Seurat", 
  "QIAN_NATURECELL_2020_combined4_raw_counts.h5Seurat", 
  "VISHWAKARMA_COMMBIOL_2021_combined6_raw_counts.h5Seurat", 
  "Zilionis_tumor.h5Seurat" )

exp <- list()
nstudy <- length(studies)

# Step 1 : loading data 

for (i in 1:nstudy) {  
  if (length(grep("h5.eurat", files[i]))) {
    data <- LoadH5Seurat(files[i])
  } else {
    data <- readRDS(files[i])
  }
  cells <- meta %>% filter(study == studies[i])
  exp[i] <- subset(data, cells = convertCellIDs(cells$cellid, data, i))
  cat(paste0("### Study : ", studies[i], " (", files[i], ")\n"))
  print (data)
  print (exp[[i]])
}

stats <- getStats(meta, exp)

names(exp) <- studies
saveRDS(exp, "myeloid_dataset.rds")

# Step 2 : adding Gene Symbols   - exp 5

exp$KRISHNA <- convertEnsemblToGeneSymbol_KeepAll2(exp$KRISHNA)

# Step 3 : Gene Symbol matching  

gs2new <- readRDS("GeneSymbol2Current.rds")
for (i in 1:nstudy) {
  obj <- doConvertGeneName(obj, gs2new, shrink=T)
}

saveRDS(exp, "myeloid_genename_updated.rds")

# *** refer to /home/murai/DC/scRNA_genelists

# Step 4 : remove low quality genes and do SCTransform

for (i in 1:nstudy) {
  study <- studies[i]
  cat ("processing ", study, "...\n")
  exp[[i]]$study <- study
  exp[[i]]$percent.mt <- PercentageFeatureSet(exp[[i]], pattern = "^MT-")
  exp[[i]] <- subset(exp[[i]], subset= nFeature_RNA > 200&nCount_RNA>300)
  exp[[i]] <- SCTransform(exp[[i]], method = "glmGamPoi", vars.to.regress = "percent.mt")
}

saveRDS(exp, "myeloid_sctransformed.rds")

