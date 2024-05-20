#### aggregate Metadata
### last updated by Jun Murai, 23 Feb, 2021
### by Jun Murai, 30 Nov, 2021

library(tidyr)
library(dplyr)
library(Seurat)

options(tibble.print_min=30)
options(width=200)

chan_exp  <- readRDS("~/DC/scRNA/chan_counts_added.rds") # raw data with annotation for all genes
chan_exp_meta <- chan_exp@meta.data

chan_meta <- readRDS("~/DC/metadata/CHAN_v1.rds")        # this is the combined table of the immune cell annotation and the myeloid cell annotation
chan_meta <- chan_meta %>% select(-n_genes_by_counts.x, -log1p_n_genes_by_counts.x, -total_counts.x, -log1p_total_counts.x, -mito_frac.x, -RBP_frac.x, -batch.x, -patient.x, -tissue.x, -treatment.x, -procedure.x, -histo.x, -n_genes_by_counts.y, -log1p_n_genes_by_counts.y, -total_counts.y, -log1p_total_counts.y, -mito_frac.y, -RBP_frac.y, -batch.y, -patient.y, -tissue.y, -treatment.y, -procedure.y, -histo.y) %>% rename(cellid = Cell)
chan_exp_meta$cellid <- rownames(chan_exp_meta)
chan_meta <- chan_exp_meta %>% full_join(chan_meta, by="cellid")

# 2: nested ifelse - ok but need to do as.char
chan_meta$celltype <- ifelse(is.na(chan_meta$cell_type.y), ifelse(is.na(chan_meta$cell_type.x), as.character(chan_meta$cell_type_fine), chan_meta$cell_type.x), chan_meta$cell_type.y)
# chan_meta$celltype  <- ifelse(is.na(chan_meta$cell_type.y), chan_meta$cell_type.x, chan_meta$cell_type.y)
# chan_meta$celltype  <- ifelse(is.na(chan_meta$celltype), chan_meta$cell_type_fine, chan_meta$celltype)

chan_meta <- chan_meta %>% mutate(celltype=gsub(celltype, pattern="Mo/M.", replacement="Mo/Mp", ignore.case=F)) # CHAN's mo/m<ph> couldn't be recognized
for (i in names(chan_meta)) {
  if (is.factor(chan_meta[[i]])) {
    chan_meta[[i]] <- as.character(chan_meta[[i]])
  }
}

saveRDS(chan_meta, "~/DC/metadata/CHAN.rds")

# chan_meta %>% group_by(celltype, celltype2, cell_type_coarse, cell_type_fine, cell_type.x, cell_type.y) %>% summarize %>% write.table("~/Chan_cell_type_processed.txt", sep="\t")

### result
#> dim(chan_meta)
#[1] 147193     25
#> head(chan_meta, 2)
#  ngenes  libsize  mito_frac  RBP_frac       batch patient           tissue
#1    255 2.594393 0.05597964 0.3155216 RU1311A_T_1 RU1311A             lung
#2   3720 4.285895 0.06010872 0.4016567      RU1215  RU1215 pleural_effusion
#         treatment     procedure histo cell_type_coarse cell_type_fine
#1 Platinum Doublet        Biopsy  SCLC         Lymphoid         T cell
#2            Naive Thoracentesis  SCLC       Epithelial         SCLC-N
#  cell_type_general clusters cell_type_med    H_knn                      cellid
#1            Immune       23        T cell 1.309847 RU1311A_T_1_165945547864806
#2        Epithelial       34        SCLC-N 0.000000      RU1215_192110488599350
#  histology_subtype cell_type.x clusters.x clusters_fine cell_type.y clusters.y
#1            SCLC-A     T cells         10            24        <NA>         NA
#2              <NA>        <NA>         NA            NA        <NA>         NA
#  histo_tissue celltype
#1         <NA>  T cells
#2         <NA>   SCLC-N
