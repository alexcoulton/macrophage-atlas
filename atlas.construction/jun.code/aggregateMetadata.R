### aggregate Metadata
### last updated by Jun Murai, 1 Mar, 2021
### by Jun Murai, 30 Nov, 2021

# memory usage: ~5GB 

rm(list=ls())
library(tidyr)
library(dplyr)

# update note:
#   21 Feb: added 4 studies
#   28 Feb: fixed duplicated metadata issue and added code to prevent celltype2cellgroup causing it

# SUBROUTINE

first_time_build_mt <- function() {
  cat ("First time build: \n");
  mt <- list()
  rdss <- c("BRAUN", "CHAN", "CHENG", "JERBY", "KIM", "KRISHNA", 
            "Pelka", "Pelka_clust", "QIAN", "VISH", "YOST_BASAL", "YOST_SCC")
  for (i in rdss) {
    cat ("loading ", i, "...\n")
    mt[[i]] <- readRDS(paste0(i,".rds"))
  }
  cat ("loading others ...\n")
  load("LI.rds") # rda
  load("ZILI.rds") #rda
  mt$LI <- li
  mt$ZILI <- zi
  
  mt$BI <- read.table("BI_Final_SCP_Metadata.txt", sep="\t", header=T)
  mt$AZIZI <- read.table("Azizi_info.csv", sep=",", quote='"')
  mt$AZIZI_clust <- read.table("Azizi_clusters.txt", sep="\t", header=T)
  
  m3 <- readRDS("new3_cell_metadata_v1.2.rds")
  rdss <- names(m3) 
  for (i in rdss) { 
    mt[[i]] <- m3[[i]]
  }
  
  rdss <- names(mt)
  for (i in rdss) {
    mt[[i]]$study <- i
  }
  rm(li, zi)
  
  # transferred to celltype2cellgroup
  mt$WU <- mt$WU %>% select (-cellgroup, -cellgroup2)
  mt$MAYNARD <- mt$MAYNARD %>% select (-cellgroup, -cellgroup2)
  
  cat ("saving combined file ... ")
  saveRDS(mt, file="metadata_combined_20.rds")
  cat ("done.\n")
  return(mt)
}

# BEGIN

cat ("::: aggregate Metadata version 0.20220222\n")

# check consistency of celltype2cellgroup
ct2cg <- read.table("celltype_cellgroup.txt", sep="\t", header=T)
if (length(ct2cg$celltype) > length(unique(ct2cg$celltype))) {
  cat ("### Duplicated celltype:\n")
  print(ct2cg[duplicated(ct2cg$celltype),]$celltype)
  stop("duplicate found in celltype_cellgroup.txt. please check!")
}

# read data (from combined or from each rds)
if (!file.exists("metadata_combined_20.rds")) {
  cat ("building metadata ...\n")
  mt <- first_time_build_mt()
} else {
  cat ("reading metadata_combined_20.rds ...\n")
  mt <- readRDS("metadata_combined_20.rds")
}

# inspect all of 
# for (i in rdss) {  cat(i, names(mt[[i]]), "\n") }
# now each table has column 'study' = name in list 'mt'
# first try - bind anything with bind_rows

mt$CHENG$cellid <- rownames(mt$CHENG)
mt$LI$cellid <- rownames(mt$LI)
mt$ZILI$cellid <- paste0(mt$ZILI$Library, "_", mt$ZILI$Barcode)
mt$Pelka_clust <- mt$Pelka_clust %>% rename(cellid = sampleID, celltype = cl295v11SubFull) %>% select(-study)

cat ("formatting and combining metadata ...\nA")
mtr <- list()
azicl <- mt$AZIZI_clust %>% select(ClusterID, final_annot, dmap, garvan, removed, annot_bulk) %>% rename(cluster = ClusterID, celltype = final_annot, major=annot_bulk)
mtr$AZIZI <- mt$AZIZI %>% left_join(azicl, by="cluster") %>% rename(cluster_id = cluster) %>% select(-replicate)
mtr$AZIZI$orig_cellid <- as.character(mtr$AZIZI$cellid)
mtr$AZIZI$cellid <- paste0("c", mtr$AZIZI$orig_cellid)
mtr$AZIZI$cancer = "BC"
mtr$AZIZI$stage  <- ifelse(mtr$AZIZI$patient == "BC2" | mtr$AZIZI$patient == "BC4", "4", "1-3")
mtr$AZIZI$sample <- paste0(mtr$AZIZI$patient, "_", mtr$AZIZI$tissue)
mtr$AZIZI <- mtr$AZIZI %>% select(-removed, -orig_cellid)

cat ("B")
mtr$BRAUN <- mt$BRAUN %>% rename(cellid = Cell_barcode, celltype = ClusterName_ImmuneCells, nFeature=nFeature_RNA, nCount = nCount_RNA, tissue=Tumor_Normal, major=ClusterName_AllCells, batch=Batch, stage=Stage, sample = Sample) %>% 
  select(-ClusterNumber_AllCells, -ClusterNumber_ImmuneCells, -ClusterNumber_TCells, -ClusterNumber_MyeloidCells)
mtr$BRAUN$cancer = "ccRCC"
mtr$BRAUN$minor <- ifelse(mtr$BRAUN$ClusterName_MyeloidCells == "-1", mtr$BRAUN$ClusterName_TCells, mtr$BRAUN$ClusterName_MyeloidCells)
mtr$BRAUN <- mtr$BRAUN %>% select (-ClusterName_MyeloidCells, -ClusterName_TCells, -Included_in_CD8_trajectory, -Included_in_myeloid_trajectory )
mtr$BRAUN$batch <- as.character(mtr$BRAUN$batch)
mtr$BRAUN$patient <- sub("_.*", "", mtr$BRAUN$sample)
mtr$BRAUN$stage <- ifelse(mtr$BRAUN$stage == "T_metast", "4", ifelse(mtr$BRAUN$stage == "T_locAdv", "3", ifelse(mtr$BRAUN$patient == "S8", "2", ifelse(mtr$BRAUN$stage == "T_early", 1, mtr$BRAUN$stage))))

cat ("C")

# About Chan's data - immune annotation and myeloid annotation have some detailed info but target cell group is limited.
#    chan_meta <- readRDS(...)   # combined immune metadata + myeloid metadata
#    chan_raw_meta <- chan@meta.data                      # metadata for all cells
#    chan_raw_meta$Cell <- rownames(chan_raw_meta)
#    chan_full_meta <- chan_raw_meta  %>% full_join(chan_meta)
#
#mtr$CHAN  <- mt$CHAN  %>% rename(cellid = Cell, nFeature=n_genes_by_counts.x, nCount = total_counts.x, percent_mito = mito_frac.x, batch = batch.x, patient = patient.x, cancer=histo.x, location=tissue.x, 
#                                 procedure=procedure.x, treatment=treatment.x, RBP_frac=RBP_frac.x) %>% 
#  select(-batch.y, -patient.y, -tissue.y, -treatment.y, -procedure.y, -histo.y, -mito_frac.y, -RBP_frac.y, -log1p_n_genes_by_counts.x, -log1p_total_counts.x, -log1p_n_genes_by_counts.y, -log1p_total_counts.y, -total_counts.y)
#mtr$CHAN$celltype <- paste0(mtr$CHAN$cell_type.x, "_", mtr$CHAN$cell_type.y)
#mtr$CHAN <- mtr$CHAN %>% mutate(celltype=gsub(celltype, pattern="_NA", replacement="", ignore.case=F))
#mtr$CHAN <- mtr$CHAN %>% mutate(celltype=gsub(celltype, pattern="Myeloid_Mo/M.", replacement="Myeloid_Mo/Mp", ignore.case=F)) # CHAN's mo/m<ph> couldn't be recognized
#mtr$CHAN <- mtr$CHAN %>% select(-cell_type.x, -cell_type.y, -libsize)
#mtr$CHAN$tissue <- ifelse(mtr$CHAN$cancer == "normal", "Normal", ifelse(mtr$CHAN$procedure == "Thoracentesis", "Effusion", "Tumor"))

mtr$CHAN <- mt$CHAN %>% rename(nFeature=ngenes, percent_mito = mito_frac, cancer=histo, location=tissue, cluster_id = clusters)
mtr$CHAN$tissue <- ifelse(mtr$CHAN$cancer == "normal", "Normal", ifelse(mtr$CHAN$procedure == "Thoracentesis", "Effusion", "Tumor"))
mtr$CHAN$cluster_id <- as.integer(mtr$CHAN$cluster_id)
mtr$CHAN <- mtr$CHAN %>% rename(major = cell_type_coarse, minor = cell_type_fine)
mtr$CHAN <- mtr$CHAN %>% select(-cell_type_general, -H_knn)
mtr$CHAN <- mtr$CHAN %>% select(-cell_type_med, -cell_type.x, -cell_type.y)  # included in celltype
mtr$CHAN$cancer <- ifelse(mtr$CHAN$cancer == "normal", "LUNG", mtr$CHAN$cancer)

mtr$CHENG <- mt$CHENG %>% rename(celltype = MajorCluster)

cat ("K")

mtr$KRISHNA <- mt$KRISHNA %>% rename(cellid = cell, celltype = cluster_name, cluster_id = cluster, tissue=type, response = Sample_name, patient = Sample, sample = Sample2, location=region)
mtr$KRISHNA$cancer = "ccRCC"

kim2tissue <- tibble(Sample_Origin = c("tLung", "tL/B", "PE", "mLN", "mBrain", "nLung", "nLN"), 
                     tissue = c("Tumor", "Tumor", "Effusion", "LymphNode", "Tumor", "Normal", "Normal"))
mtr$KIM   <- mt$KIM   %>% rename(cellid = Index) %>% select(-Barcode)  # celltype = Cell_subtype
mtr$KIM$celltype <- ifelse(is.na(mtr$KIM$Cell_subtype), mtr$KIM$Cell_type, mtr$KIM$Cell_subtype)
mtr$KIM$cancer <- "LUAD"
mtr$KIM <- mtr$KIM %>% left_join(kim2tissue, by="Sample_Origin")
#  (NG) mtr$KIM$tissue <- mtr$KIM$Sample %>% gsub(pattern="_.*", replace="")
mtr$KIM <- mtr$KIM %>% select(-Cell_subtype, -Cell_type.refined) %>% rename(major = Cell_type, sample = Sample, location = Sample_Origin)
kim_s2p <- read.delim("kim_sample2patient.txt")
mtr$KIM <- mtr$KIM %>% left_join(kim_s2p)

mtr$YOST_BASAL <- mt$YOST_BASAL %>% rename(cellid = "cell.id", celltype = cluster)
mtr$YOST_BASAL$cancer <- "SKIN_BCC"
mtr$YOST_BASAL$tissue <- "Tumor" 
mtr$YOST_SCC   <- mt$YOST_SCC   %>% rename(cellid = "cell.id", celltype = cluster)
mtr$YOST_SCC$cancer <- "SKIN_SCC"
mtr$YOST_SCC$tissue <- "Tumor"

cat ("L")

mtr$LI         <- mt$LI         %>% rename(sort = "Cell.type", batch = Seq.Batch.ID, tissue = Batch.Set.ID, major=all_mc_group, patient=PatientID, sample=Sample.Short.Name) %>% 
    select(-Date.sorted, -Date.processed, -SorterID, -SortedBy, -LibPrepBy, -NKI.Plate.ID, -FACS.Sorting.order..per.session., -Amp.Batch.ID)
mtr$LI$cancer <- "MEL"
mtr$LI$celltype <- ifelse(is.na(mtr$LI$t_nk_mc_group), mtr$LI$non_t_nk_mc_group, paste0("T/NK_", mtr$LI$t_nk_mc_group))
mtr$LI <- mtr$LI %>% select(-CD3, -CD45, -live, -Staining.panel, -TCRseq.available, -all_mc, -t_nk_mc, -t_nk_mc_group, -non_t_nk_mc, -non_t_nk_mc_group, -new.PID, -type, -Median.Net.gene.UMIs, -spike_count, -n_tumor)
mtr$LI$stage <- as.character(mtr$LI$stage)

mtr$BI <- mt$BI %>% filter(NAME != 'TYPE') %>% rename(cellid = NAME, celltype = FinalCellType, cancer = disease__ontology_label, patient=donor_id, sample=biosample_id) %>% 
    select(-library_preparation_protocol, -library_preparation_protocol__ontology_label, -organ, -species, -species__ontology_label, -disease, -ICB_Exposed, -Lineage, -InferCNV) %>%
    rename(response = ICB_Response, location=organ__ontology_label)
mtr$BI$tissue <- ifelse(mtr$BI$location=="lymph node", "LymphNode", "Tumor")

cat ("P")

mtr$Pelka <- mt$Pelka %>% rename(cellid = cellID, patient=PID, tissue=SPECIMEN_TYPE, age=Age, sex=Sex) %>% left_join(mt$Pelka_clust, by="cellid") %>% 
    select(-SizeQuantile, -NodeStatusSimple, -Ethnicity, -SINGLECELL_TYPE,
             -HistologicGradeSimple, -TissueSiteSimple, -SOURCE_HOSPITAL, -TISSUE_PROCESSING_TEAM, -MMRStatus, -MLH1Status, -MMR_IHC, -cl295v11SubShort) %>%
    rename(batch = batchID, major=clTopLevel, minor=clMidwayPr, stage=TumorStage, race=Race, sort=PROCESSING_TYPE, sample=PatientTypeID, location=TissueSite_detailed)
mtr$Pelka$cancer <- "CRC"

mtr$JERBY <- mt$JERBY %>% rename(cellid = cells, celltype = cell.types, sample=samples, nFeature=no.of.genes, nCount=no.of.reads, treated=treatment.group, cohort=Cohort)
mtr$JERBY$cancer <- "MEL"
mtr$JERBY$tissue <- "Tumor"
jerby_s <- read.delim("jerby_sample_info.txt")
mtr$JERBY <- mtr$JERBY %>% left_join(jerby_s, by="sample")

cat ("Q")

mtr$QIAN  <- mt$QIAN  %>% rename(cellid = Cell,  celltype = CellType, patient=PatientNumber, cancer=TumorType, nCount=nUMI, nFeature=nGene, location=TumorSite)
mtr$QIAN$patient <- as.character(mtr$QIAN$patient)
mtr$QIAN$tissue  <- ifelse(mtr$QIAN$CellFromTumor==1, "Tumor", "Normal")   # Normal is 'matched normal'
mtr$QIAN <- mtr$QIAN %>% select(-CellFromTumor)
mtr$QIAN$cancer[mtr$QIAN$patient == 1 | mtr$QIAN$patient == 2 | mtr$QIAN$patient == 7] <- "LUSC"
mtr$QIAN$cancer[mtr$QIAN$patient == 3 | mtr$QIAN$patient == 4 | mtr$QIAN$patient == 6] <- "LUAD"
mtr$QIAN$cancer[mtr$QIAN$patient == 5 | mtr$QIAN$patient == 8 ] <- "NSCLC"  # NSCLC other: 5:LCC, 8: Pleiomorphic

cat ("V")

mtr$VISH  <- mt$VISH  %>% rename(cellid = cells_ID, celltype = Minor, nCount = nCount_RNA, nFeature=nFeature_RNA, major=Major, tissue=type, cluster_id=Final_clusters) %>% 
  select(-integrated_snn_res.0.8, -seurat_clusters, -barcode, -orig.ident, -nCount_SCT, -nFeature_SCT)
mtr$VISH$cancer <- "ccRCC"
mtr$VISH$orig_cellid <- mtr$VISH$cellid
mtr$VISH$cellid <- paste0("VISHWAKARMA_",mtr$VISH$cellid,"-1")
mtr$VISH$patient <- mtr$VISH$sample


cat ("X")

mtr$ZILI  <- mt$ZILI  %>% rename(celltype = "Most.likely.LM22.cell.type", nCount = Total.counts, percent_mito = Percent.counts.from.mitochondrial.genes, patient=Patient, tissue=Tissue, major=Major.cell.type, minor=Minor.subset, batch=Library) %>% 
   select(-used_in_NSCLC_all_cells, -x_NSCLC_all_cells, -y_NSCLC_all_cells, -used_in_NSCLC_and_blood_immune, -x_NSCLC_and_blood_immune, -y_NSCLC_and_blood_immune, 
     -used_in_NSCLC_immune, -x_NSCLC_immune, -y_NSCLC_immune, -used_in_NSCLC_non_immune, -x_NSCLC_non_immune, -y_NSCLC_non_immune, -used_in_blood, -x_blood, -y_blood, -Barcoding.emulsion, -Barcode)
zili_cancer <- data.frame(patient=c("p1","p2","p3","p4","p5","p6","p7"), cancer=c("LUSC", "LUSC", "LUAD", "LUAD", "LUAD", "LUAD", "LUAD"))  # from paper 
mtr$ZILI <- mtr$ZILI %>% left_join(zili_cancer, by=c("patient"))

cat ("Y")

# mtr$LAMB    <- mt$LAMB
mtr$LEADER  <- mt$LEADER
mtr$LEADER <- mtr$LEADER %>% select (-cellgroup, -cellgroup2)
mtr$WU      <- mt$WU

leader_annot <- read.delim("leader_sample_info.txt")
mtr$LEADER$orig.ident <- mtr$LEADER$orig.ident %>% as.character    # as.integer don't work
leader_annot$orig.ident <- leader_annot$orig.ident %>% as.character
mtr$LEADER <- mtr$LEADER %>% left_join(leader_annot, by=c("orig.ident")) %>% rename(cluster_id=clust)
mtr$LEADER$patient <- mtr$LEADER$patient %>% as.character
mtr$LEADER$cluster_id <- mtr$LEADER$cluster_id %>% as.integer

# mtr$LEADER$cancer <- "NSCLC"

mtr$WU$tissue = "Tumor" # only tumor sample is described in the paper.

mtr$MAYNARD <- mt$MAYNARD
mtr$MAYNARD$tissue <- "Tumor"
mtr$MAYNARD$tissue[grep("_.AT", mtr$MAYNARD$patient_id)] <- "Normal"
mtr$MAYNARD <- mtr$MAYNARD %>% rename(patient = patient_id, cluster_id = main_seurat_cluster, response = best_rxn_status) 
mtr$MAYNARD <- mtr$MAYNARD %>% select(-sample_type, -DF.classifications_0.25_0.09_218, -percent.ercc, -seurat_clusters)
mtr$MAYNARD$cluster_id <- mtr$MAYNARD$cluster_id  %>% as.integer()
mtr$MAYNARD <- mtr$MAYNARD %>% rename(sample=sample_name, primary_met=primary_or_metastaic, location=biopsy_site)

cat ("Z")

cellids <- readRDS("cellids_in_rawdata.rds")
qian.cellids <- data.frame(cellid = cellids$cellid[cellids$study == "QIAN"])
qian.cellids$orig_cellid <- sub("QIAN_.+?_", "", qian.cellids$cellid)
mtr$QIAN <- mtr$QIAN %>% rename(orig_cellid = cellid) %>% left_join (qian.cellids, by="orig_cellid")

cat (" formatted.\ncombining ... ")

# bind all
rdss <- names(mtr)
res <- data.frame()
for (i in rdss) {
  cat (paste0("+",i))
  res <- res %>% bind_rows(mtr[[i]])
}

cellChecker <- function(res, cellids) {
  names <- unique(cellids$study)  # skipping new ones: they are safe
  for (i in names) {
    rawids  <- cellids$cellid[cellids$study == i]
    metaids <- res$cellid[res$study == i]
    inter <- intersect(rawids, metaids)  # setdiff(tgt,remove)
    cat (i, "::: rawdata:", length(rawids), "metadata:", length(metaids), "intersection:", length(inter), "\n")
  }
  e <- sum(is.na(res$orig_cellid)) + sum(res$cellid == "")
  cat ("Entries with blank cellid: ", e, ifelse(e==0, "OK", "NG"), " \n")
  s1 <- dim(res)[1]
  s2 <- length(unique(res$cellid))
  e <- s1-s2
  cat ("No. of entries: ", s1, "Unique cellids over studies: ", s2, "Duplicates: ", e, ifelse(e==0, "Good", "caution"), " \n")
}

cat ("fixing cellids ... ")
res$orig_cellid <- ifelse(is.na(res$orig_cellid), res$cellid, res$orig_cellid)
res$cellid      <- ifelse(is.na(res$cellid),      res$orig_cellid, res$cellid)


cat ("\nchecking intersection of cellid between raw data and metadata ...\n")
cellChecker(res, cellids)

cat ("\nA")
res <- res %>% left_join(ct2cg, by="celltype")
res <- res %>% select(-isMyeloid, -vdj_kit)   # isMyeloid used partially in pipeline; vdj_kit: no data and useless)
res <- res %>% select(study, cellid, celltype, cellgroup, cellgroup2, major, minor, patient, sample, cancer, tissue, nCount, nFeature, percent_mito, treatment, batch, response, everything())

if (dim(res)[1] != length(unique(res$cellid))) {
  warning("Duplicate cell ID exists over studies.")
  if (dim(res)[1] != dim(res %>% group_by(study, cellid) %>% summarize())[1]) {
    stop("Duplicate cell ID are detected in same study. Please check the contents of 'res'.")
  }
}

cat ("B")
res$tissue <- ifelse(res$tissue == "LYMPHNODE", "LymphNode", res$tissue)
res$tissue <- ifelse(res$tissue == "LN", "LymphNode", res$tissue)
cat ("C")
res$tissue <- ifelse(res$tissue == "NORMAL", "Normal", res$tissue)
res$tissue <- ifelse(res$tissue == "N", "Normal", res$tissue)
cat ("D")
res$tissue <- ifelse(res$tissue == "TUMOR", "Tumor", res$tissue)
res$tissue <- ifelse(res$tissue == "tumor", "Tumor", res$tissue)
cat ("E")
res$tissue <- ifelse(res$tissue == "T", "Tumor", res$tissue)
res$tissue <- ifelse(res$tissue == "PBMC", "Blood", res$tissue)
cat ("F")
res$tissue <- ifelse(res$tissue == "P", "Blood", res$tissue)
res$tissue <- ifelse(res$tissue == "BLOOD", "Blood", res$tissue)
cat ("G")
res$tissue <- ifelse(res$tissue == "blood", "Blood", res$tissue)
res$tissue <- ifelse(res$tissue == "EFFUSION", "Effusion", res$tissue)
cat ("H")
res$tissue <- ifelse(res$tissue == "pleural_effusion", "Effusion", res$tissue)
res$tissue <- ifelse(res$tissue == "LUNG", "lung", res$tissue)
cat ("I")
res$tissue <- ifelse(res$tissue == "K", "Normal", res$tissue)
res$cancer <- ifelse(res$cancer == "clear cell renal carcinoma", "ccRCC", res$cancer)
cat ("J")
res$cancer <- ifelse(res$cancer == "papillary renal cell carcinoma", "pRCC", res$cancer)
res$cancer <- ifelse(res$cancer == "KIDNEY", "RCC", res$cancer)
cat ("K")
res$cancer <- ifelse(res$cancer == "Ovarian", "OV", res$cancer)
res$cancer <- ifelse(res$cancer == "OV-FTC", "OV", res$cancer)

res$nCount   <- pmax(res$nCount, res$nCount_RNA, na.rm=T)
res$nFeature <- pmax(res$nFeature, res$nFeature_RNA, na.rm=T)
res <- res %>% select(-nCount_RNA, -nFeature_RNA)

cat (" done.\nwriting table ... ")
# ver.2   2017506 x 115

## Write everything to text file
write.table(res, file="full_metadata.txt", sep="\t", quote=F, row.names=F, col.names=T, na="")
cat ("zipping ... ")
system("gzip -f full_metadata.txt")  
cat("\nwriting RDSs ...")
saveRDS(res, file="full_metadata.rds")

cat("and smaller one ...")
# removing clonotype and study specific clusters for myeloid study
digest <- res %>% select(-CTgene, -CTnt, -CTaa, -CTstrict, -Frequency, -cloneType, -UMAP1, -UMAP2, -clusters.y, -clusters.x, -clusters_fine, -Notes)
saveRDS(digest, file="part_metadata.rds")

cat("building the myeloid subset ...")
# get myeloid subset 
myeloid <- digest %>% filter(cellgroup2 == "DC" | cellgroup2 == "Macro" | cellgroup2 == "Mono" | cellgroup2 == "Myeloid")
saveRDS(myeloid, "myeloid_metadata.rds")          # 297140 -> 437467 -> 426473/17 -> 437095

cat("limiting to tumor ...")
myeloidt <- myeloid %>% filter(tissue != "Normal" & tissue != "Blood")
saveRDS(myeloidt, "myeloid_tumor_metadata.rds")   # 212209

cat ("\nall done!\n")

# --------- all finished ---------

# How to write cell type list:
# (first time and updated time only) 
#  ctype <- data.frame(celltype=unique(res$celltype))
#  write.table(ctype, file="celltypes.txt", row.names=F, quote=F)
#  edit celltype_cellgroup.txt from celltypes.txt

