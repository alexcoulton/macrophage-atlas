### Single Cell Integration / Analysis Helper Tools
### Modified 31 Jan 2022 by Jun Murai
### 20 Jan 2022 by Jun Murai

library(Seurat)
# library(SeuratDisk)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)

#' Check the AnchorSet 
#' 
#' Check the AnchorSet (without reference) whether it could be run under limit of elements for matrix and sparse.
#' 
#' @param anchorset anchorset to inspect
#' @return list of stats
#' @examples
#' preCheckAnchorSet(anchorset)
preCheckAnchorSet <- function(anchors, verbose=T) {
    anch    <- slot(object = anchors, name = "anchors")
    offsets <- slot(object = anchors, name = "offsets")
    object.list <- slot(object = anchors, name = "object.list")
    reference.objects <- slot(object = anchors, name = "reference.objects")

    objects.ncell <- sapply(X = object.list, FUN = ncol)
    res <- list()
    res$objects.ncell <- objects.ncell 
    res$num.cells <- sum(objects.ncell) 
    res$num.study <- length(objects.ncell) 
    res$num.anchors <- dim(anchors@anchors)[1]
    if (verbose) { cat(paste0(res$num.study, " studies with ", res$num.cells, " cells, ", res$num.anchors, " anchors in total.\n")); }

    similarity.matrix <- Seurat:::CountAnchors(anchor.df = anch,
            offsets = offsets, obj.lengths = objects.ncell)
    similarity.matrix <- similarity.matrix[reference.objects,
            reference.objects]
    res$similarity.matrix <- similarity.matrix
    sample.tree <- Seurat:::BuildSampleTree(similarity.matrix = similarity.matrix)
    sample.tree <- Seurat:::AdjustSampleTree(x = sample.tree, reference.objects = reference.objects) ##
    res$sample.tree <- as.data.frame(sample.tree)
    names(res$sample.tree) <- c("dataset1", "dataset2")
    order <- Seurat:::SampleIntegrationOrder(tree = sample.tree)
    res$order <- order
    max.num <- 0
    res$sample.tree$anch.num <- NA
    res$sample.tree$msg <- NA

    names(x = object.list) <- as.character(-(1:length(x = object.list)))
    for (ii in 1:nrow(x = sample.tree)) {
        datasets   <- Seurat:::ParseMergePair(sample.tree, ii)           
        filtered.anchors <- anch[anch$dataset1 %in% datasets$object1 &
            anch$dataset2 %in% datasets$object2, ]
        anch.num <- dim(filtered.anchors)[1]
        msg <- paste0(paste(datasets$object2, collapse = " "), " into ", paste(datasets$object1,collapse = " "))
        if (max.num < anch.num) max.num <- anch.num
        res$sample.tree$anch.num[ii] <- anch.num
        res$sample.tree$msg[ii]      <- msg
    }
    res$max.num   <- max.num
    res$features  <- length(anchors@anchor.features)
    res$max.index <- as.numeric(max.num) * res$features
    res$max.gb    <- ceiling(res$max.index * 8 * 100 / 1024 / 1024 / 1024) / 100
    res$max.gb3   <- ceiling(res$max.index * 8 * 100 / 1024 / 1024 / 1024 * 3) / 100
    res$my.gb <- ceiling(object.size(anchors) %>% as.numeric() * 100 / 1024 / 1024 / 1024) / 100
    vec.limit <- 2**31-1
    if (res$max.index > vec.limit) {
      res$features.ideal <- floor(vec.limit / max.num)
      res$msg      <- paste0("Current Seurat will fail unless all routines are long vector safe. ", max.num, " anchors x ", res$features, " are out of bound for old matrix. You may be need to rebuild with about ", res$features.ideal, " or less features. Or, trying with latest R and Seurat could help you.")
    } else {
      res$features.ideal <- res$features
      res$msg      <- paste0(max.num, " anchors in one step are good with ", res$features, " features. Your dataset is ", res$my.gb, " GB and the three Matrices will use ", res$max.gb3, " GB.")
    }
    if (res$max.index <= 100000 * 3000) {
      res$gb.recommend <- "256GB or less memory"
    } else if (res$max.index <= vec.limit) {
      res$gb.recommend <- "512GB memory"
    } else {
      res$gb.recommend <- "N/A"
    }
    if (verbose) { 
      cat(res$msg, "\n")
      cat("Memory suggestion (alpha version) :", res$gb.recommend, "\n");
    }
    return(res)
}

#' compareClusters
#' 
#' Compare clusters of two Seurat objects, typically old and new version, good to trace which cluster is equivalent to each of old clusters.
#' 
#' @param old old version of cluster results from Seurat object, will be shown in rows of the result table
#' @param new another (new) version of cluster results from Seurat object, will be show in columns of the result table
#' @return Table shows relation between two Seurat object
#' @examples
#' tb <- compareClusters(Idents(old), Idents(new))
compareClusters <- function(old, new) {
  # o <- Idents(new)
  o <- data.frame(cellid = names(old), cl_o = old)
  n <- data.frame(cellid = names(new), cl_n = new)
  o <- o %>% full_join(n, by="cellid")
  tb <- o %>% group_by(cl_o, cl_n) %>% summarize(n())
  names(tb) = c("cl_o", "cl_n", "num")
  res <- pivot_wider(tb, id_cols="cl_o", names_from="cl_n", values_from="num")
  return(res)
}

# alpha version, only try Denom = min(vector[vector>0])

#' deNormalizeSeurat
#' 
#' Recover raw count data of Seurat object, in which only normalized 'data' is stored, and raw 'counts' doesn't exist or 'counts' is not correct.
#' Now it try only for raw = vec/min(vec[vec>0]), and will fail if 1 or more cell exists with min count is 2+.
#' In such case we should try raw=vec/(min/2, min/3, ...), but not yet implemented.
#' 
#' @param obj Seurat object
#' @param base2 Is it log2 normalized? TRUE: log2, FALSE: ln (FALSE for current Seurat process)
#' @return Seurat object with raw counts
#' @examples
#' deNormalizeSeurat(object)
deNormalizeSeurat <- function(obj, base2 = T) {
  # modified as.matrix used to avoid 'too large' error
  data   <- utils_big_as.matrix(obj@assays$RNA@data)
  
  # first, convert to linear scaled data. not sparse here
  if (base2) {
    scaled <- 2 ^ data -1
  } else {
    scaled <- exp(data)-1
  }
  cols <- dim(scaled)[2]
  
  for (i in 1:cols) {
    vec  <- scaled[,i]
    vpos <- vec[vec>0]   # min(vpos) may be scale factor
    vmin <- min(vpos)
    vec  <- round(vec / vmin, digits=1)   # 2 may catch errors / 6227.01
    check <- sort(unique(vec))
    for (j in 1:min(10,length(check))) {
      if (check[j] != floor(check[j])) { stop(paste("invalid numbers", check)) }  # if error, should be divided with vmin/2, vmin/3 etc.
    }
    scaled[,i] <- vec
    cat (i, "\r")
  }
  cat ("Raw count recovered.\n")
  # scaled <- round(scaled, digits=1)   # fix it after
  res <- as.sparse(scaled)
  obj@assays$RNA@counts <- res
  return(obj)
}


#' Add metadata to Seurat Object (generic version)
#' 
#' Add metadata to Seurat Object. ID should be included in df as df$id or rownames(df) and should be same as ids in the object.
#' 
#' @param obj Seurat object
#' @param df data frame including metadata
#' @return Seurat object
#' @examples
#' exp <- addMetadataToSeurat(exp, meta)
addMetadataToSeurat <- function (obj, df) {
  current <- obj@meta.data
  if (length(df$id) == 0) { 
    if (identical(rownames(wu.annot)[1:20], as.character(1:20))) {
      stop ("couldn't find row names or id for the data frame")
    } else {
      df$addAnnot_joinID = rownames(df)
      cat ("row names of df are used as the IDs.\n")
    }
  } else {
    cat ("df$id is used as the ID column.\n")
    df <- df %>% dplyr::rename(addAnnot_joinID = id)
  }
  if (length(df$addAnnot_joinID) != length(unique(df$addAnnot_joinID))) {
    stop ("The data frame has duplicate IDs (df$id).")
  }
  current$addAnnot_joinID = rownames(current)
  current <- current %>% left_join(df)
  rownames(current) <- current$addAnnot_joinID
  current <- current %>% dplyr::select(-addAnnot_joinID)
  if (!identical(rownames(current), rownames(obj@meta.data))) {
    stop("internal error: couldn't keep ids of meta data unmodified");
  }
  obj@meta.data <- current
  return(obj)
}

#' Add metadata to Seurat Object (tentative)
#' 
#' Add metadata to Seurat Object, tentative version to use old metadata (some cellid is not formatted well)
#' Designed for add metadata to seurat object of 14 studies 
#' 
#' @param obj Seurat object
#' @param meta metadata data frame (ver 0.1)
#' @return Seurat object
#' @examples
#' exp <- tmp_addMetaDataToExp(exp, meta)
tmp_addMetaDataToExp <- function(obj, meta) {
  meta$cellid[meta$study == "AZIZI"] <- paste0("c",meta$cellid[meta$study == "AZIZI"])
  meta$cellid[meta$study == "VISHWAKARMA"] <- paste0("VISHWAKARMA_",meta$cellid[meta$study == "VISHWAKARMA"],"-1")

  qiancells <- data.frame(cellid2 = obj@meta.data %>% filter(study == "QIAN") %>% rownames())
  qiancells$cellid <- sub("QIAN_[A-Z]+_", "", qiancells$cellid2)
  meta <- meta %>% left_join (qiancells)
  meta[!is.na(meta$cellid2),]$cellid <- meta[!is.na(meta$cellid2),]$cellid2
  meta <- meta %>% dplyr::select(-cellid2)
  if (length(unique(meta$cellid)) != dim(meta)[1]) { stop ("failed formatting meta.") }
 
  # update Seurat object with new metadata  # %>% select(-cancer, -tissue)
  md <- obj@meta.data 
  md$cellid <- rownames(md)
  md <- md %>% left_join(meta, by="cellid")
  md <- md %>% rename(study = study.x)
  identical(rownames(obj@meta.data), md$cellid)
  rownames(md) <- md$cellid
  # prevent breaking object
  if (!identical(rownames(obj@meta.data), rownames(md))) { stop ("failed building updated meta table.") }
  obj@meta.data <- md
  return(obj)
}

#' reportMetaDataNAs
#' 
#' Report numbers of NA's in metadata
#' 
#' @param exp Seurat Object to inspect
#' @param thr threshold of ratio of NA's to suggest removing from the table
#' @return nothing
reportMetaDataNAs <- function(exp, thr = 0.8) {
  print(dim(exp@meta.data))
  thres <- dim(exp@meta.data)[1] * thr
  suggest <- c()
  cat("Numbers of NA's in the meta data table:\n")
  for (i in names(exp@meta.data)) {
    num <- sum(is.na(exp@meta.data[[i]]))
    cat (i, ":", num, ";  ")
    if (num > thres) {
      suggest <- c(suggest, paste0("-",i))
    }
  }
  if (length(suggest) > 0) {
    suggest <- paste(suggest, collapse=", ")
    suggest <- paste0("exp@meta.data <- exp@meta.data %>% select(", suggest, ")")
  } else {
    suggest <- "none"
  }
  cat ("\nSuggestion: \n", suggest, "\n")
 # return(suggest)
}

#' getSSEfromExp
#' 
#' get Sum of Square Errors from a Seurat object
#' 
#' @param exp Seurat object
#' @param maxPC PC numbers to be tested
#' @return SSE
getSSEfromExp <- function (exp, maxPC = 50) {
  sumres <- 0
  for (i in levels(Idents(exp))) {       # Clusters
    for (j in 1:maxPC) {                 # PC numbers
      m <- mean (exp@reductions$pca@cell.embeddings[Idents(exp) == i,j])  # get mean for cluster i, PC j
      res <- exp@reductions$pca@cell.embeddings[Idents(exp) == i,j] - m   # vector of residuals: vector of [ pos - center ]
      sumres <- sumres + sum(res ** 2)   # add sum of SSE for this dimension    + sum of  (pos - center)^2
    }
  }
  sumres
}


#' saveClusterGroupPlots
#' 
#' Save a cluster/group to png file and return the plot
#' 
#' @param grouped_table meta.data table grouped by cluster and another group
#' @param label a word to add as a label on the figure and the name of the file
#' @param num.clust Flag to treat cluster names as numbers. If true, cluster ID is converted to numeric and sorted as a numeric
#' @param do.save If true(default), save the plot to a file named 'chart_LABEL.png'
#' @param fill if true(default), graphs are filled up to 100% and shows ratio, otherwise stacked graph will be drawn
#' @return list of ggplot and raw table
#' @examples
#' p <- saveClusterGroupPlots(exp@meta.data %>% group_by(integrated_snn_res.0.5, cellgroup), "cellgroup", num.clust=T)
saveClusterGroupPlots <- function(grouped_table, label="", num.clust=FALSE, do.save=TRUE, width=14, fill=TRUE) {
  fill = ifelse(fill, "fill", "stack") 
  tb <- summarize(grouped_table, n())
  if(label == "") label <- names(tb)[2]
  names(tb) <- c("cluster", "group", "cells")
  if (num.clust) tb$cluster <- as.numeric(as.character(tb$cluster))
  p <- ggplot(tb, aes(fill=group, y=cells, x=cluster)) + geom_bar(position=fill, stat="identity") + labs(fill=label)
  if (num.clust) {
    min <- min(tb$cluster)
    max <- max(tb$cluster)
    p <- p + scale_x_continuous(breaks=min:max)
  }
  if (do.save) ggsave(file=paste0("chart_",label,".png"), plot = p, width=width)
  return (list(table = tb, plot = p))
}


# swapping of ID is not allowed. For example INSX to INS and INS to INS1. original ID is always reserved for the original one.
# when doubled
#    1. original (reserved/removed from the list) -> 2. top of the list

rmBlackList <- function(ids) {
  ids <- ids[ids$ENSEMBL != "ENSG00000237541" | ids$SYMBOL != "HLA-DQA1",]  # removing ENSEMBL==ENSG00000237541 & SYMBOL==HLA-DQA1
  ids
}


doConvertGeneName <- function(obj, id2newid, shrink=F, extra.data=NA, verbose=T) {
  DefaultAssay(obj) <- "RNA"               # use assays=RNA 

  # remove extra information, only keep RNA@counts and RNA@data
  if (shrink) {
    obj <- DietSeurat(obj, assays="RNA")  
    if (verbose) { cat ("doConvertGeneName: Extra data slots were removed, except RNA@counts and RNA@data.\n") } 
  }

  # error check
  if (length(Assays(obj)) != 1) {
    msg <- paste("Other assay[s] were found:", names(obj@assays))
    if (extra.data=="ignore") { cat (msg, "\n") } else { stop (msg) }
  }
  if (dim(obj@assays[[obj@active.assay]]@scale.data)[1] != 0) stop ("Scale data exists"); 
  if (!identical(rownames(obj@assays$RNA@data), rownames(obj@assays$RNA@counts))) stop ("Row names of counts and data are not identical."); 

  id2newid<-na.omit(id2newid)                      # NA is not allowed
  id2newid<-id2newid[!duplicated(id2newid$id),]   # same id is not allowed, 2nd+ occurences will be removed

  join.df  <- data.frame(newid=rownames(obj), x=rownames(obj)) # test join 
  id2newid <- id2newid %>% left_join (join.df)
  id2newid <- id2newid[is.na(id2newid$x),]                     # remove ids with x
  id2newid <- id2newid %>% dplyr::select (-x)
  
  join.df  <- data.frame(id=rownames(obj), x=rownames(obj)) 
  id2newid <- id2newid %>% left_join (join.df)
  id2newid <- id2newid[!is.na(id2newid$x),]                     # remove ids with x
  id2newid <- id2newid %>% dplyr::select (-x)
  
  join.df  <- data.frame(id=rownames(obj), newid=rownames(obj)) 
  id2newid <- rbind(id2newid, join.df) 
  id2newid<-id2newid[!duplicated(id2newid$newid),]    # newid duplication for entry genes will be removed
  id2newid<-id2newid[!duplicated(id2newid$id),]       # finally, id=id removed with id=newid

  ids.apply <- data.frame(id=rownames(obj)) %>% left_join(id2newid)     # keep order of rownames(obj)
  rownames(obj@assays$RNA@counts) <- ids.apply$newid
  rownames(obj@assays$RNA@data)   <- ids.apply$newid
  obj
}

# test doConvertGeneName
convertEnsemblToGeneSymbol_KeepAll2 <- function(exp) {
  library(org.Hs.eg.db)
  cat(" - querying database for Gene Symbols ... \n")
  ids=select(org.Hs.eg.db, keys=rownames(exp), columns=c('ENSEMBL','SYMBOL'),keytype='ENSEMBL') %>% rmBlackList()
  names(ids) <- c("id", "newid")
  exp <- doConvertGeneName(exp, ids, shrink=T)
  return(exp)
}

convertEnsemblToGeneSymbol_KeepAll <- function(exp) {
  library(org.Hs.eg.db)
  DefaultAssay(exp) <- "RNA"
  cat(" - querying database for Gene Symbols ... \n")
  ids=select(org.Hs.eg.db, keys=rownames(exp), columns=c('ENSEMBL','SYMBOL'),keytype='ENSEMBL') %>% rmBlackList()

  cat(" - building conversion list ... \n")
  ids<-na.omit(ids)
  ids<-ids[!duplicated(ids$SYMBOL),]   # remove 2nd+ occurences
  ids<-ids[!duplicated(ids$ENSEMBL),]  # remove 2nd+ occurences  [1] 15393     2
  ids.apply <- data.frame(orig=rownames(exp))               # 16323     1
  ids.apply <- ids.apply %>% left_join(ids, by=c("orig"="ENSEMBL")) # 16323     2

  ids.apply[is.na(ids.apply$SYMBOL),]$SYMBOL = ids.apply[is.na(ids.apply$SYMBOL),]$ENSEMBL  # 16323 2
# length(unique(ids.apply$SYMBOL)) # 16323
  rownames(exp@assays$RNA@counts) <- ids.apply$SYMBOL
  rownames(exp@assays$RNA@data)  <- ids.apply$SYMBOL
  cat (" - writing Gene Symbol to the experiment - done.\n")
  return(exp)
}

convertEnsemblToGeneSymbol_RemoveUnknown <- function(exp) {
  library(org.Hs.eg.db)
  DefaultAssay(exp) <- "RNA"
  cat(" - querying database for Gene Symbols ... \n")
  ids=select(org.Hs.eg.db, keys=rownames(exp), columns=c('ENSEMBL','SYMBOL'),keytype='ENSEMBL')

  cat(" - building conversion list ... \n")
  ids<-na.omit(ids)
  ids<-ids[!duplicated(ids$SYMBOL),]   # remove 2nd+ occurences
  ids<-ids[!duplicated(ids$ENSEMBL),]  # remove 2nd+ occurences  [1] 15393     2

  pos<-match(ids$ENSEMBL, rownames(exp))    # feature numbers to use  15393: 2 3 4 5 6 7 to 16309 16310 16311 16312 16313 16314
  exp <- exp[pos,]
  exp@assays$RNA@counts@Dimnames[[1]] <- ids$SYMBOL
  exp@assays$RNA@data@Dimnames[[1]] <- ids$SYMBOL
  return(exp)
}



# Convert a large sparse matrix into a dense matrix

#' utils_big_as.matrix: Convert a large sparse matrix into a dense matrix
#' 
#' utils_big_as.matrix from CBMR-Single-Cell-Omics-Platform/SCOPfunctions: Single Cell Omics Platform Functions.
#' Avoid the following error . Error in asMethod(object) : . Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, . by slicing the matrix into submatrices, converting and cbinding them. Increases number of slices until they succeed.
#' 
#' @param sparseMat a big sparse matrix of a type coercible to dense Matrix::Matrix
#' @param n_slices_init initial number of slices. Default value 1, i.e. whole matrix
#' @param verbose print progress
#' @return a dense matrix
#' @examples
#' utils_big_as.matrix(sparseMat, n_slices_init = 1, verbose = T)
utils_big_as.matrix <- function(
  sparseMat,
  n_slices_init=1,
  verbose=T
  ) {

  n_slices <- n_slices_init-1
  while (TRUE) {
    list_densemat = list()
    n_slices = n_slices+1
    if (verbose) message(paste0("n_slices=",n_slices))
    idx_to = 0
    for (slice in 1:n_slices) {
      if (verbose) message(paste0("converting slice ",slice,"/",n_slices))
      idx_from <- idx_to+1
      idx_to <- if (slice<n_slices) as.integer(ncol(sparseMat)*slice/n_slices) else ncol(sparseMat)
      if (verbose) message(paste0("columns ", idx_from,":", idx_to))
      densemat_sub = try(
        expr = {
          as.matrix(sparseMat[,idx_from:idx_to])
        }, silent = if (verbose) FALSE else TRUE)
      if ("try-error" %in% class(densemat_sub)) {
        break # exit to while loop
      } else {
        list_densemat[[slice]] = densemat_sub
      }
    }
    if (length(list_densemat)==n_slices) break # exit while loop
  }
  if (verbose) message("cbind dense submatrices")
  densemat <- Reduce(f=cbind, x=list_densemat)
  return(densemat)
}
