# run IntegrateData with watch. 

# usage
# 
# library(Seurat)
# source("integrate_functions_with_debugcode.R")
# anch <- readRDS("integrate_anchor.rds")
# exp <- IntegrateDataW(anch, normalization.method = "SCT")
# saveRDS(exp, "Integrated_with_debugger.rds")

repSize <- function (name, obj) {
  sz <- floor(object.size(obj)/1024/1024)
  cat (":::", name, "took", sz, "MB of memory.\n")
}

IntegrateDataW <- function (anchorset, new.assay.name = "integrated", normalization.method = c("LogNormalize",
    "SCT"), features = NULL, features.to.integrate = NULL, dims = 1:30,
    k.weight = 100, weight.reduction = NULL, sd.weight = 1, sample.tree = NULL,
    preserve.order = FALSE, eps = 0, verbose = TRUE)
{
    print(gc())
    repSize("(IntegrateData) anchorset", anchorset)
    normalization.method <- match.arg(arg = normalization.method)
    repSize("(IntegrateData) normalization.method", normalization.method)
    reference.datasets <- slot(object = anchorset, name = "reference.objects")
    object.list <- slot(object = anchorset, name = "object.list")
    repSize("(IntegrateData) (shared) object.list", object.list)
    anchors <- slot(object = anchorset, name = "anchors")
    ref <- object.list[reference.datasets]
    features <- features %||% slot(object = anchorset, name = "anchor.features")
    repSize("(IntegrateData) features", features)
    unintegrated <- suppressWarnings(expr = merge(x = object.list[[1]],   # temporarily skipped
        y = object.list[2:length(x = object.list)]))
    repSize("(IntegrateData) unintegrated", unintegrated)
    if (!is.null(x = features.to.integrate)) {
        features.to.integrate <- intersect(x = features.to.integrate,
            y = Reduce(f = intersect, x = lapply(X = object.list,
                FUN = rownames)))
    }
    if (normalization.method == "SCT") {
        model.list <- list()
        for (i in 1:length(x = object.list)) {
            assay <- DefaultAssay(object = object.list[[i]])
            if (length(x = setdiff(x = features.to.integrate,
                y = features)) != 0) {
                object.list[[i]] <- GetResidual(object = object.list[[i]],
                  features = setdiff(x = features.to.integrate,
                    y = features), verbose = verbose)
            }
            model.list[[i]] <- slot(object = object.list[[i]][[assay]],
                name = "SCTModel.list")
            object.list[[i]][[assay]] <- suppressWarnings(expr = CreateSCTAssayObject(data = GetAssayData(object = object.list[[i]],
                assay = assay, slot = "scale.data")))
        }
        model.list <- unlist(x = model.list)
        slot(object = anchorset, name = "object.list") <- object.list
    }
    reference.integrated <- PairwiseIntegrateReferenceW(anchorset = anchorset,
        new.assay.name = new.assay.name, normalization.method = normalization.method,
        features = features, features.to.integrate = features.to.integrate,
        dims = dims, k.weight = k.weight, weight.reduction = weight.reduction,
        sd.weight = sd.weight, sample.tree = sample.tree, preserve.order = preserve.order,
        eps = eps, verbose = verbose)
    if (normalization.method == "SCT") {
        if (is.null(x = Tool(object = reference.integrated, slot = "Integration"))) {
            reference.sample <- slot(object = anchorset, name = "reference.objects")
        }
        else {
            reference.sample <- Seurat:::SampleIntegrationOrder(tree = slot(object = reference.integrated,
                name = "tools")$Integration@sample.tree)[1]
        }
        reference.cells <- Cells(x = object.list[[reference.sample]])
        reference.model <- NULL
        if (length(x = model.list) > 0) {
            reference.model <- sapply(X = model.list, FUN = function(model) {
                reference.check <- FALSE
                model.cells <- Cells(x = model)
                if (length(x = model.cells) > 0 & length(x = setdiff(x = model.cells,
                  y = reference.cells)) == 0) {
                  reference.check <- TRUE
                }
                return(reference.check)
            })
            reference.model <- model.list[[which(reference.model)]]
        }
    }
    if (length(x = reference.datasets) == length(x = object.list)) {
        if (normalization.method == "SCT") {
            reference.integrated[[new.assay.name]] <- CreateSCTAssayObject(data = GetAssayData(object = reference.integrated,
                assay = new.assay.name, slot = "data"), scale.data = ScaleData(object = GetAssayData(object = reference.integrated,
                assay = new.assay.name, slot = "scale.data"),
                do.scale = FALSE, do.center = TRUE, verbose = FALSE),
                SCTModel.list = reference.model)
            levels(x = reference.integrated[[new.assay.name]]) <- "refmodel"
            reference.integrated[[assay]] <- unintegrated[[assay]]
        }
        DefaultAssay(object = reference.integrated) <- new.assay.name
        VariableFeatures(object = reference.integrated) <- features
        reference.integrated[["FindIntegrationAnchors"]] <- slot(object = anchorset,
            name = "command")
        reference.integrated <- suppressWarnings(LogSeuratCommand(object = reference.integrated))
        return(reference.integrated)
    }
    else {
        active.assay <- DefaultAssay(object = ref[[1]])
        reference.integrated[[active.assay]] <- NULL
        reference.integrated[[active.assay]] <- CreateAssayObject(data = GetAssayData(object = reference.integrated[[new.assay.name]],
            slot = "data"))
        DefaultAssay(object = reference.integrated) <- active.assay
        reference.integrated[[new.assay.name]] <- NULL
        VariableFeatures(object = reference.integrated) <- features
        integrated.data <- Seurat:::MapQueryData(anchorset = anchorset,
            reference = reference.integrated, new.assay.name = new.assay.name,
            normalization.method = normalization.method, features = features,
            features.to.integrate = features.to.integrate, dims = dims,
            k.weight = k.weight, weight.reduction = weight.reduction,
            sd.weight = sd.weight, preserve.order = preserve.order,
            eps = eps, verbose = verbose)
        integrated.assay <- CreateAssayObject(data = integrated.data)
        if (normalization.method == "SCT") {
            integrated.assay <- CreateSCTAssayObject(data = integrated.data,
                scale.data = ScaleData(object = integrated.data,
                  do.scale = FALSE, do.center = TRUE, verbose = FALSE),
                SCTModel.list = reference.model)
            levels(x = integrated.assay) <- "refmodel"
        }
        unintegrated[[new.assay.name]] <- integrated.assay
        unintegrated <- SetIntegrationData(object = unintegrated,
            integration.name = "Integration", slot = "anchors",
            new.data = anchors)
        if (!is.null(x = Tool(object = reference.integrated,
            slot = "Integration"))) {
            sample.tree <- GetIntegrationData(object = reference.integrated,
                integration.name = "Integration", slot = "sample.tree")
        }
        unintegrated <- SetIntegrationData(object = unintegrated,
            integration.name = "Integration", slot = "sample.tree",
            new.data = sample.tree)
        DefaultAssay(object = unintegrated) <- new.assay.name
        VariableFeatures(object = unintegrated) <- features
        unintegrated[["FindIntegrationAnchors"]] <- slot(object = anchorset,
            name = "command")
        unintegrated <- suppressWarnings(LogSeuratCommand(object = unintegrated))
        return(unintegrated)
    }
}



PairwiseIntegrateReferenceW <- function (anchorset, new.assay.name = "integrated", normalization.method = c("LogNormalize",
    "SCT"), features = NULL, features.to.integrate = NULL, dims = 1:30,
    k.weight = 100, weight.reduction = NULL, sd.weight = 1, sample.tree = NULL,
    preserve.order = FALSE, eps = 0, verbose = TRUE)
{
    print(gc())
    object.list <- slot(object = anchorset, name = "object.list")
    reference.objects <- slot(object = anchorset, name = "reference.objects")
    features <- features %||% slot(object = anchorset, name = "anchor.features")
    features.to.integrate <- features.to.integrate %||% features
    if (length(x = reference.objects) == 1) {
        ref.obj <- object.list[[reference.objects]]
        ref.obj[[new.assay.name]] <- CreateAssayObject(data = GetAssayData(ref.obj,
            slot = "data")[features.to.integrate, ])
        DefaultAssay(object = ref.obj) <- new.assay.name
        return(ref.obj)
    }
    anchors <- slot(object = anchorset, name = "anchors")
    offsets <- slot(object = anchorset, name = "offsets")
    objects.ncell <- sapply(X = object.list, FUN = ncol)
    if (!is.null(x = weight.reduction)) {
        if (length(x = weight.reduction) == 1 | inherits(x = weight.reduction,
            what = "DimReduc")) {
            if (length(x = object.list) == 2) {
                weight.reduction <- list(NULL, weight.reduction)
            }
            else if (inherits(x = weight.reduction, what = "character")) {
                weight.reduction <- as.list(x = rep(x = weight.reduction,
                  times = length(x = object.list)))
            }
            else {
                stop("Invalid input for weight.reduction. Please specify either the names of the dimension",
                  "reduction for each object in the list or provide DimReduc objects.")
            }
        }
        if (length(x = weight.reduction) != length(x = object.list)) {
            stop("Please specify a dimension reduction for each object, or one dimension reduction to be used for all objects")
        }
        if (inherits(x = weight.reduction, what = "character")) {
            weight.reduction <- as.list(x = weight.reduction)
        }
        available.reductions <- lapply(X = object.list, FUN = FilterObjects,
            classes.keep = "DimReduc")
        for (ii in 1:length(x = weight.reduction)) {
            if (ii == 1 & is.null(x = weight.reduction[[ii]]))
                next
            if (!inherits(x = weight.reduction[[ii]], what = "DimReduc")) {
                if (!weight.reduction[[ii]] %in% available.reductions[[ii]]) {
                  stop("Requested dimension reduction (", weight.reduction[[ii]],
                    ") is not present in object ", ii)
                }
                weight.reduction[[ii]] <- object.list[[ii]][[weight.reduction[[ii]]]]
            }
        }
    }
    if (is.null(x = sample.tree)) {
        similarity.matrix <- Seurat:::CountAnchors(anchor.df = anchors,
            offsets = offsets, obj.lengths = objects.ncell)
        similarity.matrix <- similarity.matrix[reference.objects,
            reference.objects]
        sample.tree <- Seurat:::BuildSampleTree(similarity.matrix = similarity.matrix)
        sample.tree <- Seurat:::AdjustSampleTree(x = sample.tree, reference.objects = reference.objects)
    }
    cellnames.list <- list()
    for (ii in 1:length(x = object.list)) {
        cellnames.list[[ii]] <- colnames(x = object.list[[ii]])
    }
    unintegrated <- suppressWarnings(expr = merge(x = object.list[[reference.objects[[1]]]],
        y = object.list[reference.objects[2:length(x = reference.objects)]]))
    repSize("(PairwiseIntegrateReference) unintegrated", unintegrated)
    names(x = object.list) <- as.character(-(1:length(x = object.list)))
    if (!is.null(x = weight.reduction)) {
        names(x = weight.reduction) <- names(x = object.list)
    }
    if (verbose & (length(x = reference.objects) != length(x = object.list))) {
        message("Building integrated reference")
    }
    for (ii in 1:nrow(x = sample.tree)) {

        repSize("(PairwiseIntegrateReference) (current size) object.list", object.list)

        merge.pair <- as.character(x = sample.tree[ii, ])
        length1 <- ncol(x = object.list[[merge.pair[1]]])
        length2 <- ncol(x = object.list[[merge.pair[2]]])
        if (!(preserve.order) & (length2 > length1)) {
            merge.pair <- rev(x = merge.pair)
            sample.tree[ii, ] <- as.numeric(merge.pair)
        }
        if (!is.null(x = weight.reduction)) {
            weight.pair <- weight.reduction[merge.pair]
        }
        else {
            weight.pair <- NULL
        }
        object.1 <- DietSeurat(object = object.list[[merge.pair[1]]],
            assays = DefaultAssay(object = object.list[[merge.pair[1]]]),
            counts = FALSE)
        object.2 <- DietSeurat(object = object.list[[merge.pair[2]]],
            assays = DefaultAssay(object = object.list[[merge.pair[2]]]),
            counts = FALSE)
        suppressWarnings(object.1[["ToIntegrate"]] <- object.1[[DefaultAssay(object = object.1)]])
        DefaultAssay(object = object.1) <- "ToIntegrate"
        object.1 <- DietSeurat(object = object.1, assays = "ToIntegrate")
        suppressWarnings(object.2[["ToIntegrate"]] <- object.2[[DefaultAssay(object = object.2)]])
        DefaultAssay(object = object.2) <- "ToIntegrate"
        object.2 <- DietSeurat(object = object.2, assays = "ToIntegrate")
        datasets <- Seurat:::ParseMergePair(sample.tree, ii)
        if (verbose) {
            message("Merging dataset ", paste(datasets$object2,
                collapse = " "), " into ", paste(datasets$object1,
                collapse = " "))
        }
        merged.obj <- merge(x = object.1, y = object.2, merge.data = TRUE)
        if (verbose) {
            message("Extracting anchors for merged samples")
        }
        filtered.anchors <- anchors[anchors$dataset1 %in% datasets$object1 &
            anchors$dataset2 %in% datasets$object2, ]

        repSize("(PairwiseIntegrateReferenceW) object.1", object.1)
        repSize("(PairwiseIntegrateReferenceW) object.2", object.2)
        repSize("(PairwiseIntegrateReferenceW) merged.obj", merged.obj)

        integrated.matrix <- RunIntegrationW(filtered.anchors = filtered.anchors,
            normalization.method = normalization.method, reference = object.1,
            query = object.2, cellnames.list = cellnames.list,
            new.assay.name = new.assay.name, features.to.integrate = features.to.integrate,
            features = features, dims = dims, weight.reduction = weight.reduction,
            k.weight = k.weight, sd.weight = sd.weight, eps = eps,
            verbose = verbose)

        repSize("(PairwiseIntegrateReferenceW) integrated.matrix", merged.obj)

        integrated.matrix <- cbind(integrated.matrix, GetAssayData(object = object.1,
            slot = "data")[features.to.integrate, ])

        repSize("(PairwiseIntegrateReferenceW) integrated.matrix (added) ", merged.obj)

        merged.obj[[new.assay.name]] <- CreateAssayObject(data = integrated.matrix)
        DefaultAssay(object = merged.obj) <- new.assay.name

        repSize("(PairwiseIntegrateReferenceW) merged.obj to be added to object.list", merged.obj)

        object.list[[as.character(x = ii)]] <- merged.obj
        object.list[[merge.pair[[1]]]] <- NULL
        object.list[[merge.pair[[2]]]] <- NULL
        print(gc())
        invisible(x = CheckGC())
    }
    integrated.data <- GetAssayData(object = object.list[[as.character(x = ii)]],
        assay = new.assay.name, slot = "data")
    integrated.data <- integrated.data[, colnames(x = unintegrated)]
    new.assay <- new(Class = "Assay", counts = new(Class = "dgCMatrix"),
        data = integrated.data, scale.data = matrix(), var.features = vector(),
        meta.features = data.frame(row.names = rownames(x = integrated.data)),
        misc = NULL)
    unintegrated[[new.assay.name]] <- new.assay
    DefaultAssay(object = unintegrated) <- new.assay.name
    VariableFeatures(object = unintegrated) <- features
    if (normalization.method == "SCT") {
        unintegrated[[new.assay.name]] <- SetAssayData(object = unintegrated[[new.assay.name]],
            slot = "scale.data", new.data = as.matrix(x = GetAssayData(object = unintegrated[[new.assay.name]],
                slot = "data")))
    }
    unintegrated <- SetIntegrationData(object = unintegrated,
        integration.name = "Integration", slot = "anchors", new.data = anchors)
    unintegrated <- SetIntegrationData(object = unintegrated,
        integration.name = "Integration", slot = "sample.tree",
        new.data = sample.tree)
    unintegrated[["FindIntegrationAnchors"]] <- slot(object = anchorset,
        name = "command")
    suppressWarnings(expr = unintegrated <- LogSeuratCommand(object = unintegrated))
    return(unintegrated)
}

RunIntegrationW <- function (filtered.anchors, normalization.method, reference,
    query, cellnames.list, new.assay.name, features.to.integrate,
    weight.reduction, weights.matrix = NULL, no.offset = FALSE,
    features, dims, k.weight, sd.weight, eps, verbose)
{
    print(gc())
    cells1 <- colnames(x = reference)
    cells2 <- colnames(x = query)
    if (nrow(x = filtered.anchors) < k.weight) {
        warning("Number of anchors is less than k.weight. Lowering k.weight for sample pair.")
        k.weight <- nrow(x = filtered.anchors)
    }
    merged.obj <- merge(x = reference, y = query, merge.data = TRUE)
    if (no.offset) {
        cell1.offset <- filtered.anchors[, 1]
        cell2.offset <- filtered.anchors[, 2]
    }
    else {
        cell1.offset <- Seurat:::GetCellOffsets(anchors = filtered.anchors,
            dataset = 1, cell = 1, cellnames.list = cellnames.list,
            cellnames = cells1)
        cell2.offset <- Seurat:::GetCellOffsets(anchors = filtered.anchors,
            dataset = 2, cell = 2, cellnames.list = cellnames.list,
            cellnames = cells2)
    }
    filtered.anchors[, 1] <- cell1.offset
    filtered.anchors[, 2] <- cell2.offset
    integration.name <- "integrated"
    merged.obj <- Seurat:::SetIntegrationData(object = merged.obj, integration.name = integration.name,
        slot = "anchors", new.data = filtered.anchors)
    merged.obj <- Seurat:::SetIntegrationData(object = merged.obj, integration.name = integration.name,
        slot = "neighbors", new.data = list(cells1 = cells1,
            cells2 = cells2))
    repSize("(RunIntegration) merged.obj", merged.obj)

    merged.obj <- FindIntegrationMatrixW(object = merged.obj,
        integration.name = integration.name, features.integrate = features.to.integrate,
        verbose = verbose)
    print(gc())
    repSize("(RunIntegration) merged.obj (returned)", merged.obj)
    assay <- DefaultAssay(object = merged.obj)
    if (is.null(x = weights.matrix)) {
        if (is.null(x = weight.reduction) && !is.null(x = dims)) {
            if (normalization.method == "SCT") {
                centered.resids <- ScaleData(object = GetAssayData(object = merged.obj,
                  assay = assay, slot = "data"), do.scale = FALSE,
                  do.center = TRUE, verbose = FALSE)
                merged.obj[["pca"]] <- RunPCA(object = centered.resids[features,
                  ], assay = assay, npcs = max(dims), verbose = FALSE,
                  features = features)
            }
            else {
                merged.obj <- ScaleData(object = merged.obj,
                  features = features, verbose = FALSE)
                merged.obj <- RunPCA(object = merged.obj, npcs = max(dims),
                  verbose = FALSE, features = features)
            }
            dr.weights <- merged.obj[["pca"]]
        }
        else if (is.null(x = weight.reduction) && is.null(x = dims)) {
            dr.weights <- CreateDimReducObject(embeddings = as.matrix(x = t(x = GetAssayData(object = merged.obj))),
                key = "int_", assay = "ToIntegrate")
            dims <- 1:ncol(x = dr.weights)
        }
        else {
            dr <- weight.reduction[[2]]
            if (!all(cells2 %in% rownames(x = dr))) {
                stop("Query cells not present in supplied DimReduc object. Set weight.reduction to a DimReduc object containing the query cells.")
            }
            if (inherits(x = dr, what = "DimReduc")) {
                dr.weights <- dr
            }
            else {
                dr.weights <- query[[dr]]
            }
            dims <- 1:ncol(x = dr.weights)
        }
        merged.obj <- Seurat:::FindWeights(object = merged.obj, integration.name = integration.name,
            reduction = dr.weights, dims = dims, k = k.weight,
            sd.weight = sd.weight, eps = eps, verbose = verbose)
    }
    else {
        merged.obj <- SetIntegrationData(object = merged.obj,
            integration.name = "integrated", slot = "weights",
            new.data = weights.matrix)
    }
    merged.obj <- Seurat:::TransformDataMatrix(object = merged.obj, new.assay.name = new.assay.name,
        features.to.integrate = features.to.integrate, integration.name = integration.name,
        verbose = verbose)
    integrated.matrix <- GetAssayData(object = merged.obj, assay = new.assay.name,
        slot = "data")
    return(integrated.matrix[, cells2])
}

FindIntegrationMatrixW <- function (object, assay = NULL, integration.name = "integrated",
    features.integrate = NULL, verbose = TRUE)
{
    print(gc())
    assay <- assay %||% DefaultAssay(object = object)
    neighbors <- GetIntegrationData(object = object, integration.name = integration.name,
        slot = "neighbors")
    nn.cells1 <- neighbors$cells1
    nn.cells2 <- neighbors$cells2
    cat("length of neighbors: ", length(nn.cells1), ", ", length(nn.cells2), "\n")
    anchors <- GetIntegrationData(object = object, integration.name = integration.name,
        slot = "anchors")
    cat("dim anchors: ", dim(anchors), "\n")
    if (verbose) {
        message("Finding integration vectors")
    }
    features.integrate <- features.integrate %||% rownames(x = GetAssayData(object = object,
        assay = assay, slot = "data"))
    cat ("features.integrate: ", length(features.integrate), "\n")
    data.use1 <- t(x = as.matrix(GetAssayData(object = object, assay = assay, slot = "data")[features.integrate, nn.cells1]))
    repSize("(FindIntegrationMatrix) data.use1", data.use1)
    cat("dim: ", dim(data.use1), "\n")
    data.use2 <- t(x = as.matrix(GetAssayData(object = object, assay = assay, slot = "data")[features.integrate, nn.cells2]))
    repSize("(FindIntegrationMatrix) data.use2", data.use2)
    anchors1 <- nn.cells1[anchors[, "cell1"]]
    anchors2 <- nn.cells2[anchors[, "cell2"]]
    data.use1 <- data.use1[anchors1, ]
    data.use2 <- data.use2[anchors2, ]
    repSize("(anchors) data.use1", data.use1)
    repSize("(anchors) data.use2", data.use2)
    cat("dim: ", dim(data.use2), "\n")
    integration.matrix <- data.use2 - data.use1
    repSize("(FindIntegrationMatrix) integration.matrix", integration.matrix)
    object <- SetIntegrationData(object = object, integration.name = integration.name,
        slot = "integration.matrix", new.data = integration.matrix)
    repSize("(FindIntegrationMatrix) object", object)
    return(object)
}


