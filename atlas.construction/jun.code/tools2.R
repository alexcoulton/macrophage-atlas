### Single Cell Analysis Helper Tools (specific for our projects)
### Modified 8 Jul 2022 by Jun Murai
### 8 Jul 2022 by Jun Murai

### Function List

# FeaturePlot2 : updated Feature Plot to build fixed size plot (for blank map plots)
# saveFPlotR2  : updated automated save feature plot, fit for blank map plots
# saveFPlot    : old saveFPlot
# saveFPlotR   : old saveFPlot
# info         : get info and change cluster numbers for our project data 

### variables
# palette36 : color palette
# palette33 : color palette

library(cowplot)
library(patchwork)
library(Polychrome)

# changelog
#   removed title and added label title
#   removed floor for xlim/ylim

data(palette36)
palette33 <- palette36[4:36]
names(palette33) <- NULL

saveFPlotR2 <- function(feature, min.cutoff=0, max.cutoff=5, ex=exp, fn="", prefix="", cols=c("#f4f4ee", "blue")) {  #  3 colors don't work. I don't know why.
  pr <- Project(exp)
  cl <- length(levels(Idents(exp)))
  if (fn == "") fn <- feature
  fn <- paste0(prefix, fn, "_", pr, cl, ".png")
  cat ("writing", fn, "...\n")
  fp <- FeaturePlot2(ex, features=feature, max.cutoff=max.cutoff, min.cutoff=min.cutoff, raster=TRUE, raster.dpi=c(2048, 2048), cols=cols, label=TRUE, order=TRUE) + coord_fixed()
  ggsave(file=fn, fp, width=9, height=7.2)
}

FeaturePlot2 <- function (object, features, dims = c(1, 2), cells = NULL, cols = if (blend) {
    c("lightgrey", "#ff0000", "#00ff00")
} else {
    c("lightgrey", "blue")
}, pt.size = NULL, order = FALSE, min.cutoff = NA, max.cutoff = NA,
    reduction = NULL, split.by = NULL, keep.scale = "feature",
    shape.by = NULL, slot = "data", blend = FALSE, blend.threshold = 0.5,
    label = FALSE, label.size = 4, label.color = "black", repel = FALSE,
    ncol = NULL, coord.fixed = FALSE, by.col = TRUE, sort.cell = NULL,
    interactive = FALSE, combine = TRUE, raster = NULL, raster.dpi = c(512,
        512))
{
    if (!is.null(x = sort.cell)) {
        warning("The sort.cell parameter is being deprecated. Please use the order ",
            "parameter instead for equivalent functionality.",
            call. = FALSE, immediate. = TRUE)
        if (isTRUE(x = sort.cell)) {
            order <- sort.cell
        }
    }
    if (interactive) {
        return(IFeaturePlot(object = object, feature = features[1],
            dims = dims, reduction = reduction, slot = slot))
    }
    if (!(is.null(x = keep.scale)) && !(keep.scale %in% c("feature",
        "all"))) {
        stop("`keep.scale` must be set to either `feature`, `all`, or NULL")
    }
    no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold",
            size = 14, margin = margin(r = 7)))
    reduction <- reduction %||% DefaultDimReduc(object = object)
    if (length(x = dims) != 2 || !is.numeric(x = dims)) {
        stop("'dims' must be a two-length integer vector")
    }
    if (blend && length(x = features) != 2) {
        stop("Blending feature plots only works with two features")
    }
    if (blend) {
        default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
        cols <- switch(EXPR = as.character(x = length(x = cols)),
            `0` = {
                warning("No colors provided, using default colors",
                  call. = FALSE, immediate. = TRUE)
                default.colors
            }, `1` = {
                warning("Only one color provided, assuming specified is double-negative and augmenting with default colors",
                  call. = FALSE, immediate. = TRUE)
                c(cols, default.colors[2:3])
            }, `2` = {
                warning("Only two colors provided, assuming specified are for features and agumenting with '",
                  default.colors[1], "' for double-negatives",
                  call. = FALSE, immediate. = TRUE)
                c(default.colors[1], cols)
            }, `3` = cols, {
                warning("More than three colors provided, using only first three",
                  call. = FALSE, immediate. = TRUE)
                cols[1:3]
            })
    }
    if (blend && length(x = cols) != 3) {
        stop("Blending feature plots only works with three colors; first one for negative cells")
    }
    dims <- paste0(Key(object = object[[reduction]]), dims)
    cells <- cells %||% colnames(x = object)
    data <- FetchData(object = object, vars = c(dims, "ident",
        features), cells = cells, slot = slot)
    if (ncol(x = data) < 4) {
        stop("None of the requested features were found: ", paste(features,
            collapse = ", "), " in slot ", slot, call. = FALSE)
    }
    else if (!all(dims %in% colnames(x = data))) {
        stop("The dimensions requested were not found", call. = FALSE)
    }
    features <- colnames(x = data)[4:ncol(x = data)]
    min.cutoff <- mapply(FUN = function(cutoff, feature) {
        return(ifelse(test = is.na(x = cutoff), yes = min(data[,
            feature]), no = cutoff))
    }, cutoff = min.cutoff, feature = features)
    max.cutoff <- mapply(FUN = function(cutoff, feature) {
        return(ifelse(test = is.na(x = cutoff), yes = max(data[,
            feature]), no = cutoff))
    }, cutoff = max.cutoff, feature = features)
    check.lengths <- unique(x = vapply(X = list(features, min.cutoff,
        max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
    if (length(x = check.lengths) != 1) {
        stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols,
        ]$maxcolors, no = length(x = cols))
    data[, 4:ncol(x = data)] <- sapply(X = 4:ncol(x = data),
        FUN = function(index) {
            data.feature <- as.vector(x = data[, index])
            min.use <- SetQuantile(cutoff = min.cutoff[index -
                3], data.feature)
            max.use <- SetQuantile(cutoff = max.cutoff[index -
                3], data.feature)
            data.feature[data.feature < min.use] <- min.use
            data.feature[data.feature > max.use] <- max.use
            if (brewer.gran == 2) {
                return(data.feature)
            }
            data.cut <- if (all(data.feature == 0)) {
                0
            }
            else {
                as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature),
                  breaks = brewer.gran)))
            }
            return(data.cut)
        })
    colnames(x = data)[4:ncol(x = data)] <- features
    rownames(x = data) <- cells
    data$split <- if (is.null(x = split.by)) {
        RandomName()
    }
    else {
        switch(EXPR = split.by, ident = Idents(object = object)[cells,
            drop = TRUE], object[[split.by, drop = TRUE]][cells,
            drop = TRUE])
    }
    if (!is.factor(x = data$split)) {
        data$split <- factor(x = data$split)
    }
    if (!is.null(x = shape.by)) {
        data[, shape.by] <- object[[shape.by, drop = TRUE]]
    }
    plots <- vector(mode = "list", length = ifelse(test = blend,
        yes = 4, no = length(x = features) * length(x = levels(x = data$split))))
    xlims <- c(min(data[, dims[1]]), max(data[,dims[1]]))
    ylims <- c(min(data[, dims[2]]), max(data[,dims[2]]))
    if (blend) {
        ncol <- 4
        color.matrix <- BlendMatrix(two.colors = cols[2:3], col.threshold = blend.threshold,
            negative.color = cols[1])
        cols <- cols[2:3]
        colors <- list(color.matrix[, 1], color.matrix[1, ],
            as.vector(x = color.matrix))
    }
    for (i in 1:length(x = levels(x = data$split))) {
        ident <- levels(x = data$split)[i]
        data.plot <- data[as.character(x = data$split) == ident,
            , drop = FALSE]
        if (blend) {
            features <- features[1:2]
            no.expression <- features[colMeans(x = data.plot[,
                features]) == 0]
            if (length(x = no.expression) != 0) {
                stop("The following features have no value: ",
                  paste(no.expression, collapse = ", "), call. = FALSE)
            }
            data.plot <- cbind(data.plot[, c(dims, "ident")],
                BlendExpression(data = data.plot[, features[1:2]]))
            features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
        }
        for (j in 1:length(x = features)) {
            feature <- features[j]
            if (blend) {
                cols.use <- as.numeric(x = as.character(x = data.plot[,
                  feature])) + 1
                cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
            }
            else {
                cols.use <- NULL
            }
            data.single <- data.plot[, c(dims, "ident", feature,
                shape.by)]
            plot <- SingleDimPlot(data = data.single, dims = dims,
                col.by = feature, order = order, pt.size = pt.size,
                cols = cols.use, shape.by = shape.by, label = FALSE,
                raster = raster, raster.dpi = raster.dpi) + scale_x_continuous(limits = xlims) +
                scale_y_continuous(limits = ylims) + theme_cowplot() +
                CenterTitle()
            if (label) {
                plot <- LabelClusters(plot = plot, id = "ident",
                  repel = repel, size = label.size, color = label.color)
            }
            if (length(x = levels(x = data$split)) > 1) {
                plot <- plot + theme(panel.border = element_rect(fill = NA,
                  colour = "black"))
                plot <- plot + labs(title = NULL)
                if (j == length(x = features) && !blend) {
                  suppressMessages(expr = plot <- plot + scale_y_continuous(sec.axis = dup_axis(name = ident),
                    limits = ylims) + no.right)
                }
                if (j != 1) {
                  plot <- plot + theme(axis.line.y = element_blank(),
                    axis.ticks.y = element_blank(), axis.text.y = element_blank(),
                    axis.title.y.left = element_blank())
                }
                if (i != length(x = levels(x = data$split))) {
                  plot <- plot + theme(axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(), axis.text.x = element_blank(),
                    axis.title.x = element_blank())
                }
            }
            else {
                plot <- plot + labs(title = NULL)
            }
            if (!blend) {
                plot <- plot + guides(color = NULL)
                cols.grad <- cols
                if (length(x = cols) == 1) {
                  plot <- plot + scale_color_brewer(palette = cols)
                }
                else if (length(x = cols) > 1) {
                  unique.feature.exp <- unique(data.plot[, feature])
                  if (length(unique.feature.exp) == 1) {
                    warning("All cells have the same value (",
                      unique.feature.exp, ") of ", feature, ".")
                    if (unique.feature.exp == 0) {
                      cols.grad <- cols[1]
                    }
                    else {
                      cols.grad <- cols
                    }
                  }
                  plot <- (expr = plot + scale_color_gradientn(colors = cols.grad,
                    guide = "colorbar") + xlab(feature) + ylab("") + theme(axis.title.x=element_text(face="bold")))    # + guides(color = guide_legend(title = feature))
                }
            }
            if (!(is.null(x = keep.scale)) && keep.scale == "feature" &&
                !blend) {
                max.feature.value <- max(data[, feature])
                min.feature.value <- min(data[, feature])
                plot <- suppressMessages(plot & scale_color_gradientn(colors = cols,
                  limits = c(min.feature.value, max.feature.value)))
            }
            if (coord.fixed) {
                plot <- plot + coord_fixed()
            }
            plot <- plot
            plots[[(length(x = features) * (i - 1)) + j]] <- plot
        }
    }
    if (blend) {
        blend.legend <- BlendMap(color.matrix = color.matrix)
        for (ii in 1:length(x = levels(x = data$split))) {
            suppressMessages(expr = plots <- append(x = plots,
                values = list(blend.legend + scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) >
                  1, yes = levels(x = data$split)[ii], no = "")),
                  expand = c(0, 0)) + labs(x = features[1], y = features[2],
                  title = if (ii == 1) {
                    paste("Color threshold:", blend.threshold)
                  } else {
                    NULL
                  }) + no.right), after = 4 * ii - 1))
        }
    }
    plots <- Filter(f = Negate(f = is.null), x = plots)
    if (is.null(x = ncol)) {
        ncol <- 2
        if (length(x = features) == 1) {
            ncol <- 1
        }
        if (length(x = features) > 6) {
            ncol <- 3
        }
        if (length(x = features) > 9) {
            ncol <- 4
        }
    }
    ncol <- ifelse(test = is.null(x = split.by) || blend, yes = ncol,
        no = length(x = features))
    legend <- if (blend) {
        "none"
    }
    else {
        split.by %iff% "none"
    }
    if (combine) {
        if (by.col && !is.null(x = split.by) && !blend) {
            plots <- lapply(X = plots, FUN = function(x) {
                return(suppressMessages(expr = x + theme_cowplot() +
                  ggtitle("") + scale_y_continuous(sec.axis = dup_axis(name = ""),
                  limits = ylims) + no.right))
            })
            nsplits <- length(x = levels(x = data$split))
            idx <- 1
            for (i in (length(x = features) * (nsplits - 1) +
                1):(length(x = features) * nsplits)) {
                plots[[i]] <- suppressMessages(expr = plots[[i]] +
                  scale_y_continuous(sec.axis = dup_axis(name = features[[idx]]),
                    limits = ylims) + no.right)
                idx <- idx + 1
            }
            idx <- 1
            for (i in which(x = 1:length(x = plots)%%length(x = features) ==
                1)) {
                plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
                  theme(plot.title = element_text(hjust = 0.5))
                idx <- idx + 1
            }
            idx <- 1
            if (length(x = features) == 1) {
                for (i in 1:length(x = plots)) {
                  plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
                    theme(plot.title = element_text(hjust = 0.5))
                  idx <- idx + 1
                }
                ncol <- 1
                nrow <- nsplits
            }
            else {
                nrow <- split.by %iff% length(x = levels(x = data$split))
            }
            plots <- plots[c(do.call(what = rbind, args = split(x = 1:length(x = plots),
                f = ceiling(x = seq_along(along.with = 1:length(x = plots))/length(x = features)))))]
            plots <- wrap_plots(plots, ncol = nrow, nrow = ncol)
            if (!is.null(x = legend) && legend == "none") {
                plots <- plots & NoLegend()
            }
        }
        else {
            plots <- wrap_plots(plots, ncol = ncol, nrow = split.by %iff%
                length(x = levels(x = data$split)))
        }
        if (!is.null(x = legend) && legend == "none") {
            plots <- plots & NoLegend()
        }
        if (!(is.null(x = keep.scale)) && keep.scale == "all" &&
            !blend) {
            max.feature.value <- max(data[, features])
            min.feature.value <- min(data[, features])
            plots <- suppressMessages(plots & scale_color_gradientn(colors = cols,
                limits = c(min.feature.value, max.feature.value)))
        }
    }
    return(plots)
}

saveFPlot <- function(feature, min.cutoff=0, max.cutoff=5, ex=mac, fn="") {
  if (fn == "") fn <- feature
  fn <- paste0(fn, "_mp.png")
  FeaturePlot(ex, features=feature, max.cutoff=max.cutoff, min.cutoff=min.cutoff, raster=FALSE, label=TRUE) %>% ggsave(file=fn)
}

saveFPlotR <- function(feature, min.cutoff=0, max.cutoff=5, ex=mac, fn="") {
  if (fn == "") fn <- feature
  fn <- paste0(fn, "_mp.png")
  FeaturePlot(ex, features=feature, max.cutoff=max.cutoff, min.cutoff=min.cutoff, raster=TRUE, raster.dpi=c(2048, 2048), cols=c("#f4f4f0", "blue"), label=TRUE, order=TRUE) %>% ggsave(file=fn, width=7, height=7)
}

info <- function(clusters = NA, resolution = NA) {

  # modify settings
  if (!is.na(clusters)) {
    res <- prj$sses$res[prj$sses$cl == clusters]
    cat ("Resolution changed to: ", res, "\n")
    Idents(exp) <<- paste0("integrated_snn_res.",res)
  }
  if (!is.na(resolution)) {
    Idents(exp) <<- paste0("integrated_snn_res.",resolution)
  }
  
  # show info
  cat ("Project Name: ",exp@project.name,"\n")
  cat ("Active Assay: ", exp@active.assay, "\n")
  cat ("Description: ",prj$descr,"\n")
  cat ("Current Number of Clusters:", length(levels(Idents(exp))), "\n")
  cat ("\nNumber of Clusters: \n")
  print (prj$sses %>% filter(pref==1))
  cat ("\nSeurat Object Information \n")
  cat ("------------------------- \n")
  print (exp)
}

distantPal80 <- function() {
 pal <- c()
 i <- 1
 for (a in c("00", "55", "AA", "FF")) {
   for (b in c("00", "40", "80", "C0", "FF")) {
     for (c in c("00", "55", "AA", "FF")) {
       pal[i] = paste0("#", a, b, c)
       cat (pal[i], ", \n");
       i = i + 1
     }
   }
 }
 pal[80] <- "#806080";  # white shouldn't be used
 return(pal);
}

distantPal49 <- function() {
 pal <- c()
 i <- 1
 for (a in c("00", "80", "FF")) {
   for (b in c("00", "55", "AA", "FF")) {
     for (c in c("00", "55", "AA", "FF")) {
       pal[i] = paste0("#", a, b, c)
       cat (pal[i], ", \n");
       i = i + 1
     }
   }
 }
 pal[46] <- "#40D030";   # avoid clash
 pal[48] <- "#408080";  
 pal[49] <- "#C08080";  
 return(pal);
}
