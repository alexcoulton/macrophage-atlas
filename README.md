# Using a pan-cancer atlas to investigate tumour associated macrophages as regulators of immunotherapy response

This repository contains code for atlas construction and figure plotting of the macrophage atlas.

Note that as some of the count data was obtained through contact of the study authors, users will need to obtain these data themselves in order to run the atlas construction code again.

Alternatively, users can download the complete atlas as a Seurat object from the Zenodo repository linked in the paper.

# Code structure

The atlas was constructed in two phases, first as a smaller atlas by Jun Murai, whose code is in `atlas.construction/jun.code`. This was then extended by Alexander Coulton to the final version of the atlas as presented in the paper `atlas.construction/alex.code`.

# Jun atlas construction code

### get_myeloid_1.R

load 14 experiments.

### get_myeloid_2.R

load 3 more experiments.

integrate them if ready.

### get_myeloid_3.R

split myeloid dataset to 3 datasets: macrophage, DC and monocyte

### get_myeloid_4_macrophage_R4.0.R

find markers from the macrophage dataset


# Alex atlas construction code

`atlas.construction/alex.code/new.clustering`:

Extract macrophages from individual study single cell objects using author's annotation if available. If not available, define macrophages via expression of a macrophage signature per cluster after performing clustering of single study data, then extract and save macrophages.

`atlas.construction/alex.code/load.all.macrophages`:

Combine all macrophage data, including Jun's intitial atlas and further individual datasets into a single Seurat object.

`atlas.construction/alex.code/all.mac.clus.batch.script2.R`:

Perform integration of studies using SCTransform and RPCA.

`atlas.construction/alex.code/macrophage.clustering.R`:

Perform clustering of macrophages.

`atlas.construction/alex.code/find.markers2.R`:

Find markers for each cluster.

# Figure code

`figures/mac.util.R`:

General utility functions for loading the final macrophage atlas

`figures/figX.R`:

Plotting code for the respective figures in the paper.

