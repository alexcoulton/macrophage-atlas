### get myeloid 4 (R4.0 version)
### 13 Oct 2022 by Jun Murai

### Our integrated dataset is not compatible with R4.1
### So, we need to run this with R4.0, as far as we integrated dataset with R4.0

library(tidyverse)

exp  <- readRDS("macrophage_integrated_s17_20cl.rds")
meta <- readRDS("integrated_metadata_s17_220606.rds")  # "integrated_metadata_s17.rds"
prj  <- readRDS("project_info.rds")
exp@meta.data <- meta      # good to add latest metadata
Project(exp) <- "Mac"

exp <- PrepSCTFindMarkers(exp)          # doesn't work  - it means that we can't use FindMarker functions with this integrated dataset. instead we'll load markers data
markers20 <- FindAllMarkers(exp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, recorrect_umi=FALSE)
saveRDS(markers20, "markers_20.rds")
save(markers20, "markers.rda")

