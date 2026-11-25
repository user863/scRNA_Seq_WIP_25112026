
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)

library(Seurat)
library(SeuratObject)
library(scCustomize)
library(harmony)
library(SingleCellExperiment)
library(scDblFinder)
library(SoupX)
library(DropletUtils)

library(slingshot)
library(condiments)
library(tradeSeq)
library(BiocParallel)

library(CellChat)
library(ComplexHeatmap)
library(circlize)
library(magick)
library(ggwordcloud)

library(DESeq2)
library(edgeR)
library(limma)
library(speckle)
library(glmmTMB)
library(broom.mixed)
library(cluster)
library(RANN)

library(tidyverse)
library(dplyr)
library(tibble)
library(reshape2)
library(stringr)
library(readxl)
library(patchwork)
library(pheatmap)
library(ggrepel)
library(ggsci)
library(ggpubr)
library(Polychrome)
library(grid)
library(scales)
library(Matrix)
library(future)

library(qs)
library(here)
library(knitr)

`%!in%` <- Negate(`%in%`)

filter    <- dplyr::filter
select    <- dplyr::select
mutate    <- dplyr::mutate
group_by  <- dplyr::group_by
count     <- dplyr::count
rename    <- dplyr::rename
slice     <- dplyr::slice
desc      <- dplyr::desc
first     <- dplyr::first

lapply    <- base::lapply
unname    <- base::unname
intersect <- base::intersect
union     <- base::union
setdiff   <- base::setdiff
