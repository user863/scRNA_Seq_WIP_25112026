
path <- list()

path$metadata <- here("metadata")
path$data     <- here("data")
path$raw_data <- Sys.getenv("XALD_RAW_DATA", path$data)
path$sample_metadata  <- file.path(path$metadata, "sample_metadata.csv")
path$excluded_samples <- file.path(path$metadata, "excluded_samples.csv")

path$objects  <- here("objects")
path$graphs   <- here("graphs")
path$tables   <- here("tables")
path$figure_panels <- here("objects", "figure_panels")

path$umap_figures <- file.path(path$graphs, "umap_figures")
path$abundance    <- file.path(path$graphs, "abundance")
path$markers     <- file.path(path$graphs, "markers")
path$marker_umap <- file.path(path$graphs, "marker_umaps")

path$graphs_tcells <- here("graphs", "t_cells")

path$graphs_myeloid <- here("graphs", "myeloid")

path$graphs_bnk <- here("graphs", "bnk_cells")
path$tables_bnk <- here("tables", "bnk_cells")

path$graphs_cellchat           <- here("graphs", "cellchat")
path$graphs_cellchat_contrasts <- file.path(path$graphs_cellchat, "contrasts")

path$results_pseudotime <- here("results", "pseudotime_csf")
path$cache_pseudotime   <- file.path(path$results_pseudotime, "cache")

path$tables_nfl <- here("tables", "nfl")
path$graphs_nfl <- here("graphs", "nfl")

path$tables_multinichenet <- here("tables", "multinichenet")
path$graphs_multinichenet <- here("graphs", "multinichenet")

invisible(lapply(path, dir.create, showWarnings = FALSE, recursive = TRUE))

sanitize <- function(x) gsub("_+", "_", gsub("[^A-Za-z0-9]+", "_", x))
