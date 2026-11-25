
LINEAGE_DIMS        <- 1:20
LINEAGE_N_PCS       <- 50
LINEAGE_N_HVG       <- 2000
LINEAGE_BATCH       <- "patient_ID"
LINEAGE_GRAPH       <- "harmony_snn"
LINEAGE_MIN_DIST    <- 0.3
LINEAGE_GS_CONTRASTS <- c("xALD_vs_ctrl", "sev3_vs_sev12")

lineage_embed <- function(obj, batch_var = LINEAGE_BATCH, dims = LINEAGE_DIMS, n_pcs = LINEAGE_N_PCS,
                          n_hvg = LINEAGE_N_HVG, min_dist = LINEAGE_MIN_DIST, graph_name = LINEAGE_GRAPH,
                          neighbors = TRUE, umap_name = "umap_harmony") {
  obj@reductions <- list()
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = n_hvg, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = n_pcs, verbose = FALSE)
  red <- "pca"
  if (nlevels(droplevels(factor(obj@meta.data[[batch_var]]))) > 1) {
    obj <- RunHarmony(obj, group.by.vars = batch_var, verbose = FALSE)
    red <- "harmony"
  }
  obj <- RunUMAP(obj, reduction = red, dims = dims, min.dist = min_dist,
                 reduction.name = umap_name, verbose = FALSE)
  if (neighbors)
    obj <- FindNeighbors(obj, reduction = red, dims = dims,
                         graph.name = c(sub("snn$", "nn", graph_name), graph_name), verbose = FALSE)
  obj
}

lineage_res_col <- function(r) sprintf("clusters_res%s", r)

lineage_cluster <- function(so, resolutions, base_res, seed = 1, graph_name = LINEAGE_GRAPH) {
  for (r in unique(c(resolutions, base_res))) {
    so <- FindClusters(so, graph.name = graph_name, algorithm = 4, resolution = r,
                       random.seed = seed, verbose = FALSE)
    so[[lineage_res_col(r)]] <- Idents(so)
  }
  so$seurat_clusters <- so[[lineage_res_col(base_res), drop = TRUE]]
  Idents(so) <- "seurat_clusters"
  so
}

ari <- function(a, b) {
  t <- table(a, b); s <- function(x) sum(choose(x, 2))
  e <- s(rowSums(t)) * s(colSums(t)) / choose(sum(t), 2)
  (s(t) - e) / ((s(rowSums(t)) + s(colSums(t))) / 2 - e)
}

lineage_seed_stability <- function(so, base_res, seed, seeds = c(1, 2, 3, 7, 42), graph_name = LINEAGE_GRAPH) {
  ref <- so$seurat_clusters
  purrr::map_dfr(setdiff(seeds, seed), function(s) {
    cl <- FindClusters(so, graph.name = graph_name, algorithm = 4, resolution = base_res,
                       random.seed = s, verbose = FALSE)$seurat_clusters
    tibble::tibble(seed_used = seed, seed_test = s, n_clusters = dplyr::n_distinct(cl),
                   ari = round(ari(ref, cl), 3))
  })
}

knn_purity <- function(obj, cluster_col, reduction = NULL, dims = LINEAGE_DIMS, k = 20,
                       by_cluster = FALSE) {
  if (is.null(reduction)) reduction <- if ("harmony" %in% Reductions(obj)) "harmony" else "pca"
  emb  <- Embeddings(obj, reduction)
  dims <- dims[dims <= ncol(emb)]
  nn   <- RANN::nn2(emb[, dims, drop = FALSE], k = k + 1)$nn.idx[, -1, drop = FALSE]
  cl   <- as.character(obj[[cluster_col, drop = TRUE]])
  same <- matrix(cl[nn] == cl, ncol = k)
  per_cell <- rowMeans(same, na.rm = TRUE)
  if (by_cluster) tapply(per_cell, cl, mean) else mean(per_cell[!is.na(cl)])
}

LINEAGE_EVIDENCE_GENES <- list(
  T        = c("CD3D", "CD3E", "TRAC", "CD2", "IL7R"),
  myeloid  = c("LYZ", "CD14", "CST3", "C1QA", "CD68"),
  B        = c("CD79A", "MS4A1", "IGKC"),
  NK       = c("NCAM1", "KLRF1", "SH2D1B"),
  platelet = c("PPBP", "PF4", "GP9", "ITGA2B", "TUBB1", "GNG11", "CAVIN2", "GP1BB"),
  cycling  = c("MKI67", "TOP2A"))

cluster_evidence <- function(so, cluster_col = "seurat_clusters", genes = LINEAGE_EVIDENCE_GENES,
                             donor_col = "patient_ID", dims = LINEAGE_DIMS) {
  cnt <- LayerData(so, assay = "RNA", layer = "counts")
  md  <- so@meta.data
  for (g in names(genes)) {
    n <- Matrix::colSums(cnt[intersect(genes[[g]], rownames(cnt)), , drop = FALSE] > 0)
    md[[paste0(".g_", g)]] <- n >= if (g == "cycling") 1 else 2
  }
  md$.ribo <- as.numeric(unlist(PercentageFeatureSet(so, pattern = "^RP[SL][0-9]")))
  md$.phase <- tryCatch(
    CellCycleScoring(so, s.features = cc.genes.updated.2019$s.genes,
                     g2m.features = cc.genes.updated.2019$g2m.genes)$Phase,
    error = function(e) NA_character_)
  pur <- knn_purity(so, cluster_col, by_cluster = TRUE, dims = dims)
  med <- stats::median(md$nCount_RNA)
  md %>%
    dplyr::mutate(.cl = as.character(.data[[cluster_col]]), .donor = .data[[donor_col]]) %>%
    dplyr::group_by(.cl) %>%
    dplyr::summarise(n = dplyr::n(), donors = dplyr::n_distinct(.donor),
                     max_donor_share = max(table(.donor)) / dplyr::n(),
                     count_ratio = stats::median(nCount_RNA) / med, nFeature = stats::median(nFeature_RNA),
                     mt = stats::median(percent.mt), ribo = stats::median(.ribo),
                     dplyr::across(dplyr::starts_with(".g_"), mean),
                     s_g2m = mean(.phase %in% c("S", "G2M")),
                     dbl_score = stats::median(scDblFinder.score), .groups = "drop") %>%
    dplyr::rename_with(~ sub("^\\.g_", "", .x)) %>%
    dplyr::mutate(knn_purity = as.numeric(pur[.cl])) %>%
    dplyr::rename(seurat_clusters = .cl)
}

top_markers_by_cluster <- function(object, ident_col = "seurat_clusters", n = 30,
                                   min.pct = 0.25, logfc.threshold = 0.25) {
  Idents(object) <- ident_col
  FindAllMarkers(object, only.pos = TRUE, min.pct = min.pct, logfc.threshold = logfc.threshold,
                 verbose = FALSE) %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::slice_max(order_by = .data$avg_log2FC, n = n) %>%
    dplyr::ungroup() %>%
    dplyr::select(cluster, gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2)
}

plot_qc_umap <- function(obj, vars = c("nCount_RNA", "nFeature_RNA", "percent.mt", "scDblFinder.score"),
                         reduction = "umap_harmony", pt.size = 0.4, title = NULL) {
  present <- intersect(vars, colnames(obj[[]]))
  plots <- purrr::map(present, function(v)
    Seurat::FeaturePlot(obj, features = v, reduction = reduction, pt.size = pt.size, order = TRUE) +
      ggplot2::scale_colour_viridis_c(option = "magma", name = NULL) +
      ggplot2::labs(title = v) + ggplot2::coord_equal() +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                     axis.title = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(size = 11, face = "bold")))
  patchwork::wrap_plots(plots, ncol = 2) +
    patchwork::plot_annotation(title = title,
                               theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14)))
}

dotplot_modules <- function(obj, ident_col, subpops, title = NULL, cluster_order = NULL) {
  names(subpops) <- make.unique(names(subpops))
  modules    <- lapply(subpops, intersect, rownames(obj))
  modules    <- modules[lengths(modules) > 0]
  genes_all  <- unlist(modules, use.names = FALSE)
  module_all <- rep(names(modules), lengths(modules))
  keep        <- !duplicated(genes_all)
  genes       <- genes_all[keep]
  gene_module <- stats::setNames(factor(module_all[keep], levels = unique(module_all[keep])), genes)

  lv <- if (is.null(cluster_order)) sort(unique(as.character(obj[[ident_col, drop = TRUE]]))) else cluster_order
  obj[[ident_col]] <- factor(as.character(obj[[ident_col, drop = TRUE]]), levels = rev(lv))
  Idents(obj) <- ident_col

  d <- DotPlot(obj, features = genes, dot.min = 0, scale = TRUE)$data
  d$module <- gene_module[as.character(d$features.plot)]
  ggplot(d, aes(features.plot, id)) +
    geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21, colour = "grey30", stroke = 0.25) +
    scale_fill_gradient2(low = "#3B4CC0", mid = "grey95", high = "#B40426", midpoint = 0, name = "scaled\nexpr") +
    scale_size(range = c(0, 6), limits = c(0, 100), name = "% detected") +
    facet_grid(~ module, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = NULL, title = title) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
          panel.grid.major = element_line(colour = "grey93"),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(size = 7.5, angle = 90),
          panel.spacing.x = unit(1.5, "pt"))
}

check_cluster_map <- function(objs, cluster_names, removed, key = "tissue", complete = TRUE) {
  for (nm in names(objs)) {
    have <- levels(droplevels(factor(objs[[nm]]$seurat_clusters)))
    map  <- cluster_names[[nm]]
    done <- c(names(map), removed$seurat_clusters[removed[[key]] == nm])
    if (complete && length(setdiff(have, done)))
      stop(nm, ": clusters without a label or removal row: ", paste(setdiff(have, done), collapse = ", "))
    if (length(setdiff(done, have)))
      stop(nm, ": listed clusters that do not exist: ", paste(setdiff(done, have), collapse = ", "))
    if (anyDuplicated(done)) stop(nm, ": cluster both labelled and removed, or listed twice")
    bad <- names(map)[map %in% c("NA", "", "TODO") | is.na(map)]
    if (length(bad)) stop(nm, ": clusters without a real label: ", paste(bad, collapse = ", "))
  }
  invisible(TRUE)
}

apply_cluster_map <- function(so, map, label_col, min_dist = LINEAGE_MIN_DIST, ...) {
  so[[label_col]] <- unname(map[as.character(so$seurat_clusters)])
  print(table(removed = is.na(so@meta.data[[label_col]]), so$diagnosis))
  lineage_embed(subset(so, cells = colnames(so)[!is.na(so@meta.data[[label_col]])]),
                min_dist = min_dist, neighbors = FALSE, ...)
}

ab_paired_compartment <- function(objs, label_col, donor_col = "patient_ID", diagnosis_col = "diagnosis",
                                  min_total = 100, se_max = 5) {
  meta <- dplyr::bind_rows(lapply(intersect(c("CSF", "PBMC"), names(objs)), function(ti) {
    md <- objs[[ti]]@meta.data
    tibble::tibble(donor = as.character(md[[donor_col]]), diagnosis = as.character(md[[diagnosis_col]]),
                   celltype = as.character(md[[label_col]]), compartment = ti)
  })) %>%
    dplyr::filter(!is.na(celltype)) %>%
    dplyr::mutate(sample_id = paste(donor, compartment, sep = "_"),
                  compartment = factor(compartment, levels = c("CSF", "PBMC")))

  dat <- meta %>%
    dplyr::count(donor, sample_id, compartment, diagnosis, celltype, name = "n") %>%
    tidyr::complete(tidyr::nesting(donor, sample_id, compartment, diagnosis), celltype, fill = list(n = 0)) %>%
    dplyr::group_by(sample_id) %>% dplyr::mutate(total = sum(n), prop = n / total) %>% dplyr::ungroup() %>%
    dplyr::filter(total >= min_total)
  paired <- dat %>% dplyr::distinct(donor, compartment) %>% dplyr::count(donor) %>%
    dplyr::filter(n == 2) %>% dplyr::pull(donor)
  dat <- dat %>% dplyr::filter(donor %in% paired) %>% dplyr::mutate(donor = factor(donor))
  if (length(paired) < 3) {
    message("paired CSF/PBMC: ", length(paired), " donors with both compartments -- skipped")
    return(NULL)
  }

  fit <- function(df) {
    m <- try(glmmTMB::glmmTMB(cbind(n, total - n) ~ compartment + (1 | donor),
                              family = glmmTMB::betabinomial(link = "logit"), data = df), silent = TRUE)
    if (inherits(m, "try-error") || !is.finite(stats::logLik(m)))
      return(tibble::tibble(estimate = NA_real_, std.error = NA_real_, p.value = NA_real_))
    broom.mixed::tidy(m, effects = "fixed") %>%
      dplyr::filter(term == "compartmentPBMC") %>% dplyr::select(estimate, std.error, p.value)
  }
  fit_all <- function(d) d %>% dplyr::group_modify(~ fit(.x)) %>% dplyr::ungroup() %>%
    dplyr::mutate(log2OR = estimate / log(2),
                  conf.low = (estimate - 1.96 * std.error) / log(2),
                  conf.high = (estimate + 1.96 * std.error) / log(2))

  med <- dat %>% dplyr::group_by(celltype, compartment) %>%
    dplyr::summarise(median_prop = stats::median(prop), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = compartment, values_from = median_prop, names_prefix = "median_")
  res <- dat %>% dplyr::group_by(celltype) %>% fit_all() %>%
    dplyr::mutate(padj = stats::p.adjust(p.value, "BH")) %>%
    dplyr::left_join(med, by = "celltype") %>% dplyr::arrange(padj)

  n_dx <- dat %>% dplyr::group_by(diagnosis, celltype) %>%
    dplyr::summarise(n_tot = sum(n), n_nonzero = sum(n > 0), .groups = "drop")
  by_dx <- dat %>% dplyr::group_by(diagnosis, celltype) %>% fit_all() %>%
    dplyr::group_by(diagnosis) %>% dplyr::mutate(padj = stats::p.adjust(p.value, "BH")) %>% dplyr::ungroup() %>%
    dplyr::left_join(n_dx, by = c("diagnosis", "celltype")) %>%
    dplyr::mutate(
      reliable = is.finite(std.error) & std.error <= se_max & n_tot >= 20 & n_nonzero >= 4,
      status = dplyr::case_when(!reliable ~ "not estimable",
                                padj < 0.05 & log2OR > 0 ~ "PBMC-enriched",
                                padj < 0.05 & log2OR < 0 ~ "CSF-enriched",
                                TRUE ~ "n.s."),
      stars = ifelse(!reliable, "", dplyr::case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**",
                                                     padj < 0.05 ~ "*", TRUE ~ "")))

  prop <- tryCatch({
    m2 <- meta %>% dplyr::filter(donor %in% paired)
    props <- speckle::getTransformedProps(clusters = m2$celltype, sample = m2$sample_id, transform = "asin")
    si <- m2 %>% dplyr::distinct(sample_id, donor, compartment) %>%
      dplyr::slice(match(colnames(props$TransformedProps), sample_id))
    design <- stats::model.matrix(~ 0 + compartment + donor, data = si)
    colnames(design)[1:2] <- c("CSF", "PBMC")
    speckle::propeller.ttest(props, design, contrasts = c(-1, 1, rep(0, ncol(design) - 2)),
                             robust = TRUE, trend = FALSE, sort = TRUE)
  }, error = function(e) { message("propeller: ", conditionMessage(e)); NULL })

  list(dat = dat, res = res, by_dx = by_dx, propeller = prop, n_paired = length(paired), se_max = se_max)
}

plot_ab_paired <- function(x, order = NULL, what = "cells", stem, dir) {
  if (is.null(x)) return(invisible(NULL))
  ord <- if (is.null(order)) sort(unique(x$by_dx$celltype)) else
    c(intersect(order, unique(x$by_dx$celltype)), sort(setdiff(unique(x$by_dx$celltype), order)))
  ok  <- x$by_dx %>% dplyr::filter(reliable)
  lim <- if (nrow(ok)) c(floor(min(ok$conf.low)), ceiling(max(ok$conf.high))) else c(-3, 3)
  d <- x$by_dx %>%
    dplyr::mutate(conf.low = pmax(conf.low, lim[1]), conf.high = pmin(conf.high, lim[2]),
                  log2OR_c = pmin(pmax(log2OR, lim[1]), lim[2]),
                  clipped = log2OR < lim[1] | log2OR > lim[2],
                  celltype = factor(celltype, levels = rev(ord)))
  p_forest <- ggplot(d, aes(log2OR_c, celltype, colour = status)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, linewidth = 0.5) +
    geom_point(aes(shape = clipped), size = 2.6, show.legend = FALSE) +
    geom_text(aes(x = conf.high, label = stars), hjust = -0.4, size = 4, show.legend = FALSE, na.rm = TRUE) +
    facet_wrap(~ diagnosis, nrow = 1) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 62)) +
    scale_colour_manual(values = c("CSF-enriched" = "#B2182B", "PBMC-enriched" = "#2166AC",
                                   "n.s." = "grey55", "not estimable" = "grey85")) +
    coord_cartesian(xlim = lim) +
    labs(x = "log2 odds ratio   (< 0 = CSF-enriched, > 0 = PBMC-enriched)", y = NULL, colour = NULL,
         caption = paste0("Beta-binomial GLMM per diagnosis group, donor as random effect; BH within group. ",
                          "Grey = SE > ", x$se_max, " or too few cells. * <0.05, ** <0.01, *** <0.001")) +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.y = element_blank(), legend.position = "bottom")
  p_pairs <- ggplot(x$dat %>% dplyr::mutate(celltype = factor(celltype, levels = ord)),
                    aes(compartment, 100 * prop)) +
    geom_line(aes(group = donor), colour = "grey70", linewidth = 0.4) +
    geom_point(aes(colour = diagnosis), size = 2) +
    facet_wrap(~ celltype, scales = "free_y") +
    labs(x = NULL, y = paste("% of", what, "in sample")) +
    theme_bw(base_size = 10)
  ggsave(file.path(dir, paste0(stem, "_forest.png")), p_forest, width = 12, height = 7, dpi = 300)
  ggsave(file.path(dir, paste0(stem, "_pairs.png")), p_pairs, width = 10, height = 7, dpi = 300)
  list(forest = p_forest, pairs = p_pairs)
}

gs_shift_report <- function(deg_tbl, objs, comparisons, levels_use, ma_levels, celltype_col, tag, dir) {
  shift_gene <- gs_shift_diagnostics(deg_tbl, comparisons, levels_use)
  kable_save(shift_gene, paste0("geneset_", tag, "_shift_gene_level"), digits = 3,
             caption = paste0(tag, ": Wald statistic vs expression level per DEG family"))
  fam <- shift_gene %>% dplyr::filter(level %in% ma_levels) %>% dplyr::select(dplyr::all_of(GS_KEYS))
  if (nrow(fam)) {
    p_ma <- gs_plot_shift_ma(deg_tbl, fam)
    print(p_ma)
    save_umap_fig(p_ma, paste0("geneset_", tag, "_shift_ma"), dir, width_mm = FIG_WIDTH_MM[["double"]],
                  height_mm = min(FIG_MAX_DEPTH_MM, 40 + 30 * ceiling(nrow(fam) / 3)))
  }
  qc <- gs_donor_shift_qc(objs, celltype_col)
  readr::write_csv(qc$donor_values, file.path(path$tables, paste0("geneset_", tag, "_shift_donor_values.csv")))
  kable_save(qc$summary %>% dplyr::arrange(dplyr::desc(abs(rho_within_group))),
             paste0("geneset_", tag, "_shift_donor_qc"), digits = 3,
             caption = paste0(tag, ": donor shift score vs QC metrics, between and within diagnosis groups"))
  p_qc <- qc$donor_values %>%
    tidyr::pivot_longer(c(n_cells, dplyr::starts_with("median_")), names_to = "qc", values_to = "qc_value") %>%
    ggplot2::ggplot(ggplot2::aes(qc_value, score, colour = group)) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::geom_smooth(ggplot2::aes(group = group), method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.4) +
    ggplot2::facet_grid(set + tissue ~ qc, scales = "free_x") +
    ggplot2::scale_colour_manual(values = c(control = "grey55", "X-ALD" = "#B2182B"), name = NULL) +
    ggplot2::labs(x = "QC metric (donor median)", y = "shift score (mean z)",
                  title = paste0(tag, ": does the global shift follow QC within the diagnosis groups?")) +
    gs_theme()
  print(p_qc)
  save_umap_fig(p_qc, paste0("geneset_", tag, "_shift_donor_qc"), dir, width_mm = FIG_WIDTH_MM[["double"]],
                height_mm = min(FIG_MAX_DEPTH_MM, 30 + 22 * dplyr::n_distinct(paste(qc$donor_values$set,
                                                                                     qc$donor_values$tissue))))
  invisible(list(gene = shift_gene, donor = qc))
}

lineage_covariates <- function(objs) {
  md_cols <- unique(unlist(lapply(objs, function(so) colnames(so@meta.data))))
  cols <- c(
    diagnosis      = "diagnosis",
    hydrocortisone = clin_resolve(md_cols, c("Hydrocortisone"), "hydrocort|cortison", "hydrocortisone"),
    lorenzo        = clin_resolve(md_cols, c("Lorenzos_Oil", "Lorenzo_Oil"), "lorenzo", "Lorenzo's oil"),
    concurrent     = clin_resolve(md_cols, c("concurrent_diagnosis"), "concurr|comorb|addison|adrenal",
                                  "concurrent diagnosis"),
    age            = clin_resolve(md_cols, c("age", "Age"), "^age", "age"))
  num <- function(v) suppressWarnings(as.numeric(gsub(",", ".", as.character(v))))
  donor_covariates(objs, cols) %>%
    dplyr::left_join(dplyr::select(CLIN_SEVERITY, donor = patient_ID, severity, serum_nfl), by = "donor") %>%
    dplyr::mutate(age = num(age), log10_nfl = log10(serum_nfl), sev_grp = sev_group3(severity),
                  severity_active = dplyr::case_when(severity == 3 ~ "active", severity %in% 1:2 ~ "lessactive"))
}

gs_donor_associations <- function(gs_donor) {
  p_wil <- function(v, g) {
    ok <- !is.na(v) & !is.na(g); g <- droplevels(factor(g[ok]))
    if (nlevels(g) != 2 || any(table(g) < 2)) return(NA_real_)
    suppressWarnings(stats::wilcox.test(v[ok] ~ g, exact = FALSE)$p.value)
  }
  spear <- function(v, u) {
    ok <- !is.na(v) & !is.na(u)
    if (sum(ok) < 4 || stats::sd(u[ok]) == 0) return(c(rho = NA_real_, p = NA_real_, n = sum(ok)))
    ct <- suppressWarnings(stats::cor.test(v[ok], u[ok], method = "spearman", exact = FALSE))
    c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
  }
  gs_donor %>%
    dplyr::group_by(tissue, celltype, set) %>%
    dplyr::group_modify(function(x, key) {
      xa <- x[x$diagnosis %in% "X-ALD", ]
      nf <- spear(xa$score, xa$log10_nfl)
      tibble::tibble(
        n_ctrl = sum(x$diagnosis %in% "control"), n_xald = nrow(xa),
        delta_dx = stats::median(xa$score) - stats::median(x$score[x$diagnosis %in% "control"]),
        p_dx = p_wil(x$score, x$diagnosis), p_sev = p_wil(xa$score, xa$severity_active),
        p_lorenzo = p_wil(xa$score, xa$lorenzo), p_hydrocortisone = p_wil(xa$score, xa$hydrocortisone),
        rho_nfl = nf[["rho"]], p_nfl = nf[["p"]], n_nfl = nf[["n"]])
    }) %>%
    dplyr::group_by(tissue) %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with("p_"), ~ stats::p.adjust(.x, "BH"),
                                .names = "padj_{sub('^p_', '', .col)}")) %>%
    dplyr::ungroup()
}

gs_treatment_check <- function(gs_res, gs_sets, tf_res, donor_assoc, tag) {
  gc_sets <- gs_sets$provenance$set[gs_sets$provenance$panel == "confound"]
  gc_tab <- gs_res %>% dplyr::filter(set %in% gc_sets) %>%
    dplyr::select(comparison, tissue, level, subset, set, NES, padj_fgsea, padj_camera)
  if (!is.null(tf_res))
    gc_tab <- gc_tab %>% dplyr::left_join(
      tf_res %>% dplyr::filter(tf == "NR3C1") %>%
        dplyr::select(comparison, tissue, level, subset, NR3C1_score = score, NR3C1_padj = padj),
      by = c("comparison", "tissue", "level", "subset"))
  print(kable_save(gc_tab %>% dplyr::arrange(comparison, tissue, level, subset, set),
                   paste0("glucocorticoid_check_", tag), digits = 3,
                   caption = "Glucocorticoid-response sets and NR3C1 activity per contrast"))
  donor_assoc %>%
    dplyr::filter(set %in% c(gc_sets, gs_sets$provenance$set[gs_sets$provenance$panel %in% c("sterol", "lipid", "VLCFA")])) %>%
    dplyr::select(tissue, celltype, set, p_hydrocortisone, p_lorenzo, p_sev) %>%
    dplyr::arrange(tissue, p_lorenzo) %>%
    kable_save(paste0("treatment_check_", tag), digits = 3, caption = paste(
      "Within X-ALD: per-donor set scores by hydrocortisone, Lorenzo's oil and stage 3 vs 1-2.",
      "NA = fewer than 2 donors in an arm."))
}
