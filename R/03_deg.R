
DEG_RIBO_GENES <- c(
  "FAU", "RPS2", "RPS3", "RPS3A", "RPS4X", "RPS4Y1", "RPS4Y2", "RPS5", "RPS6",
  "RPS7", "RPS8", "RPS9", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15",
  "RPS15A", "RPS16", "RPS17", "RPS18", "RPS19", "RPS20", "RPS21", "RPS23",
  "RPS24", "RPS25", "RPS26", "RPS27", "RPS27A", "RPS27L", "RPS28", "RPS29", "RPSA",
  "RPL3", "RPL3L", "RPL4", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL7L1", "RPL8",
  "RPL9", "RPL10", "RPL10A", "RPL10L", "RPL11", "RPL12", "RPL13", "RPL13A",
  "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", "RPL19", "RPL21", "RPL22",
  "RPL22L1", "RPL23", "RPL23A", "RPL24", "RPL26", "RPL26L1", "RPL27", "RPL27A",
  "RPL28", "RPL29", "RPL30", "RPL31", "RPL32", "RPL34", "RPL35", "RPL35A",
  "RPL36", "RPL36A", "RPL36AL", "RPL37", "RPL37A", "RPL38", "RPL39", "RPL39L",
  "RPL41", "RPLP0", "RPLP1", "RPLP2", "UBA52")

DEG_HISTONE_PATTERN <- paste0("^(", paste(c(
  "H1-[0-9]+", "H2AC[0-9]+P?", "H2BC[0-9]+[LP]?", "H3C[0-9]+P?", "H4C[0-9]+P?", "H4-[0-9]+",
  "H2AX", "H2AZ[12]", "H2AJ", "H2AP", "H2AB[0-9]", "H2AL[0-9][A-Z0-9]*", "H2AQ1P",
  "MACROH2A[12]", "H2BW[0-9]+", "H2BL1P", "H2BK1", "H2BN1",
  "H3-[0-9]+[AB]?", "H3Y[12]", "H3P[0-9]+", "CENPA",
  "HIST[0-9]H[0-9A-Z]+", "H1F(0|X|OO|NT)", "H2AF[A-Z0-9]+", "H2BFWT", "H3F3[ABC]"),
  collapse = "|"), ")$")
DEG_MITO_PATTERN        <- "^MT-"
DEG_MITO_RIBO_PATTERN   <- "^MRP[LS][0-9]+[A-Z]?$"
DEG_LNC_UNNAMED_PATTERN <- "^[A-Z]{1,2}[0-9]{5,6}\\.[0-9]+$"

classify_genes <- function(g) {
  dplyr::case_when(
    grepl(DEG_HISTONE_PATTERN, g)     ~ "histone",
    g %in% DEG_RIBO_GENES             ~ "ribosomal",
    grepl(DEG_MITO_PATTERN, g)        ~ "mito",
    grepl(DEG_MITO_RIBO_PATTERN, g)   ~ "mito_ribosomal",
    grepl(DEG_LNC_UNNAMED_PATTERN, g) ~ "lnc_unnamed",
    TRUE                              ~ "other"
  )
}
DEG_DROP_CLASSES <- c("histone")
DEG_MASK_CLASSES <- c("ribosomal", "mito", "lnc_unnamed")

deg_reportable <- function(gene, status = NULL, mask = DEG_MASK_CLASSES) {
  ok <- !classify_genes(gene) %in% mask
  if (!is.null(status)) ok <- ok & status %in% "ok"
  ok
}

DEG_PADJ_THRESH   <- 0.05
DEG_LFC_THRESH    <- 0.5
DEG_SVALUE_THRESH <- 0.05
DEG_SHRINK        <- "ashr"
DEG_TREND_TEST    <- "LRT"

DEG_MIN_COUNT            <- 10
DEG_MIN_CELLS_PER_SAMPLE <- 10
DEG_MIN_DONORS_PER_GROUP <- 3
DEG_MIN_CELLS <- 2L * DEG_MIN_DONORS_PER_GROUP * DEG_MIN_CELLS_PER_SAMPLE

DEG_SF_MODE <- "default"

deg_is_sig <- function(padj, log2FC = NULL, padj_thresh = DEG_PADJ_THRESH, ...)
  !is.na(padj) & padj < padj_thresh

deg_is_highlight <- function(padj, log2FC, padj_thresh = DEG_PADJ_THRESH,
                             lfc_thresh = DEG_LFC_THRESH, svalue = NULL, ...) {
  sign_ok <- if (is.null(svalue)) TRUE else is.na(svalue) | svalue < DEG_SVALUE_THRESH
  deg_is_sig(padj, padj_thresh = padj_thresh) & !is.na(log2FC) & abs(log2FC) > lfc_thresh & sign_ok
}

deg_sign_uncertain <- function(svalue) !is.na(svalue) & svalue >= DEG_SVALUE_THRESH

deg_sig_label <- function(padj_thresh = DEG_PADJ_THRESH, ...) sprintf("padj < %s", padj_thresh)

deg_slim <- function(so) {
  a   <- SeuratObject::DefaultAssay(so)
  cnt <- tryCatch(SeuratObject::LayerData(so[[a]], layer = "counts"),
                  error = function(e) SeuratObject::LayerData(SeuratObject::JoinLayers(so[[a]]), layer = "counts"))
  SeuratObject::CreateSeuratObject(cnt, assay = a, meta.data = so@meta.data)
}

subset_cells <- function(so, keep, what = "filter") {
  idx <- which(keep %in% TRUE)
  if (length(idx) == 0) stop("subset_cells(): no cells left after ", what, call. = FALSE)
  so[, idx]
}

make_pseudobulk <- function(so, group_cols, min_cells = DEG_MIN_CELLS_PER_SAMPLE) {
  counts <- GetAssayData(so, layer = "counts")
  meta   <- so@meta.data
  key    <- do.call(paste, c(lapply(group_cols, function(x) as.character(meta[[x]])), sep = "__"))
  n      <- table(key)
  keep   <- names(n)[n >= min_cells]
  if (length(keep) < 2) stop("make_pseudobulk(): <2 usable samples", call. = FALSE)
  in_s <- which(key %in% keep)
  f    <- factor(key[in_s], levels = keep)
  mat  <- round(as.matrix(counts[, in_s, drop = FALSE] %*% Matrix::t(Matrix::fac2sparse(f))))
  dimnames(mat) <- list(rownames(counts), keep)
  mode(mat) <- "integer"
  sm <- as.data.frame(do.call(rbind, strsplit(keep, "__", fixed = TRUE)), stringsAsFactors = FALSE)
  colnames(sm) <- group_cols; rownames(sm) <- keep
  sm$n_cells <- as.integer(n[keep])
  list(counts = mat, meta = sm)
}

deg_size_factors <- function(dds, mode = DEG_SF_MODE) {
  switch(mode,
         default = estimateSizeFactors(dds),
         control_genes = estimateSizeFactors(dds, controlGenes = which(
           !rownames(dds) %in% unlist(gs_control_sets(rownames(dds)), use.names = FALSE))),
         tmm = {
           cnt <- counts(dds)
           sf  <- colSums(cnt) * edgeR::calcNormFactors(cnt, method = "TMM")
           sizeFactors(dds) <- sf / exp(mean(log(sf)))
           dds
         },
         stop("DEG_SF_MODE must be default | control_genes | tmm", call. = FALSE))
}

fit_dds <- function(counts_mat, sample_meta, design, sf_mode = DEG_SF_MODE, lrt_term = NULL) {
  dds  <- DESeqDataSetFromMatrix(counts_mat, sample_meta, design)
  dds  <- deg_size_factors(dds, sf_mode)
  dds  <- dds[rowSums(counts(dds)) >= DEG_MIN_COUNT &
                !(classify_genes(rownames(dds)) %in% DEG_DROP_CLASSES), ]
  S4Vectors::metadata(dds)$deg_sf_mode <- sf_mode
  if (is.null(lrt_term)) return(DESeq(dds, quiet = TRUE))
  reduced <- stats::update(design, stats::as.formula(paste("~ . -", lrt_term)))
  DESeq(dds, test = "LRT", reduced = reduced, quiet = TRUE)
}

get_res <- function(dds, contrast, shrink = DEG_SHRINK, alpha = DEG_PADJ_THRESH) {
  res  <- results(dds, contrast = contrast, alpha = alpha)
  test <- if ("LRTPvalue" %in% names(S4Vectors::mcols(dds))) "LRT" else "Wald"
  if (test == "LRT") res$stat <- sign(res$log2FoldChange) * sqrt(res$stat)
  res$log2FC_MLE <- res$log2FoldChange
  res$lfcSE_MLE  <- res$lfcSE
  res$svalue     <- NA_real_
  if (shrink != "none") {
    s <- lfcShrink(dds, contrast = contrast, res = res, type = shrink,
                   svalue = identical(shrink, "ashr"), quiet = TRUE)
    stopifnot(identical(rownames(s), rownames(res)))
    res$log2FoldChange <- s$log2FoldChange
    res$lfcSE          <- s$lfcSE
    if (!is.null(s$svalue)) res$svalue <- s$svalue
  }
  out <- as.data.frame(res)
  out$test            <- test
  out$shrink          <- shrink
  out$n_genes_tested  <- nrow(dds)
  ft <- S4Vectors::metadata(res)$filterThreshold
  out$filterThreshold <- if (is.null(ft)) NA_real_ else unname(ft[1])
  out$sf_mode         <- S4Vectors::metadata(dds)$deg_sf_mode %||% NA_character_
  out[order(out$padj), ]
}

deg_contrast_labels <- list(xALD_vs_ctrl = "X-ALD vs Control")

check_min_donors <- function(sm, group_col, min_donors = DEG_MIN_DONORS_PER_GROUP) {
  n <- table(sm[[group_col]])
  short <- n[n < min_donors]
  if (length(short))
    stop(sprintf("group(s) with <%d donors: %s", min_donors,
                 paste(sprintf("%s (n=%d)", names(short), short), collapse = ", ")), call. = FALSE)
  invisible(TRUE)
}

run_pseudobulk <- function(so, donor_col = "patient_ID", diagnosis_col = "diagnosis",
                           min_cells = DEG_MIN_CELLS_PER_SAMPLE,
                           return_dds = FALSE, sf_mode = DEG_SF_MODE) {
  pb <- make_pseudobulk(so, c(donor_col, diagnosis_col), min_cells)
  sm <- pb$meta
  names(sm)[names(sm) == diagnosis_col] <- "diagnosis"
  sm$diagnosis <- factor(sm$diagnosis, levels = c("control", "X-ALD"))
  if (anyNA(sm$diagnosis))
    warning("run_pseudobulk(): dropped samples with an unexpected diagnosis: ",
            toString(rownames(sm)[is.na(sm$diagnosis)]), call. = FALSE)
  sm <- sm[!is.na(sm$diagnosis), , drop = FALSE]
  if (nlevels(droplevels(sm$diagnosis)) < 2) stop("only one diagnosis group", call. = FALSE)
  check_min_donors(sm, "diagnosis")
  dds <- fit_dds(pb$counts[, rownames(sm), drop = FALSE], sm, ~ diagnosis, sf_mode)
  if (return_dds) return(dds)
  list(xALD_vs_ctrl = get_res(dds, c("diagnosis", "X-ALD", "control")))
}

deg_compare_results <- function(a, b, padj_thresh = DEG_PADJ_THRESH) {
  sig_of <- function(r) rownames(r)[deg_reportable(rownames(r)) & deg_is_sig(r$padj, padj_thresh = padj_thresh)]
  sa <- sig_of(a); sb <- sig_of(b); g <- intersect(rownames(a), rownames(b))
  tibble::tibble(n_sig_a = length(sa), n_sig_b = length(sb), n_shared = length(intersect(sa, sb)),
                 jaccard = length(intersect(sa, sb)) / max(1, length(union(sa, sb))),
                 cor_log2FC_MLE = stats::cor(a[g, "log2FC_MLE"], b[g, "log2FC_MLE"], use = "complete.obs"))
}

run_all_contrasts <- function(so, donor_col = "patient_ID", diagnosis_col = "diagnosis",
                              min_cells_per_sample = DEG_MIN_CELLS_PER_SAMPLE,
                              contrasts = deg_active_contrasts()) {
  failed <- function(msg, nms) setNames(as.list(rep(msg, length(nms))), nms)
  a <- if ("xALD_vs_ctrl" %in% contrasts)
    tryCatch(run_pseudobulk(so, donor_col, diagnosis_col, min_cells_per_sample),
             error = function(e) failed(conditionMessage(e), "xALD_vs_ctrl"))
  sev <- intersect(contrasts, DEG_SEVERITY_NAMES)
  s <- if (length(sev))
    tryCatch(run_pseudobulk_severity(so, donor_col, min_cells = min_cells_per_sample, which = sev),
             error = function(e) failed(conditionMessage(e), sev))
  c(a, s)
}

run_per_celltype <- function(so, celltype_col, fit_fun, fail_names, min_cells = DEG_MIN_CELLS) {
  failed <- function(msg) setNames(as.list(rep(msg, length(fail_names))), fail_names)
  ct_vals   <- as.character(so@meta.data[[celltype_col]])
  celltypes <- sort(unique(ct_vals[!is.na(ct_vals)]))
  res <- lapply(celltypes, function(ct) {
    so_sub <- subset_cells(so, ct_vals == ct, what = paste("celltype ==", ct))
    if (ncol(so_sub) < min_cells) return(failed(sprintf("skipped: %d cells < %d", ncol(so_sub), min_cells)))
    tryCatch(fit_fun(so_sub), error = function(e) failed(conditionMessage(e)))
  })
  stats::setNames(res, celltypes)
}

run_pseudobulk_per_celltype <- function(so, celltype_col, donor_col = "patient_ID",
                                        diagnosis_col = "diagnosis",
                                        min_cells = DEG_MIN_CELLS,
                                        min_cells_per_sample = DEG_MIN_CELLS_PER_SAMPLE) {
  run_per_celltype(so, celltype_col, min_cells = min_cells, fail_names = "xALD_vs_ctrl",
                   fit_fun = function(x) run_pseudobulk(x, donor_col, diagnosis_col, min_cells = min_cells_per_sample))
}

DEG_PROV_COLS <- c("test", "shrink", "n_genes_tested", "filterThreshold", "sf_mode")
DEG_TIDY_COLS <- c("baseMean", "log2FC", "lfcSE", "svalue", "log2FC_MLE", "lfcSE_MLE",
                   "stat", "pvalue", "padj", DEG_PROV_COLS)
DEG_CHR_COLS  <- c("test", "shrink", "sf_mode")

tidy_deg <- function(x, .tissue, .level, .subset, .comparison) {
  if (is.null(x)) return(NULL)
  keys <- tibble::tibble(tissue = .tissue, level = .level, subset = .subset, comparison = .comparison)
  if (is.character(x)) {
    row <- dplyr::mutate(keys, gene = NA_character_)
    for (cc in DEG_TIDY_COLS) row[[cc]] <- if (cc %in% DEG_CHR_COLS) NA_character_ else NA_real_
    row$status <- x
    return(row)
  }
  as.data.frame(x) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::rename(log2FC = log2FoldChange) %>%
    dplyr::transmute(tissue = .tissue, level = .level, subset = .subset, comparison = .comparison,
                     gene, dplyr::across(dplyr::any_of(DEG_TIDY_COLS)), status = "ok")
}

flatten_contrasts <- function(lst, .tissue, .level, .subset) {
  purrr::imap(lst, ~ tidy_deg(.x, .tissue, .level, .subset, .y)) %>% dplyr::bind_rows()
}

collect_deg_results <- function(deg_all, csv_prefix = NULL, padj_thresh = DEG_PADJ_THRESH, ...) {
  keys <- c("tissue", "level", "subset", "comparison")
  deg_all <- deg_all %>%
    dplyr::mutate(gene_class = classify_genes(gene),
                  family     = paste(tissue, level, subset, comparison, sep = " | "))
  if (is.null(csv_prefix)) return(deg_all)

  sig_now <- function(padj) deg_is_sig(padj, padj_thresh = padj_thresh)
  fam <- deg_all %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) %>%
    dplyr::summarise(status = dplyr::first(status),
                     dplyr::across(dplyr::any_of(DEG_PROV_COLS), dplyr::first),
                     n_significant = sum(deg_reportable(gene, status) & sig_now(padj)),
                     m = sum(!is.na(padj)), R = sum(sig_now(padj)), .groups = "drop") %>%
    dplyr::mutate(z_fcr = ifelse(R > 0, stats::qnorm(1 - (R * padj_thresh / m) / 2), NA_real_))
  readr::write_csv(dplyr::select(fam, -m, -R, -z_fcr), paste0(csv_prefix, "_families.csv"))
  readr::write_csv(dplyr::select(deg_all, -dplyr::any_of(DEG_PROV_COLS)), paste0(csv_prefix, "_all_genes.csv"))

  sig <- deg_all %>%
    dplyr::filter(deg_reportable(gene, status), sig_now(padj)) %>%
    dplyr::left_join(dplyr::select(fam, dplyr::all_of(keys), z_fcr), by = keys) %>%
    dplyr::mutate(highlight = deg_is_highlight(padj, log2FC, padj_thresh, svalue = svalue),
                  sign_uncertain = deg_sign_uncertain(svalue),
                  log2FC_MLE_lo_fcr = log2FC_MLE - z_fcr * lfcSE_MLE,
                  log2FC_MLE_hi_fcr = log2FC_MLE + z_fcr * lfcSE_MLE) %>%
    dplyr::select(-z_fcr)
  readr::write_csv(dplyr::select(sig, -dplyr::any_of(c(DEG_PROV_COLS, "status"))),
                   paste0(csv_prefix, "_significant.csv"))
  readr::write_csv(
    sig %>% dplyr::select(dplyr::any_of(c(
      keys, "gene", "baseMean", "log2FC", "svalue", "padj", "highlight", "sign_uncertain",
      "log2FC_MLE", "log2FC_MLE_lo_fcr", "log2FC_MLE_hi_fcr", "test",
      "loo_frac_retained", "loo_sign_flip", "loo_lost_without", "loo_status"))),
    paste0(csv_prefix, "_report.csv"))
  deg_all
}

plot_volcano <- function(deg_df, title, lfc_thresh = DEG_LFC_THRESH, padj_thresh = DEG_PADJ_THRESH,
                         lfc_col = "log2FC", padj_col = "padj", gene_col = NULL, subtitle = NULL,
                         label_subthreshold = TRUE, sub_max_overlaps = 12,
                         sub_label_size = 2.3, mask = DEG_MASK_CLASSES, ...) {
  df <- as.data.frame(deg_df)
  df$.lfc  <- df[[lfc_col]]
  df$.padj <- df[[padj_col]]
  df$.gene <- if (!is.null(gene_col)) df[[gene_col]] else rownames(df)
  df <- df[deg_reportable(df$.gene, df$status, mask), , drop = FALSE]
  pos <- df$.padj[!is.na(df$.padj) & df$.padj > 0]
  if (length(pos)) df$.padj[!is.na(df$.padj) & df$.padj == 0] <- min(pos)
  sv <- if ("svalue" %in% names(df)) df$svalue else NULL
  x_lab <- if (identical(lfc_col, "log2FC_MLE")) "log2 fold change (MLE, unshrunken)" else
    "log2 fold change (shrunken)"

  df %>%
    dplyr::mutate(.hl = deg_is_highlight(.padj, .lfc, padj_thresh, lfc_thresh, svalue = sv)) %>%
    dplyr::filter(!is.na(.padj)) %>%
    dplyr::mutate(
      .sig      = deg_is_sig(.padj, padj_thresh = padj_thresh),
      direction = dplyr::case_when(.sig & .lfc > 0 ~ "up", .sig & .lfc < 0 ~ "down", TRUE ~ "ns"),
      label     = ifelse(.hl, .gene, NA),
      label_sub = ifelse(.sig & !.hl, .gene, NA)) %>%
    ggplot(aes(x = .lfc, y = -log10(.padj), color = direction, label = label)) +
    geom_point(alpha = 0.6, size = 1.2) +
    { if (label_subthreshold)
      ggrepel::geom_text_repel(aes(label = label_sub), size = sub_label_size, colour = "grey45",
                               na.rm = TRUE, max.overlaps = sub_max_overlaps,
                               min.segment.length = 0.2, box.padding = 0.3,
                               segment.size = 0.2, segment.alpha = 0.4, show.legend = FALSE)
      else NULL } +
    ggrepel::geom_text_repel(size = 3, max.overlaps = Inf, na.rm = TRUE, min.segment.length = 0,
                             box.padding = 0.45, segment.size = 0.25, segment.alpha = 0.6) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c(up = "firebrick", down = "steelblue", ns = "grey70")) +
    labs(title = title, subtitle = subtitle, x = x_lab, y = paste0("-log10(", padj_col, ")"), color = NULL) +
    theme_classic() +
    theme(plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 7.5, face = "italic", colour = "grey35"),
          legend.position = "bottom")
}

plot_deg_summary <- function(deg_list, group_col_name = "celltype", title = NULL,
                             lfc_thresh = DEG_LFC_THRESH, padj_thresh = DEG_PADJ_THRESH,
                             lfc_col = "log2FoldChange", padj_col = "padj", subtitle = NULL,
                             gene_col = "gene", mask = DEG_MASK_CLASSES) {
  as_gene_df <- function(d) {
    d <- as.data.frame(d)
    if (!gene_col %in% names(d)) d[[gene_col]] <- rownames(d)
    d
  }
  if (!is.data.frame(deg_list))
    deg_list <- Filter(function(d) is.data.frame(d) || methods::is(d, "DataFrame"), deg_list)
  df <- if (is.data.frame(deg_list)) as_gene_df(deg_list) else
    dplyr::bind_rows(lapply(deg_list, as_gene_df), .id = group_col_name)
  df <- df[deg_reportable(df[[gene_col]], df$status, mask), , drop = FALSE]
  df$.lfc  <- df[[lfc_col]]
  df$.padj <- df[[padj_col]]

  long <- df %>%
    dplyr::group_by(.data[[group_col_name]]) %>%
    dplyr::summarise(
      up   =  sum(.lfc > 0 & deg_is_sig(.padj, padj_thresh = padj_thresh)),
      down = -sum(.lfc < 0 & deg_is_sig(.padj, padj_thresh = padj_thresh)),
      .groups = "drop") %>%
    tidyr::pivot_longer(cols = c(up, down), names_to = "direction", values_to = "n_genes") %>%
    dplyr::mutate(direction = factor(direction, levels = c("down", "up")))

  ggplot(long, aes(x = n_genes, y = reorder(.data[[group_col_name]], abs(n_genes)), fill = direction)) +
    geom_col() +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
    geom_text(data = . %>% dplyr::filter(abs(n_genes) > 0),
              aes(label = abs(n_genes), hjust = ifelse(n_genes >= 0, -0.2, 1.2)),
              size = 3, color = "black") +
    scale_fill_manual(values = c(up = "#CC0000", down = "#1A5276"),
                      labels = c(up = "Upregulated", down = "Downregulated")) +
    scale_x_continuous(labels = function(x) abs(x), expand = expansion(mult = 0.15)) +
    labs(x = sprintf("Number of DEGs (%s)", deg_sig_label(padj_thresh)),
         y = NULL, fill = NULL, title = title, subtitle = subtitle) +
    theme_classic() +
    theme(axis.text.y      = element_text(size = 9),
          axis.text.x      = element_text(size = 8),
          legend.position  = "bottom",
          legend.text      = element_text(size = 9),
          plot.title       = element_text(size = 11, face = "bold", hjust = 0.5),
          plot.subtitle    = element_text(size = 7.5, face = "italic", colour = "grey35", hjust = 0.5),
          panel.background = element_rect(fill = "white"),
          plot.background  = element_rect(fill = "white", color = NA))
}

plot_deg_summary_all <- function(df, group_col_name = "subset",
                                 comparison_col = "comparison",
                                 comparison_labels = NULL,
                                 title = NULL, subtitle = NULL,
                                 lfc_thresh = DEG_LFC_THRESH,
                                 padj_thresh = DEG_PADJ_THRESH,
                                 lfc_col = "log2FC", padj_col = "padj",
                                 gene_col = "gene", mask = DEG_MASK_CLASSES) {
  df <- as.data.frame(df)
  df <- df[deg_reportable(df[[gene_col]], df$status, mask), , drop = FALSE]
  df$.lfc  <- df[[lfc_col]]
  df$.padj <- df[[padj_col]]

  summ <- df %>%
    dplyr::group_by(.grp = .data[[group_col_name]], .cmp = .data[[comparison_col]]) %>%
    dplyr::summarise(
      up   =  sum(.lfc > 0 & deg_is_sig(.padj, padj_thresh = padj_thresh)),
      down = -sum(.lfc < 0 & deg_is_sig(.padj, padj_thresh = padj_thresh)),
      .groups = "drop")

  ord <- summ %>%
    dplyr::group_by(.grp) %>%
    dplyr::summarise(total = sum(abs(up)) + sum(abs(down)), .groups = "drop") %>%
    dplyr::arrange(total)
  summ$.grp <- factor(summ$.grp, levels = ord$.grp)

  if (!is.null(comparison_labels)) {
    lv <- intersect(names(comparison_labels), unique(as.character(summ$.cmp)))
    summ$.cmp <- factor(summ$.cmp, levels = lv, labels = unlist(comparison_labels)[lv])
  }

  long <- summ %>%
    tidyr::pivot_longer(cols = c(up, down), names_to = "direction", values_to = "n_genes") %>%
    dplyr::mutate(direction = factor(direction, levels = c("down", "up")))

  ggplot(long, aes(x = n_genes, y = .grp, fill = direction)) +
    geom_col() +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
    geom_text(data = . %>% dplyr::filter(abs(n_genes) > 0),
              aes(label = abs(n_genes), hjust = ifelse(n_genes >= 0, -0.25, 1.25)),
              size = 2.6, color = "black") +
    facet_wrap(~ .cmp, nrow = 1) +
    scale_fill_manual(values = c(up = "#CC0000", down = "#1A5276"),
                      labels = c(up = "Upregulated", down = "Downregulated")) +
    scale_x_continuous(labels = function(x) abs(x), expand = expansion(mult = 0.2)) +
    labs(x = sprintf("Number of DEGs (%s)", deg_sig_label(padj_thresh)),
         y = NULL, fill = NULL, title = title, subtitle = subtitle) +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom",
          strip.background = element_rect(fill = "grey95", color = NA),
          strip.text       = element_text(face = "bold", size = 8.5),
          plot.subtitle    = element_text(size = 7.5, face = "italic", colour = "grey35"),
          panel.background = element_rect(fill = "white"),
          plot.background  = element_rect(fill = "white", color = NA))
}

if (!exists("with_seed", mode = "function"))
  with_seed <- function(seed, expr) {
    g <- globalenv(); had <- exists(".Random.seed", envir = g, inherits = FALSE)
    if (had) old <- get(".Random.seed", envir = g)
    on.exit(if (had) assign(".Random.seed", old, envir = g) else
      if (exists(".Random.seed", envir = g, inherits = FALSE)) rm(".Random.seed", envir = g))
    set.seed(seed); force(expr)
  }

deg_null_labelings <- function(sm, mode = c("controls_split", "permute"), max_perm = 200,
                               seed = 1, min_group = DEG_MIN_DONORS_PER_GROUP) {
  mode <- match.arg(mode)
  pool <- if (mode == "controls_split") sm[sm$diagnosis == "control", , drop = FALSE] else sm
  n <- nrow(pool)
  k <- if (mode == "controls_split") floor(n / 2) else sum(sm$diagnosis == "X-ALD")
  if (min(k, n - k) < min_group)
    stop(sprintf("%d donors cannot form two groups of >= %d (mode '%s')", n, min_group, mode),
         call. = FALSE)
  idx <- with_seed(seed, {
    if (choose(n, k) <= max_perm) utils::combn(n, k, simplify = FALSE) else
      unique(replicate(max_perm, sort(sample.int(n, k)), simplify = FALSE))
  })
  if (mode == "permute") {
    real <- which(sm$diagnosis == "X-ALD")
    idx  <- Filter(function(i) !identical(sort(i), real), idx)
  } else if (2 * k == n) idx <- Filter(function(i) 1L %in% i, idx)
  lapply(idx, function(i) {
    m <- pool
    m$diagnosis <- factor(ifelse(seq_len(n) %in% i, "X-ALD", "control"),
                          levels = c("control", "X-ALD"))
    m
  })
}
