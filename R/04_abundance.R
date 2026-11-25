ab_get_meta <- function(x, compartment_col) {
  if (inherits(x, "Seurat")) {
    md <- x[[]]
    stopifnot(compartment_col %in% colnames(md))
  } else if (is.list(x)) {
    md <- purrr::imap_dfr(x, function(o, nm) {
      m <- o[[]]
      m[[compartment_col]] <- nm
      m
    })
  } else {
    stop("ab_get_meta(): x is neither a Seurat object nor a list.")
  }
  md
}

ab_pull_chr <- function(df, col) {
  v <- df[[col]]
  if (is.list(v)) v <- vapply(v, function(e) as.character(e)[1], character(1))
  as.character(v)
}

build_abundance_inputs <- function(meta_all, compartment_col, focus,
                                   donor_col, group_col, label_col,
                                   strat_col = NULL, ref_group,
                                   min_cells = 5, min_samples = 0.5,
                                   ext_totals = NULL) {
  meta <- meta_all %>% dplyr::filter(.data[[compartment_col]] == focus)

  meta <- tibble::tibble(
    donor    = ab_pull_chr(meta, donor_col),
    group    = ab_pull_chr(meta, group_col),
    strat    = if (!is.null(strat_col) && strat_col %in% colnames(meta))
                 ab_pull_chr(meta, strat_col) else NA_character_,
    celltype = ab_pull_chr(meta, label_col)
  ) %>% dplyr::filter(!is.na(celltype))
  stopifnot(!any(vapply(meta, is.list, logical(1))))

  long <- meta %>%
    dplyr::count(donor, group, strat, celltype, name = "n") %>%
    tidyr::complete(tidyr::nesting(donor, group, strat), celltype, fill = list(n = 0))

  mat <- long %>%
    dplyr::select(donor, celltype, n) %>%
    tidyr::pivot_wider(names_from = donor, values_from = n) %>%
    tibble::column_to_rownames("celltype") %>%
    as.matrix()

  sinfo <- long %>%
    dplyr::distinct(donor, group, strat) %>%
    dplyr::arrange(match(donor, colnames(mat))) %>%
    tibble::column_to_rownames("donor")
  sinfo$group <- relevel(factor(sinfo$group), ref = ref_group)
  stopifnot(identical(rownames(sinfo), colnames(mat)))

  totals <- if (is.null(ext_totals)) colSums(mat) else ext_totals[colnames(mat)]
  stopifnot(!any(is.na(totals)), all(totals > 0))

  keep <- rowSums(mat >= min_cells) >= min_samples * ncol(mat)
  if (any(!keep)) message("dropped (low abundance): ", paste(rownames(mat)[!keep], collapse = ", "))
  mat <- mat[keep, , drop = FALSE]

  list(meta = meta, long = long, mat = mat, sinfo = sinfo, totals = totals)
}

fit_ab <- function(mat, sinfo, totals, design, coef = NULL, contrast = NULL, fdr = 0.05) {
  y <- DGEList(counts = mat, samples = sinfo)
  y$samples$lib.size     <- totals
  y$samples$norm.factors <- 1
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  res <- if (is.null(contrast)) glmQLFTest(fit, coef = coef) else glmQLFTest(fit, contrast = contrast)

  df_t <- res$df.total
  if (is.null(df_t)) df_t <- fit$df.total
  if (is.null(df_t)) df_t <- fit$df.prior + fit$df.residual.zeros
  if (is.null(df_t)) df_t <- fit$df.prior + fit$df.residual
  stopifnot(is.numeric(df_t), all(is.finite(df_t) | is.infinite(df_t)))

  tt <- topTags(res, n = Inf, sort.by = "none")$table %>% rownames_to_column("celltype")
  if (length(df_t) > 1) stopifnot(length(df_t) == nrow(tt))

  tt %>%
    mutate(
      df_total = df_t,
      se    = ifelse(.data$F > 0, abs(logFC) / sqrt(.data$F), NA_real_),
      ci_lo = logFC - qt(0.975, df_total) * se,
      ci_hi = logFC + qt(0.975, df_total) * se,
      sig   = FDR < fdr
    )
}

ab_panels <- function(long, res_main, totals, pal, ct_order = names(pal),
                      ref_group = NULL, y_lab = "proportion of cells", fdr = 0.05,
                      style = c("report", "figure"), family = "", x_lab = NULL) {
  style <- match.arg(style)
  on.exit(rm(list = setdiff(ls(), "ct_order")), add = TRUE)
  fig   <- style == "figure"
  bs    <- if (fig) FIG_BASE_PT else 11
  lw    <- if (fig) FIG_LINE_MM else 0.5

  ct_order <- intersect(ct_order, unique(as.character(long$celltype)))
  grp_lv <- sort(unique(as.character(long$group)))
  if (is.null(ref_group) && "control" %in% grp_lv) ref_group <- "control"
  if (!is.null(ref_group)) grp_lv <- c(ref_group, setdiff(grp_lv, ref_group))

  props <- long %>%
    dplyr::filter(celltype %in% ct_order) %>%
    dplyr::left_join(tibble::tibble(donor = names(totals), tot = as.numeric(totals)),
                     by = "donor") %>%
    dplyr::mutate(prop = n / tot,
                  celltype = factor(celltype, levels = rev(ct_order)),
                  group = factor(group, levels = grp_lv))

  th <- ggplot2::theme_classic(base_size = bs, base_family = family) +
    ggplot2::theme(
      axis.line  = ggplot2::element_line(linewidth = lw),
      axis.ticks = ggplot2::element_line(linewidth = lw),
      axis.text  = ggplot2::element_text(colour = "black",
                                         size = if (fig) FIG_LEGEND_PT else ggplot2::rel(0.8)),
      legend.position = "top",
      legend.margin   = ggplot2::margin(0, 0, 0, 0),
      legend.key.size = grid::unit(if (fig) 3 else 5, "mm"),
      plot.tag  = ggplot2::element_text(size = if (fig) FIG_PANEL_PT else bs, face = "bold"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA))

  dw <- 0.8
  comp <- ggplot2::ggplot(props, ggplot2::aes(celltype, prop)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = celltype, alpha = group),
                          outlier.shape = NA, width = 0.68, colour = "grey30",
                          linewidth = lw, staplewidth = 0.4,
                          position = ggplot2::position_dodge(width = dw)) +
    ggplot2::geom_point(ggplot2::aes(shape = group, group = group),
                        position = ggplot2::position_jitterdodge(
                          jitter.width = 0.1, jitter.height = 0, dodge.width = dw, seed = 1),
                        size = if (fig) 0.8 else 1.9, stroke = if (fig) 0.3 else 0.6) +
    ggplot2::scale_x_discrete(limits = rev(ct_order), drop = FALSE) +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_alpha_manual(values = setNames(c(0.25, 0.8)[seq_along(grp_lv)], grp_lv),
                                name = NULL) +
    ggplot2::scale_shape_manual(values = setNames(c(1, 16)[seq_along(grp_lv)], grp_lv),
                                name = NULL) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = y_lab) +
    th + if (!fig) ggplot2::theme(panel.grid.major.y =
           ggplot2::element_line(colour = "grey92", linewidth = 0.3))

  rd <- res_main %>%
    dplyr::filter(celltype %in% ct_order) %>%
    dplyr::mutate(celltype = factor(celltype, levels = rev(ct_order)),
                  sig = factor(!is.na(FDR) & FDR < fdr, levels = c(FALSE, TRUE)))
  forest <- ggplot2::ggplot(rd, ggplot2::aes(logFC, celltype)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = lw) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
                            colour = "grey35", linewidth = lw) +
    ggplot2::geom_point(ggplot2::aes(fill = celltype, alpha = sig), shape = 21,
                        size = if (fig) 1.8 else 3.4, stroke = if (fig) 0.3 else 0.6,
                        show.legend = FALSE) +
    ggplot2::scale_y_discrete(limits = rev(ct_order), drop = FALSE) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 1), drop = FALSE) +
    ggplot2::labs(x = if (is.null(x_lab)) expression(log[2]~"fold change in proportion") else x_lab, y = NULL) +
    th + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                        axis.ticks.y = ggplot2::element_blank(),
                        axis.line.y = ggplot2::element_blank())

  list(composition = comp, forest = forest, ct_order = ct_order, groups = grp_lv)
}

plot_ab_figure <- function(long, res_main, totals, pal, ct_order = names(pal),
                           title_a = "Per-donor composition",
                           title_b_prefix = "X-ALD vs control",
                           y_lab = "proportion of cells", ref_group = NULL,
                           fdr = 0.05, out_file = NULL,
                           width = 16, height = 8, x_lab = NULL, ci_label = "CI") {
  ab <- ab_panels(long, res_main, totals, pal, ct_order = ct_order,
                  ref_group = ref_group, y_lab = y_lab, fdr = fdr, style = "report", x_lab = x_lab)
  pA <- ab$composition + ggplot2::labs(title = paste("A ", title_a))
  pB <- ab$forest + ggplot2::labs(
    title = sprintf("B %s (%.0f%% %s, filled = FDR < %.2f)", title_b_prefix, 95, ci_label, fdr))
  p <- (pA | pB) + patchwork::plot_layout(widths = c(1.25, 1))
  if (!is.null(out_file)) ggsave(out_file, p, width = width, height = height)
  p
}

run_ab_within <- function(meta_all, par, compartments, tag, within_col, within_ref, key, contrast,
                          label, disease_value = "X-ALD", min_donors = 2, ext_totals = NULL, fig = NULL,
                          quiet = FALSE, fig_dir = path$abundance, tab_dir = path$tables) {
  out <- list()
  if (!within_col %in% colnames(meta_all)) {
    message(tag, ": column '", within_col, "' absent -- skipped"); return(invisible(out))
  }

  meta_w <- meta_all %>%
    dplyr::filter(.data[[par$group_col]] == disease_value, !is.na(.data[[within_col]]))
  if (!nrow(meta_w)) { message(tag, ": no ", disease_value, " cells -- skipped")
    return(invisible(out)) }

  for (d in c(fig_dir, tab_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  out_path <- function(dir, what, comp = NULL, ext = ".csv")
    file.path(dir, paste0("abundance_", tag, "_", key, "_", what,
                          if (!is.null(comp)) paste0("_", sanitize(comp)), ext))
  written <- character(0)

  don <- meta_w %>%
    dplyr::distinct(.data[[par$donor_col]], .data[[within_col]]) %>%
    stats::setNames(c("donor", "group"))
  don_n <- dplyr::count(don, group, name = "n_donors")
  readr::write_csv(don, f <- out_path(tab_dir, "donors")); written <- c(written, f)
  if (!quiet)
    print(knitr::kable(don_n, caption = paste0(disease_value, " donors per group, ", label, " — ", tag)))

  for (comp in compartments) {
    if (!quiet) cat("\n\n#### ", comp, " (", label, ")\n\n", sep = "")

    ab <- tryCatch(
      build_abundance_inputs(
        meta_w, compartment_col = par$compartment_col, focus = comp,
        donor_col = par$donor_col, group_col = within_col,
        label_col = par$label_col, strat_col = par$group_col,
        ref_group = within_ref, min_cells = par$min_cells,
        min_samples = par$min_samples, ext_totals = ext_totals),
      error = function(e) { message(comp, ": ", conditionMessage(e)); NULL })
    if (is.null(ab)) next

    n_grp <- table(droplevels(ab$sinfo$group))
    if (length(n_grp) < 2 || any(n_grp < min_donors)) {
      message(comp, ": ", paste(sprintf("%s=%d", names(n_grp), n_grp), collapse = ", "),
              " -- below min_donors = ", min_donors, ", skipped"); next
    }
    if (nrow(ab$mat) < 2) { message(comp, ": < 2 cell types left -- skipped"); next }

    design <- model.matrix(~ group, data = ab$sinfo)
    colnames(design) <- make.names(colnames(design))
    coef <- setdiff(colnames(design), "X.Intercept.")[1]

    res <- tryCatch(fit_ab(ab$mat, ab$sinfo, ab$totals, design,
                           coef = coef, fdr = par$fdr),
                    error = function(e) { message(comp, " fit: ", conditionMessage(e)); NULL })
    if (is.null(res)) next

    test_lvl <- setdiff(names(n_grp), within_ref)[1]
    res <- res %>% dplyr::arrange(FDR) %>%
      dplyr::mutate(compartment = comp, level = tag, contrast = contrast,
                    n_ref = as.integer(n_grp[[within_ref]]),
                    n_test = as.integer(n_grp[[test_lvl]]), .before = 1)

    cnt <- tibble::rownames_to_column(as.data.frame(cbind(ab$mat, total = ab$totals)),
                                      "celltype")
    readr::write_csv(res,     f1 <- out_path(tab_dir, "edgeR",       comp))
    readr::write_csv(ab$long, f2 <- out_path(tab_dir, "counts_long", comp))
    readr::write_csv(cnt,     f3 <- out_path(tab_dir, "counts",      comp))

    pal <- umap_palette_for(fig, comp, ab$long$celltype)
    p <- plot_ab_figure(
      ab$long, res, ab$totals, pal = pal,
      title_a = paste0("Per-donor composition, ", disease_value, " only — ", comp),
      title_b_prefix = paste0("EXPLORATORY: ", label),
      y_lab = paste("proportion of all labelled", comp, "cells"),
      fdr = par$fdr, out_file = f4 <- out_path(fig_dir, "edgeR", comp, ".png"))
    ggplot2::ggsave(f5 <- out_path(fig_dir, "edgeR", comp, ".pdf"), p, width = 16, height = 8)
    written <- c(written, f1, f2, f3, f4, f5)

    if (!quiet) {
      print(knitr::kable(cbind(ab$mat, total = ab$totals),
                         caption = paste0("Cell counts per donor — ", comp,
                                          " (", disease_value, " only)")))
      print(knitr::kable(
        res %>% dplyr::select(celltype, logFC, dplyr::any_of(c("ci_lo", "ci_hi")),
                              PValue, FDR, n_ref, n_test),
        digits = 3,
        caption = paste0("EXPLORATORY: ", disease_value, " ", label, ", ", comp,
                         " — log2 FC in proportion.")))
      print(p)
    }
    out[[comp]] <- list(inputs = ab, res = res)
  }

  if (length(out)) {
    readr::write_csv(dplyr::bind_rows(lapply(out, `[[`, "res")), f <- out_path(tab_dir, "edgeR_all"))
    written <- c(written, f)
  }
  message(tag, " ", label, ": ", length(written), " files written to ", tab_dir,
          if (!identical(fig_dir, tab_dir)) paste(" and", fig_dir))
  invisible(out)
}

run_ab_sev <- function(meta_all, par, compartments, tag, ...)
  run_ab_within(meta_all, par, compartments, tag, within_col = "severity_active", within_ref = "lessactive",
            key = "sev3_vs_sev12", contrast = "sev3_vs_sev12",
            label = "stage 3 vs stages 1-2 (stage 3 confounded with hydrocortisone)", ...)

plot_ab_by_severity <- function(meta, label_col, title = NULL, umap_fig = NULL, pal = NULL,
                                donor_col = "patient_ID", compartment_col = "compartment",
                                severity_col = "severity", levels = NULL, min_cells = 20,
                                out_file = NULL) {
  d <- tibble::tibble(donor       = as.character(meta[[donor_col]]),
                      compartment = as.character(meta[[compartment_col]]),
                      severity    = meta[[severity_col]],
                      celltype    = as.character(meta[[label_col]])) %>%
    dplyr::filter(!is.na(celltype), !is.na(severity)) %>%
    dplyr::count(compartment, donor, severity, celltype, name = "n") %>%
    tidyr::complete(tidyr::nesting(compartment, donor, severity), celltype, fill = list(n = 0)) %>%
    dplyr::group_by(compartment, donor) %>%
    dplyr::mutate(total = sum(n), pct = 100 * n / total) %>%
    dplyr::ungroup() %>%
    dplyr::filter(total >= min_cells)
  lv <- if (is.null(levels)) sort(unique(d$celltype)) else
    c(intersect(levels, unique(d$celltype)), setdiff(sort(unique(d$celltype)), levels))
  stage_mean <- d %>%
    dplyr::group_by(compartment, severity, celltype) %>%
    dplyr::summarise(pct = mean(pct), n_donors = dplyr::n_distinct(donor), .groups = "drop") %>%
    dplyr::mutate(celltype = factor(celltype, levels = lv),
                  stage = factor(severity, levels = 0:3))

  comps <- c(intersect(c("CSF", "PBMC"), unique(d$compartment)),
             setdiff(sort(unique(d$compartment)), c("CSF", "PBMC")))
  panel <- function(comp) {
    x <- stage_mean[stage_mean$compartment == comp & stage_mean$pct > 0, ]
    cols <- if (!is.null(pal)) pal else umap_palette_for(umap_fig, comp, lv)
    n_st <- dplyr::distinct(stage_mean[stage_mean$compartment == comp, ], stage, n_donors)
    n_st <- stats::setNames(n_st$n_donors, n_st$stage)
    ggplot2::ggplot(x, ggplot2::aes(stage, pct, fill = celltype)) +
      ggplot2::geom_col(width = 0.75, colour = "white", linewidth = 0.2,
                        position = ggplot2::position_stack(reverse = TRUE)) +
      ggplot2::scale_fill_manual(values = cols, breaks = lv, drop = TRUE, name = NULL,
                                 guide = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::scale_x_discrete(drop = FALSE, labels = function(s)
        ifelse(s %in% names(n_st), sprintf("%s\nn=%d", s, n_st[s]), s)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::labs(title = comp, x = "Severity stage", y = "Mean share per donor (%)") +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  }
  p <- patchwork::wrap_plots(lapply(comps, panel), nrow = 1, guides = "collect") +
    patchwork::plot_annotation(
      title = title,
      caption = paste0(paste(sprintf("%s %s", names(SEVERITY_SHORT), SEVERITY_SHORT), collapse = ", "),
                       "; donors with >= ", min_cells, " cells"))
  attr(p, "ab_data") <- d
  if (!is.null(out_file))
    ggplot2::ggsave(out_file, p, width = 2.6 * length(comps) + 2.8,
                    height = max(4.5, 0.22 * length(lv) + 1.5), dpi = 300)
  p
}

ab_pick_reference <- function(mat, min_share = 0.02) {
  prop <- sweep(mat, 2, colSums(mat), "/")
  ok <- rowSums(mat > 0) == ncol(mat) & apply(prop, 1, stats::median) >= min_share
  if (!any(ok)) return(NA_character_)
  cv <- apply(prop[ok, , drop = FALSE], 1, function(p) stats::sd(p) / mean(p))
  names(which.min(cv))
}

ab_compositional <- function(ab, res_share, ref = NULL, fdr = 0.05, run_sccomp = TRUE) {
  ref <- if (is.null(ref)) ab_pick_reference(ab$mat) else ref
  design <- model.matrix(~ group, data = ab$sinfo)
  colnames(design) <- make.names(colnames(design))
  coef <- setdiff(colnames(design), "X.Intercept.")[1]

  out <- res_share %>% dplyr::select(celltype, logFC_share = logFC, FDR_share = FDR)
  if (!is.na(ref)) {
    vs_ref <- fit_ab(ab$mat[rownames(ab$mat) != ref, , drop = FALSE], ab$sinfo,
                     ab$mat[ref, ], design, coef = coef, fdr = fdr) %>%
      dplyr::select(celltype, logFC_vs_ref = logFC, FDR_vs_ref = FDR)
    out <- dplyr::left_join(out, vs_ref, by = "celltype")
  }
  if (requireNamespace("speckle", quietly = TRUE)) {
    pr <- tryCatch(speckle::propeller(clusters = ab$meta$celltype, sample = ab$meta$donor,
                                      group = ab$meta$group, transform = "logit"),
                   error = function(e) NULL)
    if (!is.null(pr)) out <- dplyr::left_join(out, tibble::rownames_to_column(pr, "celltype") %>%
                                                dplyr::select(celltype, FDR_propeller = FDR), by = "celltype")
  }
  if (isTRUE(run_sccomp) && requireNamespace("sccomp", quietly = TRUE)) {
    sc <- tryCatch(
      ab$long %>% dplyr::filter(celltype %in% rownames(ab$mat)) %>%
        dplyr::mutate(group = stats::relevel(factor(group), ref = levels(ab$sinfo$group)[1])) %>%
        sccomp::sccomp_estimate(formula_composition = ~ group, .sample = donor, .cell_group = celltype,
                                .count = n, cores = 1, verbose = FALSE) %>%
        sccomp::sccomp_test(),
      error = function(e) { message("sccomp: ", conditionMessage(e)); NULL })
    if (!is.null(sc)) {
      sc <- as.data.frame(sc)
      sc <- sc[grepl("^group", sc$parameter), , drop = FALSE]
      out <- dplyr::left_join(out, tibble::tibble(celltype = as.character(sc$celltype),
                                                  effect_sccomp = sc$c_effect, lower_sccomp = sc$c_lower,
                                                  upper_sccomp = sc$c_upper, FDR_sccomp = sc$c_FDR),
                              by = "celltype")
    }
  }
  out %>%
    dplyr::mutate(reference = ref,
                  n_methods_sig = rowSums(dplyr::across(dplyr::starts_with("FDR_"), ~ .x < fdr), na.rm = TRUE)) %>%
    dplyr::arrange(FDR_share)
}

ab_res_for_figure <- function(res_share, comp) {
  if (!all(c("effect_sccomp", "lower_sccomp", "upper_sccomp", "FDR_sccomp") %in% names(comp))) {
    message("sccomp results missing: figure shows the edgeR share fit")
    return(list(res = res_share, method = "edgeR", x_lab = NULL, ci_label = "CI"))
  }
  list(res = dplyr::transmute(comp, celltype, logFC = effect_sccomp, ci_lo = lower_sccomp, ci_hi = upper_sccomp,
                              FDR = FDR_sccomp),
       method = "sccomp", x_lab = "sccomp composition effect (logit)", ci_label = "CrI")
}
