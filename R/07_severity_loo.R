
CLIN_SEVERITY <- readr::read_csv(here::here("metadata", "clinical_severity.csv"),
                                 col_types = readr::cols(patient_ID = "c", LD_ID = "c", severity = "i",
                                                         serum_nfl = "d"),
                                 na = c("", "NA")) |>
  as.data.frame() |> tibble::as_tibble()
stopifnot(!anyDuplicated(CLIN_SEVERITY$patient_ID))
SEVERITY_SHORT <- c(`0` = "control", `1` = "AMN", `2` = "smoldering CALD", `3` = "active CALD")
SEVERITY_COLS  <- c(sev0 = "grey55", sev1 = "#56B4E9", sev2 = "#E69F00", sev3 = "#B2182B")
readr::write_csv(CLIN_SEVERITY, file.path(path$tables, "clinical_severity.csv"))

add_severity <- function(x, donor_col = "patient_ID", quiet = FALSE) {
  if (is.list(x) && !inherits(x, "Seurat"))
    return(lapply(x, add_severity, donor_col = donor_col, quiet = quiet))
  d <- as.character(x@meta.data[[donor_col]])
  i <- match(d, CLIN_SEVERITY$patient_ID)
  x$LD_ID          <- CLIN_SEVERITY$LD_ID[i]
  x$severity       <- CLIN_SEVERITY$severity[i]
  x$severity_grp   <- ifelse(is.na(x$severity), NA_character_, paste0("sev", x$severity))
  x$severity_label <- unname(SEVERITY_SHORT[as.character(x$severity)])
  x$serum_nfl      <- CLIN_SEVERITY$serum_nfl[i]
  x$severity_cald   <- dplyr::case_when(x$severity == 0 ~ "control", x$severity %in% 2:3 ~ "CALD")
  x$severity_active <- dplyr::case_when(x$severity == 3 ~ "active", x$severity %in% 1:2 ~ "lessactive")
  if (!quiet) {
    miss <- unique(d[is.na(i)])
    if (length(miss)) warning("severity missing for donor(s): ", toString(miss), call. = FALSE)
    if ("diagnosis" %in% colnames(x@meta.data)) {
      dx  <- as.character(x$diagnosis)
      bad <- unique(d[!is.na(x$severity) & !is.na(dx) & ((x$severity == 0) != (dx == "control"))])
      if (length(bad)) warning("severity 0 and diagnosis 'control' disagree for: ", toString(bad),
                               call. = FALSE)
    }
  }
  x
}

severity_report <- function(objs, name, donor_col = "patient_ID") {
  tab <- purrr::imap_dfr(objs, function(so, ti) {
    md <- so@meta.data
    tibble::tibble(tissue = ti, donor = as.character(md[[donor_col]]),
                   diagnosis = as.character(md$diagnosis),
                   CI = if ("cerebral_inflammation" %in% names(md)) as.character(md$cerebral_inflammation) else NA,
                   severity = md$severity) %>% dplyr::distinct()
  }) %>% dplyr::arrange(tissue, severity, donor) %>%
    dplyr::mutate(flag = dplyr::case_when(is.na(severity) ~ "no severity",
                                          (severity == 0) != (diagnosis == "control") ~ "severity / diagnosis disagree",
                                          TRUE ~ ""))
  readr::write_csv(tab, file.path(path$tables, paste0("clinical_severity_", name, ".csv")))
  print(table(tab$tissue, paste0("sev", tab$severity), useNA = "ifany"))
  invisible(tab)
}

SEV_CONTRASTS <- list(
  sev1_vs_ctrl = c("sev1", "sev0"), sev2_vs_ctrl = c("sev2", "sev0"), sev3_vs_ctrl = c("sev3", "sev0"),
  sev2_vs_sev1 = c("sev2", "sev1"), sev3_vs_sev1 = c("sev3", "sev1"), sev3_vs_sev2 = c("sev3", "sev2"))
SEV_POOLED <- list(
  sev23_vs_ctrl = list(test = 2:3, ref = 0L),
  sev3_vs_sev12 = list(test = 3L,  ref = 1:2))
SEV_TRENDS         <- c("sev_trend", "sev_trend_xald")
DEG_SEVERITY_NAMES <- c(SEV_TRENDS, names(SEV_POOLED), names(SEV_CONTRASTS))
deg_contrast_labels <- c(deg_contrast_labels, list(
  sev_trend      = "Severity trend 0-3 (per stage)",
  sev_trend_xald = "Severity trend within X-ALD 1-3 (per stage)",
  sev1_vs_ctrl   = "AMN (1) vs control",
  sev2_vs_ctrl   = "Smoldering CALD (2) vs control",
  sev3_vs_ctrl   = "Active CALD (3) vs control",
  sev2_vs_sev1   = "Smoldering CALD (2) vs AMN (1)",
  sev3_vs_sev1   = "Active CALD (3) vs AMN (1)",
  sev3_vs_sev2   = "Active (3) vs smoldering CALD (2)",
  sev23_vs_ctrl  = "CALD (2+3) vs control",
  sev3_vs_sev12  = "Active CALD (3) vs AMN + smoldering (1+2)"))
SEV_CONFOUND_NOTE <- "Stage 3 coincides with hydrocortisone treatment; may reflect steroid exposure."
is_sev3_comparison <- function(comparison) grepl("sev3|sev_trend", comparison)
sev_group3 <- function(severity)
  factor(dplyr::case_when(severity == 0 ~ "control", severity %in% 1:2 ~ "sev12", severity == 3 ~ "sev3"),
         levels = c("control", "sev12", "sev3"))
SEV_GROUP3_LABELS <- c(control = "Control", sev12 = "X-ALD stage 1-2", sev3 = "X-ALD stage 3")
SEV_GROUP3_COLS   <- c(control = "grey55", sev12 = "#0072B2", sev3 = "#D55E00")
DEG_SEV_REPORT <- c("sev_trend", "sev_trend_xald", "sev3_vs_ctrl", "sev3_vs_sev12")
deg_report_contrasts <- function() c("xALD_vs_ctrl", DEG_SEV_REPORT)
deg_active_contrasts <- function() c("xALD_vs_ctrl", DEG_SEVERITY_NAMES)

.sev_pb <- function(so, donor_col = "patient_ID", sev_col = "severity", min_cells = DEG_MIN_CELLS_PER_SAMPLE) {
  if (!sev_col %in% colnames(so@meta.data))
    stop("no '", sev_col, "' column -- run add_severity() on the object first", call. = FALSE)
  pb <- make_pseudobulk(so, c(donor_col, sev_col), min_cells = min_cells)
  sm <- pb$meta
  names(sm)[names(sm) == sev_col] <- "severity"
  sm$severity <- suppressWarnings(as.integer(sm$severity))
  sm <- sm[!is.na(sm$severity), , drop = FALSE]
  if (anyDuplicated(sm[[donor_col]])) stop("a donor has more than one severity", call. = FALSE)
  list(counts = pb$counts, meta = sm)
}

.sev_spec <- function(sm, nm) {
  if (nm %in% SEV_TRENDS) {
    d <- sm[if (nm == "sev_trend") rep(TRUE, nrow(sm)) else sm$severity >= 1, , drop = FALSE]
    if (length(unique(d$severity)) < 2 || nrow(d) < 4) stop("fewer than 2 stages or 4 donors", call. = FALSE)
    d$severity_c <- d$severity - mean(d$severity)
    return(list(meta = d, design = ~ severity_c, ct = "severity_c",
                lrt = if (identical(DEG_TREND_TEST, "LRT")) "severity_c"))
  }
  if (nm %in% names(SEV_POOLED)) {
    pc <- SEV_POOLED[[nm]]
    d  <- sm[sm$severity %in% c(pc$test, pc$ref), , drop = FALSE]
    d$pool <- factor(ifelse(d$severity %in% pc$test, "test", "ref"), levels = c("ref", "test"))
    check_min_donors(d, "pool")
    if (nlevels(droplevels(d$pool)) < 2) stop("one side empty", call. = FALSE)
    return(list(meta = d, design = ~ pool, ct = c("pool", "test", "ref")))
  }
  ct <- SEV_CONTRASTS[[nm]]
  d  <- sm
  d$sev_grp <- droplevels(factor(paste0("sev", d$severity), levels = paste0("sev", 0:3)))
  n <- table(d$sev_grp)
  if (nlevels(d$sev_grp) < 2) stop("fewer than 2 stages", call. = FALSE)
  if (!all(ct %in% names(n))) stop("stage(s) absent: ", toString(setdiff(ct, names(n))), call. = FALSE)
  if (any(n[ct] < DEG_MIN_DONORS_PER_GROUP))
    stop(sprintf("group(s) with <%d donors: %s", DEG_MIN_DONORS_PER_GROUP,
                 paste(sprintf("%s (n=%d)", ct, as.integer(n[ct])), collapse = ", ")), call. = FALSE)
  list(meta = d, design = ~ sev_grp, ct = c("sev_grp", ct[1], ct[2]))
}

run_pseudobulk_severity <- function(so, donor_col = "patient_ID", sev_col = "severity",
                                    min_cells = DEG_MIN_CELLS_PER_SAMPLE, which = DEG_SEVERITY_NAMES) {
  pb  <- .sev_pb(so, donor_col, sev_col, min_cells)
  cnt <- pb$counts; sm <- pb$meta
  out <- stats::setNames(vector("list", length(which)), which)
  for (nm in intersect(which, c(SEV_TRENDS, names(SEV_POOLED))))
    out[[nm]] <- tryCatch({
      s <- .sev_spec(sm, nm)
      get_res(fit_dds(cnt[, rownames(s$meta), drop = FALSE], s$meta, s$design, lrt_term = s$lrt),
              if (length(s$ct) == 1) list(s$ct) else s$ct)
    }, error = function(e) trimws(conditionMessage(e)))
  pairs <- intersect(which, names(SEV_CONTRASTS))
  if (!length(pairs)) return(out)
  sm$sev_grp <- droplevels(factor(paste0("sev", sm$severity), levels = paste0("sev", 0:3)))
  dds <- if (nlevels(sm$sev_grp) >= 2)
    tryCatch(fit_dds(cnt[, rownames(sm), drop = FALSE], sm, ~ sev_grp), error = function(e) conditionMessage(e))
  else "fewer than 2 stages"
  for (nm in pairs)
    out[[nm]] <- if (is.character(dds)) dds else
      tryCatch({ s <- .sev_spec(sm, nm); get_res(dds, s$ct) }, error = function(e) trimws(conditionMessage(e)))
  out
}

gs_registry_severity <- function(units, contrasts, donor_col = "patient_ID",
                                 min_cells = DEG_MIN_CELLS_PER_SAMPLE) {
  reg <- list()
  for (u in units) {
    so <- tryCatch(u$get(), error = function(e) NULL)
    pb <- if (is.null(so)) NULL else tryCatch(.sev_pb(so, donor_col, "severity", min_cells), error = function(e) NULL)
    if (is.null(pb)) next
    for (nm in intersect(contrasts, DEG_SEVERITY_NAMES)) {
      e <- tryCatch({
        s   <- .sev_spec(pb$meta, nm)
        dds <- DESeq2::DESeqDataSetFromMatrix(pb$counts[, rownames(s$meta), drop = FALSE], s$meta, s$design)
        keep <- rowSums(DESeq2::counts(dds)) >= DEG_MIN_COUNT & !(classify_genes(rownames(dds)) %in% DEG_DROP_CLASSES)
        list(dds = dds[keep, ], ct = s$ct)
      }, error = function(err) {
        message("gs registry ", u$tissue, " | ", u$level, " | ", u$subset, " | ", nm, ": ", conditionMessage(err)); NULL })
      if (!is.null(e)) reg <- gs_registry_add(reg, u$tissue, u$level, u$subset, e$dds, stats::setNames(list(e$ct), nm))
    }
  }
  reg
}

run_pseudobulk_severity_per_celltype <- function(so, celltype_col, donor_col = "patient_ID",
                                                 min_cells = DEG_MIN_CELLS,
                                                 min_cells_per_sample = DEG_MIN_CELLS_PER_SAMPLE) {
  celltypes <- unique(stats::na.omit(as.character(so@meta.data[[celltype_col]])))
  res <- lapply(celltypes, function(ct) {
    so_sub <- subset_cells(so, so@meta.data[[celltype_col]] == ct, what = paste("celltype ==", ct))
    if (ncol(so_sub) < min_cells) return(NULL)
    tryCatch(run_pseudobulk_severity(so_sub, donor_col, min_cells = min_cells_per_sample),
             error = function(e) NULL)
  })
  names(res) <- celltypes
  Filter(Negate(is.null), res)
}

DEG_LOO_TISSUES <- c("CSF", "PBMC")
DEG_LOO_CORES   <- as.integer(Sys.getenv("DEG_LOO_CORES", max(1L, min(8L, parallel::detectCores() - 1L))))
DEG_LOO_RESUME  <- TRUE
DEG_LOO_STABLE  <- 0.8
DEG_LOO_FRAGILE <- 0.6
DEG_LOO_CONTRASTS    <- c("xALD_vs_ctrl", "sev3_vs_ctrl", "sev3_vs_sev12")
DEG_LOO_GS_CONTRASTS <- c("xALD_vs_ctrl", "sev3_vs_sev12")
DEG_LOO_DIR     <- file.path(path$objects, "loo")
dir.create(DEG_LOO_DIR, recursive = TRUE, showWarnings = FALSE)
LOO_KEYS <- c("tissue", "level", "subset", "comparison")

deg_loo_units <- function(objs, tissues, level_all, celltype_col, lineage_col = NULL,
                          min_cells = DEG_MIN_CELLS) {
  u <- list()
  for (ti in intersect(tissues, names(objs))) {
    so <- deg_slim(objs[[ti]])
    u[[length(u) + 1]] <- list(tissue = ti, level = level_all, subset = "all", get = local({ s <- so; function() s }))
    add <- function(col, level) {
      v <- as.character(so@meta.data[[col]])
      for (lab in sort(unique(stats::na.omit(v)))) {
        idx <- which(v == lab)
        if (length(idx) < min_cells) next
        u[[length(u) + 1]] <<- list(tissue = ti, level = level, subset = lab,
                                    get = local({ s <- so; i <- idx; function() s[, i] }))
      }
    }
    if (!is.null(lineage_col)) add(lineage_col, "lineage")
    add(celltype_col, "celltype")
  }
  u
}

loo_with_min_donors <- function(f, min_d) {
  env <- environment(check_min_donors)
  old <- get("DEG_MIN_DONORS_PER_GROUP", envir = env)
  assign("DEG_MIN_DONORS_PER_GROUP", min_d, envir = env)
  on.exit(assign("DEG_MIN_DONORS_PER_GROUP", old, envir = env), add = TRUE)
  f()
}

loo_genes <- function(full, folds, n_fail, n_total,
                      min_frac_retained = DEG_LOO_STABLE, effect_ratio_min = 0.5) {
  sig_of    <- function(r) rownames(r)[deg_reportable(rownames(r)) & deg_is_sig(r$padj, r$log2FoldChange)]
  sig_genes <- sig_of(full)
  if (!length(sig_genes) || !length(folds)) return(NULL)

  as_mat <- function(col) matrix(unlist(lapply(folds, function(r) r[sig_genes, col])),
                                 nrow = length(sig_genes), dimnames = list(sig_genes, names(folds)))
  padj_m <- as_mat("padj"); lfc_m <- as_mat("log2FoldChange"); mle_m <- as_mat("log2FC_MLE")
  sig_m  <- deg_is_sig(padj_m, lfc_m)
  mle0   <- full[sig_genes, "log2FC_MLE"]
  ratio  <- mle_m / mle0
  driven <- !sig_m & !is.na(ratio) & ratio < effect_ratio_min
  safe   <- function(x, f) if (all(is.na(x))) NA_real_ else f(x, na.rm = TRUE)
  who    <- function(m) apply(m, 1, function(z) paste(colnames(m)[z %in% TRUE], collapse = ";"))

  pg <- tibble::tibble(
    gene                 = sig_genes,
    loo_frac_retained    = rowMeans(sig_m),
    loo_sign_flip        = rowSums(sign(mle_m) != sign(mle0), na.rm = TRUE) > 0,
    loo_min_effect_ratio = apply(ratio, 1, safe, min),
    loo_max_shift        = apply(abs(mle_m - mle0), 1, safe, max),
    loo_lost_without     = who(!sig_m),
    loo_driven_by        = who(driven),
    loo_n_folds          = ncol(sig_m),
    loo_n_missing        = as.integer(rowSums(is.na(padj_m))))
  pg$loo_status <- dplyr::case_when(
    n_fail > 0                                ~ sprintf("incomplete (%d/%d folds)", ncol(sig_m), n_total),
    pg$loo_sign_flip                          ~ "sign_flip",
    nzchar(pg$loo_driven_by)                  ~ "donor_driven",
    pg$loo_frac_retained >= min_frac_retained ~ "stable",
    TRUE                                      ~ "power_loss")
  pg
}

loo_gs_default <- function(stat, sets, min_size = 15, max_size = 500, seed = 1) {
  stat <- stat[is.finite(stat) & !is.na(names(stat))]
  stat <- stat[!duplicated(names(stat))]
  fg  <- with_seed(seed, suppressWarnings(fgsea::fgseaMultilevel(sets, stat, minSize = min_size, maxSize = max_size,
                                                                       BPPARAM = BiocParallel::SerialParam())))
  idx <- limma::ids2indices(sets[fg$pathway], names(stat), remove.empty = TRUE)
  cam <- limma::cameraPR(stat, idx, sort = FALSE)
  if (is.null(cam$FDR)) cam$FDR <- cam$PValue
  m <- match(fg$pathway, rownames(cam))
  tibble::tibble(set = fg$pathway, NES = fg$NES, padj_fgsea = fg$padj,
                 dir_camera = cam$Direction[m], padj_camera = cam$FDR[m])
}
loo_gs_sig <- function(g, padj_thresh = DEG_PADJ_THRESH) {
  ok <- g$padj_fgsea < padj_thresh & g$padj_camera < padj_thresh &
    g$dir_camera == ifelse(g$NES > 0, "Up", "Down")
  ok %in% TRUE
}
rank_of <- function(r) stats::setNames(r$stat, rownames(r))

loo_sets <- function(full, folds, sets, gs_fun, n_fail, n_total, min_frac_retained = DEG_LOO_STABLE) {
  g0  <- gs_fun(rank_of(full), sets)
  sig <- g0$set[loo_gs_sig(g0)]
  if (!length(sig) || !length(folds)) return(list(baseline = sig, per_set = NULL))
  nes0 <- g0$NES[match(sig, g0$set)]

  per_fold <- lapply(folds, function(r) {
    g <- gs_fun(rank_of(r), sets); g <- g[match(sig, g$set), , drop = FALSE]
    list(keep = loo_gs_sig(g) & sign(g$NES) == sign(nes0), nes = g$NES)
  })
  mat  <- function(w) matrix(unlist(lapply(per_fold, `[[`, w)), nrow = length(sig),
                             dimnames = list(sig, names(folds)))
  keep <- mat("keep"); keep[is.na(keep)] <- FALSE
  nes  <- mat("nes")
  ps <- tibble::tibble(
    set                  = sig,
    loo_frac_retained    = rowMeans(keep),
    loo_sign_flip        = rowSums(sign(nes) != sign(nes0), na.rm = TRUE) > 0,
    loo_min_effect_ratio = apply(nes / nes0, 1, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)),
    loo_lost_without     = apply(!keep, 1, function(z) paste(colnames(keep)[z], collapse = ";")),
    loo_n_folds          = ncol(keep))
  ps$loo_status <- dplyr::case_when(
    n_fail > 0                                ~ sprintf("incomplete (%d/%d folds)", ncol(keep), n_total),
    ps$loo_sign_flip                          ~ "sign_flip",
    ps$loo_frac_retained >= min_frac_retained ~ "stable",
    TRUE                                      ~ "fragile")
  list(baseline = sig, per_set = ps)
}

loo_unit <- function(so, fit_fun, comparisons, donor_col = "patient_ID",
                     fold_min_donors = max(2L, DEG_MIN_DONORS_PER_GROUP - 1L)) {
  env <- environment(get_res); old <- get("DEG_SHRINK", envir = env)
  assign("DEG_SHRINK", "none", envir = env)
  on.exit(assign("DEG_SHRINK", old, envir = env), add = TRUE)
  full  <- fit_fun(so)
  comps <- intersect(comparisons, names(full)[vapply(full, is.data.frame, logical(1))])
  if (!length(comps)) return(NULL)

  cell_donor <- as.character(so@meta.data[[donor_col]])
  n_cells    <- table(cell_donor)
  donors     <- names(n_cells)[n_cells >= DEG_MIN_CELLS_PER_SAMPLE]

  fold <- function(d) tryCatch(loo_with_min_donors(function()
    fit_fun(subset_cells(so, cell_donor != d, what = paste("drop donor", d))), fold_min_donors),
    error = function(e) conditionMessage(e))
  fold_out <- if (DEG_LOO_CORES > 1 && .Platform$OS.type == "unix")
    parallel::mclapply(donors, fold, mc.cores = DEG_LOO_CORES) else lapply(donors, fold)
  fold_out <- lapply(fold_out, function(r) if (inherits(r, "try-error")) as.character(r) else r)
  names(fold_out) <- donors

  lapply(stats::setNames(comps, comps), function(cmp) {
    f0 <- full[[cmp]]
    status <- vapply(donors, function(d) {
      r <- fold_out[[d]]
      if (is.character(r)) return(paste("failed:", r[1]))
      rc <- r[[cmp]]
      if (!is.data.frame(rc)) return(paste("failed:", if (is.character(rc)) rc else "contrast not returned"))
      same <- identical(rownames(rc), rownames(f0)) && identical(rc$padj, f0$padj)
      if (same) "not_in_fit" else "ok"
    }, character(1))
    ok <- donors[status == "ok"]
    list(full = f0,
         folds = lapply(stats::setNames(ok, ok), function(d) fold_out[[d]][[cmp]]),
         fold_tbl = tibble::tibble(donor = donors, n_cells = as.integer(n_cells[donors]), status = status))
  })
}

loo_fit_groups <- function(comps, runner = NULL, donor_col = "patient_ID") {
  if (!is.null(runner)) return(list(custom = list(fun = runner, comparisons = comps)))
  g <- list()
  if ("xALD_vs_ctrl" %in% comps)
    g$xald <- list(fun = function(so) run_pseudobulk(so, donor_col), comparisons = "xALD_vs_ctrl")
  sev <- intersect(comps, DEG_SEVERITY_NAMES)
  if (length(sev))
    g$severity <- list(fun = function(so) run_pseudobulk_severity(so, donor_col, which = sev), comparisons = sev)
  g
}

deg_loo_run <- function(units, deg_tbl, tag, sets = NULL, gs_full = NULL, runner = NULL,
                        comparisons = if (is.null(runner)) DEG_LOO_CONTRASTS,
                        gs_comparisons = if (is.null(runner)) DEG_LOO_GS_CONTRASTS else comparisons,
                        gs_fun = loo_gs_default, donor_col = "patient_ID") {
  fam <- deg_tbl %>%
    dplyr::filter(status == "ok") %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(LOO_KEYS))) %>%
    dplyr::mutate(comparison = as.character(comparison))
  if (!is.null(comparisons)) fam <- dplyr::filter(fam, comparison %in% comparisons)
  do_sets_for <- function(cmp) !is.null(sets) && (is.null(gs_comparisons) || cmp %in% gs_comparisons)

  out_g <- out_s <- out_f <- out_fam <- list(); t_all <- Sys.time()
  for (u in units) {
    fam_u <- dplyr::filter(fam, tissue == u$tissue, level == u$level, subset == u$subset)
    if (!nrow(fam_u)) next
    groups <- loo_fit_groups(fam_u$comparison, runner, donor_col)
    for (gname in names(groups)) {
      grp   <- groups[[gname]]
      comps <- intersect(grp$comparisons, fam_u$comparison)
      if (!length(comps)) next
      ck <- file.path(DEG_LOO_DIR, paste0(tag, "_v3__", sanitize(paste(u$tissue, u$level, u$subset, gname,
                                                                       paste(comps, collapse = "+"))), ".qs"))
      if (isTRUE(DEG_LOO_RESUME) && file.exists(ck)) {
        o <- qs::qread(ck)
        message(sprintf("[%s] %s | %s | %s | %s: checkpoint", tag, u$tissue, u$level, u$subset, gname))
      } else {
        t0  <- Sys.time()
        res <- tryCatch(loo_unit(u$get(), grp$fun, comps, donor_col), error = function(e) conditionMessage(e))
        secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        o <- list(genes = list(), sets = list(), folds = list(), fam = list())
        for (cmp in comps) {
          keys <- tibble::tibble(tissue = u$tissue, level = u$level, subset = u$subset, comparison = cmp)
          r    <- if (is.list(res)) res[[cmp]] else NULL
          if (is.null(r)) {
            why <- if (is.character(res)) res else "full fit did not return this contrast"
            o$fam[[cmp]] <- dplyr::mutate(keys, status = paste("failed:", why), n_donors_in_fit = NA_integer_,
                                          n_folds_ok = NA_integer_, n_folds_failed = NA_integer_,
                                          n_sig_fit = NA_integer_, n_sig_table = NA_integer_, sig_match = NA,
                                          n_sets_loo = NA_integer_, n_sets_table = NA_integer_,
                                          n_sets_shared = NA_integer_, secs = secs)
            next
          }
          ft      <- r$fold_tbl
          n_total <- sum(ft$status != "not_in_fit")
          n_fail  <- sum(grepl("^failed", ft$status))

          tbl_sig <- deg_tbl %>%
            dplyr::filter(tissue == u$tissue, level == u$level, subset == u$subset, comparison == cmp,
                          status == "ok", deg_reportable(gene), deg_is_sig(padj, log2FC)) %>%
            dplyr::pull(gene)
          fit_sig <- rownames(r$full)[deg_reportable(rownames(r$full)) &
                                        deg_is_sig(r$full$padj, r$full$log2FoldChange)]

          pg <- loo_genes(r$full, r$folds, n_fail, n_total)
          if (!is.null(pg)) o$genes[[cmp]] <- dplyr::bind_cols(keys[rep(1, nrow(pg)), ], pg)
          o$folds[[cmp]] <- dplyr::bind_cols(keys[rep(1, nrow(ft)), ], ft)

          n_sets_loo <- n_sets_table <- n_sets_shared <- NA_integer_
          if (do_sets_for(cmp)) {
            s <- tryCatch(loo_sets(r$full, r$folds, sets, gs_fun, n_fail, n_total),
                          error = function(e) { message("  gene-set LOO: ", conditionMessage(e)); NULL })
            if (!is.null(s)) {
              n_sets_loo <- length(s$baseline)
              if (!is.null(gs_full)) {
                tbl_sets <- gs_full %>%
                  dplyr::filter(tissue == u$tissue, level == u$level, subset == u$subset,
                                comparison == cmp, agree %in% TRUE) %>% dplyr::pull(set)
                n_sets_table  <- length(tbl_sets)
                n_sets_shared <- length(intersect(tbl_sets, s$baseline))
              }
              if (!is.null(s$per_set)) o$sets[[cmp]] <- dplyr::bind_cols(keys[rep(1, nrow(s$per_set)), ], s$per_set)
            }
          }
          o$fam[[cmp]] <- dplyr::mutate(keys, status = "ok", n_donors_in_fit = n_total,
                                        n_folds_ok = length(r$folds), n_folds_failed = n_fail,
                                        n_sig_fit = length(fit_sig), n_sig_table = length(tbl_sig),
                                        sig_match = setequal(tbl_sig, fit_sig), n_sets_loo = n_sets_loo,
                                        n_sets_table = n_sets_table, n_sets_shared = n_sets_shared,
                                        secs = secs)
        }
        qs::qsave(o, ck)
        message(sprintf("[%s] %s | %s | %s | %s: %d contrast(s), %.0f s", tag, u$tissue, u$level, u$subset,
                        gname, length(comps), secs))
      }
      out_g <- c(out_g, o$genes); out_s <- c(out_s, o$sets)
      out_f <- c(out_f, o$folds); out_fam <- c(out_fam, o$fam)
    }
  }
  families <- dplyr::bind_rows(out_fam)
  if (nrow(families))
    message(sprintf("[%s] %d families, %.1f min, %d failed, %d sig mismatches", tag, nrow(families),
                    as.numeric(difftime(Sys.time(), t_all, units = "mins")),
                    sum(families$status != "ok"), sum(families$sig_match %in% FALSE)))
  list(genes = dplyr::bind_rows(out_g), sets = dplyr::bind_rows(out_s),
       folds = dplyr::bind_rows(out_f), families = families)
}

deg_loo_summary <- function(loo_tbl) {
  if (is.null(loo_tbl) || !nrow(loo_tbl)) return(tibble::tibble())
  loo_tbl %>%
    dplyr::mutate(.st = sub(" \\(.*$", "", loo_status)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(LOO_KEYS))) %>%
    dplyr::summarise(n_significant      = dplyr::n(),
                     n_folds            = max(loo_n_folds),
                     mean_frac_retained = mean(loo_frac_retained, na.rm = TRUE),
                     n_below_fragile    = sum(loo_frac_retained < DEG_LOO_FRAGILE, na.rm = TRUE),
                     n_stable           = sum(.st == "stable"),
                     n_power_loss       = sum(.st == "power_loss"),
                     n_fragile          = sum(.st == "fragile"),
                     n_donor_driven     = sum(.st == "donor_driven"),
                     n_sign_flip        = sum(.st == "sign_flip"),
                     n_incomplete       = sum(.st == "incomplete"),
                     .groups = "drop")
}

loo_annotate <- function(x, loo_tbl, key = c("gene", "set")) {
  key <- match.arg(key)
  x <- dplyr::select(x, -dplyr::starts_with("loo_")) %>% dplyr::mutate(comparison = as.character(comparison))
  if (!is.null(loo_tbl) && nrow(loo_tbl))
    x <- dplyr::left_join(x, dplyr::mutate(loo_tbl, comparison = as.character(comparison)),
                          by = c(LOO_KEYS, key))
  if (!"loo_status" %in% names(x)) x$loo_status <- NA_character_
  sig  <- if (key == "gene") deg_reportable(x$gene, x$status) & deg_is_sig(x$padj, x$log2FC) else x$agree %in% TRUE
  miss <- sig & is.na(x$loo_status)
  loo_cmp <- if (key == "gene") DEG_LOO_CONTRASTS else DEG_LOO_GS_CONTRASTS
  x$loo_status[miss] <- dplyr::case_when(
    !x$tissue[miss] %in% DEG_LOO_TISSUES ~ "not run: tissue not in DEG_LOO_TISSUES",
    x$comparison[miss] %in% deg_active_contrasts() & !x$comparison[miss] %in% loo_cmp
                                         ~ "not run: contrast not in DEG_LOO_CONTRASTS",
    key == "set"                         ~ "not run: not significant under the LOO set test",
    TRUE                                 ~ "not run: no refit for this family")
  x
}

DEG_DONOR_HEATMAP_MAX <- 100
DEG_DONOR_TOP_N       <- 10
DEG_DONOR_LFC_LIMIT   <- 2
DEG_DONOR_SUBTITLE    <- "log2 vs control mean"
DEG_STAGE_GENES_FILE  <- here::here("metadata", "deg_stage_genes.csv")

deg_stage_genes <- function(compartment, file = DEG_STAGE_GENES_FILE) {
  if (!file.exists(file)) return(NULL)
  g <- readr::read_csv(file, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  g[g$compartment == compartment, c("contrast", "tissue", "level", "subset", "gene")]
}

deg_donor_lfc <- function(so, genes, donor_col = "patient_ID", min_cells = DEG_MIN_CELLS_PER_SAMPLE) {
  pb  <- make_pseudobulk(so, donor_col, min_cells)
  dds <- DESeq2::estimateSizeFactors(DESeq2::DESeqDataSetFromMatrix(pb$counts, pb$meta, ~ 1))
  lg  <- log2(DESeq2::counts(dds, normalized = TRUE)[intersect(genes, rownames(dds)), , drop = FALSE] + 1)
  sev <- CLIN_SEVERITY$severity[match(colnames(lg), CLIN_SEVERITY$patient_ID)]
  ctl <- which(sev %in% 0)
  if (length(ctl) < 2) stop("fewer than 2 control donors", call. = FALSE)
  rel <- lg - rowMeans(lg[, ctl, drop = FALSE])
  tibble::tibble(gene = rep(rownames(rel), ncol(rel)), donor = rep(colnames(rel), each = nrow(rel)),
                 severity = rep(sev, each = nrow(rel)), lfc = as.vector(rel)) %>%
    dplyr::filter(!is.na(severity)) %>%
    dplyr::mutate(stage = factor(paste0("sev", severity), levels = names(SEVERITY_COLS)),
                  group = sev_group3(severity))
}

plot_deg_donor_heatmap <- function(d, sig, title, subtitle = NULL) {
  d <- dplyr::filter(d, gene %in% sig$gene)
  donors <- dplyr::distinct(d, donor, severity) %>% dplyr::arrange(severity, donor)
  d$gene  <- factor(d$gene, levels = rev(sig$gene[order(-sig$log2FC)]))
  d$donor <- factor(d$donor, levels = donors$donor)
  d$group <- factor(SEV_GROUP3_LABELS[as.character(d$group)], levels = SEV_GROUP3_LABELS)
  lim <- DEG_DONOR_LFC_LIMIT
  ggplot2::ggplot(d, ggplot2::aes(donor, gene, fill = pmax(pmin(lfc, lim), -lim))) +
    ggplot2::geom_tile() +
    ggplot2::facet_grid(~ group, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-lim, lim),
                                  name = "log2 vs\ncontrol mean") +
    ggplot2::scale_x_discrete(labels = function(x) sub("^BIO-", "", x)) +
    ggplot2::labs(x = "donor", y = NULL, title = title, subtitle = subtitle) +
    gs_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = FIG_MIN_PT, angle = 90, vjust = 0.5, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = if (nrow(sig) > 60) FIG_MIN_PT - 1 else FIG_MIN_PT),
                   axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
}

plot_deg_top_severity <- function(d, genes, title, subtitle = NULL) {
  d <- dplyr::filter(d, gene %in% genes) %>% dplyr::mutate(gene = factor(gene, levels = genes))
  ggplot2::ggplot(d, ggplot2::aes(group, lfc)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = FIG_LINE_MM) +
    ggplot2::stat_summary(fun = stats::median, geom = "crossbar", width = 0.5, linewidth = FIG_LINE_MM,
                          colour = "grey30") +
    ggplot2::geom_jitter(ggplot2::aes(colour = stage), width = 0.12, height = 0, size = 1.1) +
    ggplot2::facet_wrap(~ gene, nrow = 2, scales = "free_y") +
    ggplot2::scale_colour_manual(values = SEVERITY_COLS, labels = SEVERITY_SHORT[sub("sev", "", names(SEVERITY_COLS))],
                                 name = "stage", drop = FALSE) +
    ggplot2::scale_x_discrete(labels = SEV_GROUP3_LABELS) +
    ggplot2::labs(x = NULL, y = "log2 vs control mean", title = title, subtitle = subtitle) +
    gs_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

deg_donor_views <- function(units, deg_tbl, out_dir, tag, contrasts = deg_report_contrasts(),
                            stage_genes = NULL, print_subsets = "all", donor_col = "patient_ID") {
  sig_all <- deg_tbl %>%
    dplyr::filter(comparison %in% contrasts, deg_reportable(gene, status), deg_is_sig(padj)) %>%
    dplyr::mutate(highlight = deg_is_highlight(padj, log2FC, svalue = svalue))
  out <- list()
  for (u in units) {
    fam <- dplyr::filter(sig_all, tissue == u$tissue, level == u$level, subset == u$subset)
    sel <- if (is.null(stage_genes)) stage_genes[0, ] else
      dplyr::filter(stage_genes, tissue == u$tissue, level == u$level, subset == u$subset)
    if (!nrow(fam) && !nrow(sel)) next
    d <- tryCatch(deg_donor_lfc(u$get(), unique(c(fam$gene, sel$gene)), donor_col), error = function(e) {
      message("deg_donor_views ", u$tissue, " | ", u$subset, ": ", conditionMessage(e)); NULL })
    if (is.null(d)) next
    out[[length(out) + 1]] <- dplyr::mutate(d, tissue = u$tissue, level = u$level, subset = u$subset, .before = 1)
    stem0 <- paste0("deg_donor_", tag, "_", u$tissue, "_", sanitize(u$subset))
    for (cn in unique(sel$contrast)) {
      g <- unique(sel$gene[sel$contrast == cn])
      miss <- setdiff(g, d$gene)
      if (length(miss)) message("deg_donor_views ", u$tissue, " | ", u$subset, ": not tested: ", toString(miss))
      g <- intersect(g, d$gene)
      if (!length(g)) next
      p2 <- plot_deg_top_severity(d, g, paste(u$tissue, "-", u$subset, "-", deg_contrast_labels[[cn]] %||% cn),
                                  DEG_DONOR_SUBTITLE)
      save_umap_fig(p2, paste0(stem0, "_", cn, "_selected_by_stage"), out_dir, width_mm = FIG_WIDTH_MM[["double"]],
                    height_mm = if (length(g) > 5) 95 else 55)
      if (u$subset %in% print_subsets) print(p2)
    }
    for (cn in intersect(contrasts, unique(fam$comparison))) {
      s <- dplyr::filter(fam, comparison == cn) %>% dplyr::arrange(padj)
      n_sig <- nrow(s); s <- utils::head(s, DEG_DONOR_HEATMAP_MAX)
      lab  <- paste(u$tissue, "-", u$subset, "-", deg_contrast_labels[[cn]] %||% cn)
      if (n_sig > nrow(s)) message("deg_donor_views ", lab, ": heatmap shows the top ", nrow(s), " of ", n_sig)
      stem <- paste0(stem0, "_", cn)
      p1 <- plot_deg_donor_heatmap(d, s, lab, DEG_DONOR_SUBTITLE)
      save_umap_fig(p1, paste0(stem, "_heatmap"), out_dir, width_mm = FIG_WIDTH_MM[["double"]],
                    height_mm = min(FIG_MAX_DEPTH_MM, round(35 + 2.1 * nrow(s))))
      if (u$subset %in% print_subsets) print(p1)
      if (!is.null(stage_genes)) next
      top <- utils::head(dplyr::arrange(s, dplyr::desc(highlight), padj)$gene, DEG_DONOR_TOP_N)
      p2 <- plot_deg_top_severity(d, top, paste0(lab, ": top ", length(top)), DEG_DONOR_SUBTITLE)
      save_umap_fig(p2, paste0(stem, "_top_by_stage"), out_dir, width_mm = FIG_WIDTH_MM[["double"]],
                    height_mm = 95)
      if (u$subset %in% print_subsets) print(p2)
    }
  }
  dplyr::bind_rows(out)
}
