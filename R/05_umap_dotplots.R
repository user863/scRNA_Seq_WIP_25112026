
FIG_WIDTH_MM     <- c(single = 89, onehalf = 120, onehalf_wide = 136, double = 183)
FIG_MAX_DEPTH_MM <- 247
FIG_BASE_PT      <- 7
FIG_MIN_PT       <- 5
FIG_LEGEND_PT    <- 6
FIG_LEGEND_WRAP  <- 24
FIG_PANEL_PT     <- 8
FIG_PANEL_CASE   <- "lower"
FIG_DPI          <- 450
FIG_LINE_MM      <- 0.18

FIG_PALETTE   <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                   "#56B4E9", "#E69F00", "#332288",
                   "#882255", "#999933", "#F0E442")
FIG_OTHER_COL <- "grey70"

FIG_UMAP_LABEL_WRAP   <- 16
FIG_UMAP_LABEL_PT     <- 5
FIG_UMAP_LABEL_COLOUR <- "darken"
FIG_UMAP_LABEL_GAP    <- 0.06
FIG_UMAP_PANEL_MM     <- 48
FIG_UMAP_LABEL_LAYOUT <- "auto"
FIG_UMAP_SIDES_ASP    <- 1.2
FIG_UMAP_LABEL_STYLE <- "box"
FIG_UMAP_BOX_PT      <- 5.5
FIG_UMAP_BOX_WRAP    <- 22

umap_label_layout <- function(asp, layout = FIG_UMAP_LABEL_LAYOUT) {
  if (!identical(layout, "auto")) return(layout)
  if (asp > FIG_UMAP_SIDES_ASP) "sides" else "bands"
}
umap_label_expand <- function(layout, style = FIG_UMAP_LABEL_STYLE) {
  if (identical(style, "box")) return(c(x = 0.05, y = 0.05))
  if (identical(layout, "bands")) c(x = 0.06, y = 0.45) else c(x = 0.55, y = 0.08)
}
umap_box_text_colour <- function(fill) {
  m <- grDevices::col2rgb(fill) / 255
  lum <- 0.2126 * m[1, ] + 0.7152 * m[2, ] + 0.0722 * m[3, ]
  ifelse(lum > 0.62, "grey10", "white")
}

umap_label_colour <- function(cols, mode = FIG_UMAP_LABEL_COLOUR, f = 0.62) {
  switch(mode,
    black      = rep("black", length(cols)),
    population = unname(cols),
    darken     = { m <- grDevices::col2rgb(cols) * f
                   grDevices::rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255) },
    stop("label colour mode must be darken, population or black", call. = FALSE))
}
FIG_OTHER_LAB <- "Other"

fig_font <- function(preferred = c("Helvetica", "Arial"), quiet = FALSE) {
  fams <- tryCatch(unique(systemfonts::system_fonts()$family),
                   error = function(e) character(0))
  hit <- preferred[preferred %in% fams]
  if (length(hit)) return(hit[1])
  if (!quiet)
    message("fig_font(): neither Helvetica nor Arial available -- using the ",
            "device default.")
  ""
}

umap_clean_labels <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA") | grepl("^[0-9]+$", x)] <- NA
  x
}

collapse_labels <- function(x, keep = NULL, max_n = length(FIG_PALETTE),
                            other = FIG_OTHER_LAB, quiet = FALSE) {
  x  <- as.character(x)
  tb <- sort(table(x[!is.na(x)]), decreasing = TRUE)
  keep <- if (is.null(keep)) {
    if (length(tb) <= max_n) names(tb) else names(tb)[seq_len(max_n - 1L)]
  } else intersect(keep, names(tb))
  merged <- setdiff(names(tb), keep)
  if (length(merged) && !quiet)
    message("collapse_labels(): ", length(merged), " label(s) -> '", other, "': ",
            toString(sprintf("%s (n = %s)", merged, format(tb[merged], big.mark = ","))))
  out <- ifelse(is.na(x), NA_character_, ifelse(x %in% keep, x, other))
  factor(out, levels = c(keep, if (length(merged)) other))
}

fig_palette <- function(levels) {
  levels <- as.character(levels); levels <- levels[!is.na(levels)]
  main <- setdiff(levels, FIG_OTHER_LAB)
  if (length(main) > length(FIG_PALETTE))
    stop("fig_palette(): ", length(main), " categories but only ",
         length(FIG_PALETTE), " colours available. Run collapse_labels() first.",
         call. = FALSE)
  p <- setNames(FIG_PALETTE[seq_along(main)], main)
  if (FIG_OTHER_LAB %in% levels) p <- c(p, setNames(FIG_OTHER_COL, FIG_OTHER_LAB))
  p[levels]
}

umap_label_shares <- function(objs, label_col) {
  fr <- lapply(objs, function(so)
    prop.table(table(umap_clean_labels(so[[label_col, drop = TRUE]]))))
  labs  <- unique(unlist(lapply(fr, names)))
  score <- vapply(labs, function(l)
    max(vapply(fr, function(f) if (l %in% names(f)) f[[l]] else 0, numeric(1))),
    numeric(1))
  list(per_panel = fr, ranked = names(sort(score, decreasing = TRUE)))
}

umap_panel_palettes <- function(objs, label_col, mode = c("auto", "shared", "per_panel"),
                                keep = NULL, n_col = length(FIG_PALETTE)) {
  mode <- match.arg(mode)
  sh   <- umap_label_shares(objs, label_col)
  if (!is.null(keep)) {
    bad <- setdiff(keep, sh$ranked)
    if (length(bad)) stop("keep label(s) not present: ", toString(bad), call. = FALSE)
    mode <- "shared"
  }
  if (mode == "auto") mode <- if (length(sh$ranked) <= n_col) "shared" else "per_panel"

  if (mode == "shared") {
    if (is.null(keep))
      keep <- if (length(sh$ranked) <= n_col) sh$ranked else sh$ranked[seq_len(n_col - 1L)]
    merged <- setdiff(sh$ranked, keep)
    pal <- fig_palette(c(keep, if (length(merged)) FIG_OTHER_LAB))
    return(list(mode = mode, shared_legend = TRUE,
                palettes = setNames(rep(list(pal), length(objs)), names(objs)),
                keep = setNames(rep(list(keep), length(objs)), names(objs)),
                merged = setNames(rep(list(merged), length(objs)), names(objs))))
  }

  n_in    <- table(unlist(lapply(sh$per_panel, names)))
  common  <- intersect(sh$ranked, names(n_in)[n_in > 1])
  common  <- utils::head(common, n_col)
  out <- lapply(names(objs), function(nm) {
    f    <- sh$per_panel[[nm]]
    own  <- setdiff(names(sort(f, decreasing = TRUE)), common)
    here <- intersect(common, names(f))
    free <- setdiff(FIG_PALETTE, FIG_PALETTE[match(here, common)])
    n_own <- length(free)
    keep_own <- if (length(own) <= n_own) own else own[seq_len(max(0L, n_own - 1L))]
    merged   <- setdiff(own, keep_own)
    pal <- c(setNames(FIG_PALETTE[match(here, common)], here),
             setNames(free[seq_along(keep_own)], keep_own))
    if (length(merged)) pal <- c(pal, setNames(FIG_OTHER_COL, FIG_OTHER_LAB))
    list(pal = pal, keep = c(here, keep_own), merged = merged)
  })
  names(out) <- names(objs)
  list(mode = mode, shared_legend = FALSE,
       palettes = lapply(out, `[[`, "pal"),
       keep     = lapply(out, `[[`, "keep"),
       merged   = lapply(out, `[[`, "merged"))
}

theme_umap_panel <- function(base_size = FIG_BASE_PT, family = fig_font(quiet = TRUE)) {
  new_gg <- utils::packageVersion("ggplot2") >= "3.5.0"
  th <- ggplot2::theme_classic(base_size = base_size, base_family = family) +
    ggplot2::theme(
      axis.text        = ggplot2::element_blank(),
      axis.ticks       = ggplot2::element_blank(),
      axis.line        = ggplot2::element_blank(),
      axis.title       = ggplot2::element_blank(),
      panel.grid       = ggplot2::element_blank(),
      legend.position  = "right",
      legend.justification = c(0, 0.5),
      legend.key.size  = grid::unit(3, "mm"),
      legend.text      = ggplot2::element_text(size = FIG_LEGEND_PT,
                                               margin = ggplot2::margin(l = 1, unit = "mm")),
      legend.title     = ggplot2::element_blank(),
      legend.margin    = ggplot2::margin(0, 0, 0, 0),
      plot.title       = ggplot2::element_text(size = base_size, hjust = 0),
      plot.tag         = ggplot2::element_text(size = FIG_PANEL_PT, face = "bold"),
      plot.tag.position = if (new_gg) "topleft" else c(0, 1),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin      = ggplot2::margin(2, 2, 2, 2, "mm"))
  if (new_gg) th + ggplot2::theme(plot.tag.location = "margin",
                                  legend.key.spacing.y = grid::unit(0.6, "mm")) else th
}

umap_axis_inset <- function(df, frac = 0.18, labs = c("UMAP 1", "UMAP 2"),
                            size_pt = FIG_MIN_PT, family = fig_font(quiet = TRUE)) {
  rx <- range(df$.x, na.rm = TRUE); ry <- range(df$.y, na.rm = TRUE)
  x0 <- rx[1] - diff(rx) * 0.04; y0 <- ry[1] - diff(ry) * 0.04
  lx <- diff(rx) * frac;         ly <- diff(ry) * frac
  arr <- grid::arrow(length = grid::unit(1.1, "mm"), type = "closed")
  list(
    ggplot2::annotate("segment", x = x0, xend = x0 + lx, y = y0, yend = y0,
                      arrow = arr, linewidth = FIG_LINE_MM),
    ggplot2::annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + ly,
                      arrow = arr, linewidth = FIG_LINE_MM),
    ggplot2::annotate("text", x = x0 + lx / 2, y = y0 - diff(ry) * 0.035,
                      label = labs[1], size = size_pt / .pt, family = family,
                      vjust = 1, colour = "black"),
    ggplot2::annotate("text", x = x0 - diff(rx) * 0.035, y = y0 + ly / 2,
                      label = labs[2], size = size_pt / .pt, family = family,
                      angle = 90, vjust = 0, colour = "black"))
}

plot_umap_panel <- function(so, label_col, reduction = "umap_harmony",
                            title = NULL, panel = NULL,
                            keep_labels = NULL, palette = NULL,
                            legend_levels = NULL,
                            max_colours = length(FIG_PALETTE),
                            pt_size = 0.18, pt_alpha = 0.9,
                            label_points = TRUE, legend = FALSE,
                            label_colour = FIG_UMAP_LABEL_COLOUR, obstacle_grid = 80,
                            label_layout = FIG_UMAP_LABEL_LAYOUT,
                            label_style = FIG_UMAP_LABEL_STYLE,
                            raster = TRUE, family = fig_font(quiet = TRUE)) {
  emb <- as.data.frame(Embeddings(so, reduction))[, 1:2]
  names(emb) <- c(".x", ".y")
  emb$.lab <- collapse_labels(umap_clean_labels(so[[label_col, drop = TRUE]]),
                              keep = keep_labels, max_n = max_colours)
  emb <- emb[!is.na(emb$.lab), , drop = FALSE]
  rm(so)
  set.seed(1)
  emb <- emb[sample(nrow(emb)), ]
  pal <- if (is.null(palette)) fig_palette(levels(emb$.lab)) else palette
  miss <- setdiff(levels(emb$.lab), names(pal))
  if (length(miss)) stop("plot_umap_panel(): palette lacks ", toString(miss), call. = FALSE)
  lv <- if (is.null(legend_levels)) intersect(names(pal), levels(emb$.lab)) else legend_levels
  wrap_n <- FIG_LEGEND_WRAP
  on.exit(rm(list = setdiff(ls(), c("wrap_n", "pt_size", "pt_alpha"))), add = TRUE)

  pts <- ggplot2::geom_point(size = pt_size, alpha = pt_alpha, shape = 16,
                             show.legend = TRUE)
  if (raster && requireNamespace("ggrastr", quietly = TRUE)) {
    dev <- if (requireNamespace("ragg", quietly = TRUE)) "ragg" else "cairo"
    pts <- tryCatch(ggrastr::rasterise(pts, dpi = FIG_DPI, dev = dev),
                    error = function(e) {
                      message("plot_umap_panel(): rasterise failed (", dev, ": ",
                              conditionMessage(e), ") -- points stay vector.")
                      ggplot2::geom_point(size = pt_size, alpha = pt_alpha, shape = 16,
                             show.legend = TRUE)
                    })
  } else if (raster) {
    message("plot_umap_panel(): ggrastr not installed -- points stay vector, ",
            "the PDF will be large. install.packages('ggrastr')")
  }

  p <- ggplot2::ggplot(emb, ggplot2::aes(.x, .y, colour = .lab)) +
    pts +
    ggplot2::scale_colour_manual(
      values = pal[lv], limits = lv, breaks = lv, drop = FALSE,
      labels = function(x) stringr::str_wrap(x, wrap_n)) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      ncol = 1, override.aes = list(size = 2, alpha = 1))) +
    ggplot2::coord_equal(clip = "off") +
    umap_axis_inset(emb, family = family) +
    theme_umap_panel(family = family)
  if (!legend) p <- p + ggplot2::theme(legend.position = "none")

  if (label_points && identical(label_style, "box")) {
    anchor <- do.call(rbind, lapply(split(emb, emb$.lab, drop = TRUE), function(g) {
      i <- which.min((g$.x - stats::median(g$.x))^2 + (g$.y - stats::median(g$.y))^2)
      data.frame(.x = g$.x[i], .y = g$.y[i], .lab = as.character(g$.lab[i]))
    }))
    anchor$.txt <- stringr::str_wrap(anchor$.lab, FIG_UMAP_BOX_WRAP)
    rx <- range(emb$.x); ry <- range(emb$.y); sx <- diff(rx); sy <- diff(ry)
    ex <- umap_label_expand(NULL, "box")
    x0 <- rx[1] - sx * 0.04; y0 <- ry[1] - sy * 0.04
    s  <- seq(0, 0.18, length.out = 10)
    obs <- data.frame(
      .x = c(x0 + s * sx, rep(x0, 10), x0 + s * sx, rep(x0 - sx * 0.05, 10)),
      .y = c(rep(y0, 10), y0 + s * sy, rep(y0 - sy * 0.05, 10), y0 + s * sy),
      .txt = "")
    lab_df <- rbind(anchor[, c(".x", ".y", ".txt")], obs)
    fill   <- unname(pal[anchor$.lab])
    p <- p + ggrepel::geom_label_repel(
      data = lab_df, inherit.aes = FALSE,
      mapping = ggplot2::aes(x = .x, y = .y, label = .txt),
      fill = c(fill, rep(NA, nrow(obs))),
      colour = c(umap_box_text_colour(fill), rep(NA, nrow(obs))),
      fontface = "bold", size = FIG_UMAP_BOX_PT / .pt, family = family, lineheight = 0.85,
      label.r = grid::unit(0.6, "mm"), label.padding = grid::unit(0.45, "mm"), label.size = 0,
      box.padding = grid::unit(0.6, "mm"), point.padding = 0, point.size = 0.5,
      force = 2, force_pull = 1.5, max.overlaps = Inf, max.time = 5, max.iter = 1e5,
      min.segment.length = grid::unit(1.5, "mm"), segment.size = FIG_LINE_MM,
      segment.colour = "grey30", seed = 1) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = ex[["x"]])) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = ex[["y"]]))
  } else if (label_points) {
    anchor <- do.call(rbind, lapply(split(emb, emb$.lab, drop = TRUE), function(g) {
      i <- which.min((g$.x - stats::median(g$.x))^2 + (g$.y - stats::median(g$.y))^2)
      data.frame(.x = g$.x[i], .y = g$.y[i], .lab = as.character(g$.lab[i]))
    }))
    anchor$.txt <- stringr::str_wrap(anchor$.lab, FIG_UMAP_LABEL_WRAP)

    rx <- range(emb$.x); ry <- range(emb$.y)
    sx <- diff(rx); sy <- diff(ry)
    lay <- umap_label_layout(sy / sx, label_layout)
    ex  <- umap_label_expand(lay, "callout")
    edge_at <- function(v0, side, along_x) {
      w <- if (along_x) abs(emb$.x - v0) < 0.06 * sx else abs(emb$.y - v0) < 0.06 * sy
      v <- if (along_x) emb$.y else emb$.x
      v <- if (any(w)) v[w] else v
      stats::quantile(v, if (side > 0) 0.995 else 0.005, names = FALSE)
    }
    if (identical(lay, "bands")) {
      half_w <- vapply(strsplit(anchor$.txt, "\n", fixed = TRUE), function(l)
        max(nchar(l)), numeric(1)) * 1.1 * FIG_UMAP_LABEL_PT / 6 / 2 / FIG_UMAP_PANEL_MM * sx * (1 + 2 * ex[["x"]])
      lo <- rx[1] - ex[["x"]] * sx; hi <- rx[2] + ex[["x"]] * sx
      target <- t(vapply(seq_len(nrow(anchor)), function(i) {
        x0  <- min(max(anchor$.x[i], lo + half_w[i]), hi - half_w[i])
        top <- edge_at(x0, 1, TRUE); bot <- edge_at(x0, -1, TRUE)
        if (top - anchor$.y[i] <= anchor$.y[i] - bot) c(x0, top + FIG_UMAP_LABEL_GAP * sy)
        else                                        c(x0, bot - FIG_UMAP_LABEL_GAP * sy)
      }, numeric(2)))
      hj <- rep(0.5, nrow(anchor))
    } else {
      target <- t(vapply(seq_len(nrow(anchor)), function(i) {
        y0 <- min(max(anchor$.y[i], ry[1] + 0.05 * sy), ry[2] - 0.05 * sy)
        rgt <- edge_at(y0, 1, FALSE); lft <- edge_at(y0, -1, FALSE)
        if (rgt - anchor$.x[i] <= anchor$.x[i] - lft) c(rgt + FIG_UMAP_LABEL_GAP * sx, y0)
        else                                        c(lft - FIG_UMAP_LABEL_GAP * sx, y0)
      }, numeric(2)))
      hj <- ifelse(target[, 1] > anchor$.x, 0, 1)
    }

    gx <- pmin(floor((emb$.x - rx[1]) / sx * obstacle_grid), obstacle_grid - 1)
    gy <- pmin(floor((emb$.y - ry[1]) / sy * obstacle_grid), obstacle_grid - 1)
    key <- gx * obstacle_grid + gy
    obs <- data.frame(.x = as.numeric(tapply(emb$.x, key, mean)),
                      .y = as.numeric(tapply(emb$.y, key, mean)))
    x0 <- rx[1] - sx * 0.04; y0 <- ry[1] - sy * 0.04
    s  <- seq(0, 0.18, length.out = 10)
    obs <- rbind(obs, data.frame(
      .x = c(x0 + s * sx, rep(x0, 10), x0 + s * sx, rep(x0 - sx * 0.05, 10)),
      .y = c(rep(y0, 10), y0 + s * sy, rep(y0 - sy * 0.05, 10), y0 + s * sy)))
    obs$.lab <- NA_character_; obs$.txt <- ""

    lab_df  <- rbind(anchor, obs)
    txt_col <- c(umap_label_colour(pal[anchor$.lab], label_colour),
                 rep("black", nrow(obs)))
    p <- p + ggrepel::geom_text_repel(
      data = lab_df, inherit.aes = FALSE, bg.color = "white", bg.r = 0.15,
      hjust = c(hj, rep(0.5, nrow(obs))),
      mapping = ggplot2::aes(x = .x, y = .y, label = .txt),
      nudge_x = c(target[, 1] - anchor$.x, rep(0, nrow(obs))),
      nudge_y = c(target[, 2] - anchor$.y, rep(0, nrow(obs))),
      colour = txt_col, size = FIG_UMAP_LABEL_PT / .pt, family = family,
      lineheight = 0.85, point.size = 0.6, point.padding = 0.1, box.padding = 0.2,
      force = 3, force_pull = 0.2, max.overlaps = Inf, max.time = 10, max.iter = 1e6,
      xlim = c(NA, NA), ylim = c(NA, NA),
      min.segment.length = 0, segment.size = FIG_LINE_MM, segment.colour = "grey50",
      seed = 1) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = ex[["x"]])) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = ex[["y"]]))
  }
  if (!is.null(title)) p <- p + ggplot2::labs(title = title)
  if (!is.null(panel))
    p <- p + ggplot2::labs(tag = if (identical(FIG_PANEL_CASE, "upper"))
      toupper(panel) else tolower(panel))
  p
}

save_umap_fig <- function(p, stem, dir, width_mm, height_mm = NULL) {
  if (is.null(height_mm)) height_mm <- round(width_mm * 0.85)
  if (height_mm > FIG_MAX_DEPTH_MM)
    stop("save_umap_fig(): ", height_mm, " mm exceeds the ",
         FIG_MAX_DEPTH_MM, " mm maximum depth.", call. = FALSE)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  w <- width_mm / 25.4; h <- height_mm / 25.4

  pdf_dev <- grDevices::pdf
  if (capabilities("cairo")) {
    ok <- tryCatch({ f <- tempfile(fileext = ".pdf")
                     grDevices::cairo_pdf(f); grDevices::dev.off()
                     unlink(f); TRUE }, error = function(e) FALSE,
                   warning = function(w) FALSE)
    if (ok) pdf_dev <- grDevices::cairo_pdf
    else message("  cairo_pdf unavailable -- using pdf(). Check font embedding.")
  }
  ggplot2::ggsave(file.path(dir, paste0(stem, ".pdf")), p,
                  width = w, height = h, units = "in", device = pdf_dev)
  png_dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else NULL
  ggplot2::ggsave(file.path(dir, paste0(stem, ".png")), p,
                  width = w, height = h, units = "in", dpi = FIG_DPI,
                  device = png_dev)
  message("  saved ", stem, ".pdf/.png  (", width_mm, " x ", height_mm, " mm)")
  invisible(file.path(dir, paste0(stem, c(".pdf", ".png"))))
}

umap_legend_mm <- function(labels) {
  w <- stringr::str_wrap(labels, FIG_LEGEND_WRAP)
  lines <- strsplit(w, "\n", fixed = TRUE)
  list(width  = 1.15 * max(nchar(unlist(lines))) + 7,
       height = sum(lengths(lines)) * 2.6 + length(labels) * 1.2 + 6)
}

umap_compartment_figure <- function(objs, label_col, stem, level_title,
                                    reduction = "umap_harmony",
                                    palette_mode = "auto", keep = NULL,
                                    label_mode = c("plot", "legend"),
                                    label_style = FIG_UMAP_LABEL_STYLE,
                                    dir = path$umap_figures,
                                    width_mm = FIG_WIDTH_MM[["double"]]) {
  objs <- objs[vapply(objs, function(so) label_col %in% colnames(so@meta.data),
                      logical(1))]
  if (!length(objs)) { message(stem, ": '", label_col, "' absent -- skipped"); return(invisible(NULL)) }

  label_mode <- match.arg(label_mode)
  on_plot <- label_mode == "plot"
  pp <- umap_panel_palettes(objs, label_col, mode = palette_mode, keep = keep)
  message(stem, ": palette mode '", pp$mode, "'")
  for (nm in names(objs)) if (length(pp$merged[[nm]]))
    message("  ", nm, ": merged into '", FIG_OTHER_LAB, "': ", toString(pp$merged[[nm]]))
  if (pp$mode == "per_panel") {
    shared_lab <- Reduce(intersect, lapply(pp$palettes, names))
    message("  same colour in every panel: ", toString(setdiff(shared_lab, FIG_OTHER_LAB)),
            "; remaining colours are panel-specific")
  }

  all_levels <- if (pp$shared_legend && !on_plot) names(pp$palettes[[1]]) else NULL
  panels <- lapply(seq_along(objs), function(i) {
    so <- objs[[i]]; comp <- names(objs)[i]
    n  <- sum(!is.na(umap_clean_labels(so[[label_col, drop = TRUE]])))
    plot_umap_panel(so, label_col, reduction = reduction,
                    keep_labels = pp$keep[[comp]], palette = pp$palettes[[comp]],
                    legend_levels = all_levels,
                    label_points = on_plot, legend = !on_plot, label_style = label_style,
                    title = sprintf("%s, %s (n = %s)", comp, level_title,
                                    format(n, big.mark = ",")),
                    panel = letters[i])
  })
  names(panels) <- names(objs)

  asp <- max(vapply(objs, function(so) {
    e <- Embeddings(so, reduction)[, 1:2]; diff(range(e[, 2])) / diff(range(e[, 1]))
  }, numeric(1)))
  lg  <- lapply(pp$palettes, function(p)
    if (on_plot) list(width = 0, height = 0) else umap_legend_mm(names(p)))
  lg_w <- max(vapply(lg, `[[`, numeric(1), "width"))
  lg_h <- max(vapply(lg, `[[`, numeric(1), "height"))
  n_p  <- length(panels)
  panel_w <- if (pp$shared_legend) (width_mm - lg_w) / n_p else width_mm / n_p - lg_w
  h <- min(FIG_MAX_DEPTH_MM, round(max((panel_w - 6) * asp + 14, lg_h + 12)))

  fig <- patchwork::wrap_plots(panels, nrow = 1)
  if (pp$shared_legend && !on_plot) fig <- fig + patchwork::plot_layout(guides = "collect")
  save_umap_fig(fig, stem, dir, width_mm = width_mm, height_mm = h)

  w1 <- FIG_WIDTH_MM[["onehalf"]]
  purrr::iwalk(panels, function(p, comp) {
    pw <- w1 - lg[[comp]]$width
    save_umap_fig(p, paste0(stem, "_", sanitize(comp)), dir, width_mm = w1,
                  height_mm = min(FIG_MAX_DEPTH_MM,
                                  round(max((pw - 6) * asp + 14, lg[[comp]]$height + 12))))
  })

  raw_asp <- vapply(objs, function(so) {
    e <- Embeddings(so, reduction)[, 1:2]; diff(range(e[, 2])) / diff(range(e[, 1]))
  }, numeric(1))
  lay_by <- vapply(raw_asp, function(a) if (on_plot) umap_label_layout(a) else "legend",
                   character(1))
  asp_by <- vapply(seq_along(raw_asp), function(i) {
    if (!on_plot) return(raw_asp[[i]])
    ex <- umap_label_expand(lay_by[[i]], label_style)
    raw_asp[[i]] * (1 + 2 * ex[["y"]]) / (1 + 2 * ex[["x"]])
  }, numeric(1))
  names(asp_by) <- names(objs)
  pal_tab <- purrr::imap_dfr(pp$palettes, function(p, comp)
    tibble::tibble(panel = comp, label = names(p), colour = unname(p)))
  readr::write_csv(pal_tab, file.path(dir, paste0(stem, "_palette.csv")))
  invisible(list(figure = fig, panels = panels, palettes = pp$palettes,
                 mode = pp$mode, merged = pp$merged, asp = asp_by,
                 label_layout = lay_by))
}

FIG_EXPR_COLS <- viridisLite::viridis(7)

umap_palette_for <- function(fig, comp, labels) {
  labels <- unique(as.character(labels[!is.na(labels)]))
  pal <- if (!is.null(fig) && !is.null(fig$palettes[[comp]])) fig$palettes[[comp]] else NULL
  if (is.null(pal)) return(make_celltype_palette(labels))
  miss <- setdiff(labels, names(pal))
  if (length(miss)) {
    free <- setdiff(FIG_PALETTE, pal)
    add  <- if (length(free) >= length(miss)) free[seq_along(miss)]
            else scales::hue_pal()(length(miss))
    message("umap_palette_for(", comp, "): not in the UMAP palette, new colours: ",
            toString(miss))
    pal <- c(pal, setNames(add, miss))
  }
  pal[intersect(names(pal), labels)]
}

plot_dot_panel <- function(so, features, label_col, palette = NULL, levels = NULL,
                           extra = NULL, title = NULL, panel = NULL, cap = 2.5,
                           colour_key = FALSE, family = fig_font(quiet = TRUE)) {
  lab <- umap_clean_labels(so[[label_col, drop = TRUE]])
  so  <- so[, !is.na(lab)]; lab <- lab[!is.na(lab)]
  present <- unique(lab)
  lv <- if (!is.null(levels)) intersect(levels, present)
        else if (!is.null(palette)) intersect(names(palette), present) else sort(present)
  lv <- c(lv, setdiff(sort(present), lv))
  so$.grp <- factor(lab, levels = lv)

  gone <- setdiff(names(features), lv)
  if (length(gone)) message("plot_dot_panel(): no cells for module(s) ", toString(gone))
  mods <- c(extra, features[intersect(lv, names(features))])
  mods <- lapply(mods, intersect, rownames(so))
  genes <- unlist(mods, use.names = FALSE)
  mod_of <- rep(names(mods), lengths(mods))
  if (any(duplicated(genes)))
    message("plot_dot_panel(): shown once only: ", toString(unique(genes[duplicated(genes)])))
  keep <- !duplicated(genes); genes <- genes[keep]; mod_of <- mod_of[keep]
  if (!length(genes)) stop("plot_dot_panel(): no marker genes present", call. = FALSE)

  d <- Seurat::DotPlot(so, features = genes, group.by = ".grp", scale = TRUE,
                       col.min = -cap, col.max = cap)$data
  d$gene   <- as.character(d$features.plot)
  d$module <- factor(mod_of[match(d$gene, genes)], levels = unique(mod_of))
  d$id     <- as.character(d$id)
  d$pct.exp[d$pct.exp == 0] <- NA
  rm(so)

  pal <- if (is.null(palette)) setNames(rep("grey60", length(lv)), lv)
         else umap_palette_for(list(palettes = list(x = palette)), "x", lv)
  wrap_n <- FIG_LEGEND_WRAP
  wrap <- function(x) stringr::str_wrap(x, wrap_n)
  on.exit(rm(list = setdiff(ls(), c("wrap_n", "wrap", "genes", "lv"))), add = TRUE)
  tag  <- if (is.null(panel)) NULL
          else if (identical(FIG_PANEL_CASE, "upper")) toupper(panel) else tolower(panel)

  base <- ggplot2::theme_classic(base_size = FIG_BASE_PT, base_family = family) +
    ggplot2::theme(
      axis.line  = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = FIG_LEGEND_PT),
      legend.text  = ggplot2::element_text(size = FIG_MIN_PT),
      legend.key.size = grid::unit(3, "mm"),
      plot.title = ggplot2::element_text(size = FIG_BASE_PT, hjust = 0),
      plot.tag   = ggplot2::element_text(size = FIG_PANEL_PT, face = "bold"),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(2, 1, 2, 1, "mm"))
  scales_expr <- list(
    ggplot2::geom_point(ggplot2::aes(size = pct.exp, colour = avg.exp.scaled), na.rm = TRUE),
    ggplot2::scale_size_area(max_size = 3, limits = c(0, 100),
                             breaks = c(25, 50, 75, 100), name = "% expressing"),
    ggplot2::scale_colour_gradientn(colours = FIG_EXPR_COLS, limits = c(-cap, cap),
                                    oob = scales::squish, name = "Scaled\nexpression"),
    ggplot2::guides(
      colour = ggplot2::guide_colourbar(order = 1, barwidth = grid::unit(2, "mm"),
                                        barheight = grid::unit(14, "mm")),
      size   = ggplot2::guide_legend(order = 2)))
  gene_txt <- ggplot2::element_text(size = FIG_LEGEND_PT, face = "italic", colour = "black")
  pop_txt  <- ggplot2::element_text(size = FIG_LEGEND_PT, colour = "black")
  key_df   <- data.frame(id = lv)

  dots <- ggplot2::ggplot(d, ggplot2::aes(factor(gene, levels = genes),
                                          factor(id, levels = rev(lv)))) +
    scales_expr +
    ggplot2::facet_grid(~ module, scales = "free_x", space = "free_x") +
    base +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = FIG_LEGEND_PT, face = "italic",
                                                       colour = "black", angle = 90,
                                                       hjust = 1, vjust = 0.5),
                   panel.spacing.x = grid::unit(1.2, "mm"))
  key <- ggplot2::ggplot(key_df, ggplot2::aes(1, factor(id, levels = rev(lv)), colour = id)) +
    ggplot2::geom_point(size = 2, shape = 16, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::scale_y_discrete(labels = wrap) +
    base + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.text.y = pop_txt) +
    ggplot2::labs(title = title, tag = tag)
  fig <- if (colour_key) {
    patchwork::wrap_plots(key, dots, nrow = 1, widths = c(3, 3.4 * length(genes)))
  } else {
    dots + ggplot2::scale_y_discrete(labels = wrap) +
      ggplot2::theme(axis.text.y = pop_txt) + ggplot2::labs(title = title, tag = tag)
  }

  dots_f <- ggplot2::ggplot(d, ggplot2::aes(factor(id, levels = lv),
                                            factor(gene, levels = rev(genes)))) +
    scales_expr +
    ggplot2::facet_grid(module ~ ., scales = "free_y", space = "free_y") +
    base +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.text.y = gene_txt,
                   panel.spacing.y = grid::unit(1.2, "mm")) +
    ggplot2::labs(title = title, tag = tag)
  key_f <- ggplot2::ggplot(key_df, ggplot2::aes(factor(id, levels = lv), 1, colour = id)) +
    ggplot2::geom_point(size = 2, shape = 16, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::scale_x_discrete(labels = wrap) +
    base + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                          axis.text.x = ggplot2::element_text(size = FIG_LEGEND_PT,
                                                              colour = "black", angle = 90,
                                                              hjust = 1, vjust = 0.5))
  fig_f <- if (colour_key) {
    patchwork::wrap_plots(dots_f, key_f, ncol = 1, heights = c(3.4 * length(genes), 3))
  } else {
    dots_f + ggplot2::scale_x_discrete(labels = wrap) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = FIG_LEGEND_PT,
                       colour = "black", angle = 90, hjust = 1, vjust = 0.5))
  }

  lab_len  <- 1.15 * max(nchar(unlist(strsplit(wrap(lv), "\n", fixed = TRUE)))) + 4
  gene_len <- 1.1 * max(nchar(genes)) + 4
  n_gap    <- length(unique(mod_of)) - 1
  key_w <- if (colour_key) 6 else 0
  dims <- list(width  = lab_len + key_w + 3.4 * length(genes) + 1.2 * n_gap + 24,
               height = max(3.8 * length(lv) + gene_len + 12, 52))
  dims_f <- list(width  = gene_len + 3.8 * length(lv) + 26,
                 height = 3.4 * length(genes) + 1.2 * n_gap + 6 + lab_len + 12)

  dims$lab_w <- lab_len; dims$gene_len <- gene_len
  dims$n_genes <- length(genes); dims$n_gap <- n_gap; dims$n_levels <- length(lv)
  list(figure = fig, dims = dims, figure_flip = fig_f, dims_flip = dims_f,
       genes = genes, levels = lv)
}

dot_compartment_figure <- function(objs, features, label_col, stem,
                                   umap_fig = NULL, extra = NULL,
                                   dir = path$umap_figures) {
  out <- list()
  for (comp in names(objs)) {
    pal <- if (!is.null(umap_fig)) umap_fig$palettes[[comp]] else NULL
    dp  <- plot_dot_panel(objs[[comp]], features, label_col, palette = pal,
                          extra = extra, title = comp)
    save_umap_fig(dp$figure, paste0(stem, "_", sanitize(comp)), dir,
                  width_mm = round(dp$dims$width),
                  height_mm = min(FIG_MAX_DEPTH_MM, round(dp$dims$height)))
    out[[comp]] <- dp
  }
  invisible(out)
}

store_figure_panels <- function(level, part, obj, dir = path$figure_panels) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  if (part == "umap") obj$figure <- NULL
  f <- file.path(dir, paste0(level, ".qs"))
  x <- if (file.exists(f)) qs::qread(f) else list()
  x[[part]] <- obj
  x$updated[[part]] <- format(Sys.time())
  qs::qsave(x, f)
  message("figure panels: ", level, "/", part, " -> ", f, " (",
          format(structure(file.size(f), class = "object_size"), units = "auto"), ")")
  invisible(f)
}

abundance_figure_panels <- function(fits, umap_fig = NULL, ref_group = NULL,
                                    contrast = "X-ALD vs control", fdr = 0.05,
                                    y_lab = "proportion of cells") {
  purrr::imap(fits, function(f, comp) {
    pal <- umap_palette_for(umap_fig, comp, f$inputs$long$celltype)
    res <- if (!is.null(f$res_fig)) f$res_fig else f$res
    c(ab_panels(f$inputs$long, res, f$inputs$totals, pal, ref_group = ref_group,
                y_lab = y_lab, fdr = fdr, style = "figure",
                family = fig_font(quiet = TRUE), x_lab = f$x_lab),
      list(contrast = contrast, fdr = fdr, method = if (is.null(f$method)) "edgeR" else f$method))
  })
}
