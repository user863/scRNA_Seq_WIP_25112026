
dotplot_gradient <- c("#1a0a4b", "#6b2fa0", "#e07b39", "#f5e642")

make_celltype_palette <- function(labels, n_base = 10) {
  labels <- sort(unique(as.character(na.omit(labels))))
  base_cols <- if (requireNamespace("ggsci", quietly = TRUE)) {
    ggsci::pal_npg("nrc")(n_base)
  } else {
    scales::hue_pal()(n_base)
  }
  cols <- if (length(labels) > n_base) {
    colorRampPalette(base_cols)(length(labels))
  } else {
    base_cols[seq_along(labels)]
  }
  setNames(cols, labels)
}

plot_labeled_umap <- function(so, group_col, reduction = "umap_harmony",
                              palette = NULL, facet_col = NULL,
                              label = TRUE, pt_size = 0.1,
                              center_fun = c("mean", "median")) {
  center_fun <- match.arg(center_fun)
  cfun <- if (center_fun == "mean") mean else median

  emb <- as.data.frame(Embeddings(so, reduction))[, 1:2]
  colnames(emb) <- c("UMAP_1", "UMAP_2")

  df <- so@meta.data %>%
    dplyr::select(dplyr::all_of(c(group_col, facet_col))) %>%
    dplyr::bind_cols(emb) %>%
    dplyr::arrange(!is.na(.data[[group_col]]))

  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = .data[[group_col]],
                   alpha = is.na(.data[[group_col]])),
               size = pt_size) +
    scale_alpha_manual(values = setNames(c(0.1, 1), c(TRUE, FALSE)), guide = "none") +
    labs(x = "UMAP1", y = "UMAP2", color = group_col) +
    theme_classic() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.5),
          legend.position = "right",
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9, face = "bold"),
          plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white", color = NA))

  if (!is.null(palette)) {
    p <- p + scale_color_manual(values = palette, na.value = "grey80", drop = TRUE)
  }

  if (label) {
    centers <- df %>%
      dplyr::filter(!is.na(.data[[group_col]])) %>%
      dplyr::group_by(.data[[group_col]]) %>%
      dplyr::summarise(x = cfun(UMAP_1, na.rm = TRUE),
                       y = cfun(UMAP_2, na.rm = TRUE), .groups = "drop")
    p <- p + geom_label(data = centers,
                        aes(x = x, y = y, label = .data[[group_col]]))
  }

  if (!is.null(facet_col)) p <- p + facet_wrap(stats::as.formula(paste("~", facet_col)))
  p
}

plot_marker_dotplot <- function(so, features, group_by, title = NULL,
                                gradient = dotplot_gradient,
                                order = NULL, order_labels = NULL) {
  p <- DotPlot(so, features = features, group.by = group_by) +
    scale_color_gradientn(colors = gradient, name = "Mean expression\nin group") +
    labs(title = title, x = NULL, y = NULL, size = "Fraction of cells\nin group (%)") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, face = "bold", size = 11),
      axis.text.y = element_text(size = 9, face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title.position = "top"
    ) +
    coord_flip() +
    scale_size(range = c(0, 8), limits = c(0, 100), breaks = c(20, 40, 60, 80))

  if (!is.null(order)) {
    p <- p + scale_y_discrete(position = "right", limits = order,
                              labels = order_labels %||% order)
  } else {
    p <- p + scale_y_discrete(position = "right")
  }
  p
}
`%||%` <- function(a, b) if (is.null(a)) b else a
