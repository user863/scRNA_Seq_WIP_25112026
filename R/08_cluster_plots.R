
medoid_centers <- function(emb, group) {
  lapply(levels(factor(group)), function(l) {
    d <- emb[group == l, 1:2, drop = FALSE]
    m <- c(stats::median(d[, 1]), stats::median(d[, 2]))
    i <- which.min((d[, 1] - m[1])^2 + (d[, 2] - m[2])^2)
    data.frame(grp = l, x = d[i, 1], y = d[i, 2])
  }) |> dplyr::bind_rows()
}

plot_clusters <- function(so, col, reduction = "umap_harmony", pt.size = 0.3) {
  emb <- as.data.frame(Embeddings(so, reduction)[, 1:2])
  colnames(emb) <- c("UMAP1", "UMAP2")
  emb$grp <- factor(as.character(so[[col, drop = TRUE]]))
  cent <- medoid_centers(emb, emb$grp)

  ggplot(emb, aes(UMAP1, UMAP2, colour = grp)) +
    geom_point(size = pt.size, alpha = 0.75, shape = 16) +
    geom_label_repel(data = cent, aes(x, y, label = grp), inherit.aes = FALSE,
                     size = 3, seed = 1, max.overlaps = Inf) +
    coord_equal() + theme_classic(base_size = 11) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          legend.position = "none") +
    labs(title = col, x = "UMAP 1", y = "UMAP 2")
}
