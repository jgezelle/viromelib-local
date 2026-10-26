#!/usr/bin/env Rscript
# 06d_plots.R - all step 06 figures, drawn from the tables written by 06, 06b and 06c.
#
# Runs locally (laptop). Get the input tables from the cluster first with 06d_fetch_inputs.sh.
#
# USAGE
#   Rscript 06d_plots.R <VIS folder>
#   where <VIS folder> holds the step 06 output folders:
#     clusters/        from 06_cluster_profiles.py
#     stall_windows/   from 06b_stall_windows.py
#     entropy/         from 06c_shannon_entropy.py   (or entropy_test/)
#   Any folder that is missing is skipped. Figures go to <folder>/figures/ as .svg and .pdf.
#
# FIGURES
#   clusters/figures/       clusters_pca, clusters_heatmap
#   stall_windows/figures/  window_vs_outside, landscape_heatmap, cluster_landscapes,
#                           cm_landscapes, taxon_landscapes
#   entropy/figures/        entropy_by_cm, entropy_vs_reactivity
#
# INPUTS (per folder; 06d_fetch_inputs.sh downloads exactly these)
#   clusters/       pca_scores.tsv pca_loadings.tsv pca_variance.tsv features_zscored.tsv
#   stall_windows/  window_features.tsv window_vs_outside.tsv window_clusters.tsv landscapes.tsv
#                   cm_similarity.tsv [taxon_similarity.tsv]
#   entropy/        entropy_windows.tsv entropy_cm_vs_none.tsv
#
# REQUIRES  R >= 4.0 with ggplot2; svglite recommended for .svg (falls back to grDevices::svg)
#           install once: install.packages(c("ggplot2", "svglite"))

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript 06d_plots.R <VIS folder>")
VIS <- args[1]

# ---- style ------------------------------------------------------------------------------------
COLORS <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
CM_COLORS <- c(Class1a = "#e87ba4", Class1b = "#008300", Class2 = "#4a3aa7", cl3 = "#e34948",
               cl4 = "#eda100", none = "#c9c8c3")
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; SURF <- "#fcfcfb"; NODATA <- "#e4e3df"

theme_dms <- function(base = 9) {
  theme_classic(base_size = base) +
    theme(plot.background = element_rect(fill = SURF, colour = NA),
          panel.background = element_rect(fill = SURF, colour = NA),
          axis.line = element_line(colour = MUTED, linewidth = 0.4),
          axis.ticks = element_line(colour = MUTED, linewidth = 0.4),
          axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2),
          plot.title = element_text(colour = INK, face = "bold", size = base + 2, hjust = 0),
          plot.title.position = "plot",
          strip.background = element_blank(),
          strip.text = element_text(colour = INK, face = "bold", hjust = 0, size = base + 1),
          legend.background = element_rect(fill = SURF, colour = NA),
          legend.key = element_rect(fill = SURF, colour = NA),
          legend.text = element_text(colour = INK2), legend.title = element_text(colour = INK2))
}
theme_set(theme_dms())
cluster_col <- function(k) COLORS[((as.integer(k) - 1) %% length(COLORS)) + 1]

save_both <- function(p, dir, name, w, h) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  svg_dev <- if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else grDevices::svg
  pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggsave(file.path(dir, paste0(name, ".svg")), p, width = w, height = h, device = svg_dev, bg = SURF)
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = pdf_dev, bg = SURF)
  message("  wrote ", file.path(dir, name), ".svg/.pdf")
}
read_tsv <- function(f) read.delim(f, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
NEEDS <- list(
  clusters = c("pca_scores.tsv", "pca_loadings.tsv", "pca_variance.tsv", "features_zscored.tsv"),
  stall_windows = c("window_features.tsv", "window_vs_outside.tsv", "window_clusters.tsv",
                    "landscapes.tsv", "cm_similarity.tsv"),
  entropy = c("entropy_windows.tsv", "entropy_cm_vs_none.tsv"))
have_inputs <- function(dir, kind) {
  miss <- NEEDS[[kind]][!file.exists(file.path(dir, NEEDS[[kind]]))]
  if (length(miss)) message("SKIPPING ", dir, ": missing ", paste(miss, collapse = ", "),
                            " (download with 06d_fetch_inputs.sh)")
  length(miss) == 0
}
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, formatC(p, format = "e", digits = 1), signif(p, 2)))

# landscapes.tsv (windows x offsets) -> long table with offset as a number
landscape_long <- function(L) {
  offs <- as.integer(sub("^o", "", colnames(L)[-1]))
  data.frame(window = rep(L[[1]], times = length(offs)),
             offset = rep(offs, each = nrow(L)),
             reactivity = unlist(L[-1], use.names = FALSE))
}
median_iqr <- function(d, by) {
  # median and IQR of reactivity per offset within groups `by`
  agg <- aggregate(reactivity ~ ., data = d[, c(by, "offset", "reactivity")],
                   FUN = function(v) c(med = median(v), q1 = quantile(v, .25, names = FALSE),
                                       q3 = quantile(v, .75, names = FALSE)))
  cbind(agg[, c(by, "offset")], as.data.frame(agg$reactivity))
}

# ---- 06: whole-RNA clusters ---------------------------------------------------------------------
plot_clusters <- function(dir) {
  message("clusters: ", dir); out <- file.path(dir, "figures")
  sc <- read_tsv(file.path(dir, "pca_scores.tsv")); ld <- read_tsv(file.path(dir, "pca_loadings.tsv"))
  vr <- read_tsv(file.path(dir, "pca_variance.tsv")); z <- read_tsv(file.path(dir, "features_zscored.tsv"))
  sc$cluster <- factor(sc$cluster)
  k_lab <- sapply(levels(sc$cluster), function(k) sprintf("cluster %s (n=%d)", k, sum(sc$cluster == k)))
  sx <- 0.75 * max(abs(sc$PC1)) / max(abs(ld$PC1)); sy <- 0.75 * max(abs(sc$PC2)) / max(abs(ld$PC2))
  ld$x <- ld$PC1 * sx; ld$y <- ld$PC2 * sy
  p <- ggplot(sc, aes(PC1, PC2)) +
    geom_segment(data = ld, aes(x = 0, y = 0, xend = x, yend = y), colour = MUTED,
                 arrow = arrow(length = unit(0.08, "inches")), linewidth = 0.4) +
    geom_point(aes(fill = cluster), shape = 21, colour = SURF, size = 2.4, stroke = 0.6) +
    geom_label(data = ld, aes(x = x * 1.12, y = y * 1.12, label = feature), size = 2.8, colour = INK2,
               fill = SURF, label.size = 0, label.padding = unit(0.1, "lines")) +
    scale_fill_manual(values = cluster_col(levels(sc$cluster)), labels = k_lab, name = NULL) +
    labs(x = sprintf("PC1 (%.0f%% of variance)", 100 * vr$var_explained[1]),
         y = sprintf("PC2 (%.0f%% of variance)", 100 * vr$var_explained[2]),
         title = sprintf("%d RNAs by A/C reactivity features (k=%d)", nrow(sc), nlevels(sc$cluster))) +
    theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1))
  save_both(p, out, "clusters_pca", 6.5, 5)

  feats <- setdiff(colnames(z), c("rna", "cluster", "tree_order"))
  zl <- data.frame(rna = rep(z$rna, length(feats)), tree_order = rep(z$tree_order, length(feats)),
                   feature = factor(rep(feats, each = nrow(z)), levels = feats),
                   value = unlist(z[feats], use.names = FALSE))
  strip <- data.frame(tree_order = z$tree_order, col = cluster_col(z$cluster))
  p <- ggplot(zl, aes(as.integer(feature), tree_order)) +
    geom_tile(aes(fill = pmax(pmin(value, 2.5), -2.5))) +
    geom_rect(data = strip, aes(xmin = -0.15, xmax = 0.35, ymin = tree_order - 0.5, ymax = tree_order + 0.5),
              fill = strip$col, inherit.aes = FALSE) +
    scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0,
                         limits = c(-2.5, 2.5), name = "z-score") +
    scale_x_continuous(breaks = seq_along(feats), labels = feats, expand = c(0, 0)) +
    scale_y_reverse(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = "Features (z-scored), RNAs in clustering-tree order; bar = cluster") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))
  save_both(p, out, "clusters_heatmap", 6.5, 8)
}

# ---- 06b: stall windows -------------------------------------------------------------------------
plot_stall_windows <- function(dir) {
  message("stall_windows: ", dir); out <- file.path(dir, "figures")
  Fe <- read_tsv(file.path(dir, "window_features.tsv")); cmp <- read_tsv(file.path(dir, "window_vs_outside.tsv"))
  C <- read_tsv(file.path(dir, "window_clusters.tsv")); L <- read_tsv(file.path(dir, "landscapes.tsv"))
  sim <- read_tsv(file.path(dir, "cm_similarity.tsv"))
  tf <- file.path(dir, "taxon_similarity.tsv"); taxsim <- if (file.exists(tf)) read_tsv(tf) else NULL
  cm_lv <- intersect(names(CM_COLORS), unique(C$cm_class))

  # window vs rest of RNA
  feats <- c("gini", "frac_low", "frac_high", "mean_rate")
  d <- do.call(rbind, lapply(feats, function(f) {
    x <- Fe[Fe$out_n_AC >= 20, ]
    data.frame(feature = f, outside = x[[paste0("out_", f)]], window = x[[paste0("win_", f)]],
               cm = ifelse(x$cm_class == "none", "no CM hit", "CM hit"))
  }))
  d <- d[complete.cases(d), ]; d$feature <- factor(d$feature, levels = feats)
  lim <- do.call(rbind, lapply(split(d, d$feature), function(x)
    data.frame(feature = x$feature[1], lo = min(x$outside, x$window), hi = max(x$outside, x$window))))
  ann <- merge(lim, cmp, by = "feature")
  ann$label <- sprintf("window higher in %.0f%%\np = %s", 100 * ann$frac_window_higher, fmt_p(ann$wilcoxon_p))
  p <- ggplot(d, aes(outside, window)) +
    geom_segment(data = lim, aes(x = lo, y = lo, xend = hi, yend = hi), colour = MUTED, linetype = 2, linewidth = 0.4) +
    geom_point(aes(fill = cm), shape = 21, colour = SURF, size = 1.8, stroke = 0.4) +
    geom_label(data = ann, aes(x = lo, y = hi, label = label), hjust = 0, vjust = 1, size = 2.6,
               colour = INK2, fill = SURF, label.size = 0) +
    facet_wrap(~feature, nrow = 1, scales = "free") +
    scale_fill_manual(values = c("CM hit" = COLORS[1], "no CM hit" = COLORS[2]), name = NULL) +
    labs(x = "rest of RNA", y = "stall window",
         title = "Stall window vs the rest of the same RNA (A/C, HQ positions); dashed = equal") +
    theme(legend.position = "bottom")
  save_both(p, out, "window_vs_outside", 12, 3.8)

  # landscape heatmap
  ll <- landscape_long(L)
  ll <- merge(ll, C[, c("window", "tree_order", "cluster", "cm_class")], by = "window")
  offs <- sort(unique(ll$offset)); x0 <- min(offs)
  rows <- unique(ll[, c("window", "tree_order", "cluster", "cm_class")])
  ks <- sort(unique(rows$cluster))
  rows$band <- ifelse(match(rows$cluster, ks) %% 2 == 1, "#d7d6d1", "#9a9893")
  rows$cmcol <- CM_COLORS[rows$cm_class]; rows$cmcol[is.na(rows$cmcol)] <- MUTED
  knum <- aggregate(tree_order ~ cluster, data = rows, FUN = mean)
  leg <- data.frame(cm_class = factor(cm_lv, levels = cm_lv), x = x0, y = 1)
  p <- ggplot(ll, aes(offset, tree_order)) +
    geom_tile(aes(fill = pmin(reactivity, 1.5))) +
    geom_rect(data = rows, aes(xmin = x0 - 9.5, xmax = x0 - 6, ymin = tree_order - 0.5, ymax = tree_order + 0.5),
              fill = rows$band, inherit.aes = FALSE) +
    geom_rect(data = rows, aes(xmin = x0 - 5.5, xmax = x0 - 2, ymin = tree_order - 0.5, ymax = tree_order + 0.5),
              fill = rows$cmcol, inherit.aes = FALSE) +
    geom_text(data = knum, aes(x = x0 - 11, y = tree_order, label = cluster), hjust = 1, size = 2.8, colour = INK2) +
    geom_point(data = leg, aes(x, y, colour = cm_class), alpha = 0, inherit.aes = FALSE) +
    scale_colour_manual(values = CM_COLORS[cm_lv], name = "CM class",
                        guide = guide_legend(override.aes = list(alpha = 1, shape = 15, size = 4))) +
    scale_fill_gradientn(colours = c("#ffffcc", "#fed976", "#fd8d3c", "#e31a1c", "#800026"),
                         limits = c(0, 1.5), na.value = NODATA, name = "normalized\nreactivity") +
    scale_y_reverse(expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(0, max(offs), 10), expand = c(0, 0), limits = c(x0 - 13, max(offs) + 0.5)) +
    coord_cartesian(clip = "off") +
    labs(x = "offset from stall site (nt)", y = NULL,
         title = sprintf("%d stall windows, smoothed A/C reactivity (gray = no data)", nrow(rows)),
         subtitle = "Left bars: cluster (numbered gray bands), then CM class. Rows in clustering-tree order.") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          plot.subtitle = element_text(colour = INK2, size = 8))
  save_both(p, out, "landscape_heatmap", 9, max(6, 0.05 * nrow(rows) + 2))

  # median landscape per cluster
  mc <- median_iqr(ll, "cluster")
  comp <- aggregate(cm_class ~ cluster, data = rows, FUN = function(v) {
    t <- sort(table(v), decreasing = TRUE); paste(sprintf("%s: %d", names(t), t), collapse = "\n") })
  comp$n <- as.integer(table(rows$cluster)[as.character(comp$cluster)])
  mc$panel <- factor(sprintf("cluster %s (n=%d)", mc$cluster, comp$n[match(mc$cluster, comp$cluster)]))
  comp$panel <- factor(sprintf("cluster %s (n=%d)", comp$cluster, comp$n), levels = levels(mc$panel))
  p <- ggplot(mc, aes(offset)) +
    geom_ribbon(aes(ymin = q1, ymax = q3, fill = factor(cluster)), alpha = 0.25) +
    geom_line(aes(y = med, colour = factor(cluster)), linewidth = 0.8) +
    geom_text(data = comp, aes(x = max(offs), y = Inf, label = cm_class), hjust = 1, vjust = 1.1,
              size = 2.3, colour = INK2, lineheight = 0.9) +
    facet_wrap(~panel, ncol = 4) +
    scale_fill_manual(values = cluster_col(ks), guide = "none") +
    scale_colour_manual(values = cluster_col(ks), guide = "none") +
    labs(x = "offset from stall (nt)", y = "A/C reactivity (median, IQR)",
         title = "Median stall-window landscape per cluster")
  save_both(p, out, "cluster_landscapes", 12, 2.8 * ceiling(length(ks) / 4) + 0.6)

  # median landscape per CM class
  mcm <- median_iqr(ll, "cm_class"); mcm$cm_class <- factor(mcm$cm_class, levels = cm_lv)
  lab <- data.frame(cm_class = factor(cm_lv, levels = cm_lv))
  lab$n <- as.integer(table(rows$cm_class)[cm_lv])
  lab$p <- sim$p_value[match(cm_lv, sim$cm_class)]
  lab$label <- ifelse(is.na(lab$p), sprintf("n=%d", lab$n), sprintf("n=%d\nshared pattern p=%s", lab$n, fmt_p(lab$p)))
  p <- ggplot(mcm, aes(offset)) +
    geom_ribbon(aes(ymin = q1, ymax = q3, fill = cm_class), alpha = 0.25) +
    geom_line(aes(y = med, colour = cm_class), linewidth = 0.8) +
    geom_text(data = lab, aes(x = max(offs), y = Inf, label = label), hjust = 1, vjust = 1.1, size = 2.5, colour = INK2) +
    facet_wrap(~cm_class, nrow = 1) +
    scale_fill_manual(values = CM_COLORS, guide = "none") + scale_colour_manual(values = CM_COLORS, guide = "none") +
    labs(x = "offset from stall (nt)", y = "A/C reactivity (median, IQR)",
         title = "Stall-window reactivity by Infernal CM class")
  save_both(p, out, "cm_landscapes", 3 * length(cm_lv), 3)

  # taxa whose no-CM windows share a pattern (p < 0.05); skip a taxon whose windows repeat one shown
  if (!is.null(taxsim) && nrow(taxsim)) {
    ranks <- c("genus", "family")
    h <- taxsim[taxsim$subset == "no_CM" & taxsim$p_value < 0.05, ]
    h <- h[order(h$p_value, match(h$rank, ranks)), ]
    seen <- character(0); keep <- logical(nrow(h))
    for (i in seq_len(nrow(h))) {
      w <- sort(C$window[C$cm_class == "none" & C[[h$rank[i]]] %in% h$group[i]])
      key <- paste(w, collapse = ";"); keep[i] <- !(key %in% seen); seen <- c(seen, key)
    }
    h <- h[keep, ]
    if (nrow(h)) {
      h$panel <- factor(sprintf("%s (%s)", h$group, h$rank), levels = sprintf("%s (%s)", h$group, h$rank))
      tl <- do.call(rbind, lapply(seq_len(nrow(h)), function(i) {
        w <- C$window[C$cm_class == "none" & C[[h$rank[i]]] %in% h$group[i]]
        x <- ll[ll$window %in% w, c("window", "offset", "reactivity")]; x$panel <- h$panel[i]; x }))
      tm <- aggregate(reactivity ~ panel + offset, data = tl, FUN = median)
      h$label <- sprintf("%d windows, %d RNAs\np=%s  q=%s", h$n_windows, h$n_rnas, fmt_p(h$p_value), fmt_p(h$q_value))
      p <- ggplot(tl, aes(offset, reactivity)) +
        geom_line(aes(group = window), colour = MUTED, linewidth = 0.3, alpha = 0.7, na.rm = TRUE) +
        geom_line(data = tm, colour = COLORS[1], linewidth = 0.9) +
        geom_text(data = h, aes(x = max(offs), y = Inf, label = label), hjust = 1, vjust = 1.1, size = 2.4, colour = INK2) +
        facet_wrap(~panel, ncol = 4) +
        labs(x = "offset from stall (nt)", y = "A/C reactivity",
             title = "No-CM windows: taxa whose windows share a reactivity pattern (gray = each window, blue = median)")
      save_both(p, out, "taxon_landscapes", 12, 2.8 * ceiling(nrow(h) / 4) + 0.6)
    }
  }
}

# ---- 06c: entropy ---------------------------------------------------------------------------------
plot_entropy <- function(dir) {
  message("entropy: ", dir); out <- file.path(dir, "figures")
  S <- read_tsv(file.path(dir, "entropy_windows.tsv")); Tt <- read_tsv(file.path(dir, "entropy_cm_vs_none.tsv"))
  cm_lv <- intersect(names(CM_COLORS), unique(S$cm_class))
  d <- rbind(data.frame(fold = "seq", cm_class = S$cm_class, entropy = S$mean_entropy_seq),
             data.frame(fold = "dms", cm_class = S$cm_class, entropy = S$mean_entropy_dms))
  d <- d[!is.na(d$entropy), ]
  nlab <- sapply(cm_lv, function(c) sprintf("%s\n(n=%d)", c, sum(S$cm_class == c)))
  d$cm_class <- factor(d$cm_class, levels = cm_lv)
  ftitle <- c(seq = "sequence only", dms = "with DMS restraints")
  d$panel <- factor(sprintf("%s: CM hit vs none p = %s", ftitle[d$fold], fmt_p(Tt$mannwhitney_p[match(d$fold, Tt$fold)])))
  d$panel <- factor(d$panel, levels = unique(d$panel[order(match(d$fold, c("seq", "dms")))]))
  med <- aggregate(entropy ~ panel + cm_class, data = d, FUN = median)
  p <- ggplot(d, aes(cm_class, entropy)) +
    geom_jitter(aes(fill = cm_class), shape = 21, colour = SURF, size = 1.8, stroke = 0.3, width = 0.18, height = 0) +
    geom_errorbar(data = med, aes(ymin = entropy, ymax = entropy), width = 0.6, colour = INK, linewidth = 0.6) +
    facet_wrap(~panel, nrow = 1) +
    scale_fill_manual(values = CM_COLORS, guide = "none") +
    scale_x_discrete(labels = nlab) +
    labs(x = NULL, y = "window mean Shannon entropy", title = "Shannon entropy of the stall windows")
  save_both(p, out, "entropy_by_cm", 10, 3.8)

  S$cm_class <- factor(S$cm_class, levels = cm_lv)
  p <- ggplot(S, aes(mean_reactivity_AC, mean_entropy_dms)) +
    geom_point(aes(fill = cm_class), shape = 21, colour = SURF, size = 2, stroke = 0.3, na.rm = TRUE) +
    scale_fill_manual(values = CM_COLORS, name = "CM class") +
    labs(x = "window mean A/C reactivity", y = "window mean entropy (DMS fold)",
         title = "Low reactivity + low entropy = well-defined fold")
  save_both(p, out, "entropy_vs_reactivity", 5.5, 4.2)
}

# ---- run ------------------------------------------------------------------------------------------
if (!dir.exists(VIS)) stop("Folder not found: ", VIS)
did <- 0
for (kind in c("clusters", "stall_windows")) {
  d <- file.path(VIS, kind)
  if (dir.exists(d) && have_inputs(d, kind)) {
    if (kind == "clusters") plot_clusters(d) else plot_stall_windows(d); did <- did + 1 }
}
for (e in c("entropy", "entropy_test")) {
  d <- file.path(VIS, e)
  if (dir.exists(d) && have_inputs(d, "entropy")) { plot_entropy(d); did <- did + 1 }
}
if (did == 0) stop("Nothing plotted: no complete step 06 output folders in ", VIS)
message("done")
