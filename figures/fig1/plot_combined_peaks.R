# plot_combined_peaks.R ─────────────────────────────────────────────────────────
#
# Visualisations of cross-replicate Xrn1 stall peaks.
#
# Input  : `post_process_peaks.py` output (raw, per-peak CSV) where the
#          source was a `combine_replicates.py` cross-rep CSV.  Expected
#          columns include the meta_* statistics, n_reps_passing, per-rep
#          `treated_sf_R*` / `tile_reads_treated_R*`, taxonomy columns
#          (family / genus / species), and one `*_hit` column per Infernal
#          .tbl that was annotated.
#
# Output : Four plots, each bound to a named object so they can be displayed
#          interactively or saved with ggsave().
#          p1   — fine-grained class ranking: -log10(meta-p) vs fractional
#                 rank, coloured by xrRNA class (1a / 1b / 2 / 3 / None).
#          p1b  — same ranking with broad class grouping (Flavivirus /
#                 Plant infecting / None), faded for unclassified, with
#                 inline labels for highlighted virus exemplars.  This is
#                 the publication-style version of p1.  Saved to rank_plot.svg
#                 in the DISPLAY section.
#          p2   — stall fraction with per-rep min/max as error bars and
#                 coverage as size, restricted to high-SF peaks.
#          p3   — per-rep vs median SF, consistency diagnostic over all
#                 peaks.
#
# Adjust the CONFIG section to change the input file, thresholds, the
# highlighted-virus label set, or the CM column mapping.

# ── CONFIG ─────────────────────────────────────────────────────────────────────

input_csv         <- "../work/10_postproc/123_peaks_annotated.csv"
min_reps_passing  <- 2     # cross-rep call must pass in at least N replicates
high_sf_threshold <- 0.5   # threshold (inclusive) for "strong stall" plots p1/p1b/p2

# Highlighted peaks: annotated inline on p1b.  One row per peak you want
# labelled.  Multi-line labels are written with literal \n.  Peaks not
# in this table get no label.
highlight_labels <- tibble::tribble(
  ~accession,                                 ~stall_position, ~label,
  "NC_015843.2|10322-10990_inferred_3UTR",    142,             "TMUV\nNC_015843.2",
  "NC_012532.1|genome",                       10379,           "ZIKV\nNC_012532.1",
  "NC_001564.2|10129-10682_3UTR",             131,             "CFAV\nNC_001564.2",
  "NC_029906.1|3265-3557_intergenic",         63,              "BVG\nNC_029906.1"
)

# Infernal CM-hit columns produced by post_process_peaks (one per .tbl file).
# Listed in priority order: case_when below assigns the first matching label
# to each peak's `cm_class`.
CM_HIT_1A <- "CL_Class_1a_xrRNA_final.cm_hit"
CM_HIT_1B <- "CL_Class_1b_xrRNA_final.cm_hit"
CM_HIT_2  <- "CL_Class_2_xrRNA_final.cm_hit"
CM_HIT_3  <- "cl3_NAR_2025.align.mpkadded.edited.confident.dedup.IDcurated.SS_cons.cm_hit"
NO_CM_HIT <- "None found"

# Broad grouping of the fine-grained xrRNA classes for p1b.
# Helper returns a character vector of broad labels for a vector of cm_class
# values; the explicit Inf-relevel later keeps "None found" last in the legend.
broad_cm_class_for <- function(cm) {
  dplyr::case_when(
    cm %in% c("1a", "1b", "2") ~ "Flavivirus",
    cm == "3"                  ~ "Plant infecting",
    .default = NO_CM_HIT
  )
}


# ── LOAD ───────────────────────────────────────────────────────────────────────

df <- readr::read_csv(input_csv, show_col_types = FALSE)


# ── DERIVE ─────────────────────────────────────────────────────────────────────
#
# Build the visualisation frame.  Two dedup steps come first:
#   (1) require multi-replicate calls (n_reps_passing >= min_reps_passing);
#   (2) when two rows share an (accession, stall_position) — which can happen
#       when multiple overlapping tiles detect the same site — keep only the
#       most-significant row, to avoid one site being plotted twice.
# Then add derived columns used by the plots: mlog10p, an id string, the
# fine-grained xrRNA class assignment, the broad xrRNA class for p1b, an
# accession base name (sans the `|...` suffix), and `point_alpha` (used by
# p1b to fade unclassified peaks into the background).  Inline labels for
# the highlighted-virus set are joined from `highlight_labels`.

vis_df <- df |>
  dplyr::filter(n_reps_passing >= min_reps_passing) |>
  dplyr::filter(meta_pvalue == min(meta_pvalue),
                .by = c(accession, stall_position)) |>
  dplyr::mutate(
    mlog10p  = -log10(meta_pvalue),
    rank     = rank(-mlog10p),
    id       = glue::glue("{accession}-{tile_start}-{tile_end}-{stall_position}"),
    cm_class = dplyr::case_when(
      !is.na(.data[[CM_HIT_1A]]) ~ "1a",
      !is.na(.data[[CM_HIT_1B]]) ~ "1b",
      !is.na(.data[[CM_HIT_2]])  ~ "2",
      !is.na(.data[[CM_HIT_3]])  ~ "3",
      .default = NO_CM_HIT
    ),
    broad_cm_class = broad_cm_class_for(cm_class) |>
      forcats::fct_relevel(NO_CM_HIT, after = Inf),
    # Unclassified peaks fade to half opacity on p1b so the classified
    # points (the story) stand out.
    point_alpha    = ifelse(cm_class == NO_CM_HIT, 0.5, 1.0),
    accession_base = stringr::str_remove(accession, "[|].*$"),
    .before = 1
  ) |>
  dplyr::left_join(highlight_labels, by = c("accession", "stall_position")) |>
  dplyr::mutate(id = forcats::fct_reorder(id, rank, .desc = TRUE))


# ── PLOT HELPERS ───────────────────────────────────────────────────────────────

# Palettes.  cm_class_colours: colourblind-friendly Paul-Tol set for the four
# fine-grained classes (1a/1b/2/3) plus grey for unclassified.
# broad_cm_class_colours: hand-picked accessible pair for the two broad groups
# plus grey, keyed by the levels of `broad_cm_class` (Flavivirus /
# Plant infecting / None found).
cm_class_colours       <- c(rev(ggthemes::ptol_pal()(3)), "darkgray")
broad_cm_class_colours <- c("#CC6677", "#117733", "gray")

# Shared theme for the rank plots: ggplot2's theme_bw with the vertical
# gridlines suppressed — the x-axis is a dense rank and the gridlines just
# add visual noise.
rank_axis_theme <- ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank()
  )

# Helper: add a `fractional_rank` column (∈ [0, 1]) sorted by `metric`.
# Used by p1/p1b/p2 to put the most-significant point on the right of the
# x-axis with the rest distributed proportionally.
add_fractional_rank <- function(d, metric, descending = TRUE) {
  sign <- if (descending) -1 else 1
  d |>
    dplyr::mutate(
      subset_rank     = rank(sign * .data[[metric]]),
      fractional_rank = (subset_rank - 1) / max(subset_rank),
      .before = 1
    )
}


# ── PLOT 1: Significance ranking by fine-grained xrRNA class ──────────────────
# Restricted to high-SF peaks: fractional rank of -log10(meta-p) on the
# x-axis, value on y, coloured by 1a/1b/2/3/None.  Useful as a sanity check
# / detail view; p1b is the polished publication version.

p1 <- vis_df |>
  dplyr::filter(median_treated_sf >= high_sf_threshold) |>
  add_fractional_rank("mlog10p") |>
  ggplot2::ggplot(ggplot2::aes(fractional_rank, mlog10p, color = cm_class)) +
  ggplot2::geom_point(size = 0.7) +
  ggplot2::scale_color_manual(values = cm_class_colours, na.translate = FALSE) +
  ggplot2::guides(color = ggplot2::guide_legend(title = "xrRNA Class")) +
  ggplot2::labs(y = expression(-log[10](P)), x = "Fractional Rank") +
  rank_axis_theme


# ── PLOT 1b: Broad-class ranking with highlighted virus labels ────────────────
# Same data and axes as p1.  Differences:
#   - Coloured by broad_cm_class (Flavivirus / Plant infecting / None).
#   - Unclassified peaks faded via point_alpha = 0.5 to emphasise the story.
#   - Inline ggrepel labels for the highlighted-virus set defined at the top
#     of this script (highlight_labels tibble).
#   - "P-value" spelled out in the y-axis label and an annotation on the x-axis.

p1b <- vis_df |>
  dplyr::filter(median_treated_sf >= high_sf_threshold) |>
  add_fractional_rank("mlog10p") |>
  ggplot2::ggplot(ggplot2::aes(fractional_rank, mlog10p,
                               color = broad_cm_class, alpha = point_alpha)) +
  ggplot2::geom_point(size = 0.7) +
  ggrepel::geom_label_repel(
    ggplot2::aes(label = label),
    nudge_x = 0.05, nudge_y = 1, hjust = 0, show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(values = broad_cm_class_colours, na.translate = FALSE) +
  ggplot2::scale_alpha_identity() +
  ggplot2::guides(color = ggplot2::guide_legend(title = "xrRNA Class")) +
  ggplot2::labs(
    y = expression(-log[10](`P-value`)),
    x = "Fractional Rank\n(0 = highest signal)"
  ) +
  rank_axis_theme


# ── PLOT 2: Stall fraction with replicate spread ──────────────────────────────
# Same ranking as p1/p1b but the y-axis shows median treated_sf with
# per-replicate min/max as error bars and median tile read count (log scale)
# as point size.  Combines magnitude and across-replicate consistency in one
# plot.

p2 <- vis_df |>
  dplyr::filter(median_treated_sf >= high_sf_threshold) |>
  add_fractional_rank("mlog10p") |>
  # Pivot the per-rep sf columns to long form so we can compute per-id
  # min/max for the error bars.
  tidyr::pivot_longer(dplyr::starts_with("treated_sf_"),
                      names_to = "rep", values_to = "rep_sf") |>
  dplyr::mutate(
    sf_lower = min(rep_sf, na.rm = TRUE),
    sf_upper = max(rep_sf, na.rm = TRUE),
    .by = id
  ) |>
  ggplot2::ggplot(ggplot2::aes(
    x = fractional_rank, y = median_treated_sf,
    ymin = sf_lower, ymax = sf_upper,
    size = median_tile_reads_treated, color = cm_class
  )) +
  ggplot2::geom_pointrange() +
  ggplot2::scale_size(range = c(0, 0.75), transform = "log10") +
  ggplot2::scale_color_manual(values = cm_class_colours, na.translate = FALSE) +
  ggplot2::guides(
    color = ggplot2::guide_legend(title = "xrRNA Class"),
    size  = ggplot2::guide_legend(title = "Median Tile Reads")
  ) +
  ggplot2::labs(y = "Stall Fraction", x = "Fractional Rank") +
  rank_axis_theme


# ── PLOT 3: Per-rep vs median SF (consistency diagnostic) ─────────────────────
# For each peak, plots each replicate's stall fraction against the median
# across replicates.  Vertical scatter at a given median tells you the
# magnitude of per-rep disagreement.

p3 <- vis_df |>
  add_fractional_rank("median_treated_sf", descending = FALSE) |>
  tidyr::pivot_longer(dplyr::starts_with("treated_sf_"),
                      names_to = "replicate", values_to = "rep_sf") |>
  ggplot2::ggplot(ggplot2::aes(median_treated_sf, rep_sf, color = replicate)) +
  ggplot2::geom_point(size = 0.5, alpha = 0.7) +
  ggplot2::labs(x = "Median treated stall fraction",
                y = "Per-replicate stall fraction") +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.title = ggplot2::element_blank())


# ── DISPLAY / SAVE ────────────────────────────────────────────────────────────
# Plots are bound to p* names above so you can save them individually.  The
# main published artefact is p1b (broad classes + virus labels), saved to
# rank_plot.svg with a sensible aspect ratio for a figure panel.

ggplot2::ggsave("rank_plot_plant.svg", plot = p1b,
                width = 7, height = 5.7, units = "in")

# Uncomment to save the other plots:
#   ggplot2::ggsave("plot_1_rank.png",         p1,  width = 6, height = 4, dpi = 150)
#   ggplot2::ggsave("plot_2_sf_spread.png",    p2,  width = 6, height = 4, dpi = 150)
#   ggplot2::ggsave("plot_3_consistency.png",  p3,  width = 6, height = 4, dpi = 150)

p1
p1b
p2
p3
