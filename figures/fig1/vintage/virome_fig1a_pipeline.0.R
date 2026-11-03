# =============================================================================
# Virome library — Figure 1A-style plot (styled to match Seo et al. 2023)
#
# Left panel:
#   GREY bar   = total accessions per family in library
#   COLOUR bar = species with annotated 3′UTR ≥ MIN_UTR_LEN (overlaid)
# Right panel: 3′ UTR length boxplot per family
#
# Install once:
#   install.packages(c("tidyverse","patchwork","scales","ggtext","glue"))
# =============================================================================

library(tidyverse)
library(patchwork)
library(scales)
library(ggtext)
library(glue)

# ── 0. Config ────────────────────────────────────────────────────────────────
FASTA_FILE        <- "pilot_library_20260319.uniqueIDs.fa"
CSV_FILE          <- "virus_stats.csv"
NCBI_COUNTS_FILE  <- "ncbi_family_counts.csv"
OUT_FIGURE        <- "fig1a_style.pdf"
OUT_CSV           <- "virome_summary.csv"
FIGURE_W          <- 14
POINT_ALPHA       <- 0.45
MIN_UTR_LEN       <- 100
BAR_WIDTH         <- 0.55   # same width for both bars so they align perfectly
LABEL_FONT_SIZE   <- 3.2    # family label font size (in ggplot units; ~9pt)

# ── 1. Parse FASTA headers ───────────────────────────────────────────────────
message("Parsing FASTA headers...")

headers   <- readLines(FASTA_FILE)
headers   <- sub("^>", "", headers[startsWith(headers, ">")])
rx        <- "^([^|]+)\\|(\\d+)-(\\d+)_(.+)$"
matched   <- regmatches(headers, regexec(rx, headers))

parse_one <- function(m) {
  if (length(m) < 5) return(NULL)
  data.frame(
    acc_clean  = sub("\\.\\d+$", "", m[2]),
    utr_length = as.integer(m[4]) - as.integer(m[3]),
    feature    = m[5],
    stringsAsFactors = FALSE
  )
}

fasta_df <- bind_rows(lapply(matched, parse_one))
all_accs <- unique(fasta_df$acc_clean)

utr_raw  <- filter(fasta_df, grepl("3UTR", feature, ignore.case = TRUE))
n_stub   <- sum(utr_raw$utr_length < MIN_UTR_LEN)
utr_df   <- filter(utr_raw, utr_length >= MIN_UTR_LEN) %>%
            select(acc_clean, utr_length)

message(glue("  {length(all_accs)} unique accessions in library"))
message(glue("  {nrow(utr_df)} with annotated 3\u2032UTR \u2265 {MIN_UTR_LEN} nt"))
message(glue("  {n_stub} stubs excluded"))

# ── 2. Taxonomy ──────────────────────────────────────────────────────────────
message("\nLoading taxonomy CSV...")

artifact_phyla <- c("Macronovirus", "Tobacco albetovirus 3")

tax_df <- read_csv(CSV_FILE, show_col_types = FALSE) %>%
  rename_with(tolower) %>%
  mutate(
    acc_clean = sub("\\.\\d+$", "", accession),
    phylum    = if_else(phylum %in% artifact_phyla, "Unannotated", phylum),
    family    = if_else(is.na(family) | family == "Unannotated", "Unannotated", family),
    phylum    = if_else(is.na(phylum) | phylum == "Unannotated", "Unannotated", phylum),
    order     = if_else(is.na(order)  | order  == "Unannotated", "Unannotated", order)
  )

# ── 3. Counts ────────────────────────────────────────────────────────────────
all_df <- tibble(acc_clean = all_accs) %>%
  inner_join(tax_df, by = "acc_clean")

lib_counts <- all_df %>%
  group_by(family, phylum, order) %>%
  summarise(n_lib = n_distinct(acc_clean), .groups = "drop")

utr_joined <- inner_join(utr_df, tax_df, by = "acc_clean") %>%
  filter(!is.na(family))

colour_counts <- utr_joined %>%
  group_by(family) %>%
  summarise(n_utr = n_distinct(acc_clean), .groups = "drop")

# NCBI totals capped at library count
ncbi <- read_csv(NCBI_COUNTS_FILE, show_col_types = FALSE)
message("NCBI counts loaded")

grey_counts <- lib_counts %>%
  left_join(ncbi, by = "family") %>%
  mutate(
    ncbi_total = replace_na(as.integer(ncbi_total), 0L),
    n_total    = pmin(ncbi_total, n_lib)   # cap NCBI at library count
  ) %>%
  mutate(n_total = if_else(n_total == 0, n_lib, n_total))  # fallback if no NCBI match

family_sum <- grey_counts %>%
  left_join(colour_counts, by = "family") %>%
  replace_na(list(n_utr = 0)) %>%
  left_join(
    utr_joined %>%
      group_by(family) %>%
      summarise(
        utr_median = median(utr_length),
        utr_mean   = mean(utr_length),
        utr_min    = min(utr_length),
        utr_max    = max(utr_length),
        utr_q25    = quantile(utr_length, 0.25),
        utr_q75    = quantile(utr_length, 0.75),
        .groups    = "drop"),
    by = "family"
  )

write_csv(family_sum, OUT_CSV)
message(glue("{sum(family_sum$n_total)} library total | ",
             "{sum(family_sum$n_utr)} with annotated UTR | ",
             "{n_distinct(family_sum$family)} families"))

# ── 4. Colour palette ────────────────────────────────────────────────────────
phylum_colors <- c(
  "Pisuviricota"    = "#4472C4",
  "Kitrinoviricota" = "#C9542A",
  "Lenarviricota"   = "#70AD47",
  "Unannotated"     = "#A9A9A9"
)
extra_phyla <- setdiff(unique(family_sum$phylum), names(phylum_colors))
if (length(extra_phyla) > 0)
  phylum_colors <- c(phylum_colors,
                     setNames(hue_pal()(length(extra_phyla)), extra_phyla))

# ── 5. Y-axis order ─────────────────────────────────────────────────────────
phylum_rank <- c("Pisuviricota" = 1, "Kitrinoviricota" = 2,
                 "Lenarviricota" = 3, "Unannotated" = 99)

fam_plot <- family_sum %>%
  mutate(
    phylum_idx = coalesce(phylum_rank[phylum], 50L),
    family     = factor(family)
  ) %>%
  arrange(phylum_idx, n_total) %>%
  mutate(family = fct_inorder(family))

fam_levels <- levels(fam_plot$family)

obs_df <- utr_joined %>%
  mutate(family = factor(family, levels = fam_levels)) %>%
  filter(!is.na(family))

# ── 6. Left panel ────────────────────────────────────────────────────────────
max_count <- max(fam_plot$n_total)
label_w   <- max_count * 0.55     # label box width
label_x   <- -(label_w / 2)       # label box centre

label_df  <- fam_plot %>% mutate(fam_idx = as.integer(family))

p_left <- ggplot(fam_plot, aes(y = family)) +

  # coloured label boxes
  geom_tile(
    data        = label_df,
    aes(x       = label_x,
        width   = label_w,
        height  = BAR_WIDTH + 0.15,   # box slightly taller than bar
        fill    = phylum),
    alpha       = 0.20,
    show.legend = FALSE
  ) +

  # family name — italic, phylum colour, centred in box
  geom_text(
    data        = label_df,
    aes(x       = label_x,
        label   = family,
        colour  = phylum),
    size        = LABEL_FONT_SIZE,
    hjust       = 0.5,
    fontface    = "italic",
    show.legend = FALSE
  ) +

  # grey bar: total in library (drawn first, underneath)
  geom_col(
    aes(x = n_total),
    fill  = "grey76",
    width = BAR_WIDTH     # exact same width as coloured bar
  ) +

  # coloured bar: annotated UTR subset (drawn on top, same width)
  geom_col(
    aes(x = n_utr, fill = phylum),
    width       = BAR_WIDTH,
    show.legend = FALSE
  ) +

  scale_fill_manual(values   = phylum_colors) +
  scale_colour_manual(values = phylum_colors) +
  scale_x_continuous(
    limits = c(label_x - label_w / 2 - 1, max_count * 1.07),
    breaks = pretty_breaks(4)(0:max_count),
    labels = label_number(accuracy = 1),
    expand = expansion(0)
  ) +
  labs(x = "Species count", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x  = element_text(size = 7.5, colour = "grey40"),
    axis.title.x = element_text(size = 8.5, colour = "grey40",
                                margin = margin(t = 4)),
    panel.grid   = element_blank(),
    plot.margin  = margin(5, 0, 5, 4)
  )

# phylum group labels (rotated, right margin)
group_pos <- fam_plot %>%
  mutate(fam_idx = as.integer(family)) %>%
  group_by(phylum) %>%
  summarise(y_mid = mean(range(fam_idx)), .groups = "drop")

p_left <- p_left +
  geom_text(
    data        = group_pos,
    aes(x       = max_count * 1.05,
        y       = y_mid,
        label   = phylum,
        colour  = phylum),
    inherit.aes = FALSE,
    angle = 90, size = 2.4, hjust = 0.5, fontface = "bold"
  )

# ── 7. Right panel ───────────────────────────────────────────────────────────
x_max <- max(obs_df$utr_length, na.rm = TRUE)

p_right <- ggplot(obs_df, aes(y = family, x = utr_length, fill = phylum)) +

  # faint row shading for families with UTR data
  geom_tile(
    data = fam_plot %>%
      filter(n_utr > 0) %>%
      mutate(family = factor(family, levels = fam_levels)),
    aes(y      = family,
        x      = x_max / 2,
        width  = x_max * 1.15,
        height = BAR_WIDTH + 0.15,
        fill   = phylum),
    alpha       = 0.07,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +

  geom_jitter(
    aes(colour = phylum),
    shape  = 16, height = 0.18, size = 0.85,
    alpha  = POINT_ALPHA, show.legend = FALSE
  ) +

  geom_boxplot(
    outlier.shape = NA,
    width         = BAR_WIDTH,     # match left panel bar width
    linewidth     = 0.32,
    colour        = "grey25",
    show.legend   = TRUE
  ) +

  scale_fill_manual(
    values = phylum_colors,
    name   = "Phylum",
    guide  = guide_legend(
      title.position = "top",
      override.aes   = list(alpha = 1, colour = NA, shape = NA,
                            linetype = 0, linewidth = 0)
    )
  ) +
  scale_colour_manual(values = phylum_colors, guide = "none") +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.04)),
    breaks = pretty_breaks(5),
    labels = label_number(suffix = " nt", big.mark = ",")
  ) +
  labs(x = "3\u2032 UTR length (nt)", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.text.x        = element_text(size = 7.5, colour = "grey40"),
    axis.title.x       = element_text(size = 8.5, colour = "grey40",
                                      margin = margin(t = 4)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    legend.position    = "right",
    legend.text        = element_text(size = 8),
    legend.title       = element_text(size = 9, face = "bold"),
    legend.key.size    = unit(0.7, "lines"),
    plot.margin        = margin(5, 8, 5, 0)
  )

# ── 8. Assemble ──────────────────────────────────────────────────────────────
n_fam <- length(fam_levels)
fig_h <- max(8, n_fam * 0.34 + 2.5)

fig <- p_left + p_right +
  plot_layout(widths = c(1.1, 1.4)) +
  plot_annotation(
    title    = "+ssRNA virome library \u2014 family distribution and 3\u2032 UTR lengths",
    subtitle = glue(
      "{sum(family_sum$n_total)} accessions in library  \u00b7  ",
      "{sum(family_sum$n_utr)} with 3\u2032UTR \u2265 {MIN_UTR_LEN} nt  \u00b7  ",
      "{n_distinct(family_sum$family)} families  \u00b7  ",
      "{n_distinct(family_sum$phylum[family_sum$phylum != 'Unannotated'])} phyla"
    ),
    caption = glue(
      "Left: grey bar = total accessions per family in library; ",
      "coloured bar = species with 3\u2032UTR \u2265 {MIN_UTR_LEN} nt. ",
      "Sequences < {MIN_UTR_LEN} nt excluded (n = {n_stub}). ",
      "Stubs < {MIN_UTR_LEN} nt and inferred 3\u2032UTRs excluded (n = {n_stub}). ",
      "Right: annotated 3\u2032UTR length \u2014 box = IQR, line = median, ",
      "points = individual sequences."
    ),
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold", colour = "grey10"),
      plot.subtitle = element_text(size = 8.5, colour = "grey45"),
      plot.caption  = element_text(size = 7,   colour = "grey55",
                                   margin = margin(t = 5)),
      plot.margin   = margin(10, 10, 8, 10)
    )
  )

ggsave(OUT_FIGURE, fig, width = FIGURE_W, height = fig_h)
message("\nFigure saved \u2192 ", OUT_FIGURE)
message("Summary CSV  \u2192 ", OUT_CSV)
