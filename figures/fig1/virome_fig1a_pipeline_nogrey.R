# =============================================================================
# Virome library — Figure 1A-style plot
# Seo et al. 2023 Cell Fig. 1A layout for +ssRNA libraries
#
# Panels:  [phylum label] | [species count ←] | [3′ UTR length boxplot →]
#
# Inputs:
#   FASTA_FILE : FASTA with headers like >NC_001564.2|10129-10682_3UTR
#   CSV_FILE   : virus_stats.csv  (accession, family, phylum, order, ...)
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
FASTA_FILE  <- "pilot_library_20260319.uniqueIDs.fa"   # <── your FASTA
CSV_FILE    <- "virus_stats.csv"                     # <── your taxonomy CSV
OUT_FIGURE  <- "fig1a_style.pdf"
OUT_CSV     <- "virome_summary.csv"
FIGURE_W    <- 13
POINT_ALPHA <- 0.40
MIN_UTR_LEN <- 100    # exclude UTRs shorter than this (annotation stubs)

# ── 1. Parse FASTA headers ───────────────────────────────────────────────────
message("Parsing FASTA headers...")

headers <- readLines(FASTA_FILE)
headers <- sub("^>", "", headers[startsWith(headers, ">")])

rx      <- "^([^|]+)\\|(\\d+)-(\\d+)_(.+)$"
matched <- regmatches(headers, regexec(rx, headers))

parse_one <- function(m) {
  if (length(m) < 5) return(NULL)
  data.frame(
    accession  = m[2],
    acc_clean  = sub("\\.\\d+$", "", m[2]),
    start      = as.integer(m[3]),
    end        = as.integer(m[4]),
    utr_length = as.integer(m[4]) - as.integer(m[3]),
    feature    = m[5],
    stringsAsFactors = FALSE
  )
}

fasta_df <- bind_rows(lapply(matched, parse_one))

# annotated 3UTR only — drop inferred_3UTR
utr_df <- fasta_df %>%
  filter(feature == "3UTR") %>%
  select(acc_clean, utr_length, feature) %>%
  mutate(stub = utr_length < MIN_UTR_LEN)

n_stub <- sum(utr_df$stub)

utr_df <- filter(utr_df, !stub) %>%
  select(acc_clean, utr_length)

message(glue("  {nrow(utr_df)} annotated 3\u2032UTR entries retained"))
message(glue("  {n_stub} stub entries excluded (< {MIN_UTR_LEN} nt)"))
message(glue("  Length range: {min(utr_df$utr_length)}\u2013{max(utr_df$utr_length)} nt",
             "  |  median: {median(utr_df$utr_length)} nt"))

# ── 2. Load taxonomy CSV ─────────────────────────────────────────────────────
message("\nLoading taxonomy CSV...")

tax_df <- read_csv(CSV_FILE, show_col_types = FALSE) %>%
  rename_with(tolower) %>%
  mutate(
    acc_clean = sub("\\.\\d+$", "", accession),
    family    = if_else(is.na(family) | family == "Unannotated", "Unannotated", family),
    phylum    = if_else(is.na(phylum) | phylum == "Unannotated", "Unannotated", phylum),
    order     = if_else(is.na(order)  | order  == "Unannotated", "Unannotated", order)
  )

# ── 3. Join ──────────────────────────────────────────────────────────────────
joined <- utr_df %>%
  inner_join(tax_df, by = "acc_clean") %>%
  filter(!is.na(family))

message(glue("  {nrow(joined)} / {nrow(utr_df)} UTR entries matched to taxonomy"))

# ── 4. Per-family summaries ──────────────────────────────────────────────────
family_sum <- joined %>%
  group_by(phylum, order, family) %>%
  summarise(
    n_species  = n_distinct(acc_clean),
    utr_median = median(utr_length),
    utr_mean   = mean(utr_length),
    utr_min    = min(utr_length),
    utr_max    = max(utr_length),
    utr_q25    = quantile(utr_length, 0.25),
    utr_q75    = quantile(utr_length, 0.75),
    .groups    = "drop"
  )

write_csv(family_sum, OUT_CSV)
message(glue("\nSummary: {sum(family_sum$n_species)} accessions | ",
             "{n_distinct(family_sum$family)} families | ",
             "{n_distinct(joined$order[joined$order != 'Unannotated'])} orders"))

# ── 5. Colours ───────────────────────────────────────────────────────────────
phylum_colors <- c(
  "Pisuviricota"    = "#185FA5",
  "Kitrinoviricota" = "#C9542A",
  "Lenarviricota"   = "#3B8C3B",
  "Unannotated"     = "#999999"
)
extra_phyla <- setdiff(unique(joined$phylum), names(phylum_colors))
if (length(extra_phyla) > 0)
  phylum_colors <- c(phylum_colors,
                     setNames(hue_pal()(length(extra_phyla)), extra_phyla))

# ── 6. Y-axis order ─────────────────────────────────────────────────────────
phylum_rank <- c("Pisuviricota" = 1, "Kitrinoviricota" = 2,
                 "Lenarviricota" = 3, "Unannotated" = 99)

fam_order_df <- family_sum %>%
  group_by(family) %>%
  slice_max(n_species, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(phylum_idx = coalesce(phylum_rank[phylum], 50L)) %>%
  arrange(phylum_idx, n_species)

fam_levels <- fam_order_df$family

# ── 7. Plot data ─────────────────────────────────────────────────────────────
fam_total <- fam_order_df %>%
  mutate(family = factor(family, levels = fam_levels))

obs_df <- joined %>%
  mutate(family = factor(family, levels = fam_levels))

# ── 8. Left panel — species count bars ──────────────────────────────────────
max_count <- max(fam_total$n_species)

p_left <- ggplot(fam_total, aes(y = family, x = n_species, fill = phylum)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  scale_fill_manual(values = phylum_colors) +
  scale_x_reverse(
    limits = c(max_count * 1.05, 0),
    expand = expansion(mult = c(0.01, 0)),
    breaks = pretty_breaks(4),
    labels = label_number(accuracy = 1)
  ) +
  labs(x = "Species count", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y        = element_text(size = 7.5, hjust = 1,
                                      colour = "grey20", face = "italic"),
    axis.text.x        = element_text(size = 8, colour = "grey40"),
    axis.title.x       = element_text(size = 8.5, colour = "grey40"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    plot.margin        = margin(5, 2, 5, 4)
  )

# Phylum group labels
group_pos <- fam_total %>%
  mutate(fam_idx = as.integer(family)) %>%
  group_by(phylum) %>%
  summarise(mid_idx = mean(range(fam_idx)), .groups = "drop")

p_left <- p_left +
  geom_text(
    data        = group_pos,
    aes(x       = max_count * 1.04,
        y       = mid_idx,
        label   = phylum,
        colour  = phylum),
    inherit.aes = FALSE,
    angle = 90, size = 2.5, hjust = 0.5, fontface = "bold"
  ) +
  scale_colour_manual(values = phylum_colors, guide = "none")

# ── 9. Right panel — 3′ UTR boxplots ────────────────────────────────────────
p_right <- ggplot(obs_df, aes(y = family, x = utr_length, fill = phylum)) +
  geom_jitter(
    aes(colour = phylum),
    shape = 16, height = 0.25, size = 0.9,
    alpha = POINT_ALPHA, show.legend = FALSE
  ) +
  geom_boxplot(
    outlier.shape = NA, width = 0.55, linewidth = 0.35,
    show.legend   = TRUE
  ) +
  scale_fill_manual(
    values = phylum_colors, name = "Phylum",
    guide  = guide_legend(title.position = "top",
                          override.aes   = list(alpha = 1, colour = NA))
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
    axis.text.x        = element_text(size = 8, colour = "grey40"),
    axis.title.x       = element_text(size = 8.5, colour = "grey40"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    legend.position    = "right",
    legend.text        = element_text(size = 8),
    legend.title       = element_text(size = 9, face = "bold"),
    legend.key.size    = unit(0.65, "lines"),
    plot.margin        = margin(5, 8, 5, 2)
  )

# ── 10. Assemble & save ──────────────────────────────────────────────────────
n_fam <- length(fam_levels)
fig_h <- max(8, n_fam * 0.32 + 2.5)

fig <- p_left + p_right +
  plot_layout(widths = c(1, 1.6)) +
  plot_annotation(
    title    = "+ssRNA virome library \u2014 family distribution and 3\u2032 UTR lengths",
    subtitle = glue(
      "{sum(family_sum$n_species)} sequences  \u00b7  ",
      "{n_distinct(family_sum$family)} families  \u00b7  ",
      "{n_distinct(joined$order[joined$order != 'Unannotated'])} orders  \u00b7  ",
      "{n_distinct(joined$phylum[joined$phylum != 'Unannotated'])} phyla  \u00b7  ",
      "annotated 3\u2032UTRs only  \u00b7  stubs < {MIN_UTR_LEN} nt excluded (n = {n_stub})"
    ),
    caption = paste(
      "Left: species count per family, grouped by phylum.",
      "Right: annotated 3\u2032 UTR length \u2014 box = IQR, line = median, points = individual sequences.",
      "Inferred 3\u2032UTRs and stubs excluded. Lengths from GenBank feature coordinates."
    ),
    theme = theme(
      plot.title    = element_text(size = 13, face = "bold", colour = "grey10"),
      plot.subtitle = element_text(size = 9,  colour = "grey40"),
      plot.caption  = element_text(size = 7.5, colour = "grey55",
                                   margin = margin(t = 6)),
      plot.margin   = margin(10, 10, 8, 10)
    )
  )

ggsave(OUT_FIGURE, fig, width = FIGURE_W, height = fig_h)

message("\nFigure saved \u2192 ", OUT_FIGURE)
message("Summary CSV  \u2192 ", OUT_CSV)
