# =============================================================================
# Virome library — Figure 1A-style plot (compact, Seo et al. 2023)
#
# Three panels flush together:
#   [A] Phylum label  (rotated text, coloured)
#   [B] Species count (grey = library total, colour = has 3′UTR ≥ MIN_UTR_LEN)
#   [C] 3′UTR length  (boxplot)
#
# install.packages(c("tidyverse","patchwork","scales","glue"))
# =============================================================================

library(tidyverse)
library(patchwork)
library(scales)
library(glue)

# ── 0. Config ────────────────────────────────────────────────────────────────
FASTA_FILE       <- "pilot_library_20260319.uniqueIDs.fa"
CSV_FILE         <- "virus_stats.csv"
OUT_FIGURE       <- "fig1a_style.pdf"
OUT_CSV          <- "virome_summary.csv"
MIN_UTR_LEN      <- 100
BAR_WIDTH        <- 0.65
POINT_ALPHA      <- 0.40
LABEL_SIZE       <- 2.8    # family name font (ggplot units)
FIGURE_W         <- 10     # inches — compact for a figure panel
FIGURE_H         <- 9      # fixed height — adjust if labels overlap

# ── 1. Parse FASTA ───────────────────────────────────────────────────────────
message("Parsing FASTA...")

headers  <- sub("^>", "", readLines(FASTA_FILE)[startsWith(readLines(FASTA_FILE), ">")])
rx       <- "^([^|]+)\\|(\\d+)-(\\d+)_(.+)$"
matched  <- regmatches(headers, regexec(rx, headers))

parse_one <- function(m) {
  if (length(m) < 5) return(NULL)
  data.frame(acc_clean  = sub("\\.\\d+$", "", m[2]),
             utr_length = as.integer(m[4]) - as.integer(m[3]),
             feature    = m[5], stringsAsFactors = FALSE)
}
fasta_df <- bind_rows(lapply(matched, parse_one))
all_accs <- unique(fasta_df$acc_clean)

utr_raw  <- filter(fasta_df, grepl("3UTR", feature, ignore.case = TRUE))
n_short  <- sum(utr_raw$utr_length < MIN_UTR_LEN)
utr_df   <- filter(utr_raw, utr_length >= MIN_UTR_LEN) %>% select(acc_clean, utr_length)

message(glue("  {length(all_accs)} accessions | {nrow(utr_df)} with 3'UTR >= {MIN_UTR_LEN} nt"))

# ── 2. Taxonomy ──────────────────────────────────────────────────────────────
artifact_phyla <- c("Macronovirus", "Tobacco albetovirus 3")

tax_df <- read_csv(CSV_FILE, show_col_types = FALSE) %>%
  rename_with(tolower) %>%
  mutate(
    acc_clean = sub("\\.\\d+$", "", accession),
    phylum    = if_else(phylum %in% artifact_phyla | is.na(phylum) |
                          phylum == "Unannotated", "Unannotated", phylum),
    family    = if_else(is.na(family) | family == "Unannotated", "Unannotated", family),
    order     = if_else(is.na(order)  | order  == "Unannotated", "Unannotated", order)
  )

# ── 3. Summaries ─────────────────────────────────────────────────────────────
all_df <- tibble(acc_clean = all_accs) %>% inner_join(tax_df, by = "acc_clean")

lib_counts <- all_df %>%
  group_by(family) %>%
  summarise(n_lib  = n_distinct(acc_clean),
            phylum = names(sort(table(phylum), decreasing = TRUE))[1],
            .groups = "drop")

utr_joined <- inner_join(utr_df, tax_df, by = "acc_clean") %>% filter(!is.na(family))

colour_counts <- utr_joined %>%
  group_by(family) %>% summarise(n_utr = n_distinct(acc_clean), .groups = "drop")

utr_stats <- utr_joined %>%
  group_by(family) %>%
  summarise(utr_med = median(utr_length), utr_q25 = quantile(utr_length, .25),
            utr_q75 = quantile(utr_length, .75), utr_min = min(utr_length),
            utr_max = max(utr_length), .groups = "drop")

family_sum <- lib_counts %>%
  left_join(colour_counts, by = "family") %>%
  replace_na(list(n_utr = 0)) %>%
  left_join(utr_stats, by = "family")

write_csv(family_sum, OUT_CSV)

# ── 4. Colours ───────────────────────────────────────────────────────────────
phylum_colors <- c(
  "Pisuviricota"    = "#4472C4",
  "Kitrinoviricota" = "#C0D6DF",
  "Lenarviricota"   = "#28965A",
  "Unannotated"     = "#999999"
)

# ── 5. Y-axis order ─────────────────────────────────────────────────────────
phylum_rank <- c("Pisuviricota"=1,"Kitrinoviricota"=2,"Lenarviricota"=3,"Unannotated"=99)

fam_plot <- family_sum %>%
  mutate(phylum_idx = coalesce(phylum_rank[phylum], 50L)) %>%
  arrange(phylum_idx, n_lib) %>%          # ascending → top of group = most species
  mutate(family = fct_inorder(factor(family)))

fam_levels <- levels(fam_plot$family)
n_fam      <- length(fam_levels)

obs_df <- utr_joined %>%
  mutate(family = factor(family, levels = fam_levels)) %>%
  filter(!is.na(family))

# ── 6. Shared theme ──────────────────────────────────────────────────────────
theme_panel <- function() {
  theme_minimal(base_size = 8) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.text.x        = element_text(size = 6.5, colour = "grey40"),
    axis.title.x       = element_text(size = 7.5, colour = "grey40",
                                      margin = margin(t = 3)),
    panel.grid         = element_blank(),
    plot.margin        = margin(2, 0, 2, 0)
  )
}

# ── 7. Panel A — phylum labels ───────────────────────────────────────────────
# One coloured rotated label per phylum group, positioned at group midpoint

phylum_pos <- fam_plot %>%
  mutate(y = as.integer(family)) %>%
  group_by(phylum) %>%
  summarise(y_mid = mean(range(y)), y_min = min(y), y_max = max(y), .groups = "drop")

p_phylum <- ggplot(phylum_pos, aes(x = 0, y = y_mid, label = phylum, colour = phylum)) +
  # bracket lines showing group extent
  geom_segment(aes(x = 0.3, xend = 0.3, y = y_min - 0.3, yend = y_max + 0.3),
               linewidth = 0.5, show.legend = FALSE) +
  geom_text(aes(x = 0), angle = 90, hjust = 0.5, fontface = "bold",
            size = 2.2, show.legend = FALSE) +
  scale_colour_manual(values = phylum_colors) +
  scale_y_continuous(limits = c(0.5, n_fam + 0.5), expand = expansion(0)) +
  scale_x_continuous(limits = c(-0.5, 0.5), expand = expansion(0)) +
  labs(x = NULL, y = NULL) +
  theme_panel() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        plot.margin = margin(2, 1, 2, 2))

# ── 8. Panel B — family labels + species count bars ─────────────────────────
max_count <- max(fam_plot$n_lib)

p_count <- ggplot(fam_plot, aes(y = family)) +

  # family name — italic, coloured by phylum, right-justified against bar
  geom_text(aes(x = 0, label = family, colour = phylum),
            hjust = 1, size = LABEL_SIZE, fontface = "italic",
            show.legend = FALSE) +

  # grey bar (library total)
  geom_col(aes(x = n_lib),  fill = "grey76", width = BAR_WIDTH) +

  # coloured bar (UTR subset)
  geom_col(aes(x = n_utr, fill = phylum), width = BAR_WIDTH, show.legend = FALSE) +

  scale_fill_manual(values = phylum_colors) +
  scale_colour_manual(values = phylum_colors) +
  scale_x_continuous(
    limits = c(-max_count * 0.75, max_count * 1.05),
    breaks = pretty_breaks(3)(0:max_count),
    labels = label_number(accuracy = 1),
    expand = expansion(0)
  ) +
  labs(x = "Species", y = NULL) +
  theme_panel() +
  theme(plot.margin = margin(2, 0, 2, 0))

# ── 9. Panel C — 3′UTR boxplot ───────────────────────────────────────────────
x_max <- max(obs_df$utr_length, na.rm = TRUE)

p_utr <- ggplot(obs_df, aes(y = family, x = utr_length, fill = phylum)) +

  geom_jitter(aes(colour = phylum), shape = 16, height = 0.2,
              size = 0.6, alpha = POINT_ALPHA, show.legend = FALSE) +

  geom_boxplot(outlier.shape = NA, width = BAR_WIDTH,
               linewidth = 0.28, colour = "grey30", show.legend = TRUE) +

  scale_fill_manual(
    values = phylum_colors, name = "Phylum",
    guide  = guide_legend(
      title.position = "top",
      override.aes   = list(colour = NA, linetype = 0, alpha = 1, shape = NA)
    )
  ) +
  scale_colour_manual(values = phylum_colors, guide = "none") +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.04)),
    breaks = pretty_breaks(4),
    labels = label_number(suffix = " nt", big.mark = ",")
  ) +
  labs(x = "3\u2032 UTR length (nt)", y = NULL) +
  theme_panel() +
  theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
        legend.position    = "right",
        legend.text        = element_text(size = 7),
        legend.title       = element_text(size = 8, face = "bold"),
        legend.key.size    = unit(0.6, "lines"),
        plot.margin        = margin(2, 4, 2, 0))

# ── 10. Assemble three panels flush ─────────────────────────────────────────
fig <- p_phylum + p_count + p_utr +
  plot_layout(widths = c(0.08, 0.55, 1), nrow = 1) +
  plot_annotation(
    title    = "+ssRNA virome library \u2014 family distribution and 3\u2032 UTR lengths",
    subtitle = glue("{sum(family_sum$n_lib)} accessions  \u00b7  ",
                    "{sum(family_sum$n_utr)} with 3\u2032UTR \u2265 {MIN_UTR_LEN} nt  \u00b7  ",
                    "{n_fam} families  \u00b7  ",
                    "{sum(phylum_rank[family_sum$phylum] <= 3, na.rm=TRUE) > 0 && 3} phyla"),
    caption  = glue(
      "Grey bar = total accessions per family; coloured bar = with 3\u2032UTR \u2265 {MIN_UTR_LEN} nt ",
      "(annotated + inferred; {n_short} sequences < {MIN_UTR_LEN} nt excluded). ",
      "Box = IQR, line = median, points = individual sequences."),
    theme = theme(
      plot.title    = element_text(size = 10, face = "bold", colour = "grey10"),
      plot.subtitle = element_text(size = 7.5, colour = "grey45"),
      plot.caption  = element_text(size = 6.5, colour = "grey55", margin = margin(t = 4)),
      plot.margin   = margin(8, 6, 6, 6))
  )

ggsave(OUT_FIGURE, fig, width = FIGURE_W, height = FIGURE_H)
message("Figure saved \u2192 ", OUT_FIGURE)
message("Summary CSV  \u2192 ", OUT_CSV)
