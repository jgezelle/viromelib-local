# =============================================================================
# Virome library — Figure 1A-style
# install.packages(c("tidyverse","scales","glue","cowplot"))
# =============================================================================

library(tidyverse)
library(scales)
library(glue)
library(cowplot)

# ── 0. Config ────────────────────────────────────────────────────────────────
FASTA_FILE  <- "pilot_library_20260319.uniqueIDs.fa"
CSV_FILE    <- "virus_stats.csv"
OUT_FIGURE  <- "fig1a_style.pdf"
OUT_CSV     <- "virome_summary.csv"
MIN_UTR_LEN <- 100
BAR_WIDTH   <- 0.7
POINT_ALPHA <- 0.4
FIGURE_W    <- 9
FIGURE_H    <- 8

# ── 1. Parse FASTA ───────────────────────────────────────────────────────────
message("Parsing FASTA...")
lines   <- readLines(FASTA_FILE)
headers <- sub("^>", "", lines[startsWith(lines, ">")])
rx      <- "^([^|]+)\\|(\\d+)-(\\d+)_(.+)$"
matched <- regmatches(headers, regexec(rx, headers))

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

# ── 2. Taxonomy ──────────────────────────────────────────────────────────────
artifact_phyla <- c("Macronovirus", "Tobacco albetovirus 3")
tax_df <- read_csv(CSV_FILE, show_col_types = FALSE) %>%
  rename_with(tolower) %>%
  mutate(
    acc_clean = sub("\\.\\d+$", "", accession),
    phylum = if_else(phylum %in% artifact_phyla | is.na(phylum) |
                       phylum == "Unannotated", "Unannotated", phylum),
    family = if_else(is.na(family) | family == "Unannotated", "Unannotated", family)
  )

# ── 3. Summaries ─────────────────────────────────────────────────────────────
all_df <- tibble(acc_clean = all_accs) %>% inner_join(tax_df, by = "acc_clean")

lib_counts <- all_df %>%
  group_by(family) %>%
  summarise(n_lib  = n_distinct(acc_clean),
            phylum = names(sort(table(phylum), decreasing = TRUE))[1],
            .groups = "drop")

utr_joined <- utr_df %>%
  inner_join(tax_df %>% select(acc_clean, family), by = "acc_clean") %>%
  filter(!is.na(family))

colour_counts <- utr_joined %>%
  group_by(family) %>% summarise(n_utr = n_distinct(acc_clean), .groups = "drop")

family_sum <- lib_counts %>%
  left_join(colour_counts, by = "family") %>%
  replace_na(list(n_utr = 0))

write_csv(family_sum, OUT_CSV)
n_phyla <- n_distinct(family_sum$phylum[family_sum$phylum != "Unannotated"])
message(glue("{nrow(family_sum)} families | {sum(family_sum$n_lib)} accessions | {n_phyla} phyla"))

# ── 4. Colours ───────────────────────────────────────────────────────────────
phylum_colors <- c(
  "Pisuviricota"    = "#4472C4",
  "Kitrinoviricota" = "#C0D6DF",
  "Lenarviricota"   = "#28965A",
  "Unannotated"     = "#999999"
)

# ── 5. Factor order ───────────────────────────────────────────────────────────
phylum_rank <- c("Pisuviricota"=1,"Kitrinoviricota"=2,"Lenarviricota"=3,"Unannotated"=99)

fam_plot <- family_sum %>%
  mutate(phylum_idx = coalesce(phylum_rank[phylum], 50L)) %>%
  arrange(phylum_idx, n_lib) %>%
  mutate(family = fct_inorder(factor(family)))

fam_levels <- levels(fam_plot$family)
n_fam      <- length(fam_levels)

obs_df <- utr_joined %>%
  mutate(family = factor(family, levels = fam_levels)) %>%
  left_join(fam_plot %>% select(family, phylum), by = "family") %>%
  filter(!is.na(family))

# ── 6. Left panel: species count bars ─────────────────────────────────────────
# Clean horizontal bars, family labels as y axis text coloured by phylum group
# Use geom_col with discrete y — simplest correct approach

p_count <- ggplot(fam_plot, aes(y = family)) +
  geom_col(aes(x = n_lib), fill = "grey78", width = BAR_WIDTH) +
  geom_col(aes(x = n_utr, fill = phylum), width = BAR_WIDTH, show.legend = FALSE) +
  scale_fill_manual(values = phylum_colors) +
  scale_y_discrete(limits = fam_levels) +
  scale_x_continuous(
    limits   = c(0, max(fam_plot$n_lib) * 1.05),
    breaks   = c(0, 25, 50, 75, 100, 125, 150),
    expand   = expansion(mult = c(0, 0))
  ) +
  # phylum group bracket: coloured segment on right margin
  geom_segment(
    data = fam_plot %>%
      mutate(yi = as.integer(family)) %>%
      group_by(phylum) %>%
      summarise(y_lo = min(yi) - 0.45, y_hi = max(yi) + 0.45,
                y_mid = mean(range(yi)), .groups = "drop"),
    aes(x = max(fam_plot$n_lib) * 1.04,
        xend = max(fam_plot$n_lib) * 1.04,
        y = y_lo, yend = y_hi, colour = phylum),
    linewidth = 1.2, show.legend = FALSE
  ) +
  scale_colour_manual(values = phylum_colors) +
  labs(x = "Species count", y = NULL) +
  theme_minimal(base_size = 7.5) +
  theme(
    axis.text.y        = element_text(size = 6.5, face = "italic",
                                      colour = "grey20", hjust = 1),
    axis.ticks.y       = element_blank(),
    axis.text.x        = element_text(size = 6, colour = "grey40"),
    axis.title.x       = element_text(size = 7, colour = "grey40",
                                      margin = margin(t = 3)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.2),
    plot.margin        = margin(4, 8, 4, 4)
  )

# ── 7. Right panel: UTR boxplot ───────────────────────────────────────────────
p_utr <- ggplot(obs_df, aes(y = family, x = utr_length, fill = phylum)) +
  geom_jitter(aes(colour = phylum), shape = 16, height = 0.2,
              size = 0.45, alpha = POINT_ALPHA, show.legend = FALSE) +
  geom_boxplot(outlier.shape = NA, width = BAR_WIDTH, linewidth = 0.25,
               colour = "grey25", show.legend = TRUE) +
  scale_fill_manual(
    values = phylum_colors, name = "Phylum",
    guide  = guide_legend(title.position = "top",
                          override.aes   = list(colour = NA, linetype = 0,
                                                alpha = 1, shape = NA))
  ) +
  scale_colour_manual(values = phylum_colors, guide = "none") +
  scale_y_discrete(limits = fam_levels) +
  scale_x_continuous(
    limits = c(0, NA),
    breaks = c(0, 500, 1000, 1500),
    labels = c("0", "500", "1,000", "1,500 nt"),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(x = "3\u2032 UTR length (nt)", y = NULL) +
  theme_minimal(base_size = 7.5) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.text.x        = element_text(size = 6, colour = "grey40"),
    axis.title.x       = element_text(size = 7, colour = "grey40",
                                      margin = margin(t = 3)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.2),
    legend.position    = "right",
    legend.text        = element_text(size = 6.5),
    legend.title       = element_text(size = 7, face = "bold"),
    legend.key.size    = unit(0.5, "lines"),
    plot.margin        = margin(4, 4, 4, 2)
  )

# ── 8. Combine with cowplot ────────────────────────────────────────────────────
# align="v" + axis="lr" syncs the discrete y axis heights across panels
plots <- cowplot::align_plots(p_count, p_utr, align = "v", axis = "lr")

fig <- cowplot::plot_grid(
  plots[[1]], plots[[2]],
  nrow       = 1,
  rel_widths = c(1, 1.1)
)

final <- cowplot::plot_grid(
  cowplot::ggdraw() + cowplot::draw_label(
    glue("+ssRNA virome library \u2014 family distribution and 3\u2032 UTR lengths"),
    fontface = "bold", size = 9, x = 0.01, hjust = 0, colour = "grey10"),
  cowplot::ggdraw() + cowplot::draw_label(
    glue("{sum(family_sum$n_lib)} accessions  \u00b7  ",
         "{sum(family_sum$n_utr)} with 3\u2032UTR \u2265 {MIN_UTR_LEN} nt  \u00b7  ",
         "{nrow(family_sum)} families  \u00b7  {n_phyla} phyla"),
    size = 7, x = 0.01, hjust = 0, colour = "grey45"),
  fig,
  cowplot::ggdraw() + cowplot::draw_label(
    glue("Grey = total per family; colour = with 3\u2032UTR \u2265 {MIN_UTR_LEN} nt ",
         "({n_short} < {MIN_UTR_LEN} nt excluded). ",
         "Box = IQR, line = median, points = individual sequences."),
    size = 6, x = 0.01, hjust = 0, colour = "grey55"),
  ncol        = 1,
  rel_heights = c(0.04, 0.03, 1, 0.04)
)

ggsave(OUT_FIGURE, final, width = FIGURE_W, height = FIGURE_H, device = cairo_pdf)
message("Saved \u2192 ", OUT_FIGURE)
