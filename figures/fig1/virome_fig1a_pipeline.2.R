# =============================================================================
# Virome library — Figure 1A-style, 3 panels
#
# Panel 1: Species count bars (grey = library total, colour = with 3'UTR >= MIN)
# Panel 2: 3'UTR + intergenic length distributions (boxplots, embedded legend)
# Panel 3: Detected accessions per family across 3 replicates (stacked bars)
#
# install.packages(c("tidyverse","scales","glue","cowplot"))
# =============================================================================

library(tidyverse)
library(scales)
library(glue)
library(cowplot)

# ── 0. Config ────────────────────────────────────────────────────────────────
FASTA_FILE    <- "pilot_library_20260319.uniqueIDs.fa"
CSV_FILE      <- "virus_stats.csv"
REP1_FILE     <- "1_peaks_per_accession.csv"
REP2_FILE     <- "2_peaks_per_accession.csv"
REP3_FILE     <- "3_peaks_per_accession.csv"
OUT_FIGURE    <- "fig1a_style.pdf"
OUT_CSV       <- "virome_summary.csv"
MIN_UTR_LEN   <- 100
MIN_IG_LEN    <- 100
BAR_WIDTH     <- 0.7
POINT_ALPHA   <- 0.35
FIGURE_W      <- 14
FIGURE_H      <- 8

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

utr_raw <- filter(fasta_df, grepl("3UTR", feature, ignore.case = TRUE))
n_short_utr <- sum(utr_raw$utr_length < MIN_UTR_LEN)
utr_df  <- filter(utr_raw, utr_length >= MIN_UTR_LEN) %>% select(acc_clean, utr_length)

ig_raw  <- filter(fasta_df, grepl("intergenic", feature, ignore.case = TRUE))
n_short_ig <- sum(ig_raw$utr_length < MIN_IG_LEN)
ig_df   <- filter(ig_raw, utr_length >= MIN_IG_LEN) %>%
           select(acc_clean, ig_length = utr_length)

message(glue("  {nrow(utr_df)} 3'UTR >= {MIN_UTR_LEN} nt | {nrow(ig_df)} intergenic >= {MIN_IG_LEN} nt"))

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

# ── 3. Library summaries ──────────────────────────────────────────────────────
all_df <- tibble(acc_clean = all_accs) %>% inner_join(tax_df, by = "acc_clean")

lib_counts <- all_df %>%
  group_by(family) %>%
  summarise(n_lib  = n_distinct(acc_clean),
            phylum = names(sort(table(phylum), decreasing = TRUE))[1],
            .groups = "drop")

utr_joined <- utr_df %>%
  inner_join(tax_df %>% select(acc_clean, family), by = "acc_clean") %>%
  filter(!is.na(family))

ig_joined <- ig_df %>%
  inner_join(tax_df %>% select(acc_clean, family), by = "acc_clean") %>%
  filter(!is.na(family))

colour_counts <- utr_joined %>%
  group_by(family) %>% summarise(n_utr = n_distinct(acc_clean), .groups = "drop")

family_sum <- lib_counts %>%
  left_join(colour_counts, by = "family") %>%
  replace_na(list(n_utr = 0))

write_csv(family_sum, OUT_CSV)
n_phyla <- n_distinct(family_sum$phylum[family_sum$phylum != "Unannotated"])

# ── 4. Detection data ─────────────────────────────────────────────────────────
read_rep <- function(path, rep_id) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(replicate = rep_id,
           acc_clean = sub("\\.\\d+$", "", sub("\\|.*", "", accession)),
           family    = if_else(is.na(family) | family == "Unannotated",
                               "Unannotated", family))
}

det_all <- bind_rows(
  read_rep(REP1_FILE, 1),
  read_rep(REP2_FILE, 2),
  read_rep(REP3_FILE, 3)
)

# Count replicates each accession appears in, then summarise per family
det_summary <- det_all %>%
  group_by(family, accession) %>%
  summarise(n_reps = n_distinct(replicate), .groups = "drop") %>%
  group_by(family, n_reps) %>%
  summarise(n_accs = n(), .groups = "drop") %>%
  mutate(rep_label = factor(n_reps, levels = c(1, 2, 3),
                            labels = c("1 replicate", "2 replicates", "3 replicates")))

# Write detection summary CSV
det_wide <- det_summary %>%
  pivot_wider(names_from = n_reps, values_from = n_accs,
              names_prefix = "n_reps_", values_fill = 0) %>%
  mutate(total_detected = rowSums(across(starts_with("n_reps_"))))
write_csv(det_wide, "detection_summary.csv")
message(glue("Detection summary written: {nrow(det_wide)} families detected"))

# ── 5. Colours ───────────────────────────────────────────────────────────────
phylum_colors <- c(
  "Pisuviricota"    = "#4472C4",
  "Kitrinoviricota" = "#C9542A",
  "Lenarviricota"   = "#70AD47",
  "Unannotated"     = "#999999"
)

rep_colors <- c(
  "1 replicate"  = "#CCCCCC",
  "2 replicates" = "#888888",
  "3 replicates" = "#333333"
)

utr_ig_colors <- c(
  "3\u2032 UTR"     = "#E07B54",
  "Intergenic" = "#6BAED6"
)

# ── 6. Factor order ───────────────────────────────────────────────────────────
phylum_rank <- c("Pisuviricota"=1,"Kitrinoviricota"=2,"Lenarviricota"=3,"Unannotated"=99)

fam_plot <- family_sum %>%
  mutate(phylum_idx = coalesce(phylum_rank[phylum], 50L)) %>%
  arrange(phylum_idx, n_lib) %>%
  mutate(family = fct_inorder(factor(family)))

fam_levels <- levels(fam_plot$family)
n_fam      <- length(fam_levels)

# Add phylum to obs data frames
obs_utr <- utr_joined %>%
  mutate(family = factor(family, levels = fam_levels)) %>%
  left_join(fam_plot %>% select(family, phylum), by = "family") %>%
  filter(!is.na(family)) %>%
  mutate(region = "3\u2032 UTR", length = utr_length)

obs_ig <- ig_joined %>%
  mutate(family = factor(family, levels = fam_levels)) %>%
  left_join(fam_plot %>% select(family, phylum), by = "family") %>%
  filter(!is.na(family)) %>%
  mutate(region = "Intergenic", length = ig_length)

obs_both <- bind_rows(obs_utr, obs_ig) %>%
  mutate(region = factor(region, levels = c("3\u2032 UTR", "Intergenic")))

det_plot <- det_summary %>%
  mutate(family = factor(family, levels = fam_levels)) %>%
  left_join(fam_plot %>% select(family, phylum), by = "family") %>%
  filter(!is.na(family))

# ── 7. Shared theme ───────────────────────────────────────────────────────────
theme_fig <- function(show_y = FALSE) {
  theme_minimal(base_size = 7.5) +
  theme(
    axis.text.y        = if (show_y) element_text(size = 6.5, face = "italic",
                                                   colour = "grey20", hjust = 1)
                         else element_blank(),
    axis.ticks.y       = element_blank(),
    axis.text.x        = element_text(size = 6, colour = "grey40"),
    axis.title.x       = element_text(size = 7, colour = "grey40",
                                      margin = margin(t = 3)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.2),
    plot.margin        = margin(4, 4, 4, 4)
  )
}

# ── 8. Panel 1: species count ─────────────────────────────────────────────────
phylum_brackets <- fam_plot %>%
  mutate(yi = as.integer(family)) %>%
  group_by(phylum) %>%
  summarise(y_lo = min(yi) - 0.45, y_hi = max(yi) + 0.45,
            y_mid = mean(range(yi)), .groups = "drop")

p_count <- ggplot(fam_plot, aes(y = family)) +
  geom_col(aes(x = n_lib), fill = "grey78", width = BAR_WIDTH) +
  geom_col(aes(x = n_utr, fill = phylum), width = BAR_WIDTH, show.legend = FALSE) +
  geom_segment(data = phylum_brackets,
               aes(x = max(fam_plot$n_lib) * 1.06,
                   xend = max(fam_plot$n_lib) * 1.06,
                   y = y_lo, yend = y_hi, colour = phylum),
               linewidth = 1.5, show.legend = FALSE) +
  geom_text(data = phylum_brackets,
            aes(x = max(fam_plot$n_lib) * 1.18, y = y_mid,
                label = phylum, colour = phylum),
            angle = 90, size = 1.8, hjust = 0.5, fontface = "bold",
            show.legend = FALSE) +
  scale_fill_manual(values = phylum_colors) +
  scale_colour_manual(values = phylum_colors) +
  scale_y_discrete(limits = fam_levels) +
  scale_x_continuous(limits = c(0, max(fam_plot$n_lib) * 1.25),
                     breaks = c(0, 50, 100, 150),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Species count", y = NULL) +
  theme_fig(show_y = TRUE)

# ── 9. Panel 2: UTR + intergenic lengths ──────────────────────────────────────
p_lengths <- ggplot(obs_both, aes(y = family, x = length, fill = region)) +
  geom_boxplot(aes(colour = region), outlier.shape = NA,
               width = BAR_WIDTH * 0.45, linewidth = 0.25,
               position = position_dodge(width = BAR_WIDTH),
               show.legend = TRUE) +
  geom_point(aes(colour = region),
             shape = 16, size = 0.4, alpha = POINT_ALPHA,
             position = position_jitterdodge(jitter.height = 0.15,
                                             dodge.width = BAR_WIDTH),
             show.legend = FALSE) +
  scale_fill_manual(values = utr_ig_colors, name = NULL) +
  scale_colour_manual(values = utr_ig_colors, name = NULL) +
  scale_y_discrete(limits = fam_levels) +
  scale_x_continuous(limits = c(0, NA),
                     breaks = c(0, 500, 1000, 1500),
                     labels = c("0", "500", "1,000", "1,500 nt"),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Length (nt)", y = NULL) +
  theme_fig(show_y = FALSE) +
  theme(
    legend.position      = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.text          = element_text(size = 6),
    legend.key.size      = unit(0.45, "lines"),
    legend.background    = element_rect(fill = "white", colour = NA),
    legend.margin        = margin(2, 4, 2, 4)
  )

# ── 10. Panel 3: detection across replicates ──────────────────────────────────
p_detect <- ggplot(det_plot,
                   aes(y = family, x = n_accs, fill = rep_label)) +
  geom_col(width = BAR_WIDTH, position = "stack") +
  scale_fill_manual(values = rep_colors, name = "Detected in") +
  scale_y_discrete(limits = fam_levels) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = pretty_breaks(4)) +
  labs(x = "Accessions detected", y = NULL) +
  theme_fig(show_y = FALSE) +
  theme(
    legend.position      = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.text          = element_text(size = 6),
    legend.key.size      = unit(0.45, "lines"),
    legend.background    = element_rect(fill = "white", colour = NA),
    legend.margin        = margin(2, 4, 2, 4)
  )

# ── 11. Align and combine ─────────────────────────────────────────────────────
plots <- cowplot::align_plots(p_count, p_lengths, p_detect,
                              align = "v", axis = "lr")

fig <- cowplot::plot_grid(
  plots[[1]], plots[[2]], plots[[3]],
  nrow       = 1,
  rel_widths = c(1, 1.1, 0.8)
)

final <- cowplot::plot_grid(
  cowplot::ggdraw() + cowplot::draw_label(
    glue("+ssRNA virome library \u2014 family distribution, feature lengths, and detection"),
    fontface = "bold", size = 9, x = 0.01, hjust = 0, colour = "grey10"),
  cowplot::ggdraw() + cowplot::draw_label(
    glue("{sum(family_sum$n_lib)} accessions  \u00b7  {nrow(family_sum)} families  \u00b7  ",
         "{n_phyla} phyla  \u00b7  ",
         "{n_distinct(det_all$accession)} unique accessions detected across 3 replicates"),
    size = 7, x = 0.01, hjust = 0, colour = "grey45"),
  fig,
  cowplot::ggdraw() + cowplot::draw_label(
    glue("Panel 1: grey = total in library, colour = with 3\u2032UTR \u2265 {MIN_UTR_LEN} nt. ",
         "Panel 2: box = IQR, line = median ({n_short_utr} UTR / {n_short_ig} intergenic < {MIN_UTR_LEN} nt excluded). ",
         "Panel 3: stacked by number of replicates detected."),
    size = 6, x = 0.01, hjust = 0, colour = "grey55"),
  ncol        = 1,
  rel_heights = c(0.04, 0.03, 1, 0.04)
)

ggsave(OUT_FIGURE, final, width = FIGURE_W, height = FIGURE_H, device = cairo_pdf)
message("Saved \u2192 ", OUT_FIGURE)
message("Detection summary \u2192 detection_summary.csv")
