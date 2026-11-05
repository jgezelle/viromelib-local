## -----------------------------------------------------------------
## Norm. Perc. Resistance vs Mean SF, by organism, with trendline +
## shaded SD/SE ribbon and correlation
## -----------------------------------------------------------------

## install.packages(c("ggplot2","dplyr","ggrepel"))  # run once if needed
library(ggplot2)
library(dplyr)
library(ggrepel)

## 1. Load & tidy -----------------------------------------------------
df <- read.csv("cy5_bandquant_plot_-_plot.csv", stringsAsFactors = FALSE)

names(df) <- c("replicate", "Organism", "Accession",
               "mean_sf", "perc_resistance", "norm_perc_resistance")

df$replicate <- factor(df$replicate)
df$Organism  <- factor(df$Organism)

## 2. Correlation (across all replicate points) -----------------------
ct <- cor.test(df$mean_sf, df$norm_perc_resistance, method = "pearson")
r_val <- round(unname(ct$estimate), 3)
p_val <- signif(ct$p.value, 3)
cor_label <- paste0("Pearson r = ", r_val, "\np = ", p_val)

## 3. Per-organism mean +/- SD (from the 2 replicates) -----------------
summary_df <- df %>%
  group_by(Organism, mean_sf) %>%
  summarise(
    mean_resist = mean(norm_perc_resistance),
    sd_resist   = sd(norm_perc_resistance),
    .groups = "drop"
  )

## 4. Plot ---------------------------------------------------------------
p <- ggplot(df, aes(x = mean_sf, y = norm_perc_resistance)) +

  # trendline + shaded standard-error ribbon (fit across all replicate points)
  geom_smooth(
    method = "lm",              # swap to "loess" for a curved fit
    se = TRUE,
    color = "grey20",
    fill = "steelblue",
    alpha = 0.2,                # faded/semi-transparent band
    linewidth = 0.8
  ) +

  # optional: per-organism SD error bars around the replicate mean
  geom_errorbar(
    data = summary_df,
    aes(x = mean_sf, y = mean_resist,
        ymin = mean_resist - sd_resist, ymax = mean_resist + sd_resist),
    width = 0.015, color = "grey40", inherit.aes = FALSE
  ) +

  # raw replicate points, colored by organism
  geom_point(aes(color = Organism), size = 3, alpha = 0.9) +

  # organism labels (repelled so the two replicate points don't overlap text)
  geom_text_repel(
    data = summary_df,
    aes(x = mean_sf, y = mean_resist, label = Organism),
    size = 3.2, color = "grey20", max.overlaps = Inf
  ) +

  # correlation annotation
  annotate(
    "text",
    x = min(df$mean_sf) + 0.02,
    y = max(df$norm_perc_resistance) * 0.98,
    label = cor_label,
    hjust = 0, vjust = 1, size = 3.6, fontface = "italic"
  ) +

  labs(
    x = "Mean SF",
    y = "Normalized % Resistance",
    color = "Organism",
    title = "Normalized % Resistance vs Mean SF",
    subtitle = "Points = individual replicates; error bars = SD of the 2 replicates; band = SE of trendline"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(p)

## 5. Save to file ---------------------------------------------------
ggsave("resistance_vs_sf.png", p, width = 8, height = 6, dpi = 300)
