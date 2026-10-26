#!/usr/bin/env Rscript
# 05_plot_5p_sbcnv5p.R: figure of RNA 5' ends (A3480) on the SbCNV-5 genome. See SbCNV_5prime.md §6.
#
# Usage (any machine with R; needs ggplot2, patchwork, ggrepel, svglite):
#   Rscript 05_plot_5p_sbcnv5p.R \
#     --counts 260930-scnv-blue-5end-A3480.5p_counts.tsv \   (or the plus-strand 5'-end .bedGraph)
#     --subs   substitutions.tsv \
#     --ann    sbcnv5_annotations.tsv \
#     --out    SbCNV5_5prime_ends.svg
# Also writes <out>_linear_break.svg/.pdf: linear read axis with one scale break, overview + element 1/2 zooms.
# All arguments are optional (defaults below). From R / R Markdown, set them with options() before source():
#   options(sbcnv5.counts = "...", sbcnv5.subs = "...", sbcnv5.ann = "...", sbcnv5.out = "plots/x.svg")
# Writes the .svg, a .pdf, and two tables:
#   <out>_peaks.tsv   every 5' end with reads, %, product length, nearest feature
#   <out>_stalls.tsv  each assay stall site: reads at the site and at +/-1, +/-2, nearest major 5' end
#
# What is plotted:
#   on-target reads only (R2 = RT primer), + strand, 5' ends <= 18760 (ends past that are mostly primer),
#   counts from the patched-reference alignment (MAPQ >= 10). Read counts on a log10 axis.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(ggrepel) })

# ---------------- settings ----------------
args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(k, d) {                       # command line (--k value) > options(sbcnv5.k = ...) > default
  i <- match(paste0("--", k), args); if (!is.na(i)) return(args[i + 1])
  v <- getOption(paste0("sbcnv5.", k)); if (is.null(v)) d else v
}
counts_f   <- getarg("counts", "260930-scnv-blue-5end-A3480.5p_counts.tsv")
subs_f     <- getarg("subs",   "substitutions.tsv")
ann_f      <- getarg("ann",    "sbcnv5_annotations.tsv")
out_svg    <- getarg("out",    "SbCNV5_5prime_ends.svg")
GENOME_LEN <- as.numeric(getarg("genome_len", 19199))
PRIMER_END <- 18798            # 3' end of the RT primer site: product length = PRIMER_END - end5 + 1
MAX_END    <- 18760            # drop ends whose R1 would be mostly primer sequence
WIN_B      <- c(17900, 19199)  # zoom 1: all products through the genome end (past the primer = not measurable)
WIN_C      <- c(18555, 18700)  # zoom 2: elements 1 and 2 at single-nt resolution
LABEL_B    <- 0.8              # label peaks >= this % of on-target reads in zoom 1 (one label per 10-nt cluster)
LABEL_C    <- 1.0              # label every position >= this % in zoom 2

# colors (validated: CVD-safe, labels carry identity; element fill is a tint behind labeled boxes)
col <- list(ink = "#0b0b0b", ink2 = "#52514e", muted = "#8a8984", grid = "#e7e6e2",
            stem = "#2a78d6", element = "#1baf7a", stall = "#eb6834", primer = "#4a3aa7", cds = "#b9b8b2", nm = "#f0efeb")
base_family <- "Helvetica"

# ---------------- data ----------------
if (grepl("\\.bedgraph$", counts_f, ignore.case = TRUE)) {
  # 5'-end bedGraph (0-based start, end = nt): one line per nt, value = on-target 5' ends there
  bg <- read.table(counts_f, sep = "\t", comment.char = "", skip = sum(grepl("^(track|browser|#)", readLines(counts_f, 5))))
  cn <- data.frame(chrom = bg[[1]], pos = bg[[3]], strand = "+", n_all = bg[[4]], n_ontarget = bg[[4]])
} else cn <- read.delim(counts_f)      # 5p_counts.tsv from step 04 (chrom, pos, strand, n_all, n_ontarget)
cn <- subset(cn, strand == "+" & n_ontarget > 0 & pos <= MAX_END)
total <- sum(cn$n_ontarget)
cn$pct     <- 100 * cn$n_ontarget / total
cn$product <- PRIMER_END - cn$pos + 1

subs <- if (file.exists(subs_f)) read.delim(subs_f) else data.frame(pos = numeric(0))
ann  <- read.delim(ann_f, comment.char = "#", stringsAsFactors = FALSE)
ann$name[is.na(ann$name)] <- ""
ann$type <- tolower(trimws(ann$type))
el  <- subset(ann, type == "element"); cds <- subset(ann, type == "cds"); utr <- subset(ann, type == "utr")
pr  <- subset(ann, type == "primer");  st  <- subset(ann, type == "stall")
st$high <- tolower(st$name) == "high"
gen <- rbind(transform(cds, kind = "CDS"), transform(utr, kind = "UTR"))

# ---------------- peak labels ----------------
pick_labels <- function(d, win, thr, cluster = 10) {
  d <- d[d$pos >= win[1] & d$pos <= win[2] & d$pct >= thr, ]
  d <- d[order(-d$n_ontarget), ]
  keep <- integer(0)
  for (i in seq_len(nrow(d))) if (cluster == 0 || all(abs(d$pos[i] - d$pos[keep]) > cluster)) keep <- c(keep, i)
  d[keep, ]
}
lab_B <- pick_labels(cn, WIN_B, LABEL_B, cluster = 10)
lab_C <- pick_labels(cn, WIN_C, LABEL_C, cluster = 0)

# ---------------- themes / helpers ----------------
th <- theme_classic(base_size = 8, base_family = base_family) +
  theme(axis.line = element_line(linewidth = 0.3, colour = col$ink2),
        axis.ticks = element_line(linewidth = 0.3, colour = col$ink2),
        axis.text = element_text(colour = col$ink2, size = 7),
        axis.title = element_text(colour = col$ink, size = 8),
        panel.grid.major.y = element_line(linewidth = 0.25, colour = col$grid),
        plot.title = element_text(size = 9, face = "bold", colour = col$ink, hjust = 0),
        plot.subtitle = element_text(size = 7, colour = col$ink2, hjust = 0),
        plot.margin = margin(4, 8, 2, 4))
th_track <- th + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                       panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = col$ink, size = 7))
kfmt  <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
rfmt  <- function(x) ifelse(x >= 1000, paste0(x / 1000, "k"), as.character(x))
Y0    <- 0.6                     # stem baseline on the log axis, so 1-read positions stay visible
ymax  <- 10^ceiling(log10(max(cn$n_ontarget) * 1.5))

stem_panel <- function(win, labs, lw, title, subtitle = NULL, show_x = TRUE) {
  d <- cn[cn$pos >= win[1] & cn$pos <= win[2], ]
  e <- el[el$end >= win[1] & el$start <= win[2], ]
  s <- st[st$start >= win[1] & st$start <= win[2], , drop = FALSE]
  p <- ggplot() +
    { if (win[2] > MAX_END) annotate("rect", xmin = MAX_END + 0.5, xmax = win[2] + 0.5, ymin = Y0, ymax = ymax,
                                     fill = col$nm) } +
    { if (win[2] > PRIMER_END + 150) annotate("text", x = (PRIMER_END + win[2]) / 2, y = ymax / 4, size = 2.4,
                                              family = base_family, colour = col$ink2, lineheight = 0.9,
                                              label = "not measured\n(RT primer and downstream)") } +
    { if (nrow(s)) geom_vline(data = s, aes(xintercept = start), colour = col$stall, linewidth = 0.3,
                              linetype = "22") } +
    { if (nrow(e)) geom_rect(data = e, aes(xmin = pmax(start, win[1]) - 0.5, xmax = pmin(end, win[2]) + 0.5,
                                          ymin = Y0, ymax = ymax), fill = col$element, alpha = 0.10) } +
    geom_segment(data = d, aes(x = pos, xend = pos, y = Y0, yend = n_ontarget),
                 colour = col$stem, linewidth = lw, lineend = "butt") +
    geom_point(data = labs, aes(pos, n_ontarget), colour = col$stem, size = 0.9) +
    geom_text_repel(data = labs, aes(pos, n_ontarget, label = kfmt(pos)), family = base_family,
                    size = 2.5, colour = col$ink, nudge_y = 0.25, direction = "both", min.segment.length = 0.1,
                    segment.size = 0.2, segment.colour = col$muted, box.padding = 0.35, max.overlaps = Inf, seed = 1) +
    scale_y_log10(limits = c(Y0, ymax), breaks = 10^(0:6), labels = rfmt, expand = c(0, 0)) +
    scale_x_continuous(limits = c(win[1] - 0.5, win[2] + 0.5), breaks = scales::breaks_pretty(8), labels = kfmt, expand = c(0, 0)) +
    labs(x = NULL, y = "5' ends (reads)", title = title, subtitle = subtitle) + th
  if (!show_x) p <- p + theme(axis.text.x = element_blank())
  p
}

track_panel <- function(win, show_axis = TRUE) {
  lanes <- c("CDS / UTR", "Element", "Stall (assay)", "RT primer", "Strain variant")
  has <- c(nrow(gen) > 0, nrow(el) > 0, nrow(st) > 0, nrow(pr) > 0, nrow(subs) > 0)
  inwin <- function(d) d[d$end >= win[1] & d$start <= win[2], , drop = FALSE]
  fit <- function(name, x0, x1) {                    # shorten labels that would overflow their box
    w <- (x1 - x0) / diff(win); short <- sub(" \\(.*", "", name)
    abbr <- ifelse(grepl("^Element", name), paste0("E", gsub("\\D", "", name)), "")
    ifelse(w >= nchar(name) * 0.009, name, ifelse(w >= nchar(short) * 0.009, short, ifelse(w >= nchar(abbr) * 0.009, abbr, "")))
  }
  has <- has & c(nrow(inwin(gen)) > 0, nrow(inwin(el)) > 0, nrow(inwin(st)) > 0, nrow(inwin(pr)) > 0,
                 sum(subs$pos >= win[1] & subs$pos <= win[2]) > 0)
  lanes <- lanes[has]; y <- setNames(rev(seq_along(lanes)), lanes)
  p <- ggplot() + scale_y_continuous(breaks = y, labels = names(y), limits = c(0.4, length(y) + 0.6), expand = c(0, 0))
  if ("CDS / UTR" %in% lanes) { d <- inwin(gen); yy <- y["CDS / UTR"]
    d$x0 <- pmax(d$start, win[1]) - 0.5; d$x1 <- pmin(d$end, win[2]) + 0.5
    p <- p + geom_rect(data = d, aes(xmin = x0, xmax = x1, ymin = yy - ifelse(kind == "CDS", 0.3, 0.18),
                                     ymax = yy + ifelse(kind == "CDS", 0.3, 0.18)),
                       fill = ifelse(d$kind == "CDS", col$cds, "white"), colour = col$cds, linewidth = 0.4) +
      geom_text(data = d, aes(x = (x0 + x1) / 2, y = yy, label = fit(name, x0, x1)), size = 2.4, family = base_family, colour = col$ink) }
  if ("Element" %in% lanes && nrow(inwin(el))) { d <- inwin(el)
    d$x0 <- pmax(d$start, win[1]) - 0.5; d$x1 <- pmin(d$end, win[2]) + 0.5
    p <- p + geom_rect(data = d, aes(xmin = x0, xmax = x1,
                                     ymin = y["Element"] - 0.3, ymax = y["Element"] + 0.3), fill = col$element, alpha = 0.55) +
      geom_text(data = d, aes(x = (x0 + x1) / 2, y = y["Element"], label = fit(name, x0, x1)),
                size = 2.4, family = base_family, colour = col$ink) }
  if ("Stall (assay)" %in% lanes && nrow(inwin(st))) { d <- inwin(st); d <- d[order(d$high), ]
    p <- p + geom_segment(data = d, aes(x = start, xend = end, y = y["Stall (assay)"], yend = y["Stall (assay)"]),
                          colour = col$stall, linewidth = 0.6) +
      geom_point(data = d, aes(x = (start + end) / 2, y = y["Stall (assay)"]), shape = 25, size = 1.9,
                 fill = ifelse(d$high, col$stall, "white"), colour = col$stall, stroke = 0.5) }
  if ("RT primer" %in% lanes && nrow(inwin(pr))) { d <- inwin(pr)
    p <- p + geom_segment(data = d, aes(x = end, xend = start, y = y["RT primer"], yend = y["RT primer"]),
                          colour = col$primer, linewidth = 0.8, arrow = arrow(length = unit(1.2, "mm"), type = "closed")) }
  if ("Strain variant" %in% lanes) { d <- subs[subs$pos >= win[1] & subs$pos <= win[2], , drop = FALSE]
    if (nrow(d)) p <- p + geom_segment(data = d, aes(x = pos, xend = pos, y = y["Strain variant"] - 0.25,
                                                   yend = y["Strain variant"] + 0.25), colour = col$muted, linewidth = 0.4) }
  p <- p + scale_x_continuous(limits = c(win[1] - 0.5, win[2] + 0.5), breaks = scales::breaks_pretty(8), labels = kfmt, expand = c(0, 0)) +
    labs(x = if (show_axis) "Position in SbCNV-5 genome (nt)" else NULL, y = NULL) + th_track
  p
}

# ---------------- panel A: genome map ----------------
pA <- ggplot() +
  geom_segment(aes(x = 1, xend = GENOME_LEN, y = 0, yend = 0), colour = col$ink2, linewidth = 0.6) +
  { if (nrow(cds)) geom_rect(data = cds, aes(xmin = start, xmax = end, ymin = -0.35, ymax = 0.35), fill = col$cds) } +
  { if (nrow(cds)) geom_text(data = cds, aes(x = (start + end) / 2, y = 0, label = name), size = 2.4,
                             family = base_family, colour = col$ink) } +
  { if (nrow(utr)) geom_rect(data = utr, aes(xmin = start, xmax = end, ymin = -0.2, ymax = 0.2), fill = "white",
                             colour = col$cds, linewidth = 0.4) } +
  { if (nrow(utr)) geom_text(data = utr, aes(x = (start + end) / 2, y = -0.72, label = name), size = 2.3,
                             family = base_family, colour = col$ink2) } +
  annotate("rect", xmin = WIN_B[1], xmax = WIN_B[2], ymin = -0.45, ymax = 0.6, fill = NA, colour = col$stem, linewidth = 0.5) +
  annotate("text", x = WIN_B[1] - 150, y = 0.95, label = "B", hjust = 1, size = 2.6, colour = col$stem,
           family = base_family, fontface = "bold") +
  scale_x_continuous(limits = c(0, GENOME_LEN * 1.005), breaks = seq(0, GENOME_LEN, 2000),
                     labels = function(x) paste0(x / 1000, " kb"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.85, 1.3)) +
  labs(x = NULL, y = NULL, title = "A  SbCNV-5 genome (NC_024077, 19,199 nt)") +
  th + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
             panel.grid.major.y = element_blank())

# ---------------- panels B, C ----------------
sub_B <- sprintf("%s on-target reads; RT primer at %s-%s; product length = %s - 5' end + 1",
                 kfmt(total), kfmt(min(pr$start, PRIMER_END - 18)), kfmt(PRIMER_END), kfmt(PRIMER_END))
pB  <- stem_panel(WIN_B, lab_B, 0.25, "B  RNA 5' ends, 3' end of the genome", sub_B, show_x = FALSE)
tB  <- track_panel(WIN_B, show_axis = TRUE)
pC  <- stem_panel(WIN_C, lab_C, 0.9, "C  Elements 1 and 2, single-nucleotide resolution", NULL, show_x = FALSE)
tC  <- track_panel(WIN_C, show_axis = TRUE)

fig <- pA / pB / tB / pC / tC + plot_layout(heights = c(0.55, 3, 1, 2.4, 1))
w_in <- 7.2; h_in <- 8.4
dir.create(dirname(out_svg), recursive = TRUE, showWarnings = FALSE)
ggsave(out_svg, fig, width = w_in, height = h_in, device = svglite::svglite)
ggsave(sub("\\.svg$", ".pdf", out_svg), fig, width = w_in, height = h_in, device = cairo_pdf)

# ---------------- figure 2: linear read axis with a scale break ----------------
# Same counts, linear scale, so the size differences between the dominant 5' ends read directly.
# The y axis is broken once: a lower segment up to just above the 2nd-largest end, and an upper segment around
# the largest (18681). Both segments use the SAME reads-per-mm scale (panel heights are proportional to their ranges).
WIN_L  <- c(17900, MAX_END)                  # region for the overview
ZOOM_L <- list(c(18562, 18586), c(18668, 18690))  # element 1 and element 2 starts, single-nt bars
LAB_L  <- 1000                               # label every 5' end with >= this many reads
dl <- cn[cn$pos >= WIN_L[1] & cn$pos <= WIN_L[2], ]
srt <- sort(dl$n_ontarget, decreasing = TRUE)
lo_max <- ceiling(srt[2] * 1.12 / 1000) * 1000                 # lower segment: 0 .. just above 2nd-largest
hi_rng <- c(floor(srt[1] * 0.93 / 1000) * 1000, ceiling(srt[1] * 1.10 / 1000) * 1000)  # upper segment
brk_lo <- pretty(c(0, lo_max), 5); brk_lo <- brk_lo[brk_lo <= lo_max]
brk_hi <- pretty(hi_rng, 3);       brk_hi <- brk_hi[brk_hi >= hi_rng[1] & brk_hi <= hi_rng[2]]

broken <- function(win, bars, labs, title, lab_fmt) {
  d <- cn[cn$pos >= win[1] & cn$pos <= win[2], ]
  e <- el[el$end >= win[1] & el$start <= win[2], ]
  s <- st[st$start >= win[1] & st$start <= win[2], , drop = FALSE]
  xs <- scale_x_continuous(limits = c(win[1] - 0.5, win[2] + 0.5), breaks = scales::breaks_pretty(6),
                           labels = kfmt, expand = c(0, 0))
  geom_data <- if (bars) geom_col(data = d, aes(pos, n_ontarget), fill = col$stem, width = 0.75) else
    geom_segment(data = d, aes(x = pos, xend = pos, y = 0, yend = n_ontarget), colour = col$stem, linewidth = 0.35)
  bg <- list(if (nrow(e)) geom_rect(data = e, aes(xmin = pmax(start, win[1]) - 0.5, xmax = pmin(end, win[2]) + 0.5,
                                                  ymin = -Inf, ymax = Inf), fill = col$element, alpha = 0.10),
             if (nrow(s)) geom_vline(data = s, aes(xintercept = start), colour = col$stall, linewidth = 0.3, linetype = "22"))
  lab_layer <- function(dd) if (bars)
      geom_text(data = dd, aes(pos, n_ontarget, label = kfmt(n_ontarget)), angle = 90, hjust = -0.15,
                family = base_family, size = 2.3, colour = col$ink) else
      geom_text_repel(data = dd, aes(pos, n_ontarget, label = lab_fmt(pos, n_ontarget)),
      family = base_family, size = 2.4, colour = col$ink, lineheight = 0.9, direction = "x", nudge_y = 0,
      min.segment.length = 0.1, segment.size = 0.2, segment.colour = col$muted, box.padding = 0.3,
      max.overlaps = Inf, seed = 1, vjust = -0.4)
  if (bars) xs <- scale_x_continuous(limits = c(win[1] - 0.5, win[2] + 0.5), expand = c(0, 0),
                                     breaks = sort(labs$pos), labels = kfmt)   # tick + label at each labeled end
  # slanted break marks drawn across the y axis line
  mark <- function(yv, yr) { dx <- 0.015 * (diff(win) + 1); dy <- 0.08 * diff(yr)   # slash just inside the axis
    annotate("segment", x = win[1] - 0.5, xend = win[1] - 0.5 + dx, y = yv - dy, yend = yv + dy,
             colour = col$ink2, linewidth = 0.5) }
  top <- ggplot() + bg + geom_data + lab_layer(labs[labs$n_ontarget >= hi_rng[1], ]) + xs +
    scale_y_continuous(breaks = brk_hi, labels = kfmt, expand = c(0, 0)) +
    coord_cartesian(ylim = hi_rng) + mark(hi_rng[1], hi_rng) +
    labs(x = NULL, y = NULL, title = title) + th +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
          plot.margin = margin(4, 8, 0, 4))
  bot <- ggplot() + bg + geom_data + lab_layer(labs[labs$n_ontarget < hi_rng[1], ]) + xs +
    scale_y_continuous(breaks = brk_lo, labels = kfmt, expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, lo_max)) + mark(lo_max, c(0, lo_max)) +
    labs(x = "Position in SbCNV-5 genome (nt)", y = "5' ends (reads)") + th +
    theme(plot.margin = margin(2, 8, 4, 4), axis.title.y = element_text(hjust = 0.75)) +
    { if (bars) theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5)) }
  list(top = top, bot = bot)
}
fmt_pos  <- function(p, n) kfmt(p)
fmt_both <- function(p, n) paste0(kfmt(p), "\n", kfmt(n))
labs_ov <- pick_labels(transform(dl, pct = n_ontarget), WIN_L, LAB_L, cluster = 10)   # one label per 10-nt cluster
pO <- broken(WIN_L, FALSE, labs_ov, "A  5' ends, linear scale with axis break (17,900-18,760)", fmt_pos)
zl <- lapply(seq_along(ZOOM_L), function(i) {
  w <- ZOOM_L[[i]]; dz <- cn[cn$pos >= w[1] & cn$pos <= w[2] & cn$n_ontarget >= LAB_L, ]
  broken(w, TRUE, dz, paste0(c("B", "C")[i], "  ", c("Element 1 start", "Element 2 start")[i], " (reads above bars)"), fmt_both)
})
hr <- c(diff(hi_rng), lo_max) / lo_max               # equal reads-per-mm in both segments
fig2 <- wrap_plots(A = pO$top, B = pO$bot, C = zl[[1]]$top, D = zl[[2]]$top, E = zl[[1]]$bot, F = zl[[2]]$bot,
                   design = "AA\nBB\nCD\nEF") +
        plot_layout(heights = c(hr[1], hr[2], hr[1], hr[2]))
out2 <- sub("\\.svg$", "_linear_break.svg", out_svg)
ggsave(out2, fig2, width = 7.2, height = 7.4, device = svglite::svglite)
ggsave(sub("\\.svg$", ".pdf", out2), fig2, width = 7.2, height = 7.4, device = cairo_pdf)

# ---------------- tables ----------------
nearest_feature <- function(p) {
  f <- rbind(el, cds, utr, pr)
  if (!nrow(f)) return("")
  inside <- f[p >= f$start & p <= f$end, ]
  if (nrow(inside)) return(paste0(inside$name[1], " +", p - inside$start[1]))
  d <- ifelse(p < f$start, f$start - p, p - f$end); i <- which.min(d)
  paste0(d[i], " nt ", ifelse(p < f$start[i], "upstream of ", "downstream of "), f$name[i])
}
pk <- cn[order(-cn$n_ontarget), c("pos", "n_ontarget", "pct", "product")]
pk$pct <- round(pk$pct, 3); pk$feature <- vapply(pk$pos, nearest_feature, "")
base <- sub("\\.svg$", "", out_svg)
write.table(pk, paste0(base, "_peaks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (nrow(st)) {
  at <- function(p) { r <- cn$n_ontarget[cn$pos == p]; if (length(r)) r else 0 }
  so <- do.call(rbind, lapply(seq_len(nrow(st)), function(i) {
    p <- st$start[i]; cov <- p <= MAX_END
    w <- cn[abs(cn$pos - p) <= 5, ]                       # strongest 5' end within +/-5 nt
    top <- if (cov && nrow(w)) w[which.max(w$n_ontarget), ] else NULL
    data.frame(stall_pos = p, high_score = st$high[i], measurable = cov,
               reads_at = at(p), reads_pm1 = sum(vapply((p - 1):(p + 1), at, 0)),
               reads_pm2 = sum(vapply((p - 2):(p + 2), at, 0)),
               top_5p_within5 = if (is.null(top)) NA else top$pos,
               top_reads = if (is.null(top)) NA else top$n_ontarget,
               offset = if (is.null(top)) NA else top$pos - p)
  }))
  write.table(so, paste0(base, "_stalls.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
cat(sprintf("on-target 5' ends plotted: %s reads at %d positions\nwrote %s (+ .pdf, _linear_break.svg/.pdf, _peaks.tsv%s)\n",
            kfmt(total), nrow(cn), out_svg, if (nrow(st)) ", _stalls.tsv" else ""))
