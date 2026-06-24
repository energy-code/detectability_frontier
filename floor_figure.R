# ============================================================================
# floor_figure.R  --  figure only (not the rest of the analysis)
#
# Draws the "possibility floor" from the same closed form:
#   (A) ggplot2 2-D filled-contour heatmap   -> floor_heatmap.png
#
# Floor:  Z_floor(beta, N) = K / (beta * sqrt(N * p))   [Gray], *1000 -> mGy
# This is the minimum RMS cumulative dose a study of N people can resolve for
# an effect of size beta. Below it, no study of that size can detect the
# effect, whatever the source.
# ============================================================================

## ---- packages -------------------------------------------------------------
# install.packages(c("ggplot2", "geomtextpath", "scales"))   # run once if needed
library(ggplot2)
library(geomtextpath)   # labels written into the contour lines

## ---- constants (identical to detectability_frontier.py) -------------------
ALPHA_1SIDED <- 0.05
POWER        <- 0.80
K            <- qnorm(1 - ALPHA_1SIDED) + qnorm(POWER)   # ~2.486
P_BASE       <- 0.20
N_US         <- 3.35e8
N_WORLD      <- 8e9
BETA_LNT     <- 0.5

## ---- the floor, in mGy ----------------------------------------------------
zfloor_mGy <- function(beta, N, p = P_BASE, k = K) {
  k / (beta * sqrt(N * p)) * 1e3
}

# sanity check -> should print 0.61 and 0.12
cat(sprintf("US floor    (beta=0.5): %.3f mGy\n", zfloor_mGy(BETA_LNT, N_US)))
cat(sprintf("world floor (beta=0.5): %.3f mGy\n", zfloor_mGy(BETA_LNT, N_WORLD)))

## ---- shared grid ----------------------------------------------------------
beta_seq <- 10^seq(log10(0.1),  log10(5.0),  length.out = 200)   # ERR / Gy
N_seq    <- 10^seq(log10(1e4),  log10(1e10), length.out = 200)   # people

# headline points to annotate
N_NPP <- 1e6   # ~ US population within 10 miles of a nuclear power plant

pts <- data.frame(
  N     = c(N_US, N_WORLD, N_NPP),
  beta  = c(BETA_LNT, BETA_LNT, BETA_LNT),
  floor = zfloor_mGy(BETA_LNT, c(N_US, N_WORLD, N_NPP)),
  label = c("US: 0.61 mGy", "world: 0.12 mGy", "10-mi NPP: 11 mGy"),
  hjust = c(-0.08, -0.08, -0.08),
  vjust = c(-0.8, 1.8, 1.8),
  col   = c("navy", "darkred", "deepskyblue3")
)

# ============================================================================
# 2-D filled-contour heatmap
# ============================================================================
# x = RMS dose, y = beta, fill = people required N (inverse of the floor form)
n_required <- function(dose_mGy, beta, p = P_BASE, k = K) {
  k^2 / (p * (dose_mGy / 1e3)^2 * beta^2)
}

dose_seq <- 10^seq(log10(0.01), log10(200), length.out = 200)   # mGy
grid <- expand.grid(dose = dose_seq, beta = beta_seq)
grid$N <- n_required(grid$dose, grid$beta)
grid$x <- log10(grid$dose)    # work in log space so geom_tile cells stay even
grid$y <- log10(grid$beta)

fill_breaks <- 10^c(4, 6, 8, 10, 12)   # colourbar levels
ref_breaks  <- 10^c(5, 7, 9, 11)       # white reference contours (clear of the pops)
log_lab     <- scales::trans_format("log10", scales::math_format(10^.x))  # legend

# unicode "10^k" for the inline labels: geomtextpath gaps the line around plain
# text but NOT around parsed plotmath, so use real superscript characters
sup <- function(k) {
  d <- c("0"="⁰","1"="¹","2"="²","3"="³","4"="⁴",
         "5"="⁵","6"="⁶","7"="⁷","8"="⁸","9"="⁹")
  vapply(k, function(v) paste0("10", paste0(d[strsplit(as.character(v), "")[[1]]],
                                            collapse = "")), character(1))
}

# the contours are straight lines of slope -1 in (log dose, log beta); build
# each line's in-frame segment so geomtextpath can draw the line AND write the
# label *into* it. geomtextpath rotates the text along the path automatically,
# so no aspect-ratio pinning is needed.
xlim <- range(grid$x); ylim <- range(grid$y)
line_seg <- function(cc, n = 80) {                 # the line x + y = cc
  x_lo <- max(xlim[1], cc - ylim[2])
  x_hi <- min(xlim[2], cc - ylim[1])
  xs   <- seq(x_lo, x_hi, length.out = n)
  data.frame(x = xs, y = cc - xs)
}
c_of_N <- function(N) log10(zfloor_mGy(1, N))      # x where the line crosses beta = 1

# place each label below the beta = 0.5 line (toward the plot bottom), clear of
# its point marker; hjust is the fraction along the path landing at that height
y_lab    <- log10(BETA_LNT) - 0.35
hjust_at <- function(cc) {
  s <- line_seg(cc)
  (cc - y_lab - s$x[1]) / (s$x[nrow(s)] - s$x[1])
}

# coloured population contours with the name written into the line
pop_paths <- lapply(seq_len(nrow(pts)), function(i) {
  cc <- c_of_N(pts$N[i])
  geom_textpath(data = line_seg(cc), aes(x, y), label = pts$label[i],
                inherit.aes = FALSE, colour = pts$col[i], linewidth = 0.8,
                size = 3, straight = TRUE, gap = TRUE,
                hjust = hjust_at(cc), vjust = 0.5)
})

# reference contours drawn as straight segments too, so the line-gap works;
# placed in the upper portion (hjust 0.28) to stay clear of the population labels
ref_paths <- lapply(seq_along(ref_breaks), function(i) {
  cc <- c_of_N(ref_breaks[i])
  geom_textpath(data = line_seg(cc), aes(x, y),
                label = sup(round(log10(ref_breaks[i]))), inherit.aes = FALSE,
                colour = "grey90", linewidth = 0.4, size = 2.6,
                straight = TRUE, gap = TRUE, vjust = 0.5, hjust = 0.28)
})

p_heat <- ggplot(grid, aes(x, y)) +
  geom_tile(aes(fill = N)) +
  # reference contours (straight segments) with the N written into the line
  ref_paths +
  pop_paths +
  # the LNT-scale slice (beta = 0.5), where the headline floors live
  geom_hline(yintercept = log10(BETA_LNT), linetype = "dashed",
             colour = "grey15", linewidth = 0.5) +
  # keep the floor points (markers); the names ride inside their contour lines
  geom_point(data = pts, aes(log10(floor), log10(beta)),
             colour = pts$col, size = 2.2) +
  scale_x_continuous(breaks = -2:2,
                     labels = c("0.01", "0.1", "1", "10", "100")) +
  scale_y_continuous(breaks = log10(c(0.1, 0.5, 1, 5)),
                     labels = c("0.1", "0.5", "1", "5")) +
  scale_fill_viridis_c(trans = "log10", name = "people\nrequired  N",
                       option = "D",
                       breaks = fill_breaks, labels = log_lab) +
  labs(
    title = "The possibility floor",
    x = "RMS cumulative dose (mGy)",
    y = expression("effect size  " * beta * "  (ERR/Gy)")
  ) +
  coord_cartesian(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(colour = "grey40", linewidth = 0.3),
        plot.title.position = "plot",  
        plot.title = element_text(hjust = 0.5),     
        axis.ticks.length = grid::unit(3, "pt"))

ggsave("floor_heatmap.png", p_heat, width = 8, height = 6, dpi = 600)
cat("wrote floor_heatmap.png\n")
