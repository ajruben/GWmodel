repo <- Sys.getenv("GWMODEL_REPO", unset = getwd())
if (!file.exists(file.path(repo, "R", "gtwr.R")) &&
    file.exists(file.path(getwd(), "R", "gtwr.R"))) repo <- getwd()

results_rds <- file.path(repo, "bench", "results.rds")
if (!file.exists(results_rds)) {
  stop("Run bench/bench_st_dist.R first — no results.rds found.")
}
df <- readRDS(results_rds)

plot_dir <- file.path(repo, "bench", "plots")
dir.create(plot_dir, showWarnings = FALSE)

st_df <- df[df$mode %in% c("rp.given", "symmetric") & !is.na(df$n), ]
st_df <- st_df[order(st_df$mode, st_df$n), ]

png(file.path(plot_dir, "st_dist_speedup.png"),
    width = 900, height = 560, res = 140)
op <- par(mar = c(4.5, 4.5, 3, 1), mgp = c(2.6, 0.8, 0))
modes <- unique(st_df$mode)
cols  <- c("rp.given" = "#1f77b4", "symmetric" = "#d62728")
ylim  <- c(1, max(st_df$speedup_median) * 1.15)
plot(NA, xlim = range(st_df$n), ylim = ylim, log = "x",
     xlab = "n (matrix side)", ylab = "speedup vs master (×)",
     main = "st.dist: vectorized vs nested-loop master")
grid(col = "grey85", lty = 1)
abline(h = 1, lty = 2, col = "grey50")
for (m in modes) {
  sub <- st_df[st_df$mode == m, ]
  lines(sub$n, sub$speedup_median, type = "b", pch = 19, lwd = 2, col = cols[m])
  text(sub$n, sub$speedup_median,
       labels = sprintf("%.1f×", sub$speedup_median),
       pos = 3, cex = 0.8, col = cols[m])
}
legend("bottomright", legend = names(cols), col = cols, lty = 1, lwd = 2,
       pch = 19, bty = "n")
par(op); dev.off()

png(file.path(plot_dir, "st_dist_wall_time.png"),
    width = 900, height = 560, res = 140)
op <- par(mar = c(4.5, 4.5, 3, 1), mgp = c(2.6, 0.8, 0))
ylim <- range(c(st_df$old_median_s, st_df$new_median_s))
plot(NA, xlim = range(st_df$n), ylim = ylim, log = "xy",
     xlab = "n (matrix side)", ylab = "median wall time (s)",
     main = "st.dist wall time")
grid(col = "grey85", lty = 1)
for (m in modes) {
  sub <- st_df[st_df$mode == m, ]
  lines(sub$n, sub$old_median_s, type = "b", pch = 1, lty = 2, lwd = 2, col = cols[m])
  lines(sub$n, sub$new_median_s, type = "b", pch = 19, lty = 1, lwd = 2, col = cols[m])
}
legend("topleft",
       legend = c("rp.given master", "rp.given vectorized",
                  "symmetric master", "symmetric vectorized"),
       col = c(cols["rp.given"], cols["rp.given"], cols["symmetric"], cols["symmetric"]),
       pch = c(1, 19, 1, 19), lty = c(2, 1, 2, 1), lwd = 2, bty = "n", cex = 0.85)
par(op); dev.off()

help_df <- df[!(df$mode %in% c("rp.given", "symmetric")), ]
if (nrow(help_df) > 0) {
  png(file.path(plot_dir, "helpers_speedup.png"),
      width = 900, height = 480, res = 140)
  op <- par(mar = c(4.5, 8, 3, 1), mgp = c(2.6, 0.8, 0))
  ord <- order(help_df$speedup_median)
  h <- help_df[ord, ]
  h$speedup_capped <- pmin(h$speedup_median, 1000)
  bp <- barplot(h$speedup_capped, names.arg = h$mode, horiz = TRUE, las = 1,
                col = "#2ca02c", border = NA, xlim = c(0, max(h$speedup_capped) * 1.2),
                xlab = "speedup vs master (×)",
                main = "Helper vectorizations (capped at 1000×)")
  text(h$speedup_capped, bp,
       labels = sprintf("%.0f×%s", h$speedup_median,
                        ifelse(h$speedup_median > 1000, " (capped)", "")),
       pos = 4, cex = 0.85)
  par(op); dev.off()
}

par_rds <- file.path(repo, "bench", "parallel_results.rds")
if (file.exists(par_rds)) {
  pdf_df <- readRDS(par_rds)
  png(file.path(plot_dir, "parallel_speedup.png"),
      width = 900, height = 560, res = 140)
  op <- par(mar = c(4.5, 4.5, 3, 1), mgp = c(2.6, 0.8, 0))
  ns <- sort(unique(pdf_df$n))
  cols_n <- setNames(rainbow(length(ns), start = 0, end = 0.7), as.character(ns))
  ylim <- c(0, max(pdf_df$speedup) * 1.15)
  plot(NA, xlim = range(pdf_df$cores), ylim = ylim,
       xlab = "cores", ylab = "speedup over serial (×)",
       main = "gtwr_parallel scaling")
  grid(col = "grey85", lty = 1)
  abline(h = 1, lty = 2, col = "grey50")
  abline(0, 1, lty = 3, col = "grey60")
  for (nn in ns) {
    sub <- pdf_df[pdf_df$n == nn, ]
    sub <- sub[order(sub$cores), ]
    lines(sub$cores, sub$speedup, type = "b", pch = 19, lwd = 2,
          col = cols_n[as.character(nn)])
  }
  legend("topleft", legend = sprintf("n=%d", ns), col = cols_n[as.character(ns)],
         lty = 1, lwd = 2, pch = 19, bty = "n")
  par(op); dev.off()
}

cat("Wrote plots to", plot_dir, "\n")
list.files(plot_dir, full.names = TRUE)
