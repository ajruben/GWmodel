repo <- Sys.getenv("GWMODEL_REPO", unset = getwd())
if (!file.exists(file.path(repo, "R", "gtwr.r")) &&
    file.exists(file.path(getwd(), "R", "gtwr.r"))) repo <- getwd()

sys.source(file.path(repo, "R", "gtwr_parallel.R"), envir = globalenv())

time_reps <- function(expr, reps) {
  gc(verbose = FALSE)
  t <- numeric(reps)
  for (k in seq_len(reps)) {
    t0 <- Sys.time()
    force(expr())
    t[k] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  t
}

make_workload <- function(n_ctx, p) {
  X <- cbind(1, matrix(rnorm(n_ctx * (p - 1)), n_ctx))
  y <- rnorm(n_ctx)
  d <- runif(n_ctx, 0.1, 10)
  function(i) {
    w  <- exp(-(d + rnorm(1, sd = 0.01))^2 / 4)
    Xw <- X * w
    A  <- crossprod(X, Xw)
    Ainv <- solve(A)
    beta <- Ainv %*% crossprod(Xw, y)
    Ci   <- Ainv %*% t(Xw)
    S_row <- X[i %% n_ctx + 1L, , drop = FALSE] %*% Ci
    se_sq <- diag(Ci %*% t(Ci))
    list(beta = beta, S_row = S_row, se_sq = se_sq)
  }
}

sizes      <- c(300L, 800L, 1500L)
cores_grid <- c(1L, 2L, 4L, max(1L, parallel::detectCores(logical = FALSE)))
cores_grid <- sort(unique(cores_grid))
reps       <- 3L
n_ctx      <- 1500L
p          <- 6L

set.seed(42)
worker <- make_workload(n_ctx, p)

results <- list()
for (n in sizes) {
  cat(sprintf("[synthetic] n=%d serial baseline ...\n", n))
  t_ser <- time_reps(function() .gtwr_dispatch(worker, n, 1L, FALSE), reps)
  ser_med <- median(t_ser)
  for (co in cores_grid) {
    cat(sprintf("[synthetic] n=%d cores=%d ...\n", n, co))
    t_par <- time_reps(function() .gtwr_dispatch(worker, n, co, FALSE), reps)
    par_med <- median(t_par)
    results[[length(results) + 1L]] <- data.frame(
      n = n, cores = co,
      serial_median_s = ser_med,
      parallel_median_s = par_med,
      speedup = ser_med / par_med,
      reps = reps
    )
  }
}

pdf_df <- do.call(rbind, results)
print(pdf_df, row.names = FALSE)
saveRDS(pdf_df, file.path(repo, "bench", "parallel_results.rds"))

md <- file.path(repo, "bench", "BENCHMARK_PARALLEL.md")
sysinfo <- Sys.info()
lines <- c(
  "# .gtwr_dispatch parallel-scaling benchmark (synthetic)",
  "",
  sprintf("R %s.%s on %s %s (%s)", R.version$major, R.version$minor,
          sysinfo["sysname"], sysinfo["release"], sysinfo["machine"]),
  sprintf("Cores detected (physical): %d", parallel::detectCores(logical = FALSE)),
  sprintf("Cluster type: %s", if (.Platform$OS.type == "unix") "fork (mclapply)" else "psock (parLapply)"),
  "",
  paste("Workload per point: small `crossprod` + `solve` at p =", p,
        "and", n_ctx, "context observations, ~equivalent to the arithmetic that `gw_reg` does at each regression point without requiring the compiled GWmodel package.",
        "Compares `.gtwr_dispatch(..., cores = 1)` (serial) vs `cores = k`."),
  "",
  "| n | cores | serial (s) | parallel (s) | speedup | reps |",
  "|---:|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(pdf_df))) {
  r <- pdf_df[i, ]
  lines <- c(lines, sprintf("| %d | %d | %.3f | %.3f | %.2f× | %d |",
                            r$n, r$cores, r$serial_median_s, r$parallel_median_s,
                            r$speedup, r$reps))
}
lines <- c(lines,
  "",
  "> The end-to-end `gtwr()` bench (`bench/bench_parallel.R`) needs the compiled `GWmodel` package installed (which in turn needs system GDAL/PROJ/GEOS). If you have those, run it and it will overwrite `parallel_results.rds` with real end-to-end numbers.",
  "")
writeLines(lines, md)
cat("Wrote", md, "\n")
