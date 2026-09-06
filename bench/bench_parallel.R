repo <- Sys.getenv("GWMODEL_REPO", unset = getwd())
if (!file.exists(file.path(repo, "R", "gtwr.R")) &&
    file.exists(file.path(getwd(), "R", "gtwr.R"))) repo <- getwd()

if (!requireNamespace("GWmodel", quietly = TRUE) ||
    !requireNamespace("sp", quietly = TRUE)) {
  cat("SKIP: GWmodel/sp not installed; parallel bench needs the compiled package\n")
  quit(status = 0)
}
suppressPackageStartupMessages({ library(sp); library(GWmodel) })

sys.source(file.path(repo, "R", "gtwr_parallel.R"), envir = globalenv())
sys.source(file.path(repo, "R", "gtwr.R"),          envir = globalenv())

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

make_case <- function(n, seed = 1L) {
  set.seed(seed)
  coords <- matrix(runif(2 * n, 0, 100), ncol = 2)
  obs    <- as.numeric(sample.int(1e6, n))
  x1 <- runif(n); x2 <- runif(n)
  y  <- 1 + 0.5 * x1 - 0.3 * x2 + rnorm(n, sd = 0.2)
  df <- data.frame(y = y, x1 = x1, x2 = x2)
  spdf <- SpatialPointsDataFrame(coords = coords, data = df,
                                 proj4string = CRS("+proj=longlat +datum=WGS84"))
  st.dMat <- st.dist(coords, obs.tv = obs, lamda = 0.3, ksi = 0.1)
  list(spdf = spdf, obs = obs, st.dMat = st.dMat)
}

sizes <- c(150, 400, 800)
cores_grid <- c(1L, 2L, 4L)
reps <- c(3, 3, 2)
results <- list()

for (i in seq_along(sizes)) {
  n <- sizes[i]; r <- reps[i]
  cs <- make_case(n)
  cat(sprintf("[bench parallel] n=%d serial baseline ...\n", n))
  t_ser <- time_reps(function() {
    gtwr(y ~ x1 + x2, data = cs$spdf, obs.tv = cs$obs, st.bw = 40, adaptive = TRUE,
         lamda = 0.3, ksi = 0.1, st.dMat = cs$st.dMat,
         cores = 1L, verbose = FALSE)
  }, r)
  ser_med <- median(t_ser)
  for (co in cores_grid) {
    cat(sprintf("[bench parallel] n=%d cores=%d ...\n", n, co))
    t_par <- time_reps(function() {
      gtwr(y ~ x1 + x2, data = cs$spdf, obs.tv = cs$obs, st.bw = 40, adaptive = TRUE,
           lamda = 0.3, ksi = 0.1, st.dMat = cs$st.dMat,
           cores = co, verbose = FALSE)
    }, r)
    par_med <- median(t_par)
    results[[length(results) + 1L]] <- data.frame(
      n = n, cores = co,
      serial_median_s = ser_med,
      parallel_median_s = par_med,
      speedup = ser_med / par_med,
      reps = r
    )
  }
}

pdf_df <- do.call(rbind, results)
print(pdf_df, row.names = FALSE)
saveRDS(pdf_df, file.path(repo, "bench", "parallel_results.rds"))
