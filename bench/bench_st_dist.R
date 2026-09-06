repo <- Sys.getenv("GWMODEL_REPO", unset = getwd())
if (!file.exists(file.path(repo, "R", "gtwr.r")) &&
    file.exists(file.path(getwd(), "R", "gtwr.r"))) repo <- getwd()

master_path <- file.path(repo, "bench", "gtwr_master.R")
new_path    <- file.path(repo, "R",     "gtwr.r")

stopifnot(file.exists(master_path), file.exists(new_path))

env_old <- new.env(parent = globalenv())
env_new <- new.env(parent = globalenv())
sys.source(master_path, envir = env_old)
sys.source(new_path,    envir = env_new)

st.dist_old <- env_old$st.dist
st.dist_new <- env_new$st.dist
get.uloat   <- env_new$get.uloat
get.ts      <- env_new$get.ts
ti.distm    <- env_new$ti.distm
ti.distv    <- env_new$ti.distv
ti.dist     <- env_new$ti.dist

environment(st.dist_old) <- environment()
environment(st.dist_new) <- environment()

make_case <- function(n_dp, n_rp, seed = 1L) {
  set.seed(seed)
  dp <- matrix(runif(2 * n_dp), ncol = 2)
  rp <- matrix(runif(2 * n_rp), ncol = 2)
  obs.tv <- as.numeric(sample.int(1e6, n_dp))
  reg.tv <- as.numeric(sample.int(1e6, n_rp))
  ucoord_dp <- get.uloat(dp)[[1]]
  ucoord_rp <- get.uloat(rp)[[1]]
  # direct cross-distance: dist(rbind(dp, rp)) would materialise a
  # (n_dp + n_rp)^2 matrix just to discard three quarters of it, which
  # becomes the binding constraint well before st.dist itself does.
  s.dMat <- sqrt(pmax(outer(rowSums(ucoord_dp^2), rowSums(ucoord_rp^2), "+") -
                      2 * tcrossprod(ucoord_dp, ucoord_rp), 0))
  uts_obs <- get.ts(obs.tv)[[1]]
  uts_reg <- get.ts(reg.tv)[[1]]
  t.dMat  <- ti.distm(uts_obs, uts_reg, units = "auto")
  list(dp = dp, rp = rp, obs = obs.tv, reg = reg.tv,
       s.dMat = s.dMat, t.dMat = t.dMat)
}

make_case_sym <- function(n, seed = 1L) {
  set.seed(seed)
  dp <- matrix(runif(2 * n), ncol = 2)
  obs.tv <- as.numeric(sample.int(1e6, n))
  ucoord_dp <- get.uloat(dp)[[1]]
  s.dMat <- sqrt(pmax(outer(rowSums(ucoord_dp^2), rowSums(ucoord_dp^2), "+") -
                      2 * tcrossprod(ucoord_dp), 0)) 
  uts_obs <- get.ts(obs.tv)[[1]]
  t.dMat  <- ti.distm(uts_obs, units = "auto")
  list(dp = dp, obs = obs.tv, s.dMat = s.dMat, t.dMat = t.dMat)
}

# Returns times plus the peak heap R reported while running them, so the
# memory cost of each implementation is visible alongside the speed.
time_reps <- function(expr, reps) {
  gc(verbose = FALSE, reset = TRUE)
  t <- numeric(reps)
  for (k in seq_len(reps)) {
    t0 <- Sys.time()
    force(expr())
    t[k] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  attr(t, "peak_mb") <- sum(gc()[, 6])
  t
}

bench_rp_given <- function(n, reps, lamda = 0.3, ksi = 0.2) {
  cs <- make_case(n, n)
  f_old <- function() st.dist_old(cs$dp, cs$rp, cs$obs, cs$reg,
                                  focus = 0, lamda = lamda, ksi = ksi,
                                  s.dMat = cs$s.dMat, t.dMat = cs$t.dMat)
  f_new <- function() st.dist_new(cs$dp, cs$rp, cs$obs, cs$reg,
                                  focus = 0, lamda = lamda, ksi = ksi,
                                  s.dMat = cs$s.dMat, t.dMat = cs$t.dMat)
  a <- f_old(); b <- f_new()
  ok <- isTRUE(all.equal(a, b, tolerance = 1e-10))
  reps_old <- if (n <= 500) reps else max(1L, reps %/% 3L)
  told <- time_reps(f_old, reps_old)
  tnew <- time_reps(f_new, reps)
  list(n = n, mode = "rp.given", equal = ok,
       old_peak_mb = attr(told, "peak_mb"), new_peak_mb = attr(tnew, "peak_mb"),
       old_median = median(told), new_median = median(tnew),
       old_min = min(told), new_min = min(tnew),
       reps_old = reps_old, reps_new = reps)
}

bench_sym <- function(n, reps, lamda = 0.3, ksi = 0.2) {
  cs <- make_case_sym(n)
  f_old <- function() st.dist_old(cs$dp, obs.tv = cs$obs,
                                  focus = 0, lamda = lamda, ksi = ksi,
                                  s.dMat = cs$s.dMat, t.dMat = cs$t.dMat)
  f_new <- function() st.dist_new(cs$dp, obs.tv = cs$obs,
                                  focus = 0, lamda = lamda, ksi = ksi,
                                  s.dMat = cs$s.dMat, t.dMat = cs$t.dMat)
  a <- f_old(); b <- f_new()
  ok <- isTRUE(all.equal(a, b, tolerance = 1e-10))
  reps_old <- if (n <= 500) reps else max(1L, reps %/% 3L)
  told <- time_reps(f_old, reps_old)
  tnew <- time_reps(f_new, reps)
  list(n = n, mode = "symmetric", equal = ok,
       old_peak_mb = attr(told, "peak_mb"), new_peak_mb = attr(tnew, "peak_mb"),
       old_median = median(told), new_median = median(tnew),
       old_min = min(told), new_min = min(tnew),
       reps_old = reps_old, reps_new = reps)
}

fmt_time <- function(s) {
  if (s < 1e-3) sprintf("%.0f µs", s * 1e6)
  else if (s < 1) sprintf("%.1f ms", s * 1e3)
  else sprintf("%.2f s", s)
}

# Sizes are capped by the nested-loop side: at n=8000 one master call is
# several minutes, so it gets a single rep while the vectorised side keeps 2-3.
sizes  <- c(100, 300, 1000, 2000, 4000, 8000, 12000, 16000)
reps_n <- c(20,  10,  5,    3,    3,    2,    1,     1)

results <- list()
for (i in seq_along(sizes)) {
  cat(sprintf("[rp.given ] n=%d ...\n", sizes[i]))
  results[[length(results) + 1L]] <- bench_rp_given(sizes[i], reps_n[i])
  cat(sprintf("[symmetric] n=%d ...\n", sizes[i]))
  results[[length(results) + 1L]] <- bench_sym(sizes[i], reps_n[i])
}

bench_helper <- function(label, f_old, f_new, reps) {
  gc(verbose = FALSE)
  told <- time_reps(f_old, reps)
  tnew <- time_reps(f_new, reps)
  a <- f_old(); b <- f_new()
  list(mode = label, n = NA_integer_,
       old_peak_mb = attr(told, "peak_mb"), new_peak_mb = attr(tnew, "peak_mb"),
       equal = isTRUE(all.equal(a, b, tolerance = 1e-10)),
       old_median = median(told), new_median = median(tnew),
       old_min = min(told), new_min = min(tnew),
       reps_old = reps, reps_new = reps)
}

cat("[helper ti.distv]\n")
{
  set.seed(9); n <- 5000
  obs <- as.numeric(sample.int(1e7, n))
  focal <- as.numeric(sample.int(1e7, 1))
  results[[length(results) + 1L]] <- bench_helper(
    "ti.distv n=5000",
    function() env_old$ti.distv(focal, obs, "auto"),
    function() env_new$ti.distv(focal, obs, "auto"),
    reps = 5)
}
cat("[helper get.ts]\n")
{
  set.seed(11); n <- 5000
  tv <- as.numeric(sample.int(1e7, n))
  results[[length(results) + 1L]] <- bench_helper(
    "get.ts n=5000",
    function() env_old$get.ts(tv),
    function() env_new$get.ts(tv),
    reps = 5)
}
cat("[helper get.uloat]\n")
{
  set.seed(13); n <- 5000
  coords <- matrix(runif(2 * n), ncol = 2)
  coords <- rbind(coords, coords[sample.int(n, n %/% 10), ])
  results[[length(results) + 1L]] <- bench_helper(
    "get.uloat n=5500",
    function() env_old$get.uloat(coords),
    function() env_new$get.uloat(coords),
    reps = 5)
}

df <- do.call(rbind, lapply(results, function(r) {
  data.frame(mode = r$mode, n = r$n,
             equal = r$equal,
             old_peak_mb = r$old_peak_mb, new_peak_mb = r$new_peak_mb,
             old_median_s = r$old_median,
             new_median_s = r$new_median,
             old_min_s    = r$old_min,
             new_min_s    = r$new_min,
             speedup_median = r$old_median / r$new_median,
             reps_old = r$reps_old, reps_new = r$reps_new,
             stringsAsFactors = FALSE)
}))

df <- df[order(df$mode, df$n), ]

cat("\n=== RESULTS ===\n")
print(df, row.names = FALSE)

saveRDS(df, file.path(repo, "bench", "results.rds"))
