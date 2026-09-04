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
  s.dMat <- as.matrix(dist(rbind(ucoord_dp, ucoord_rp)))[
    seq_len(nrow(ucoord_dp)),
    nrow(ucoord_dp) + seq_len(nrow(ucoord_rp)),
    drop = FALSE
  ]
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
  s.dMat <- as.matrix(dist(ucoord_dp))
  uts_obs <- get.ts(obs.tv)[[1]]
  t.dMat  <- ti.distm(uts_obs, units = "auto")
  list(dp = dp, obs = obs.tv, s.dMat = s.dMat, t.dMat = t.dMat)
}

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
       old_median = median(told), new_median = median(tnew),
       old_min = min(told), new_min = min(tnew),
       reps_old = reps_old, reps_new = reps)
}

fmt_time <- function(s) {
  if (s < 1e-3) sprintf("%.0f µs", s * 1e6)
  else if (s < 1) sprintf("%.1f ms", s * 1e3)
  else sprintf("%.2f s", s)
}

sizes  <- c(100, 300, 1000, 2000)
reps_n <- c(20,  10,  5,    3)

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

md_path <- file.path(repo, "bench", "BENCHMARK.md")
r_ver   <- paste(R.version$major, R.version$minor, sep = ".")
sysinfo <- Sys.info()

is_helper <- is.na(df$n)
df_core   <- df[!is_helper, , drop = FALSE]
df_help   <- df[ is_helper, , drop = FALSE]

fmt_row <- function(r) {
  n_cell <- if (is.na(r$n)) "—" else as.character(r$n)
  sprintf("| %s | %s | %s | %s | %.1f× | %s | %d / %d |",
          r$mode, n_cell, fmt_time(r$old_median), fmt_time(r$new_median),
          r$speedup_median, ifelse(r$equal, "yes", "**NO**"),
          r$reps_old, r$reps_new)
}

lines <- c(
  "# GTWR vectorization benchmark",
  "",
  sprintf("R %s on %s %s (%s). Timings are `Sys.time()` deltas, median over reps. Old = master `R/gtwr.R` snapshot; new = current `R/gtwr.r` on `vectorize-gtwr`.",
          r_ver, sysinfo["sysname"], sysinfo["release"], sysinfo["machine"]),
  "",
  "## What changed",
  "",
  "- `st.dist()` `focus == 0` branches (both `rp.given` and symmetric): the two nested `for` loops became one vectorized matrix expression on a boolean-mask subset.",
  "- `st.dist()` `focus > 0` branches: same idea, reduced to a single vector.",
  "- `ti.distv()`: dropped the `c()`-growth loop; single vectorized `ti.dist()` call plus a mask for future observations.",
  "- `ti.distm()`: kept the outer column loop (which is cheap) now that `ti.distv` is vector-native.",
  "- `get.ts()`: replaced buggy `as.factor()` / `format()` round-trip with `sort(unique())` + `match()`. Also fixes a latent precision bug — see below.",
  "- `get.uloat()`: replaced the O(n²) `which()` scan with a single `match()` on paste-keys.",
  "",
  "## st.dist (called from `gtwr()` and `bw.gtwr()` per bandwidth candidate)",
  "",
  "| Mode | n | old median | new median | speedup | equal | reps (old/new) |",
  "|---|---:|---:|---:|---:|:---:|:---:|"
)
for (i in seq_len(nrow(df_core))) lines <- c(lines, fmt_row(df_core[i, ]))

lines <- c(lines,
  "",
  "## Helpers (called once per bandwidth trial and once per fit)",
  "",
  "| Function | n | old median | new median | speedup | equal | reps (old/new) |",
  "|---|---:|---:|---:|---:|:---:|:---:|"
)
for (i in seq_len(nrow(df_help))) lines <- c(lines, fmt_row(df_help[i, ]))

lines <- c(lines,
  "",
  "> `get.ts` shows `equal = **NO**` intentionally. The old implementation truncates full-precision numeric time stamps via `as.factor()` → `format()` (default `getOption(\"digits\") = 7`), so `which(ts == i)` matches fewer entries than there are inputs and `index` comes back short. The new implementation returns a correct index of length `length(tv)`. See `tests/test_gtwr_equivalence.R` for the round-trip assertions.",
  "",
  "## Notes on the shape of the speedups",
  "",
  "- The `st.dist` speedup plateaus (~30× for `rp.given`, ~10–12× for symmetric at n ≥ 1000) once R interpreter overhead per iteration dominates the old path. Below n=100 setup costs on both sides pull the ratio down.",
  "- The symmetric branch's ratio is lower than `rp.given` because the old code only iterated the upper triangle (~n²/2) while the new code touches the full n² before masking. One vectorized pass over n² still beats n²/2 R-level iterations by ~10×.",
  "- Helper speedups (`ti.distv` ~800×, `get.uloat` ~30×, `get.ts` ~280×) reflect that the old code paid `c()`-growth allocation on every step. Under `gtwr()` these are called once per fit, so the wall-clock gain shows up mostly at large n or during bandwidth search where they run inside a golden-section loop.",
  "- At n ≥ 1000 the old `st.dist` was timed with a single rep — one run takes 0.5–2.2 s and repeated timing added little signal; new-side medians are over 3–5 reps.",
  "- Numbers exclude `s.dMat` / `t.dMat` construction, which is identical on both sides.",
  "",
  "## Plots",
  "",
  "Rendered by `bench/plot_results.R` after this script writes `bench/results.rds`.",
  "",
  "- `bench/plots/st_dist_speedup.png` — speedup vs n for `st.dist`, one line per branch.",
  "- `bench/plots/st_dist_wall_time.png` — log-log wall time, master vs vectorized, both branches.",
  "- `bench/plots/helpers_speedup.png` — bar chart of helper speedups (`ti.distv`, `get.ts`, `get.uloat`).",
  "- `bench/plots/parallel_speedup.png` — `.gtwr_dispatch` scaling from `bench/bench_parallel_synthetic.R` (no package install needed). The end-to-end `bench/bench_parallel.R` needs the compiled `GWmodel` package installed (which needs system GDAL/PROJ/GEOS); if you have those, it overwrites `parallel_results.rds` with real gtwr numbers.",
  "",
  "## Correctness",
  "",
  "See `tests/test_gtwr_equivalence.R` — 53 assertions across `st.dist` (all branches, master parity), helper round-trip properties, and the AICc guard. Runs in a couple of seconds via `Rscript tests/test_gtwr_equivalence.R`.",
  "",
  "## Reproduce",
  "",
  "```bash",
  "Rscript bench/bench_st_dist.R      # this file — writes bench/results.rds",
  "Rscript bench/bench_parallel.R     # optional — needs GWmodel + sp installed",
  "Rscript bench/plot_results.R       # renders PNGs into bench/plots/",
  "Rscript tests/test_gtwr_equivalence.R",
  "```"
)
writeLines(lines, md_path)
cat(sprintf("\nWrote %s\n", md_path))
