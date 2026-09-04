repo <- Sys.getenv("GWMODEL_REPO", unset = getwd())
if (!file.exists(file.path(repo, "R", "gtwr.R")) &&
    file.exists(file.path(getwd(), "R", "gtwr.R"))) repo <- getwd()

if (!requireNamespace("GWmodel", quietly = TRUE) ||
    !requireNamespace("sp", quietly = TRUE)) {
  cat("SKIP: GWmodel or sp not installed; parallel test needs the compiled package\n")
  quit(status = 0)
}
suppressPackageStartupMessages({ library(sp); library(GWmodel) })

sys.source(file.path(repo, "R", "gtwr_parallel.R"), envir = globalenv())
sys.source(file.path(repo, "R", "gtwr.R"),          envir = globalenv())

pass <- 0L; fail <- 0L; failures <- character(0)
assert_equal <- function(a, b, name, tol = 1e-8) {
  eq <- all.equal(a, b, tolerance = tol)
  if (isTRUE(eq)) {
    pass <<- pass + 1L
    cat(sprintf("  ok   %s\n", name))
  } else {
    fail <<- fail + 1L
    failures <<- c(failures, name)
    cat(sprintf("  FAIL %s [%s]\n", name, paste(eq, collapse = "; ")))
  }
}

set.seed(11); n <- 120
coords <- matrix(runif(2 * n, 0, 100), ncol = 2)
obs.tv <- as.numeric(sample.int(300, n))
x1 <- runif(n); x2 <- runif(n)
y  <- 1 + 0.5 * x1 - 0.3 * x2 + rnorm(n, sd = 0.2)
df <- data.frame(y = y, x1 = x1, x2 = x2)
spdf <- SpatialPointsDataFrame(coords = coords, data = df,
                               proj4string = CRS("+proj=longlat +datum=WGS84"))
st.dMat <- st.dist(coords, obs.tv = obs.tv, lamda = 0.9, ksi = 0.1)

get_betas <- function(fit) as.matrix(fit$SDF@data[, c("Intercept", "x1", "x2")])

fit_ser <- gtwr(y ~ x1 + x2, data = spdf, obs.tv = obs.tv, st.bw = 300, adaptive = FALSE,
                lamda = 0.9, ksi = 0.1, st.dMat = st.dMat,
                cores = 1L, verbose = FALSE)

if (.Platform$OS.type == "unix") {
  cat("== gtwr cores=2 (fork) matches serial ==\n")
  fit_p2 <- gtwr(y ~ x1 + x2, data = spdf, obs.tv = obs.tv, st.bw = 300, adaptive = FALSE,
                 lamda = 0.9, ksi = 0.1, st.dMat = st.dMat,
                 cores = 2L, verbose = FALSE)
  assert_equal(get_betas(fit_p2), get_betas(fit_ser), "cores=2 fork betas match")
  assert_equal(fit_p2$GTW.diagnostic$RSS.gw, fit_ser$GTW.diagnostic$RSS.gw,
               "cores=2 fork RSS matches")
  assert_equal(fit_p2$GTW.diagnostic$AICc, fit_ser$GTW.diagnostic$AICc,
               "cores=2 fork AICc matches")
  assert_equal(as.matrix(fit_p2$SDF@data), as.matrix(fit_ser$SDF@data),
               "cores=2 fork full SDF data matches")

  cat("== gtwr cores=4 (fork) matches serial ==\n")
  fit_p4 <- gtwr(y ~ x1 + x2, data = spdf, obs.tv = obs.tv, st.bw = 300, adaptive = FALSE,
                 lamda = 0.9, ksi = 0.1, st.dMat = st.dMat,
                 cores = 4L, verbose = FALSE)
  assert_equal(get_betas(fit_p4), get_betas(fit_ser), "cores=4 fork betas match")
  assert_equal(fit_p4$GTW.diagnostic$RSS.gw, fit_ser$GTW.diagnostic$RSS.gw,
               "cores=4 fork RSS matches")
} else {
  cat("SKIP fork tests (non-unix); PSOCK path not exercised here\n")
}

cat("== gtwr reporting guard ==\n")
{
  set.seed(21); n2 <- 80
  cds <- matrix(runif(2 * n2, 0, 100), ncol = 2)
  tv  <- as.numeric(sample.int(200, n2))
  xa  <- runif(n2); xb <- runif(n2)
  ya  <- xa + rnorm(n2, sd = 0.1)
  sm  <- SpatialPointsDataFrame(coords = cds,
                                data = data.frame(y = ya, xa = xa, xb = xb),
                                proj4string = CRS("+proj=longlat +datum=WGS84"))
  dm  <- st.dist(cds, obs.tv = tv, lamda = 0.9, ksi = 0.1)

  fit_ok <- suppressWarnings(
    gtwr(y ~ xa + xb, data = sm, obs.tv = tv, st.bw = 200, adaptive = FALSE,
         lamda = 0.9, ksi = 0.1, st.dMat = dm, cores = 1L, verbose = FALSE))
  assert_equal(is.finite(fit_ok$GTW.diagnostic$AICc), TRUE,
               "reasonable bw: AICc finite")
  assert_equal(is.finite(fit_ok$GTW.diagnostic$gwR2.adj), TRUE,
               "reasonable bw: gwR2.adj finite")

  ww <- character(0)
  fit_bad <- withCallingHandlers(
    gtwr(y ~ xa + xb, data = sm, obs.tv = tv, st.bw = 1, adaptive = FALSE,
         lamda = 0.9, ksi = 0.1, st.dMat = dm, cores = 1L, verbose = FALSE),
    warning = function(w) { ww <<- c(ww, conditionMessage(w)); invokeRestart("muffleWarning") })
  guard_tripped <- is.na(fit_bad$GTW.diagnostic$AICc) ||
                   any(grepl("unreliable|RSS|degrees of freedom", ww))
  assert_equal(guard_tripped, TRUE, "over-fit bw: guard trips (warning or NA)")
}

cat(sprintf("\nResults: %d passed, %d failed\n", pass, fail))
if (fail > 0) {
  for (f in failures) cat("  -", f, "\n")
  quit(status = 1)
}
