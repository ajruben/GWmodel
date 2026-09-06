repo <- Sys.getenv("GWMODEL_REPO", unset = getwd())
if (!file.exists(file.path(repo, "R", "gtwr.r")) &&
    file.exists(file.path(getwd(), "R", "gtwr.r"))) repo <- getwd()

master_path <- file.path(repo, "bench", "gtwr_master.r")
new_path    <- file.path(repo, "R",     "gtwr.r")
stopifnot(file.exists(master_path), file.exists(new_path))

env_old <- new.env(parent = globalenv()); sys.source(master_path, envir = env_old)
env_new <- new.env(parent = globalenv()); sys.source(new_path,    envir = env_new)

get.uloat_new <- env_new$get.uloat
get.ts_new    <- env_new$get.ts
ti.distm_new  <- env_new$ti.distm

pass <- 0L; fail <- 0L; failures <- character(0)
assert <- function(cond, name) {
  if (isTRUE(cond)) {
    pass <<- pass + 1L
    cat(sprintf("  ok   %s\n", name))
  } else {
    fail <<- fail + 1L
    failures <<- c(failures, name)
    cat(sprintf("  FAIL %s\n", name))
  }
}
assert_equal <- function(a, b, name, tol = 1e-10) {
  eq <- all.equal(a, b, tolerance = tol)
  assert(isTRUE(eq), if (isTRUE(eq)) name else sprintf("%s [%s]", name, paste(eq, collapse="; ")))
}

call_st_dist <- function(env, ...) {
  f <- env$st.dist
  environment(f) <- env_new
  f(...)
}

build_case <- function(n_dp, n_rp, seed = 1L) {
  set.seed(seed)
  dp <- matrix(runif(2 * n_dp), ncol = 2)
  rp <- matrix(runif(2 * n_rp), ncol = 2)
  obs <- as.numeric(sample.int(1e6, n_dp))
  reg <- as.numeric(sample.int(1e6, n_rp))
  ucoord_dp <- get.uloat_new(dp)[[1]]
  ucoord_rp <- get.uloat_new(rp)[[1]]
  s.dMat <- as.matrix(dist(rbind(ucoord_dp, ucoord_rp)))[
    seq_len(nrow(ucoord_dp)),
    nrow(ucoord_dp) + seq_len(nrow(ucoord_rp)),
    drop = FALSE
  ]
  uts_obs <- get.ts_new(obs)[[1]]
  uts_reg <- get.ts_new(reg)[[1]]
  t.dMat  <- ti.distm_new(uts_obs, uts_reg, units = "auto")
  list(dp = dp, rp = rp, obs = obs, reg = reg, s.dMat = s.dMat, t.dMat = t.dMat)
}

build_case_sym <- function(n, seed = 1L) {
  set.seed(seed)
  dp <- matrix(runif(2 * n), ncol = 2)
  obs <- as.numeric(sample.int(1e6, n))
  ucoord_dp <- get.uloat_new(dp)[[1]]
  s.dMat <- as.matrix(dist(ucoord_dp))
  uts_obs <- get.ts_new(obs)[[1]]
  t.dMat  <- ti.distm_new(uts_obs, units = "auto")
  list(dp = dp, obs = obs, s.dMat = s.dMat, t.dMat = t.dMat)
}

cat("== st.dist focus=0, rp.given ==\n")
for (n in c(50, 200, 500)) {
  cs <- build_case(n, n)
  a <- call_st_dist(env_old, cs$dp, cs$rp, cs$obs, cs$reg,
                    focus = 0, lamda = 0.3, ksi = 0.2,
                    s.dMat = cs$s.dMat, t.dMat = cs$t.dMat)
  b <- call_st_dist(env_new, cs$dp, cs$rp, cs$obs, cs$reg,
                    focus = 0, lamda = 0.3, ksi = 0.2,
                    s.dMat = cs$s.dMat, t.dMat = cs$t.dMat)
  assert_equal(a, b, sprintf("rp.given n=%d matches master", n))
  assert(is.matrix(b) && all(dim(b) == c(n, n)),
         sprintf("rp.given n=%d shape %dx%d", n, n, n))
}

cat("== st.dist focus=0, symmetric ==\n")
for (n in c(50, 200, 500)) {
  cs <- build_case_sym(n)
  a <- call_st_dist(env_old, cs$dp, obs.tv = cs$obs,
                    focus = 0, lamda = 0.3, ksi = 0.2,
                    s.dMat = cs$s.dMat, t.dMat = cs$t.dMat)
  b <- call_st_dist(env_new, cs$dp, obs.tv = cs$obs,
                    focus = 0, lamda = 0.3, ksi = 0.2,
                    s.dMat = cs$s.dMat, t.dMat = cs$t.dMat)
  assert_equal(a, b, sprintf("symmetric n=%d matches master", n))
  assert_equal(b, t(b), sprintf("symmetric n=%d is symmetric", n))
  assert(all(diag(b) == 0), sprintf("symmetric n=%d diagonal is zero", n))
}

cat("== st.dist focus>0, rp.given ==\n")
for (n in c(50, 200)) {
  cs <- build_case(n, n)
  full_old <- call_st_dist(env_old, cs$dp, cs$rp, cs$obs, cs$reg,
                           focus = 0, lamda = 0.3, ksi = 0.2,
                           s.dMat = cs$s.dMat, t.dMat = cs$t.dMat)
  full_new <- call_st_dist(env_new, cs$dp, cs$rp, cs$obs, cs$reg,
                           focus = 0, lamda = 0.3, ksi = 0.2,
                           s.dMat = cs$s.dMat, t.dMat = cs$t.dMat)
  for (i in c(1, ceiling(n/3), n)) {
    a <- call_st_dist(env_old, cs$dp, cs$rp, cs$obs, cs$reg,
                      focus = i, lamda = 0.3, ksi = 0.2,
                      s.dMat = cs$s.dMat, t.dMat = cs$t.dMat)
    b <- call_st_dist(env_new, cs$dp, cs$rp, cs$obs, cs$reg,
                      focus = i, lamda = 0.3, ksi = 0.2,
                      s.dMat = cs$s.dMat, t.dMat = cs$t.dMat)
    assert_equal(a, b, sprintf("focus=%d rp.given n=%d master parity", i, n))
    assert_equal(as.numeric(b), full_new[, i],
                 sprintf("focus=%d rp.given n=%d equals column of full matrix", i, n))
  }
}

cat("== st.dist focus>0, symmetric (master parity only; cross-path consistency intentionally not asserted) ==\n")
for (n in c(50, 200)) {
  cs <- build_case_sym(n)
  for (i in c(1, ceiling(n/3), n)) {
    a <- call_st_dist(env_old, cs$dp, obs.tv = cs$obs,
                      focus = i, lamda = 0.3, ksi = 0.2,
                      s.dMat = cs$s.dMat, t.dMat = cs$t.dMat)
    b <- call_st_dist(env_new, cs$dp, obs.tv = cs$obs,
                      focus = i, lamda = 0.3, ksi = 0.2,
                      s.dMat = cs$s.dMat, t.dMat = cs$t.dMat)
    assert_equal(a, b, sprintf("focus=%d symmetric n=%d master parity", i, n))
  }
}

cat("== helpers: get.ts round-trip ==\n")
{
  set.seed(42)
  tv_int  <- as.numeric(sample.int(1e6, 100))
  tv_frac <- runif(100, 0, 100)
  for (label in c("int", "frac")) {
    tv <- if (label == "int") tv_int else tv_frac
    r  <- get.ts_new(tv)
    ts_new <- r[[1]]; idx_new <- r[[2]]
    assert(length(idx_new) == length(tv),
           sprintf("get.ts %s: index length = input length", label))
    assert(!anyNA(idx_new),
           sprintf("get.ts %s: index has no NA", label))
    assert(all(idx_new >= 1 & idx_new <= length(ts_new)),
           sprintf("get.ts %s: indices within range", label))
    assert(isTRUE(all.equal(ts_new[idx_new], tv)),
           sprintf("get.ts %s: ts[index] reconstructs input", label))
    assert(!is.unsorted(ts_new),
           sprintf("get.ts %s: ts is sorted", label))
    assert(!anyDuplicated(ts_new),
           sprintf("get.ts %s: ts is unique", label))
  }
}

cat("== helpers: get.uloat round-trip ==\n")
{
  set.seed(7)
  coords <- matrix(runif(200), ncol = 2)
  coords <- rbind(coords, coords[1:5, , drop = FALSE])
  r <- get.uloat_new(coords)
  ucoords <- r[[1]]; idx <- r[[2]]
  assert(length(idx) == nrow(coords),
         "get.uloat: index length equals nrow(coords)")
  assert(nrow(ucoords) < nrow(coords),
         "get.uloat: unique coords fewer than input (duplicates removed)")
  assert(isTRUE(all.equal(unname(ucoords[idx, ]), unname(coords))),
         "get.uloat: ucoords[index,] reconstructs input")
}

cat("== helpers: ti.distv basic ==\n")
{
  obs <- as.numeric(c(1, 2, 3, 4, 5))
  d <- env_new$ti.distv(3, obs, units = "auto")
  assert_equal(d, c(2, 1, 0, Inf, Inf), "ti.distv future obs -> Inf")
  d2 <- env_new$ti.distv(5, obs, units = "auto")
  assert_equal(d2, c(4, 3, 2, 1, 0), "ti.distv all past = correct distances")
}

cat("== gtwr.aic guard: well-conditioned matches unguarded formula ==\n")
{
  aicc_formula <- function(y, X, betas, S,
                           rss.floor = 0, enp.margin = 0) {
    n <- length(y); tr.S <- sum(diag(S))
    if (tr.S >= n - 2 - enp.margin) return(Inf)
    yhat <- rowSums(X * betas)
    rss  <- sum((y - yhat)^2)
    gTSS <- sum((y - mean(y))^2)
    rss.eff <- max(rss, rss.floor * gTSS)
    n * log(rss.eff / n) + n * log(2 * pi) +
      n * ((n + tr.S) / (n - 2 - tr.S))
  }
  set.seed(3); n <- 40; p <- 3
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n))
  beta_true <- runif(p, -1, 1)
  Y <- X %*% beta_true + rnorm(n, sd = 0.5)
  betas <- matrix(rep(beta_true, each = n), n, p)
  S <- diag(n) * 0.2
  guarded   <- aicc_formula(Y, X, betas, S, rss.floor = 1e-8, enp.margin = 1)
  unguarded <- aicc_formula(Y, X, betas, S, rss.floor = 0,    enp.margin = 0)
  assert_equal(guarded, unguarded, "AICc guard: well-conditioned unchanged")

  S_over  <- diag(n) * ((n - 2 - 0.5) / n)
  aicc_over <- aicc_formula(Y, X, betas, S_over, rss.floor = 1e-8, enp.margin = 1)
  assert(is.infinite(aicc_over) && aicc_over > 0,
         "AICc guard: tr(S) too large -> Inf")

  betas_perfect <- solve(t(X) %*% X) %*% t(X) %*% Y
  betas_perfect <- matrix(rep(betas_perfect, each = n), n, p)
  Y_perfect <- X %*% betas_perfect[1, ]
  aicc_floored <- aicc_formula(Y_perfect, X, betas_perfect,
                               diag(n) * 0.2, rss.floor = 1e-8, enp.margin = 1)
  assert(is.finite(aicc_floored),
         "AICc guard: near-zero RSS floored to finite AICc")
}

cat(sprintf("\nResults: %d passed, %d failed\n", pass, fail))
if (fail > 0) {
  cat("Failures:\n")
  for (f in failures) cat("  -", f, "\n")
  quit(status = 1)
}
