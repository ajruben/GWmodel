.gtwr_dispatch <- function(fit_fn, n, cores, verbose)
{
  cores <- as.integer(cores)
  if (is.na(cores) || cores < 1L) cores <- 1L
  if (cores == 1L) {
    pb <- if (isTRUE(verbose) && n > 1L)
            utils::txtProgressBar(min = 0, max = n, style = 3) else NULL
    out <- vector("list", n)
    for (i in seq_len(n)) {
      out[[i]] <- fit_fn(i)
      if (!is.null(pb)) utils::setTxtProgressBar(pb, i)
    }
    if (!is.null(pb)) { close(pb); cat("\n") }
    return(out)
  }
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makePSOCKcluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, { requireNamespace("GWmodel", quietly = TRUE) })
    fn_env <- environment(fit_fn)
    parallel::clusterExport(cl, varlist = ls(envir = fn_env), envir = fn_env)
    return(parallel::parLapply(cl, seq_len(n), fit_fn))
  }
  parallel::mclapply(seq_len(n), fit_fn, mc.cores = cores)
}
