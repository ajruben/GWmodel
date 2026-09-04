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
    if (!is.null(pb)) { close(pb); cat("
") }
    return(out)
  }
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makePSOCKcluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # library(), not requireNamespace(): the workers need GWmodel's exported
    # functions (gw.weight, gw.dist, ...) on their search path, and loading a
    # namespace without attaching it does not put them there.
    parallel::clusterEvalQ(cl, { suppressPackageStartupMessages(library(GWmodel)) })
    # parLapply serialises fit_fn together with its enclosing environment, so
    # x, y, st.dMat etc. travel to the workers on their own. What does NOT
    # travel is an argument still held as an unforced promise: it points at
    # the caller's environment (often the global env), which a worker does
    # not have -- hence "object 'use_adaptive_bw' not found". Forcing each
    # binding here caches the value in the frame so serialisation carries it.
    # Forcing in place rather than clusterExport()ing the frame keeps a
    # single copy per worker; st.dMat alone is 8*n^2 bytes, so shipping it
    # twice is what exhausts memory on large fits. try() skips formals that
    # were never supplied (e.g. regression.points, which gtwr() only ever
    # tests with missing()) -- those have no value to send.
    fn_env <- environment(fit_fn)
    for (.v in ls(fn_env, all.names = TRUE))
      try(get(.v, envir = fn_env), silent = TRUE)
    # The patched GWmodel R files are sys.source()d into the master's global
    # env rather than loaded from the package, and a worker's global env
    # starts empty, so those helpers must be shipped explicitly.
    def_env <- parent.env(fn_env)
    helpers <- ls(def_env, all.names = TRUE)
    helpers <- helpers[startsWith(helpers, ".gtwr_")]
    if (length(helpers))
      parallel::clusterExport(cl, varlist = helpers, envir = def_env)
    return(parallel::parLapply(cl, seq_len(n), fit_fn))
  }
  parallel::mclapply(seq_len(n), fit_fn, mc.cores = cores)
}

# Chunked variant of the above. chunk_fn(idx, dcols) fits a contiguous block
# of regression points; slice_fn(idx) produces just the columns of the
# distance matrix that block needs, so a worker never holds the whole 8*n^2
# matrix. chunk_fn is expected to carry an environment free of large objects
# (see gtwr()), otherwise the closure would drag them along regardless.
.gtwr_dispatch_chunks <- function(chunk_fn, n, cores, verbose, slice_fn)
{
  cores <- as.integer(cores)
  if (is.na(cores) || cores < 1L) cores <- 1L
  cores  <- max(1L, min(cores, n))
  nchunk <- cores
  chunks <- unname(split(seq_len(n),
                         as.integer(cut(seq_len(n), breaks = nchunk, labels = FALSE))))
  chunks <- chunks[lengths(chunks) > 0L]
  if (cores == 1L || length(chunks) == 1L) {
    pb <- if (isTRUE(verbose) && length(chunks) > 1L)
            utils::txtProgressBar(min = 0, max = length(chunks), style = 3) else NULL
    out <- vector("list", length(chunks))
    for (j in seq_along(chunks)) {
      out[[j]] <- chunk_fn(chunks[[j]], slice_fn(chunks[[j]]))
      if (!is.null(pb)) utils::setTxtProgressBar(pb, j)
    }
    if (!is.null(pb)) { close(pb); cat("
") }
    return(out)
  }
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makePSOCKcluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, { suppressPackageStartupMessages(library(GWmodel)) })
    def_env <- parent.env(environment(chunk_fn))
    helpers <- ls(def_env, all.names = TRUE)
    helpers <- helpers[startsWith(helpers, ".gtwr_")]
    if (length(helpers))
      parallel::clusterExport(cl, varlist = helpers, envir = def_env)
    # Ship each worker its slice one at a time and let the master drop its
    # copy immediately. Building the full list of slices up front (as
    # clusterMap requires) would hold a second complete n x n matrix on the
    # master -- 9.3 GB at n=34013, which is what tipped the full-size fit
    # over. This keeps the master to st.dMat plus a single chunk.
    # These setters are re-parented to the global env on purpose. Defined
    # inline they would close over THIS frame, which holds slice_fn, whose own
    # environment is gtwr()'s frame containing st.dMat -- so every send would
    # serialise the entire matrix again, exactly what the chunking avoids.
    put_fn    <- function(f) assign(".gtwr_chunk_fn", f, envir = globalenv())
    put_chunk <- function(ix, sl) {
      assign(".gtwr_idx",   ix, envir = globalenv())
      assign(".gtwr_slice", sl, envir = globalenv())
      invisible(NULL)
    }
    environment(put_fn)    <- globalenv()
    environment(put_chunk) <- globalenv()
    parallel::clusterCall(cl, put_fn, chunk_fn)
    for (w in seq_along(chunks)) {
      cw <- cl[w]; class(cw) <- class(cl)
      parallel::clusterCall(cw, put_chunk, chunks[[w]], slice_fn(chunks[[w]]))
    }
    return(parallel::clusterEvalQ(cl,
             .gtwr_chunk_fn(.gtwr_idx, .gtwr_slice)))
  }
  parallel::mcmapply(chunk_fn, chunks, lapply(chunks, slice_fn),
                     SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.cores = cores)
}
