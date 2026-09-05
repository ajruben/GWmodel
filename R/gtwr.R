#Geographically and temporally weighted regression
###########################################
#Bo Huang , Bo Wu & Michael Barry (2010) Geographically and temporally weighted regression for modeling spatio-temporal 
#variation in house prices, International Journal of Geographical Information Science, 24:3, 383-401
# Bo Wu, Rongrong Li & Bo Huang (2014) A geographically and temporally
# weighted autoregressive model with application to housing prices, International Journal of Geographical Information Science, 28:5, 1186-1204

#Then an optimal spatial bandwidth is specified for each time period based on a goodnessof-
#  fit criterion such as cross-validation (CV) or Akaike information criterion (AIC). Using the
# optimal spatial bandwidth, the optimal temporal bandwidth is determined, again based on CV or
# AIC. Once both optimal spatial and temporal bandwidths are derived, they can be used to
# construct the spatiotemporal weight matrix W, which allows local parameters to be estimated
# using equation (2).

###########################################
#Calibrate the GTWR model
gtwr<- function(formula, data, regression.points, obs.tv, reg.tv, st.bw, kernel="bisquare",
                 adaptive=FALSE, p=2, theta=0, longlat=F,lamda=0.05,t.units = "auto",ksi=0, st.dMat,
                 cores = 1L, verbose = interactive(),
                 aicc.rss.floor = 1e-8, aicc.enp.margin = 1)
{
  ##Record the start time
  timings <- list()
  timings[["start"]] <- Sys.time()
  ###################################macth the variables
  this.call <- match.call()
  p4s <- as.character(NA)
  polygons <- NULL
  sfdf <- TRUE
  ##Data points{
  if (inherits(data, "Spatial"))
  {
    p4s <- proj4string(data)
    dp.locat<-coordinates(data)
    if(is(data, "SpatialPolygonsDataFrame"))
       polygons <- polygons(data)
    data <- as(data, "data.frame")
    sfdf <- FALSE
  }
  else if(inherits(data, "sf")) {
    if(any((st_geometry_type(data)=="POLYGON")) | any(st_geometry_type(data)=="MULTIPOLYGON"))
       dp.locat <- st_coordinates(st_centroid(st_geometry(data)))
    else
       dp.locat <- st_coordinates(st_geometry(data))
  }
  else
  {
    stop("Given regression data must be a Spatial*DataFrame or sf object")
  }
  dp.n <- nrow(dp.locat)
  ####Check the time stamps given for the data
  if(missing(obs.tv))
  {
    stop("Please provide the corresponding time stamps for the observations!")
  }
  else
  {
    if(is(obs.tv, "Date")||is(obs.tv, "POSIXlt")||is(obs.tv, "POSIXct")||is(obs.tv, "numeric")||is(obs.tv, "yearmon")||is(obs.tv, "yearqtr"))
    {
      if(length(obs.tv)!=dp.n)
        stop("The given time stamps must correspond strictly to the observation data")
    }
    else
    {
      stop("Please provide the time stamps in accepted format: numeric, Date, POSIXlt, POSIXct, yearmon or yearqtr")
    }
  }
  #####Check the given data frame and regression points
  #####Regression points
  
  if (missing(regression.points))
  {
    rp.given <- FALSE
    rp.locat<-dp.locat
    hatmatrix<-T
    reg.tv <- obs.tv
    if(sfdf){
      regression.points <- data
    }
    else{
      if(is.null(polygons))
         regression.points <- SpatialPointsDataFrame(coords=rp.locat, data=data)
      else
        regression.points <-SpatialPolygonsDataFrame(Sr=polygons, data=data,match.ID=F)
    }
  }
  else 
  {
    rp.given <- TRUE
    hatmatrix<-F
     if(inherits(regression.points, "Spatial")) 
    {
      rp.locat<-coordinates(regression.points)
      if (is(regression.points, "SpatialPolygonsDataFrame"))
         polygons<-polygons(regression.points)
    }
    else if (inherits(regression.points, "sf"))
    {
      if (any((st_geometry_type(regression.points)=="POLYGON")) | any(st_geometry_type(regression.points)=="MULTIPOLYGON"))
         rp.locat <- st_coordinates(st_centroid(st_geometry(regression.points)))
      else
         rp.locat<- st_coordinates(st_centroid(st_geometry(regression.points)))
    }
    else if (is.numeric(regression.points) && dim(regression.points)[2] == 2)
      rp.locat<-regression.points
    else
    {
      warning("Output loactions are not packed in a Spatial object,and it has to be a two-column numeric vector")
      rp.given <- F
      rp.locat<-dp.locat
      reg.tv <- obs.tv
    }
    rp.n <- nrow(rp.locat)
    ###time stamps for regression locations
    if(missing(reg.tv))
      stop("Please provide the corresponding time stamps for the regression points!")
    else
    {
      if(is(reg.tv, "Date")||is(reg.tv, "POSIXlt")||is(reg.tv, "POSIXct")||is(reg.tv, "numeric")||is(reg.tv, "yearmon")||is(reg.tv, "yearqtr"))
      {
        if(length(reg.tv)!=rp.n)
          stop("The given time stamps must correspond strictly to the regression data")
      }
      else
      {
        stop("Please provide the time stamps in accepted format: numeric, Date, POSIXlt, POSIXct, yearmon or yearqtr")
      }
    }
  }
 #  if(class(obs.tv)==class(reg.tv))
 # {
 #   if(class(obs.tv)=="numeric" && t.units != "years")
#      t.units <- "auto"
#    if (class(obs.tv)=="yearmon")
#      t.units <- "months"
#    if (class(obs.tv)=="yearqtr")
#      t.units <- "quarters"
#  }
  ###
  ####################
  ######Extract the data frame
  ####Refer to the function lm
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  var.n<-ncol(x)
  rp.n<-nrow(rp.locat)
  dp.n<-nrow(data)
  betas <-matrix(nrow=rp.n, ncol=var.n)
  betas.SE <-matrix(nrow=rp.n, ncol=var.n)
  betas.TV <-matrix(nrow=rp.n, ncol=var.n)
  ##S: hatmatrix -- never materialised. Every use of S below is a
  ##   reduction over its rows, accumulated by the chunked dispatch.
  #C.M<-matrix(nrow=dp.n,ncol=dp.n)
  idx1 <- match("(Intercept)", colnames(x))
  if(!is.na(idx1))
    colnames(x)[idx1]<-"Intercept" 
  colnames(betas) <- colnames(x)
  #colnames(betas)[1]<-"Intercept"
  
  ##################################################
  #####Linear regression
  lm.res <- lm(formula,data=data)
  lm.res$x <- x
  lm.res$y <- y
  gTSS <- c(cov.wt(matrix(y, ncol=1), wt=rep(as.numeric(1), dp.n), method="ML")$cov*dp.n)
  #####GTWR
  ######### Spatial Distance matrix is given or not
  if (missing(st.dMat))
  {
    DM.given<-F
    if(dp.n+rp.n <10000)
    {
      if(rp.given)
      {
          st.dMat <- st.dist(dp.locat, rp.locat, obs.tv, reg.tv, p=p, theta=theta, longlat=longlat,lamda=lamda,t.units = t.units,ksi=ksi)
      }
      else
      {
          st.dMat <- st.dist(dp.locat, obs.tv=obs.tv, p=p, theta=theta, longlat=longlat, lamda=lamda,t.units = t.units,ksi=ksi)
      }
        
      DM.given <- T
    }
  }
  else
  {
    DM.given<-T
    dim.stdMat<-dim(st.dMat)
    if (dim.stdMat[1]!=dp.n||dim.stdMat[2]!=rp.n)
      stop("Dimensions of spatio-temporal distance matrix sdMat are not correct")
  }
  # Local fits are dispatched in contiguous chunks rather than one point at a
  # time. Two things fall out of that:
  #   * each worker is handed only its own columns of st.dMat, so per-worker
  #     memory is 8*n^2/cores rather than 8*n^2 -- broadcasting the whole
  #     matrix to every worker is what exhausted RAM on large fits;
  #   * the hat matrix is never materialised. Every downstream use of S is a
  #     reduction over its rows -- diag(S), sum(S^2), S %*% y, colSums(S^2) --
  #     so each chunk accumulates those and returns O(n) numbers instead of an
  #     n x n block travelling back over a socket.
  # fit_chunk is given an environment holding only the small objects it needs;
  # st.dMat is deliberately absent so it cannot ride along inside the closure.
  # mget() also forces the formals, which a PSOCK worker could not resolve.
  w_env <- list2env(mget(c("x", "y", "st.bw", "kernel", "adaptive", "hatmatrix",
                           "dp.n", "DM.given", "dp.locat", "rp.locat", "obs.tv",
                           "reg.tv", "p", "theta", "longlat", "lamda", "t.units", "ksi")),
                    parent = globalenv())
  fit_chunk <- function(idx, dcols) {
    m <- length(idx); nb <- ncol(x)
    o <- list(idx   = idx,
              betas = matrix(NA_real_, m, nb),
              se    = if (hatmatrix) matrix(NA_real_, m, nb) else NULL,
              sdiag = if (hatmatrix) numeric(m) else NULL,
              yhat  = if (hatmatrix) numeric(m) else NULL,
              colS2 = if (hatmatrix) numeric(dp.n) else NULL,
              sumS2 = 0, nfail = 0L, emsg = NULL)
    for (k in seq_len(m)) {
      i <- idx[k]
      st.disti <- if (DM.given) dcols[, k]
                  else st.dist(dp.locat, rp.locat, obs.tv, reg.tv, focus = i,
                               p = p, theta = theta, longlat = longlat, lamda = lamda,
                               t.units = t.units, ksi = ksi)
      f <- .gtwr_point_fit(i, x, y, st.disti, st.bw, kernel, adaptive, hatmatrix)
      if (!is.null(f$.failed)) {
        o$nfail <- o$nfail + 1L
        if (is.null(o$emsg)) o$emsg <- f$.failed
      }
      o$betas[k, ] <- f$beta
      if (hatmatrix) {
        sr <- f$S_row
        o$se[k, ]  <- f$se_sq
        o$sdiag[k] <- sr[i]
        o$yhat[k]  <- sum(sr * y)
        sr2        <- sr * sr
        o$sumS2    <- o$sumS2 + sum(sr2)
        o$colS2    <- o$colS2 + sr2
      }
    }
    o
  }
  environment(fit_chunk) <- w_env
  slice_fn <- if (DM.given) function(ix) st.dMat[, ix, drop = FALSE]
              else function(ix) NULL
  chunk_res <- .gtwr_dispatch_chunks(fit_chunk, rp.n, cores, verbose, slice_fn)
  nfail <- sum(vapply(chunk_res, function(r) r$nfail, integer(1)))
  if (nfail > 0L) {
    ex <- Filter(Negate(is.null), lapply(chunk_res, function(r) r$emsg))[[1]]
    warning(sprintf("Local fit failed at %d / %d regression points (bandwidth too small or design singular). First error: %s",
                    nfail, rp.n, ex))
  }
  if (hatmatrix) {
    s.diag <- numeric(dp.n); yhat.v <- numeric(dp.n)
    colS2  <- numeric(dp.n); sumS2  <- 0
  }
  for (r in chunk_res) {
    betas[r$idx, ] <- r$betas
    if (hatmatrix) {
      betas.SE[r$idx, ] <- r$se
      s.diag[r$idx]     <- r$sdiag
      yhat.v[r$idx]     <- r$yhat
      colS2             <- colS2 + r$colS2
      sumS2             <- sumS2 + r$sumS2
    }
  }
  ########################Diagnostic information
  
  GTW.diagnostic<-NA
  if (hatmatrix)
  {
    # S was never formed; the chunked dispatch accumulated its reductions.
    # With A = I - S,  y'A'A y == ||A y||^2 == sum(residual^2)  and
    # diag(A'A) == colSums(S^2) - 2*diag(S) + 1, so RSS and the leverage
    # terms cost O(n) here. The original built Q = A'A explicitly: an O(n^3)
    # matmul plus three n x n matrices. RSS.gw is kept as a 1x1 matrix so
    # downstream code indexing it as such still works.
    tr.S   <- sum(s.diag)
    tr.StS <- sumS2
    yhat<-yhat.v
    residual<-y-yhat
    RSS.gw<-matrix(sum(residual^2), 1L, 1L)
    edf.raw <- dp.n - 2 * tr.S + tr.StS
    if (!is.finite(edf.raw) || edf.raw <= aicc.enp.margin) {
      warning(sprintf("Effective degrees of freedom (%.2f) at or below %g; standard errors set to NA.",
                      edf.raw, aicc.enp.margin))
      sigma.hat1    <- NA_real_
      Stud_residual <- rep(NA_real_, dp.n)
      betas.SE[]    <- NA_real_
      betas.TV[]    <- NA_real_
    } else {
      sigma.hat1    <- RSS.gw / edf.raw
      q.diag        <- colS2 - 2 * s.diag + 1
      Stud_residual <- as.numeric(residual) / sqrt(as.numeric(sigma.hat1) * q.diag)
      betas.SE      <- sqrt(as.numeric(sigma.hat1) * betas.SE)
      betas.TV      <- betas / betas.SE
    }
    yss.g <- sum((y - mean(y))^2)
    RSS.num <- as.numeric(RSS.gw)
    if (!is.finite(RSS.num)) {
      RSS.eff <- NA_real_
    } else {
      RSS.eff <- max(RSS.num, aicc.rss.floor * yss.g)
      if (RSS.eff > RSS.num)
        warning(sprintf("Reported RSS (%.3g) is below %g * TSS; using floor %.3g for AIC/AICc.",
                        RSS.num, aicc.rss.floor, RSS.eff))
    }
    sigma.hat2 <- if (is.finite(RSS.eff)) RSS.eff / dp.n else NA_real_
    enp.max <- dp.n - 2 - aicc.enp.margin
    if (!is.finite(tr.S) || tr.S >= enp.max || !is.finite(RSS.eff)) {
      warning(sprintf("AIC/AICc unreliable at this bandwidth (tr(S)=%.2f, n - 2 - %g = %.2f, RSS.eff=%.3g); set to NA.",
                      as.numeric(tr.S), aicc.enp.margin, enp.max, as.numeric(RSS.eff)))
      AIC <- NA_real_; AICc <- NA_real_
    } else {
      AIC  <- dp.n * log(sigma.hat2) + dp.n * log(2 * pi) + dp.n + tr.S
      AICc <- dp.n * log(sigma.hat2) + dp.n * log(2 * pi) +
              dp.n * ((dp.n + tr.S) / (dp.n - 2 - tr.S))
    }
    edf <- dp.n - 2 * tr.S + tr.StS
    enp <- 2 * tr.S - tr.StS
    gw.R2 <- if (is.finite(RSS.num) && yss.g > 0) 1 - RSS.num / yss.g else NA_real_
    gwR2.adj <- if (is.finite(edf) && edf > 1 && is.finite(gw.R2))
                  1 - (1 - gw.R2) * (dp.n - 1) / (edf - 1)
                else NA_real_
    GTW.diagnostic <- list(RSS.gw = RSS.gw, RSS.eff = RSS.eff,
                           AIC = AIC, AICc = AICc, enp = enp, edf = edf,
                           gw.R2 = gw.R2, gwR2.adj = gwR2.adj)
  }
  
  ####encapsulate the GWR results
  GTW.arguments<-list(formula=formula,rp.given=rp.given,hatmatrix=hatmatrix,st.bw=st.bw,lamda=lamda, ksi=ksi,kernel=kernel,
                      adaptive=adaptive,p=p, theta=theta, longlat=longlat,DM.given=DM.given,units=units)
  
  #observed y: y
  #fitted y : yhat 
  #residual : y-yhat
  #Studentised residual: Stud_residual=residual.i/(sigma.hat*sqrt(q.ii))
  if (hatmatrix)                                         
  {
    gtwres.df<-data.frame(betas,y,yhat,residual,obs.tv, Stud_residual,betas.SE,betas.TV)
    colnames(gtwres.df)<-c(c(c(colnames(betas),c("y","yhat","residual","time_stamp","Stud_residual")),
                             paste(colnames(betas), "SE", sep="_")),paste(colnames(betas), "TV", sep="_"))
    
  }
  else
  {
    
    gtwres.df<-data.frame(betas, reg.tv)
    colnames(gtwres.df)<- c(colnames(betas),"time_stamp") 
  }
  rownames(rp.locat)<-rownames(gtwres.df)
  
  if(inherits(regression.points, "Spatial")) 
  {
     if (!is.null(polygons))
    {
       rownames(gtwres.df) <- sapply(slot(polygons, "polygons"),
                                  function(i) slot(i, "ID"))
       SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=gtwres.df,match.ID=F)
    }
    else
    {
      SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gtwres.df, proj4string=CRS(p4s), match.ID=F)
    }
  }
  else if(inherits(regression.points, "sf"))
  {
     SDF <- st_sf(gtwres.df, geometry = st_geometry(regression.points))
  }
  else
     SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gtwres.df, proj4string=CRS(p4s), match.ID=F)
  
  timings[["stop"]] <- Sys.time()
  ##############
  res<-list(GTW.arguments=GTW.arguments,GTW.diagnostic=GTW.diagnostic,lm=lm.res,SDF=SDF,
            timings=timings,this.call=this.call)
  class(res) <-"gtwrm"
  invisible(res)
}

.gtwr_point_fit <- function(i, x, y, st.disti, st.bw, kernel, adaptive, hatmatrix)
{
  W.i <- gw.weight(st.disti, st.bw, kernel, adaptive)
  res <- tryCatch(gw_reg(x, y, W.i, hatmatrix, i), error = function(e) e)
  if (inherits(res, "error")) {
    p <- ncol(x); n <- nrow(x)
    out <- list(beta = rep(NA_real_, p), .failed = conditionMessage(res))
    if (hatmatrix) { out$S_row <- rep(NA_real_, n); out$se_sq <- rep(NA_real_, p) }
    return(out)
  }
  out <- list(beta = res[[1]])
  if (hatmatrix) {
    out$S_row <- res[[2]]
    Ci        <- res[[3]]
    out$se_sq <- diag(Ci %*% t(Ci))
  }
  out
}

############################Layout function for outputing the GWR results
##Author: BL
print.gtwrm<-function(x, ...)
{
  if(!inherits(x, "gtwrm")) stop("It's not a gwm object")
  cat("   ***********************************************************************\n")
  cat("   *                       Package   GWmodel                             *\n")
  cat("   ***********************************************************************\n")
  cat("   Program starts at:", as.character(x$timings$start), "\n")
  cat("   Call:\n")
  cat("   ")
  print(x$this.call)
  vars<-all.vars(x$GTW.arguments$formula)
  var.n<-length(x$lm$coefficients)
  cat("\n   Dependent (y) variable: ",vars[1])
  cat("\n   Independent variables: ",vars[-1])
  dp.n<-length(x$lm$residuals)
  cat("\n   Number of data points:",dp.n)
  ################################################################ Print Linear
  cat("\n   ***********************************************************************\n")
  cat("   *                    Results of Global Regression                     *\n")
  cat("   ***********************************************************************\n")
  print(summary.lm(x$lm))
  cat("   ***Extra Diagnostic information\n")
  lm_RSS<-sum(x$lm$residuals^2)
  lm_Rank<-x$lm$rank     
  cat("   Residual sum of squares:", lm_RSS)
  #lm_sigma<-sqrt(lm_RSS/(dp.n-lm_Rank-2))
  lm_sigma<-sqrt(lm_RSS/(dp.n-2))
  cat("\n   Sigma(hat):", lm_sigma)
  lm_AIC<-dp.n*log(lm_RSS/dp.n)+dp.n*log(2*pi)+dp.n+2*(var.n + 1)
  #AIC = dev + 2.0 * (double)(MGlobal + 1.0);
  cat("\n   AIC: ", lm_AIC)
  ##AICc = 	dev + 2.0 * (double)N * ( (double)MGlobal + 1.0) / ((double)N - (double)MGlobal - 2.0);
  lm_AICc= dp.n*log(lm_RSS/dp.n)+dp.n*log(2*pi)+dp.n+2*dp.n*(var.n+1)/(dp.n-var.n-2)
  cat("\n   AICc: ", lm_AICc)
  #lm_rdf <- x$dfsidual
  
  #########################################################################
  cat("\n   ***********************************************************************\n")
    cat("   *    Results of Geographically and Temporally Weighted Regression     *\n")
  cat("   ***********************************************************************\n")
  cat("\n   *********************Model calibration information*********************\n")
  cat("   Kernel function for geographically and temporally weighting:", x$GTW.arguments$kernel, "\n")
  if(x$GTW.arguments$adaptive)
    cat("   Adaptive bandwidth for geographically and temporally  weighting: ", x$GTW.arguments$st.bw, " (number of nearest neighbours)\n", sep="")
  else
    cat("   Fixed bandwidth for geographically and temporally weighting: ", x$GTW.arguments$st.bw, "\n")
  if(x$GTW.arguments$rp.given) 
    cat("   Regression points: A seperate set of regression points is used.\n")
  else
    cat("   Regression points: the same locations as observations are used.\n")
  if (x$GTW.arguments$DM.given)
    cat("   Distance metric for geographically and temporally  weighting: A distance matrix is specified for this model calibration.\n")
  else
  {
    if (x$GTW.arguments$longlat)
      cat("   Distance metric for geographically weighting: Great Circle distance metric is used.\n")
    else if (x$GTW.arguments$p==2)
      cat("   Distance metric for geographically weighting: Euclidean distance metric is used.\n")
    else if (x$GTW.arguments$p==1)
      cat("   Distance metric for geographically weighting: Manhattan distance metric is used.\n") 
    else if (is.infinite(x$GTW.arguments$p))
      cat("   Distance metric for geographically weighting: Chebyshev distance metric is used.\n")
    else 
      cat("   Distance metric for geographically weighting: A generalized Minkowski distance metric is used with p=",x$GTW.arguments$p,".\n")
    if (x$GTW.arguments$theta!=0&&x$GTW.arguments$p!=2&&!x$GTW.arguments$longlat)
      cat("   Coordinate rotation: The coordinate system is rotated by an angle", x$GTW.arguments$theta, "in radian.\n")
    cat("   The temporal distance is calculated in ", x$GTW.arguments$units, ".\n")
    cat("   The adjustment parameter for calculating spatio-temporal distances lamda is: ", x$GTW.arguments$lamda, ".\n")
    cat("   The adjustment parameter for calculating spatio-temporal distances ksi is: ", x$GTW.arguments$ksi, "in radian.\n")    
  } 
  
  cat("\n   ****************Summary of GTWR coefficient estimates:*****************\n")       
  if(inherits(x$SDF, "Spatial"))
       df0 <- as(x$SDF, "data.frame")[,1:var.n, drop=FALSE]
    else
       df0 <- st_drop_geometry(x$SDF)[,1:var.n, drop=FALSE]
  if (any(is.na(df0))) {
    df0 <- na.omit(df0)
    warning("NAs in coefficients dropped")
  }
  CM <- t(apply(df0, 2, summary))[,c(1:3,5,6)]
  if(var.n==1) 
  { 
    CM <- matrix(CM, nrow=1)
    colnames(CM) <- c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.")
    rownames(CM) <- names(x$SDF)[1]
  }
  rnames<-rownames(CM)
  for (i in 1:length(rnames))
    rnames[i]<-paste("   ",rnames[i],sep="")
  rownames(CM) <-rnames 
  printCoefmat(CM)
  if (x$GTW.arguments$hatmatrix) 
  {	
    cat("   ************************Diagnostic information*************************\n")
    cat("   Number of data points:", dp.n, "\n")
    cat("   Effective number of parameters (2trace(S) - trace(S'S)):", x$GTW.diagnostic$enp, "\n")
    cat("   Effective degrees of freedom (n-2trace(S) + trace(S'S)):", x$GTW.diagnostic$edf, "\n")
    cat("   AICc (GWR book, Fotheringham, et al. 2002, p. 61, eq 2.33):",
        x$GTW.diagnostic$AICc, "\n")
    cat("   AIC (GWR book, Fotheringham, et al. 2002,GWR p. 96, eq. 4.22):", x$GTW.diagnostic$AIC, "\n")
    cat("   Residual sum of squares:", x$GTW.diagnostic$RSS.gw, "\n")
    cat("   R-square value: ",x$GTW.diagnostic$gw.R2,"\n")
    cat("   Adjusted R-square value: ",x$GTW.diagnostic$gwR2.adj,"\n")	
  }
  cat("\n   ***********************************************************************\n")
  cat("   Program stops at:", as.character(x$timings$stop), "\n")
  invisible(x)
}

ti.distv <- function(focal.t, obs.tv, units="auto")
{
  d <- ti.dist(obs.tv, focal.t, units = units)
  d[obs.tv > focal.t] <- Inf
  d
}
ti.distm <- function(obs.tv, reg.tv, units="auto")
{
  if (missing(reg.tv)) reg.tv <- obs.tv
  n <- length(obs.tv); m <- length(reg.tv)
  dist.tm <- matrix(numeric(m * n), nrow = n)
  for (i in seq_len(m))
    dist.tm[, i] <- ti.distv(reg.tv[i], obs.tv, units)
  dist.tm
}
#calculate the time distance
#units can be "auto", "secs", "mins", "hours","days", "weeks","months","years"
ti.dist <- function(t1,t2,units="auto")
{
  tcl <- class(t1)[1]
  switch(tcl,
         Date = abs(as.numeric(difftime(t1,t2,units = units))),
         POSIXlt = abs(as.numeric(difftime(t1,t2,units = units))),
         POSIXct = abs(as.numeric(difftime(t1,t2,units = units))),
         numeric = abs(t2-t1),
         integer = abs(t2-t1),
         yearmon = abs((t2-t1)*12),
         yearqtr = abs((t2-t1)*4)
         )
}
#Get the forward time stamps before a specific time stamp
get.before.ti <- function(ti, ts, time.lag,units)
{
  ts.ti <- c()
  idx <- which(ts==ti)
  idx <- idx-1
  while(idx>0 && ti.dist(ti, ts[idx], units)<time.lag)
  {
    ts.ti <- c(ts.ti,ts[idx])
    idx <- idx-1
  }
  ts.ti
}
get.ts <- function(tv)
{
  ts <- sort(unique(tv))
  list(ts = ts, index = match(tv, ts))
}
get.uloat <- function(coords)
{
  TF.dup   <- duplicated(coords)
  ucoords  <- coords[!TF.dup, , drop = FALSE]
  keys.all <- paste(coords[, 1],  coords[, 2],  sep = "\r")
  keys.u   <- paste(ucoords[, 1], ucoords[, 2], sep = "\r")
  list(ucoords, match(keys.all, keys.u))
}
##Calculate the spatial distance matrix with duplicated locations removed
sdist.mat <- function(dp.locat, rp.locat, p=2, theta=0, longlat=F)
{
   if (missing(dp.locat)||!is.numeric(dp.locat)||dim(dp.locat)[2]!=2)
      stop("Please input correct coordinates of data points")
  
   if (!missing(rp.locat)) 
   {
       rp.given<-T
   }
   else
   {
       rp.given<-F 
       rp.locat<- dp.locat
   } 
   if (!is.numeric(rp.locat))
      stop("Please input correct coordinates of regression points")
   else
      rp.locat <- matrix(rp.locat, ncol=2)
   ###Remove the duplicated locations
   ucoord.dp <- get.uloat(dp.locat) 
   coord.dp.idx <- ucoord.dp[[2]]
   ucoord.dp <- ucoord.dp[[1]]
   ucoord.rp <- get.uloat(rp.locat) 
   coord.rp.idx <- ucoord.dp[[2]]
   ucoord.rp <- ucoord.rp[[1]]
   if(rp.given)
   {
     dists <- gw.dist(ucoord.dp, ucoord.rp, p=p, theta=theta, longlat=longlat) 
   }
   else
   {
     dists <- gw.dist(ucoord.dp, p=p, theta=theta, longlat=longlat) 
   }
   dists
}
tdist.mat <- function(obs.tv, reg.tv, t.units = "auto")
{
   uts.obv <- get.ts(obs.tv)
   uts.obv.idx <- uts.obv[[2]]
   uts.obv <- uts.obv[[1]]
   if(missing(reg.tv))
      dists <- ti.distm(uts.obv, units=t.units)
   else
   {
     uts.reg <- get.ts(reg.tv)
     uts.reg.idx <- uts.reg[[2]]
     uts.reg <- uts.reg[[1]]
     dists <- ti.distm(uts.obv,uts.reg, units=t.units)
   }
   dists
}
#Calculate the ST distance matrix
st.dist <- function(dp.locat, rp.locat, obs.tv, reg.tv,focus=0, p=2, theta=0, longlat=F,lamda=0.05,t.units = "auto",ksi=0, s.dMat,t.dMat)
{
   if (missing(dp.locat)||!is.numeric(dp.locat)||dim(dp.locat)[2]!=2)
      stop("Please input correct coordinates of data points")
  
   if (!missing(rp.locat)) 
   {
       rp.given<-T
   }
   else
   {
       rp.given<-F 
       rp.locat<- dp.locat
   } 
   if (!is.numeric(rp.locat))
      stop("Please input correct coordinates of regression points")
   else
      rp.locat <- matrix(rp.locat, ncol=2)
   if (focus<0||focus>length(rp.locat[,1]))
      stop("No regression point is fixed")

   n.rp<-length(rp.locat[,1])
   n.dp<-length(dp.locat[,1])
   
   if (focus>0)
       dists<-numeric(n.dp) 
   else
       dists<-matrix(numeric(n.rp*n.dp),nrow=n.dp)
   ###Remove the duplicated locations
   ucoord.dp <- get.uloat(dp.locat) 
   coord.dp.idx <- ucoord.dp[[2]]
   ucoord.dp <- ucoord.dp[[1]]
   ucoord.rp <- get.uloat(rp.locat) 
   coord.rp.idx <- ucoord.rp[[2]]
   ucoord.rp <- ucoord.rp[[1]]
   ######Calculate the spatial distance matrix
   if(missing(s.dMat))
   {
      if(rp.given)
        s.dMat <- gw.dist(ucoord.dp, ucoord.rp, p=p, theta=theta, longlat=longlat)
      else
        s.dMat <- gw.dist(ucoord.dp, p=p, theta=theta, longlat=longlat) 
   }
   else
   {
      if(!(dim(s.dMat)[1]==nrow(ucoord.dp)&&dim(s.dMat)[2]==nrow(ucoord.rp)))
        stop("s.dMat is of dimnensions with duplicated locations removed")
   }
   
   ####Calculate the temporal distance matrix
   uts.obv <- get.ts(obs.tv)
   uts.obv.idx <- uts.obv[[2]]
   uts.obv <- uts.obv[[1]]
   if(rp.given)
   {
     uts.reg <- get.ts(reg.tv)
     uts.reg.idx <- uts.reg[[2]]
     uts.reg <- uts.reg[[1]]
   }
   if(missing(t.dMat))
   {
      if(rp.given)
        t.dMat <- ti.distm(uts.obv,uts.reg, units=t.units)
      else
        t.dMat <- ti.distm(uts.obv, units=t.units)
   }
   ####Calculate the distance matrix
   if(focus>0)
   {
     if(rp.given)
     {
       s_vec <- s.dMat[coord.dp.idx, coord.rp.idx[focus]]
       t_vec <- t.dMat[uts.obv.idx, uts.reg.idx[focus]]
     }
     else
     {
       s_vec <- s.dMat[coord.dp.idx, coord.dp.idx[focus]]
       t_vec <- t.dMat[uts.obv.idx, uts.obv.idx[focus]]
     }
     finite <- is.finite(t_vec)
     dists[] <- Inf
     sf <- s_vec[finite]; tf <- t_vec[finite]
     dists[finite] <- lamda*sf + (1-lamda)*tf +
                      2*sqrt(lamda*(1-lamda)*sf*tf)*cos(ksi)
   }
   else
   {
     # Built one column block at a time. Holding S, Tm, t(Tm), the two logical
     # masks and the arithmetic temporaries all at n x n peaked at ~57 GB to
     # produce an 8.6 GB result at n=34013 -- the reason full-size runs could
     # not start. s.dMat and t.dMat only ever hold unique coordinates and
     # unique time stamps (2918 and 12 respectively for a London LSOA month
     # panel), so a block re-expands a small matrix rather than keeping a
     # second copy of the large one. The arithmetic is left verbatim so the
     # result is bit-identical to the unblocked version.
     blk  <- max(1L, min(n.rp, as.integer(ceiling(1e7 / max(1L, n.dp)))))
     rows <- seq_len(n.dp)
     for (.start in seq(1L, n.rp, by = blk))
     {
       cols <- .start:min(.start + blk - 1L, n.rp)
       m    <- length(cols)
       if(rp.given)
       {
         S  <- s.dMat[coord.dp.idx, coord.rp.idx[cols], drop=FALSE]
         Tm <- t.dMat[uts.obv.idx, uts.reg.idx[cols], drop=FALSE]
       }
       else
       {
         S  <- s.dMat[coord.dp.idx, coord.dp.idx[cols], drop=FALSE]
         # Tm[i,j] <- t.dMat[a[min(i,j)], a[max(i,j)]] reproduces the
         # lower.tri(Tm) <- t(Tm)[lower.tri(Tm)] symmetrisation exactly,
         # without ever forming the n x n transpose.
         I  <- matrix(rows, nrow=n.dp, ncol=m)
         Cc <- matrix(cols, nrow=n.dp, ncol=m, byrow=TRUE)
         a  <- uts.obv.idx
         Tm <- matrix(t.dMat[cbind(a[pmin(I,Cc)], a[pmax(I,Cc)])], n.dp, m)
       }
       finite <- is.finite(Tm)
       db <- matrix(Inf, n.dp, m)
       sf <- S[finite]; tf <- Tm[finite]
       db[finite] <- lamda*sf + (1-lamda)*tf +
                     2*sqrt(lamda*(1-lamda)*sf*tf)*cos(ksi)
       dists[, cols] <- db
     }
   }
   dists
}

