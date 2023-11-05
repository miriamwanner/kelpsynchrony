# my function for plotting the spline correlogram

my_spline_plot <- function (x, ylim = c(-1, 1), add = FALSE, ...) 
{
  args.default <- list(xlab = "Distance (km)", ylab = "Correlation")
  args.input <- list(...)
  args <- c(args.default[!names(args.default) %in% names(args.input)], 
            args.input)
  cbar <- x$real$cbar
  if (!add) {
    do.call(plot, c(list(x = x$real$predicted$x, y = x$real$predicted$y, 
                         ylim = ylim, type = "l"), args))
  }
  if (!is.null(x$boot$boot.summary)) {
    polygon(c(x$boot$boot.summary$predicted$x, rev(x$boot$boot.summary$predicted$x)), 
            c(x$boot$boot.summary$predicted$y["0.025", ], rev(x$boot$boot.summary$predicted$y["0.975", 
            ])), col = gray(0.8), lty = 0)
  }
  lines(x$real$predicted$x, x$real$predicted$y)
  lines(c(0, max(x$real$predicted$x)), c(0, 0))
  # lines(c(0, max(x$real$predicted$x)), c(cbar, cbar))
}







my_sncf <- function (x, z, w = NULL, df = NULL, type = "boot", resamp = 1000, 
          npoints = 300, save = FALSE, filter = FALSE, fw = 0, max.it = 25, 
          xmax = FALSE, na.rm = FALSE, latlon = FALSE, circ = FALSE, 
          quiet = FALSE) 
{
  NAO <- FALSE
  if (any(!is.finite(unlist(z)))) {
    if (na.rm) {
      warning("Missing values exist; Pairwise deletion will be used")
      NAO <- TRUE
    }
    else {
      stop("Missing values exist; use na.rm = TRUE for pairwise deletion")
    }
  }
  if (is.null(w)) {
    n <- dim(z)[1]
    p <- dim(z)[2]
    z <- as.matrix(z) + 0
    moran <- cor2(t(z), circ = circ)
  }
  else {
    n <- dim(z)[1]
    p <- dim(z)[2]
    z <- as.matrix(z) + 0
    w <- as.matrix(w) + 0
    moran <- cor2(t(z), t(w), circ = circ)
  }
  if (is.null(df)) {
    df <- sqrt(n)
  }
  if (latlon) {
    xdist <- x
    print(xdist)
    # xdist <- gcdist(x, y)
  }
  else {
    xdist <- x
    # xdist <- sqrt(outer(x, x, "-")^2 + outer(y, y, "-")^2)
  }
  maxdist <- ifelse(!xmax, max(na.omit(xdist)), xmax)
  if (is.null(w)) {
    triang <- lower.tri(xdist)
  }
  else {
    triang <- is.finite(xdist)
  }
  u <- xdist[triang]
  v <- moran[triang]
  sel <- is.finite(v) & is.finite(u)
  u <- u[sel]
  v <- v[sel]
  v <- v[u <= maxdist]
  u <- u[u <= maxdist]
  xpoints <- seq(0, maxdist, length = npoints)
  out <- gather(u = u, v = v, w = w, moran = moran, df = df, 
                xpoints = xpoints, filter = filter, fw = fw)
  real <- list(cbar = out$cbar, x.intercept = out$xint, e.intercept = out$eint, 
               y.intercept = out$yint, cbar.intercept = out$cint, predicted = list(x = matrix(out$x, 
                                                                                              nrow = 1), y = matrix(out$y, nrow = 1)))
  boot <- list(NULL)
  boot$boot.summary <- list(NULL)
  if (resamp != 0) {
    boot$boot.summary$x.intercept <- matrix(NA, nrow = resamp, 
                                            ncol = 1)
    boot$boot.summary$y.intercept <- matrix(NA, nrow = resamp, 
                                            ncol = 1)
    boot$boot.summary$e.intercept <- matrix(NA, nrow = resamp, 
                                            ncol = 1)
    boot$boot.summary$cbar.intercept <- matrix(NA, nrow = resamp, 
                                               ncol = 1)
    boot$boot.summary$cbar <- matrix(NA, nrow = resamp, ncol = 1)
    predicted <- list(x = matrix(NA, nrow = 1, ncol = npoints), 
                      y = matrix(NA, nrow = resamp, ncol = npoints))
    type <- charmatch(type, c("boot", "perm"), nomatch = NA)
    if (is.na(type)) 
      stop("method should be \"boot\", or \"perm\"")
    for (i in 1:resamp) {
      whn <- pretty(c(1, resamp), n = 10)
      if (!quiet & any(i == whn)) {
        cat(i, " of ", resamp, "\r")
        flush.console()
      }
      if (type == 1) {
        trekkx <- sample(1:n, replace = TRUE)
        trekky <- trekkx
      }
      if (type == 2) {
        trekky <- sample(1:n, replace = FALSE)
        trekkx <- 1:n
      }
      xdistb <- xdist[trekkx, trekkx]
      if (is.null(w)) {
        triang <- lower.tri(xdistb)
      }
      else {
        triang <- is.finite(xdistb)
      }
      xdistb <- xdistb[triang]
      moranb <- moran[trekky, trekky][triang]
      if (type == 1 & is.null(w)) {
        moranb <- moranb[!(xdistb == 0)]
        xdistb <- xdistb[!(xdistb == 0)]
      }
      u <- xdistb
      v <- moranb
      sel <- is.finite(v) & is.finite(u)
      u <- u[sel]
      v <- v[sel]
      v <- v[u <= maxdist]
      u <- u[u <= maxdist]
      out <- gather(u = u, v = v, w = w, moran = moranb, 
                    df = df, xpoints = xpoints, filter = filter, 
                    fw = fw)
      boot$boot.summary$cbar[i, 1] <- out$cbar
      boot$boot.summary$y.intercept[i, 1] <- out$yint
      boot$boot.summary$x.intercept[i, 1] <- out$xint
      boot$boot.summary$e.intercept[i, 1] <- out$eint
      boot$boot.summary$cbar.intercept[i, 1] <- out$cint
      predicted$x[1, ] <- out$x
      predicted$y[i, ] <- out$y
    }
    if (save == TRUE) {
      boot$boot <- list(predicted = predicted)
    }
    else {
      boot$boot <- NULL
    }
    ty <- apply(predicted$y, 2, quantile, probs = c(0, 0.025, 
                                                    0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1), 
                na.rm = TRUE)
    dimnames(ty) <- list(c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 
                           0.75, 0.9, 0.95, 0.975, 1), NULL)
    tx <- predicted$x
    boot$boot.summary$predicted <- list(x = tx, y = ty)
  }
  else {
    boot <- NULL
    boot.summary <- NULL
  }
  res <- list(real = real, boot = boot, max.distance = maxdist, 
              call = deparse(match.call()))
  class(res) <- "Sncf"
  res
}


