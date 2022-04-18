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