### poor person's kde

#' @import stats
#' @import graphics
#' @import utils
NULL

#' Kernel Density Estimate
#'
#' A pure R implementation of a one-dimensional KDE with a
#' gaussian kernel.  Mainly of programming interest.
#'
#' @param y Numeric vector of observed values
#' @param bw Numeric value: bandwidth, standard deviation of the kernel density
#' @param N  integer: number of equally spaced points for the kde
#' @param lower Numeric value: lower range limit for the kde
#' @param upper Numeric value: upper range limit for the kde
#' @param an obeject of class "kde"
#' @param xlab,ylab,col base graphics parameters
#' @param restore_graphics_state logical: should the prior graphics
#'        state be restored?
#' @param ...  Extra arguments passed on to methods.
#'
#' @return a list object with results of the kde estimation.
#' @export
#'
#' @examples
#' set.seed(1234)
#' .x <- rchisq(1000, 2)
#' kx <- kde(.x, lower = 0)
#' kx
#' plot(kx, restore_graphics_state = FALSE, ylim = c(0,0.5))
#' curve(dchisq(x, 2), add = TRUE, col = "grey", n = 500)
#' rug(.x, col = 2)
#' rm(.x)
kde <- function(y, bw = bw.nrd0(y), N = 512,
                lower = ry[1] - 3*bw, upper = ry[2] + 3*bw) {
  data_name <- deparse(substitute(y))
  y <- na.omit(y)
  ry <- range(y)
  dx <- (upper - lower)/N
  fr <- tabulate(pmax(1, pmin(N, 1+floor((y - lower)/dx))), nbins = N)
  x <- lower - dx/2 + dx*(1:N)
  yo <- numeric(N)
  for(i in 1:N) {
    w <- dnorm(x, x[i], bw)
    yo[i] <- sum(fr*w)/sum(w)
  }
  structure(list(x = x, y = yo/(sum(yo)*dx), bw = bw, n = length(y),
                 lower = lower, upper = upper,
                 data_name = data_name),
            class = "kde")
}

#' @rdname kde
#' @export
plot.kde <- function(x, ..., col = "steel blue",
                     xlab = bquote(x == italic(.(x$data_name))),
                     ylab = expression(kde(italic(x))),
                     restore_graphics_state = TRUE) {
  if(restore_graphics_state) {
    oldPar <- par(no.readonly = TRUE)
    on.exit(par(oldPar))
  }
  par(mar = c(4, 5.5, 1, 1) + 0.1, las = 1,
      cex.axis = 0.7, cex.lab = 1)
  with(x, {
    plot(x = c(lower, x, upper), y = c(0, y, 0),
         axes = FALSE, type = "l", xlab = xlab, ylab = ylab,
         col = col, panel.first = grid(), ...)
    axis(1)
    axis(2)
  })
  invisible(x)
}
# .S3method("plot", "kde")

#' @rdname kde
#' @export
print.kde <- function(x, ...) {
  cat("A Kernel Density Estimate\n\n")
  labs <- c("Response", "Sample size", "Range", "Bandwidth")
  values <- with(x, c(data_name, n,
                      paste0(signif(lower, 4), " < x < ",
                             signif(upper, 4)),
                      signif(bw, 4)))
  cat(paste0(format(labs), ": ", values), sep = "\n")
  generics <- attr(methods(class = "kde"), "info")$generic
  cat("\nMethods exist for generics:",
      paste(generics, collapse = ", "), "\n")
  invisible(x)
}
# .S3method("print", "kde")

#' Two Dimensional Kernel Density Estimate
#'
#' A pure R implementation of two direct algorithms.
#'
#' @param x,y Numeric vectors of the response Variables
#' @param N Integer: number of equally spaced points for the kde.  May be one or two for x and y respectively
#' @param bw Numeric: bandwidth(s) for x and y kernels respectively.  A single value will be duplicated.
#' @param lower_x,upper_x,lower_y,upper_y Numeric: kde range limits
#' @param algorithm1 Logical: should the smaller sample algorithm be used?
#' @param k Ineger: number of bandwidths by which to increase the range for the default limits.
#' @param x An object of class kde_2d
#' @param xlab,ylab,col base graphics parameters
#' @param restore_graphics_state logical: should the prior graphics
#'        state be restored?
#' @param ...  Extra arguments passed on to methods.
#'
#' @return An object of class "kde_2d"
#' @export
#'
#' @examples
#' set.seed(1234)
#' .time <- rnorm(10000, 10)
#' .money <- runif(.time) + .time/10
#' kxy <- kde_2d(.time, .money, algorithm1 = TRUE)
#' kxy
#' plot(kxy, restore_graphics_state = FALSE)
#' contour(kxy, add = TRUE)
#' rm(.time, .money)
kde_2d <- function(x, y, N = 128, bw = c(x = bw.nrd0(x), y = bw.nrd0(y)),
                   lower_x = limits$x$lower, upper_x = limits$x$upper,
                   lower_y = limits$y$lower, upper_y = limits$y$upper,
                   algorithm1 = length(x) < 2*Nxy, k = 1) {
  data_name_x <- deparse(substitute(x))
  data_name_y <- deparse(substitute(y))
  N <- setNames(rep(N, length.out = 2), c("x", "y"))
  if(any(nas <- is.na(x) | is.na(y))) {
    x <- x[!nas]
    y <- y[!nas]
  }
  bw <- setNames(rep(bw, length.out = 2), c("x", "y"))
  rx <- range(x)
  ry <- range(y)
  limits <- list(x = list(lower = rx[1] - k*bw[["x"]],
                          upper = rx[2] + k*bw[["x"]]),
                 y = list(lower = ry[1] - k*bw[["y"]],
                          upper = ry[2] + k*bw[["y"]]))
  Nxy <- N[["x"]]*N[["y"]]
  dx <- (upper_x - lower_x)/N[["x"]]
  dy <- (upper_y - lower_y)/N[["y"]]
  xo <- lower_x - dx/2 + dx*(1:N[["x"]])
  yo <- lower_y - dy/2 + dy*(1:N[["y"]])
  if(algorithm1) {
    z <- matrix(0, N[["x"]], N[["y"]])
    for(i in seq_along(x)) {
      z <- z + outer(dnorm(xo, x[i], bw[["x"]]),
                     dnorm(yo, y[i], bw[["y"]]))
    }
  } else {
    fr <- tabulate(1 + floor((x - lower_x)/dx) +
                     N[["x"]]*floor((y - lower_y)/dy),
                   nbins = Nxy)
    z <- with(expand.grid(x = xo, y = yo), {
      z <- numeric(Nxy)
      for(i in 1:Nxy) {
        w <- dnorm(x, x[i], bw[["x"]])*dnorm(y, y[i], bw[["y"]])
        z[i] <- sum(w*fr)/sum(w)
      }
      z
    })
    z <- matrix(z/(sum(z)*dx*dy), N[["x"]], N[["y"]])
  }
  structure(list(x = xo, y = yo, z = z, bw = bw, n = length(x),
                 lower = c(x = lower_x, y = lower_y),
                 upper = c(x = upper_x, y = upper_y),
                 data_name = c(x = data_name_x, y = data_name_y)),
            class = "kde_2d")
}

#' @rdname kde_2d
#' @export
print.kde_2d <- function(x, ...) {
  cat("A Two-dimensional Kernel Density Estimate\n\n")
  labs <- c("Responses", "Sample size", "Ranges", "Bandwidths")
  values <- with(x, c(paste(data_name, collapse = ", "),
                      n,
                      paste(paste0(signif(lower["x"], 4), " < x < ",
                                   signif(upper["x"], 4)),
                            paste0(signif(lower["y"], 4), " < y < ",
                                   signif(upper["y"], 4)),
                            sep = ", "),
                      paste(signif(bw, 4), collapse = ", ")))
  cat(paste0(format(labs), ": ", values), sep = "\n")
  generics <- attr(methods(class = "kde_2d"), "info")$generic
  cat("\nMethods exist for generics:",
      paste(generics, collapse = ", "), "\n")
  invisible(x)
}
# .S3method("print", "kde_2d")

#' @rdname kde_2d
#' @export
plot.kde_2d <- function(x, ...,
                        xlab = bquote(italic(.(x$data_name[["x"]]))),
                        ylab = bquote(italic(.(x$data_name[["y"]]))),
                        col = hcl.colors(25, rev = TRUE),
                        restore_graphics_state = TRUE) {
  if(restore_graphics_state) {
    oldPar <- par(no.readonly = TRUE)
    on.exit(par(oldPar))
  }
  par(mar = c(4, 4, 1, 1) + 0.1, las = 1, cex.axis = 0.7)
  with(x, image(x, y, z, xlab = xlab, ylab = ylab, col = col, ...))
  invisible(x)
}
# .S3method("plot", "kde_2d")
