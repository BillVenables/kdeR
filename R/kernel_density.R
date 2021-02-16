### poor person's kde
###

#' @useDynLib kdeR
#' @import Rcpp
#' @import stats
#' @import graphics
#' @importFrom grDevices hcl.colors
#' @importFrom utils methods
NULL

#' Kernel Density Estimate
#'
#' A pure R implementation of a one-dimensional KDE with a
#' gaussian kernel.  Mainly of programming interest.
#'
#' @param y Numeric vector of observed values
#' @param bw Numeric value: bandwidth, standard deviation of
#'        the kernel density
#' @param N  integer: number of equally spaced points for the kde
#' @param lower Numeric value: lower range limit for the kde
#' @param upper Numeric value: upper range limit for the kde
#' @param x an obeject of class "kde"
#' @param xlab,ylab,col base graphics parameters
#' @param row.names,optional As for the as.data.frame generic
#' @param responseName character string: name for the response
#'         variable
#' @param ...  Extra arguments passed on to methods.
#'
#' @return a list object with results of the kde estimation.
#' @export
#'
#' @examples
#' set.seed(1234)
#' x_value <- rchisq(10000, 2)
#' kx <- kde(x_value, lower = 0)
#' kx
#' if("package:ggplot2" %in% search()) {
#'   ggplot(as.data.frame(kx)) +
#'     aes(x = x_value, y = density) +
#'     geom_line(colour = "steel blue") +
#'     geom_function(fun = function(x) dchisq(x, 2),
#'                   colour = "firebrick3") +
#'     theme_bw() +
#'     xlab(expression(x == italic(x_value))) +
#'     ylab(expression(kde(italic(x))))
#' } else {
#'   oldPar <- par(mar = c(3, 3, 1, 1) + 0.1, las = 1, mgp = c(1.5, 0.5, 0),
#'                 cex.axis = 0.7, cex.lab = 0.8, tck = -0.015)
#'   plot(kx, ylim = c(0,0.5))
#'   curve(dchisq(x, 2), add = TRUE, col = 2, n = 1000)
#'   rug(x_value, col = "dark green")
#'   box()
#'   par(oldPar)
#' }
#' rm(oldPar)
kde <- function(y, bw = bw.nrd0(y), N = 512,
                lower = ry[1] - 3*bw, upper = ry[2] + 3*bw) {
  data_name <- deparse(substitute(y))
  y <- na.omit(y)
  n <- length(y)
  ry <- range(y)
  dx <- (upper - lower)/N
  x <- lower - dx/2 + dx*(1:N)
  yo <- numeric(N)
  # if(!missing(lower) || !missing(upper)) {
  #   lims <- c(lower, upper)
  #   for(i in 1:n) {
  #     yo <- yo + dnorm(x, y[i], bw)/diff(pnorm(lims, y[i], bw))
  #   }
  # } else {
    fr <- tabulate(pmax(1, pmin(N, 1+floor((y - lower)/dx))), nbins = N)
    for(i in 1:N) {
      w <- dnorm(x, x[i], bw)
      yo[i] <- sum(fr*w)/sum(w)
    }
  # }
  structure(list(x = x, y = yo/(sum(yo)*dx), bw = bw, n = n,
                 lower = lower, upper = upper,
                 data_name = data_name),
            class = "kde")
}

#' @rdname kde
#' @export
plot.kde <- function(x, ..., col = "steel blue",
                     xlab = bquote(x == italic(.(x$data_name))),
                     ylab = expression(kde(italic(x)))) {
  with(x, {
    plot(x = c(lower, x, upper), y = c(0, y, 0),
         type = "l", xlab = xlab, ylab = ylab,
         col = col, panel.first = grid(), ...)
  })
  invisible(x)
}
# .S3method("plot", "kde")

#' @rdname kde
#' @export
print.kde <- function(x, ...) {
  cat("A Kernel Density Estimate\n\n")
  labs <- c("Response", "Sample size", "Range", "Bandwidth", "Resolution")
  values <- with(x, c(data_name, n,
                      paste0(signif(lower, 4), " < x < ",
                             signif(upper, 4)),
                      signif(bw, 4),
                      length(x)))
  cat(paste0(format(labs), ": ", values), sep = "\n")
  generics <- attr(methods(class = "kde"), "info")$generic
  cat("\nMethods exist for generics:",
      paste(generics, collapse = ", "), "\n")
  invisible(x)
}
# .S3method("print", "kde")

#' @rdname kde
#' @export
as.data.frame.kde <- function(x, row.names, optional, ...,
                              responseName = "density") {
  with(unclass(x),
       setNames(data.frame(X = c(lower, x, upper),
                           Y = c(0, y, 0), ...),
                c(data_name, responseName)))
}

#' @importFrom MASS kde2d bandwidth.nrd
NULL

#' Two Dimensional Kernel Density Estimate
#'
#' A pure R implementation of two direct algorithms.
#'
#' @param x Numeric vector or, in generics, a kde object
#' @param y Numeric vector of the response Variables
#' @param N Integer: number of equally spaced points for the
#'        kde.  May be one or two for x and y respectively
#' @param bw Numeric: bandwidth(s) for x and y kernels
#'        respectively.  A single value will be duplicated.
#' @param lower_x,upper_x,lower_y,upper_y Numeric: kde range limits
#' @param k Ineger: number of bandwidths by which to increase
#'        the range for the default limits.
#' @param xlab,ylab,col base graphics parameters
#' @param row.names,optional As for the as.data.frame generic
#' @param borders logical: should zero borders be included around
#'        the density when the object is coerced to a data frame.
#' @param responseName character string: name for the density vector.
#' @param ...  Extra arguments passed on to methods.
#'
#' @return An object of class "kde_2d"
#' @export
#'
#' @examples
#' set.seed(1234)
#' x_value <- abs(rnorm(5000))
#' y_value <- rnorm(x_value)
#' kxy <- kde_2d(x_value, y_value, lower_x = 0)
#' kxy
#' if("package:ggplot2" %in% search()) {
#'     ggplot(within(as.data.frame(kxy),
#'                   pDens <- 2*dnorm(x_value)*dnorm(y_value))) +
#'       aes(x = x_value, y = y_value) +
#'       geom_raster(aes(fill = density)) +
#'       geom_contour(aes(z = pDens), colour = "black") +
#'       scale_fill_viridis_c(direction = -1) +
#'       theme_bw() + theme(legend.position = "none")
#' } else {
#'   with(unclass(kxy), {
#'     oldPar <- par(mar = c(3, 3, 1, 1) + 0.1, las = 1, mgp = c(1.5, 0.5, 0),
#'                  cex.axis = 0.7, cex.lab = 0.8, tck = -0.015)
#'     plot(kxy)
#'     contour(x, y, z = outer(x, y, function(x, y) 2*dnorm(x)*dnorm(y)),
#'             nlevels = 6, add = TRUE)
#'     box()
#'     par(oldPar)
#'   })
#' }
#'
kde_2d <- function(x, y, N = 256,
                   bw = c(x = bw.nrd0(x), y = bw.nrd0(y)),
                   lower_x = limits$x$lower, upper_x = limits$x$upper,
                   lower_y = limits$y$lower, upper_y = limits$y$upper,
                   k = 1) {
  use_mass <- !any(missing(lower_x), missing(upper_x),
                   missing(lower_y), missing(upper_y))
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
  # z <- matrix(0, N[["x"]], N[["y"]])
  # for(i in seq_along(x)) {
  #   z <- z + outer(dnorm(xo, x[i], bw[["x"]]),
  #                  dnorm(yo, y[i], bw[["y"]]))
  # }
  # z <- z/sum(z*dx*dy)
  z <- if(use_mass) {
    MASS::kde2d(x, y, h = bw*212/45, n = N, lims = c(range(xo), range(yo)))
  } else {
    kde2dCpp(x, y, xo, yo,
             c(lower_x, upper_x),
             c(lower_y, upper_y), bw, dx*dy)
  }
  structure(c(z, list(bw = bw, n = length(x),
                      lower = c(x = lower_x, y = lower_y),
                      upper = c(x = upper_x, y = upper_y),
                      data_name = c(x = data_name_x, y = data_name_y))),
              class = "kde_2d")
}

#' @rdname kde_2d
#' @export
print.kde_2d <- function(x, ...) {
  cat("A Two-dimensional Kernel Density Estimate\n\n")
  labs <- c("Responses", "Sample size", "Ranges", "Bandwidths", "Resolution")
  values <- with(x, c(paste(data_name, collapse = ", "),
                      n,
                      paste(paste0(signif(lower["x"], 4), " < x < ",
                                   signif(upper["x"], 4)),
                            paste0(signif(lower["y"], 4), " < y < ",
                                   signif(upper["y"], 4)),
                            sep = ", "),
                      paste(signif(bw, 4), collapse = ", "),
                      paste(length(x), length(y), sep = ", ")))
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
                        col = hcl.colors(25, rev = TRUE)) {
  with(x, {
    image(x, y, z, xlab = xlab, ylab = ylab, col = col, ...)
  })
  invisible(x)
}
# .S3method("plot", "kde_2d")


#' @rdname kde_2d
#' @export
as.data.frame.kde_2d <- function(x, row.names, optional, ...,
                                 borders = FALSE,
                                 responseName = "density") {
  with(unclass(x),
       setNames(if(borders) {
         cbind(expand.grid(x = c(lower[["x"]], x, upper[["x"]]),
                           y = c(lower[["y"]], y, upper[["y"]])),
               z = as.vector(cbind(0, rbind(0, z, 0), 0)))
       } else {
         cbind(expand.grid(x = x, y = y), z = as.vector(z))
       }, c(data_name, responseName)))
}

