kde2d_mass <- function (x, y, h, n = 256, lims = c(range(x), range(y)), B = 5e6) {
  nx <- length(x)
  if (length(y) != nx)
    stop("data vectors must be the same length")
  if (any(!is.finite(x)) || any(!is.finite(y)))
    stop("missing or infinite values in the data are not allowed")
  if (any(!is.finite(lims)))
    stop("only finite values are allowed in 'lims'")
  n <- rep(n, length.out = 2L)
  xo <- seq.int(lims[1L], lims[2L], length.out = n[1L])
  yo <- seq.int(lims[3L], lims[4L], length.out = n[2L])
  h <- if (missing(h))
    c(bandwidth.nrd(x), bandwidth.nrd(y))/4
  # c(bw.nrd0(x), bw.nrd0(y))
  else rep(h, length.out = 2L)
  if (any(h <= 0))
    stop("bandwidths must be strictly positive")
  # h <- h/4
  # ax <- outer(xo, x, "-")/h[1L]
  # ay <- outer(yo, y, "-")/h[2L]
  # z <- tcrossprod(matrix(dnorm(ax), , nx),
  #                 matrix(dnorm(ay), , nx))/(nx * h[1L] * h[2L])
  z <- matrix(0, n[1], n[2])
  kinc <- floor(B/(n[1] * n[2]))
  k2 <- 0
  while(k2 < nx) {
    k1 <- k2 + 1
    k2 <- min(nx, k2 + kinc)
    x_kernels <- outer(xo, x[k1:k2], function(xo, x) dnorm((xo - x)/h[1]))
    y_kernels <- outer(yo, y[k1:k2], function(yo, y) dnorm((yo - y)/h[2]))
    z <- z + tcrossprod(x_kernels, y_kernels)
  }
  # list(x = xo, y = yo, z = z/(nx*h[1]*h[2]))
  list(x = xo, y = yo, z = z/(sum(z)*mean(diff(x0))*mean(diff(yo))))
}
