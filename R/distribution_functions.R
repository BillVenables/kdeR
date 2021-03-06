#' Kernel Distribution Functions
#'
#' A set of functions to generate p-, d-, q- or r-prefix functions,
#' as are available for standard probability distributions, but for
#' distributions defined by a kernel density estimate.
#'
#' @param x Either a kde object or a numeric vector from which a kde
#'          object can be generated.
#' @param ...  Additional arguments forwarded to the \code{kde} function.
#'
#' @return  An R function for probabilities, densities, quantiles or
#'          simulated samples, as indicated by the prefix.
#' @export
#'
#' @examples
#' set.seed(1234)
#' x <- rnorm(10000)
#' kx <- kde(x)
#' pkx <- p_fun(kx)
#' dkx <- d_fun(kx)
#' qkx <- q_fun(kx)
#' rkx <- r_fun(kx)
#' rbind(quantile(x, 0:10/10), qkx(0:10/10))
#' pkx(qkx(1:9/10))
#'     par(mar = c(3, 3, 1, 1) + 0.1, las = 1, mgp = c(1.5, 0.5, 0),
#'         cex.axis = 0.7, cex.lab = 0.8, tck = -0.015)
#' curve(dkx, xlim = qkx(0:1), col = 2)
#' curve(dnorm, add = TRUE)
#' grid()
#' rx <- rkx(x)
#' qqplot(x, rx, pch = ".")
#' abline(0, 1, col=2)
#' grid()
#' rm(x, rx)
p_fun <- function(x, ...) {
  UseMethod("p_fun")
}

#' @rdname p_fun
#' @export
p_fun.kde <- function(x, ...) {
  with(x, {
    dx <- diff(c(lower, x, upper))
    dy <- (c(0, y) + c(y, 0))/2
    pr <- c(0, cumsum(dx*dy), 1)
    xo <- c(lower, (c(lower, x) + c(x, upper))/2, upper)
    approxfun(x = xo, y = pr, rule = 2, yleft = 0, yright = 1)
  })
}


#' @rdname p_fun
#' @export
p_fun.numeric <- function(x, ...) {
  p_fun.kde(kde(x, ...))
}


#' @rdname p_fun
#' @export
d_fun <- function(x, ...) {
  UseMethod("d_fun")
}


#' @rdname p_fun
#' @export
d_fun.kde <- function(x, ...) {
  with(x, {
    dy <- (c(0, y) + c(y, 0))/2
    pr <- c(0, dy, 0)
    xo <- c(lower, (c(lower, x) + c(x, upper))/2, upper)
    approxfun(x = xo, y = pr, rule = 2, yleft = 0, yright = 0)
  })
}


#' @rdname p_fun
#' @export
d_fun.numeric <- function(x, ...) {
  d_fun.kde(kde(x, ...))
}


#' @rdname p_fun
#' @export
q_fun <- function(x, ...) {
  UseMethod("q_fun")
}


#' @rdname p_fun
#' @export
q_fun.kde <- function(x, ...) {
  with(x, {
    dx <- diff(c(lower, x, upper))
    dy <- (c(0, y) + c(y, 0))/2
    pr <- c(0, cumsum(dx*dy), 1)
    xo <- c(lower, (c(lower, x) + c(x, upper))/2, upper)
    approxfun(y = xo, x = pr, rule = 2, yleft = lower, yright = upper)
  })
}


#' @rdname p_fun
#' @export
q_fun.numeric <- function(x, ...) {
  q_fun.kde(kde(x, ...))
}


#' @rdname p_fun
#' @export
r_fun <- function(x, ...) {
  UseMethod("r_fun")
}


#' @rdname p_fun
#' @export
r_fun.kde <- function(x, ...) {
  .qf <- q_fun(x)
  function(n) {
    .qf(stats::runif(n))
  }
}


#' @rdname p_fun
#' @export
r_fun.numeric <- function(x, ...) {
  r_fun.kde(kde(x, ...))
}

#' Coercion to class kde
#'
#' A coercion function to class \code{"kde"}, currently the
#' only method is for objects of class \code{"density"} from
#' the \code{stats} package.
#'
#' @param x an object
#' @param ... extra arguments for methods.  Currently not used.
#'
#' @return an object of class \code{"kde"}
#' @export
#'
#' @examples
#' set.seed(1234)
#' .x <- rnorm(5000)
#' dx <- density(.x)
#' oldPar <- par(mfrow = c(1,2), pty = "s", cex.main = 0.8,
#'   mar = c(3, 3, 1, 1) + 0.1, las = 1, mgp = c(1.5, 0.5, 0),
#'   cex.axis = 0.7, cex.lab = 0.8, tck = -0.015)
#' plot(dx)
#' plot(as_kde(dx))
#' par(oldPar)
#' rm(.x, oldPar)
as_kde <- function(x, ...) {
  UseMethod("as_kde")
}

#' @rdname as_kde
#' @export
as_kde.density <- function(x, ...) {
  with(unclass(x), {
    dx <- mean(diff(x))/2
    structure(list(x = x, y = y, bw = bw, n = n,
                   lower = min(x) - dx, upper = max(x) + dx,
                   data_name = data.name),
              class = "kde")
  })
}

#' @rdname p_fun
#' @export
p_fun.density <- function(x, ...) {
  p_fun.kde(as_kde(x))
}

#' @rdname p_fun
#' @export
d_fun.density <- function(x, ...) {
  d_fun.kde(as_kde(x))
}

#' @rdname p_fun
#' @export
q_fun.density <- function(x, ...) {
  q_fun.kde(as_kde(x))
}

#' @rdname p_fun
#' @export
r_fun.density <- function(x, ...) {
  r_fun.kde(as_kde(x))
}
