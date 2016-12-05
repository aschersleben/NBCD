#' Variance calculation (updatable) by Youngs and Cramer
#'
#' @param x \code{[numeric]}\cr
#'   vector (finite, no missing values)
#'
#' @param object \code{[ycVar | NULL]}\cr
#'   list (object from updateVarYC) or NULL
#'
#' @return \code{[ycVar]}\cr
#'   \itemize{
#'     \item \code{$var}: (updated) variance via Youngs-Cramer-algorithm
#'     \item \code{$n}: number of included elements
#'   }
#'
#' @section Alternatives:
#' \itemize{
#'   \item Cotton, E.W. (1975), Remark on "Stably Updating Mean and Standard
#'    Deviation of Data," Communications of the Association for Computing
#'    Machinery, 18, 458.
#'   \item Hanson, R.J. (1975), "Stably Updating Mean and Standard Deviation of
#'    Data," Communications of the Association for Computing Machinery, 18,
#'    57-58.
#'   \item Welford, B.P. (1962), "Note on a Method for Calculating Corrected Sums
#'    of Squares and Products," Technometrics, 4, 419-420.
#'   \item West, D.H.D. (1979), "Updating Mean and Variance Estimates: An Improved
#'    Method," Communications of the Association for Computing Machinery,
#'    22, 532-535.
#' }
#'
#' @references
#' Youngs, Edward A., and Elliot M. Cramer. "Some results relevant to choice of
#' sum and sum-of-product algorithms." Technometrics 13.3 (1971): 657--665,
#' \href{http://dx.doi.org/10.2307/1267176}{DOI:10.2307/1267176}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' yc.fun = function(x) {
#'   yc = NULL
#'   for (i in seq_along(x)) {
#'     yc = NBCD:::updateVarYC(x[i], yc)
#'   }
#'   return(yc$var)
#' }
#' yc.fun(1:10)
#' }
#'

updateVarYC = function(x, object = NULL) {

  # no argument checks for speed reasons
  # assertNumeric(x, finite = TRUE, any.missing = FALSE)
  # assert(checkNull(object), checkClass(object, "ycVar"))

  if (length(x) == 1) {
    if (is.null(object)) {
      j = n = 1
      T_1j = x
      S_1n = 0
    } else {
      j = n = object$n + 1
      T_1j = object$ycT + x
      S_1n = object$ycS + (j * x - T_1j)^2 / (j * (j - 1))
    }
  } else {
    if (is.null(object)) {
      n = length(x)
      j = 2:n
      T_1j = cumsum(x)[-1]
      S_1n = sum((j * x[-1] - T_1j)^2 / (j * (j - 1)))
    } else {
      i = object$n
      n = length(x) + i
      j = (i + 1):n
      T_1j = object$ycT + cumsum(x)
      S_1n = object$ycS + sum((j * x - T_1j)^2 / (j * (j - 1)))
    }
  }

  return(structure(list(var = if (n > 1) S_1n / (n - 1) else 0,
                        ycS = S_1n, ycT = T_1j[length(T_1j)], n = n),
                   class = "ycVar"))
}



# "Normalize" objects.
#
# @param x \code{[numeric]}\cr
# Object
#
# @return
# Object divided by sum of its values.
#
# @examples
# probnorm(1:5)
#

probnorm <- function(x) {
  x / sum(x)
}



# #' 0-1-"normalize" objects.
# #'
# #' Values > 1 or < 0 will be subsetted to 1 or 0 resp. Then probnorm...
# #'
# #' @param x \code{[any]}\cr
# #' Object
# #'
# #' @return
# #' Object with values in [0, 1] and divided by its sum.
# #'
# #' @examples
# #' contnorm(c(0, 0.2, 0.4, 1.2))
# #'
# #' contnorm(c(-2, 0.3, 0.7, 10))
# #'
#
# contnorm <- function(x) {
#   x[x < 0] <- 0
#   x[x > 1] <- 1
#   probnorm(x)
# }



#' Waiting Time Calculation
#'
#' Calculates the mean waiting time for a special absorbing markov chain
#' (directed acyclic graph) that is based on an n-point-distribution. The value
#' of this function answers the question:\cr\cr
#' "How long do I have to wait (on average) until every
#' value of a discrete distribution (with a finite set of values) appears
#' \emph{at least} once?"\cr\cr
#' Uses \code{Rcpp} with memoization for faster calculations, but the problem
#' is still \eqn{O(2^n)}: 0.2s for \code{length(x) = 15}, and 11s
#' for \code{length(x) = 20}. An approximation approach is given, see argument
#' \code{"max.len"}.
#'
#'
#' @param x \code{[numeric]}\cr
#' Probabilities for which to calculate the mean waiting time.
#'
#' @param max.len \code{[numeric(1)]}\cr
#' Controls calculation time. If x is longer than max.len, the waiting time will
#' be approximated: x then just contains the \eqn{max.len-1} smallest values and
#' the sum of all others. Default is 10.
#'
#' @return \code{[numeric(1)]} Mean waiting time.
#'
#' @export
#'
#' @examples
#' # How long until every side of a dice appears at least once (on average)?
#' getWaitingTime(rep(1/6, 6))
#'
#' # Approximation:
#' \dontrun{
#' set.seed(1909)
#' tmp <- probnorm(runif(18))
#' system.time(print(getWaitingTime(tmp, max.len = 18)))
#' system.time(print(getWaitingTime(tmp, max.len = 10)))
#' system.time(print(getWaitingTime(tmp, max.len = 5)))
#' }
#'

getWaitingTime <- function(x, max.len = 10) {
  x <- as.numeric(x)
  x <- x[!(x < sqrt(.Machine$double.eps))]
  checkmate::assertNumeric(x, lower = 0, finite = TRUE, any.missing = FALSE,
                           min.len = 1)
  if (!isTRUE(all.equal(sum(x), 1))) x <- probnorm(x)
  n <- length(x)
  if (n > max.len) {
    x <- sort(x, partial = max.len)
    x <- c(x[1:(max.len - 1)], sum(x[max.len:n]))
  }
  return(kMcpp(x))
}



#' Create hyperplane data with independant normal distributions
#'
#' @param n \code{[numeric(1)]}\cr
#'
#' @param mean \code{[numeric(1) | numeric(d)]}\cr
#'
#' @param sd \code{[numeric(1) | numeric(d) | matrix(d, d)]}\cr
#'   If matrix, it's a VCV with variances on the diagonal.
#'   Else, interpreted as standard deviations.
#'
#' @param d \code{[numeric(1)]}\cr
#'   Dim.
#'
#' @export
#'
#'

rnhyper <- function(n, mean = 0, sd = 1, d) {

  if (missing(d) && length(mean) < 3) d <- 2 else d <- length(mean)
  checkmate::assert_integerish(d, lower = 1, len = 1)

  checkmate::assert(
    checkmate::test_numeric(mean, finite = TRUE, any.missing = FALSE, len = 1),
    checkmate::test_numeric(mean, finite = TRUE, any.missing = FALSE, len = d)
  )
  if (length(mean) == 1) mean <- rep(mean, d)
  if (is.matrix(sd)) {
    out <- mvtnorm::rmvnorm(n, mean, sd)
  } else {
    checkmate::assert(
      checkmate::test_numeric(sd, lower = 0, finite = TRUE, any.missing = FALSE, len = 1),
      checkmate::test_numeric(sd, lower = 0, finite = TRUE, any.missing = FALSE, len = d)
    )
    if (length(sd) == 1) sd <- rep(sd, d)
    out <- mvtnorm::rmvnorm(n, mean, diag(sd^2))
  }

  return(out)
}



#' Generates data of rotating normal (two-dimensional) distributions.
#'
#' @param angle \code{[numeric(1)]}\cr
#'   Rotation angle of the means in relation to the y axis (clockwise).
#'
#' @param n.obs \code{[numeric(1)]}\cr
#'   Number of observations to draw from the "concept" (distribution).
#'
#' @param sd \code{[numeric(1)]}\cr
#'   Standard deviation of the normal distributions.
#'
#' @param both \code{[logical(1)]}\cr
#'   If TRUE, observations are equally distributed to both classes.
#'   Default is FALSE.
#'
#' @param time \code{[numeric]}\cr
#'   The time stamp to attach to the data.
#'   If missing, it's the "angle" argument.
#'
#' @return \code{[list]}\cr
#'
#' @export
#'
#' @examples
#' set.seed(1909)
#' dat <- getDataRotNorm2(45, 100)
#' plot(dat$x, col = dat$class)
#'

getDataRotNorm2 <- function(angle, n.obs = 1, sd = sqrt(0.2), both = FALSE,
                            time) {
  out <- list(x = as.data.frame(rnhyper(n.obs, sd = diag(rep(sd^2, 2)))))
  colnames(out$x) <- c("x1", "x2")
  if (both) {
    out$class <- sample(rep(1:2, ceiling(n.obs / 2)), n.obs)
  } else {
    out$class <- sample(1:2, n.obs, TRUE)
  }
  out$class <- factor(out$class, 1:2)
  out$x <- out$x + ifelse(out$class == 1, c(1, 1), c(-1, -1)) %*%
    pol2cart(angle, measure = "d", start = "N", direction = "CW")
  out$time <- if (missing(time)) angle else time
  return(out)
}



#' Generates data of "rotating hyperplane" example.
#'
#' @param angle \code{[numeric(1)]}\cr
#'   Rotation angle of the hyperplane in relation to the x axis (clockwise).
#'
#' @param n.obs \code{[numeric(1)]}\cr
#'   Number of observations to draw from the "concept" (distribution).
#'
#' @param time \code{[numeric]}\cr
#'   The time stamp to attach to the data.
#'   If missing, it's the "angle" argument.
#'
#' @return \code{[list]}\cr
#'
#' @export
#'
#' @examples
#' set.seed(1909)
#' dat <- getDataRotHyp2(45, 100)
#' plot(dat$x, col = dat$class)
#'

getDataRotHyp2 <- function(angle, n.obs = 1, time) {
  slopes <- function(x) {
    if (x < 181) {
      circ <- (c(x, x + 180) / 180 * pi) - pi/2
      tmp <- cbind(sin(circ), cos(circ))
      b <- (tmp[2, 2] - tmp[1, 2]) / (tmp[2, 1] - tmp[1, 1])
    } else {
      circ <- (c(360 - x, 540 - x) / 180 * pi) - pi/2
      tmp <- cbind(sin(circ), cos(circ))
      b <- (tmp[2, 2] - tmp[1, 2]) / (tmp[1, 1] - tmp[2, 1])
    }
    return(b)
  }

  out <- list()
  out$x <- as.data.frame(cbind(runif(n.obs, -1, 1), runif(n.obs, -1, 1)))
  colnames(out$x) <- c("x1", "x2")
  out$class <- ifelse(slopes(angle) * out$x$x1 - out$x$x2 < 0,
                      if (angle <= 90 | angle >= 270) 1 else 2,
                      if (angle <= 90 | angle >= 270) 2 else 1)
  out$time <- if (missing(time)) angle else time
  return(out)
}



#' Polar to Cartesian Coordinates
#'
#' @param phi \code{[numeric]}\cr
#'
#' @param radius \code{[numeric]}\cr
#'
#' @param measure \code{[character(1)]}\cr
#'   Measure of "phi" (degrees, radians, or turns)
#'
#' @param start \code{[character(1)]}\cr
#'   Start position, default is east ("E").
#'
#' @param direction \code{[character(1)]}\cr
#'   Anti-clockwise (default) or clockwise.
#'
#'
#' @return \code{[matrix]}\cr
#'   Cartesian coordinates.
#'
#'
#' @export
#'
#' @examples
#' all.equal(pol2cart(pi / 2), pol2cart(90, measure = "d"))
#' all.equal(pol2cart(pi / 2), pol2cart(1/4, measure = "t"))
#'
#' oldpar <- par(pty = "s")
#' plot(pol2cart(1:360, measure = "d"), col = rainbow(360),
#'      xlim = c(-1, 1), ylim = c(-1, 1))
#' par(oldpar)
#'
#' oldpar <- par(pty = "s")
#' plot(pol2cart(seq(0, 2, len = 500), seq(0, 1, len = 500), measure = "t"),
#'      col = rev(grey.colors(500)), pch = 16, xlim = c(-1, 1), ylim = c(-1, 1))
#' par(oldpar)
#'
#' oldpar <- par(pty = "s")
#' plot(NA, xlim = c(-3, 3), ylim = c(-3, 3), xlab = "x", ylab = "y")
#' for (i in 1:360)
#'   points(matrix(rnorm(2, as.numeric(pol2cart(i, 2, m = "d")), 0.2), 1))
#' par(oldpar)
#'
#' oldpar <- par(pty = "s")
#' plot(NA, xlim = c(-3, 3), ylim = c(-3, 3), xlab = "x", ylab = "y")
#' points(rnhyper(100, pol2cart(45, 2, m = "d"), sd = 0.5), col = 2)
#' points(rnhyper(100, pol2cart(45 + 225, 2, m = "d"), sd = 1), col = 3)
#' points(rnhyper(100, pol2cart(45 + 90, 2, m = "d"), sd = 0.5), col = 4)
#' points(pol2cart(45 + c(0, 90, 225), 2, m = "d"), pch = 4, lwd = 3)
#' lines(pol2cart(seq(0, 1, len = 400), 2, m = "t"))
#' par(oldpar)
#'

pol2cart <- function(phi, radius = 1L, measure = c("rad", "deg", "turn"),
                     start = c("E", "N", "W", "S"),
                     direction = c("ACW", "CW")) {

  measure <- match.arg(measure)
  start <- match.arg(start)
  direction <- match.arg(direction)
  checkmate::assertNumeric(phi, finite = TRUE)
  checkmate::assert(
    checkmate::checkNumber(radius, lower = 0, finite = TRUE),
    checkmate::checkNumeric(radius, lower = 0, finite = TRUE,
                            any.missing = FALSE, len = length(phi))
  )

  phi <- switch(measure, rad = phi %% (2L * pi),
                deg = phi %% 360L * pi / 180L,
                turn = phi %% 1L * 2L * pi)

  if (direction == "CW")
    phi <- 2L * pi - phi

  if (start != "E")
    phi <- switch(start, N = phi + pi / 2,
                  W = phi + pi,
                  S = phi + 3 / 2 * pi)

  return(cbind(x = radius * cos(phi), y = radius * sin(phi)))

}



#' Paste head of object to a string with a maximal length.
#'
#' @param x \code{[numeric | character]}\cr
#' Object to paste.
#'
#' @param len \code{[numeric]}\cr
#' Maximal length.
#'
#' @param coll \code{[character]}\cr
#' Separation string.
#'
#' @param end \code{[character]}\cr
#' End string.
#'
#' @return \code{[character(1)]}\cr
#'

pastehead <- function(x, len = 3, coll = " / ", end = "...") {
  if (is.null(x) | length(x) == 0) return("none") else
    paste(c(head(x, len - 1), if (length(x) > len)
      end else if (length(x) == len)
        x[len]), collapse = coll)
}

