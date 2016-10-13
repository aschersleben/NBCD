#' nb2: Extended Naive Bayes - Updateable Version
#'
#' Modified and extended copy of \link[e1071]{naiveBayes} from \code{e1071}.
#'
#' @param x \code{[matrix | data.frame]}\cr
#' A numeric matrix, or a data frame of categorical and/or numeric variables.
#'
#' @param y \code{[numeric | factor | character]}\cr
#' Class vector.
#'
#' @param laplace \code{[numeric(1)]}\cr
#' positive double controlling Laplace smoothing.
#' The default (0) disables Laplace smoothing.
#'
#' @param discretize \code{[character(1)]}\cr
#' If missing: No discretization will be performed and discParams argument is
#'   ignored.
#' If "fixed": number of intervals or there boundaries have to be given in the
#'   discParam argument
#'
#' @param discParams \code{[numeric | list]}\cr
#' If discretize = "fixed": Has to be a named list giving the boundaries for
#'   the variables to be discretized.
#'
#' @param ...
#' Currently ignored.
#'
#' @return \code{[nb2]}\cr
#' additional element: "yc" contains Youngs-Cramer variance parameters
#'
#' @author
#' \itemize{
#'   \item David Meyer \email{David.Meyer@@R-project.org}: original \link[e1071]{naiveBayes}
#'   \item Philipp Aschersleben \email{aschersleben@@statistik.tu-dortmund.de}:
#'   modifications and extensions for the NBCD package
#' }
#'
#' @export
#'
#' @seealso \code{\link{update.nb2}}, \code{\link{plot.nb2}},
#' \code{\link{print.nb2}}, \code{\link{predict.nb2}}
#'
#' @examples
#' mod <- nb2(iris[, 1:4], iris[, 5])
#' mod
#'
#' # Example for easy discretization:
#' discParam <- list(Sepal.L = c(5, 6, 7), Sepal.W = c(2.5, 3.5))
#' mod2 <- nb2(iris[, 1:4], iris[, 5], discretize = "fixed", discParams = discParam)
#' mod2
#'

nb2 <- function(x, y, laplace = 0, discretize, discParams, ...) {
  call <- match.call()
  Yname <- deparse(substitute(y))
  x <- as.data.frame(x)
  y <- factor(y)
  if (nlevels(y) == 1) y <- factor(y, c(levels(y), "ZZZ_TEMP"))

  if (!missing(discretize)) {
    discretize <- match.arg(discretize, c("fixed"))
    if (discretize == "fixed") {
      discNames <- names(discParams) <- match.arg(names(discParams), names(x), several.ok = TRUE)
      out <- sapply(discNames, simplify = FALSE, function(dN) {
        breaks <- unique(c(Inf, -Inf, discParams[[dN]]))
        return(structure(cut(x[, dN], breaks = breaks), breaks = breaks))
      })
      x[discNames] <- out
      params <- list(breaks = lapply(x[discNames], attr, "breaks"))
    # } else if (discretize == "PiD") {
    #
    }
    disc <- list(discretize = discretize, params = params)
  } else {
    disc <- NULL
  }

  est <- function(var) {
    if (is.numeric(var)) {
      sd.yc <- tapply(var, y, updateVarYC)
      out <- cbind(tapply(var, y, mean, na.rm = TRUE),
                   sqrt(sapply(sd.yc, function(z)
                     if (is.null(z)) NA else z[["var"]])))
    } else {
      tab <- table(y, var)
      out <- (tab + laplace)/(rowSums(tab) + laplace * nlevels(var))
      sd.yc <- list()
    }
    return(list(out, sd.yc))
  }
  apriori <- table(y)
  apriori.list <- sapply(names(x), function(y) apriori, simplify = FALSE)
  tablest <- lapply(x, est)
  tables <- lapply(tablest, "[[", 1)
  yc <- lapply(tablest, "[[", 2)

  for (i in 1:length(tables))
    names(dimnames(tables[[i]])) <- c(Yname, colnames(x)[i])
  names(dimnames(apriori)) <- Yname
  apriori.list <- sapply(apriori.list, function(z) { names(dimnames(z)) <- Yname; z },
                         simplify = FALSE)
  structure(list(apriori = apriori, apriori.list = apriori.list,
                 tables = tables, levels = levels(y),
                 call = call, yc = yc, disc = disc),
            class = c("nb2", "naiveBayes"))
}



#' Update Naive Bayes Classifier with New Data
#'
#' @param object \code{[nb2]}\cr
#' Naive Bayes model to be updated.
#'
#' @param newdata \code{[data.frame]}\cr
#' Data.frame with same column names as in object.
#'
#' @param y \code{[factor | numeric | character]}\cr
#' Class vector.
#'
#' @param ...
#' Arguments passed to \code{\link{nb2}()}.
#'
#' @export
#'
#' @seealso \code{\link{nb2}}, \code{\link{plot.nb2}},
#' \code{\link{print.nb2}}, \code{\link{predict.nb2}}
#'
#' @examples
#' data = data.frame(x1 = 1:10, x2 = 1:10)
#' class = factor(rep(1:2, each = 5))
#' (object = nb2(data, class))
#' newdata = data.frame(x1 = 11:14, x2 = 11:14)
#' y = factor(rep(1:2, each = 2))
#' update(object, newdata, y)
#'
#' data = data.frame(x1 = 1:10, x2 = 1:10)
#' class = factor(rep(1, 10))
#' (object = nb2(data, class))
#' newdata = data.frame(x1 = 11:14, x2 = 11:14)
#' y = factor(rep(1:2, each = 2))
#' update(object, newdata, y)
#'
#' \dontrun{
#' set.seed(1909)
#' dat = mlbench::mlbench.smiley(500)
#' obj = nb2(x = dat$x[1:100, ], y = dat$classes[1:100])
#' plot(obj, xlim = c(-2, 2), ylim = c(-2, 2), gg = FALSE)
#' i = j = 101
#' while (i < nrow(dat$x)) {
#'   j = min(j + max(rpois(1, 10), 2), nrow(dat$x))
#'   newd = dat$x[i:j, ]
#'   newc = dat$classes[i:j]
#'   obj = update(obj, newd, newc)
#'   cat(i, j)
#'   print(all.equal(obj, nb2(dat$x[1:j, ], dat$classes[1:j]),
#'                   check.attributes = FALSE, check.names = FALSE))
#'   plot(obj, xlim = c(-2, 2), ylim = c(-2, 2), gridsize = 50, gg = FALSE)
#'   points(dat$x[1:j, ], col = dat$classes[1:j], pch = 4, lwd = 4)
#'   if (requireNamespace("BBmisc", quietly = TRUE)) BBmisc::pause()
#'   i = j + 1
#' }
#' all.equal(obj, nb2(dat$x, dat$classes), check.attributes = FALSE, check.names = FALSE)
#' }
#'
#'
#' ### mixed attributes:
#' set.seed(1909)
#' daten <- mlbench::mlbench.2dnormals(1000)
#' daten$x <- data.frame(daten$x, X3 = factor(as.numeric(daten$classes) + rbinom(1000, 1, 0.2)))
#' plot(daten)
#'
#' mod1 <- nb2(daten$x[, 1:2], daten$classes)
#' preds1 <- predict(mod1, daten$x[, 1:2])
#' mean(preds1 != daten$classes)
#'
#' mod2 <- nb2(daten$x, daten$classes)
#' preds2 <- predict(mod2, daten$x)
#' mean(preds2 != daten$classes)
#'
#' k <- sample(1:1000, 100)
#' newdata <- daten$x[k, ]
#' y <- daten$classes[k]
#'
#' all.equal(update(mod2, newdata, y), nb2(rbind(daten$x, newdata), c(daten$classes, y)),
#'           check.attributes = FALSE, check.names = FALSE)
#'
#'
#' iris.disc <- iris[, 1:4]
#' discParam <- list(Sepal.Length = c(5, 6, 7), Sepal.Width = c(2.5, 3.5))
#' iris.disc[, 1:2] <- sapply(names(discParam), simplify = FALSE, function(x)
#'   cut(iris[, x], breaks = unique(c(Inf, -Inf, discParam[[x]]))))
#' mod.disc1 <- nb2(iris.disc[, 1:4], iris[, 5])
#'
#' mod.disc2 <- nb2(iris[, 1:4], iris[, 5], discretize = "fixed", discParam = discParam)
#'
#' all.equal(mod.disc1, mod.disc2)
#'
#' plot(mod.disc1)
#' plot(mod.disc2, xlim = c(4, 8), ylim = c(1.5, 4.5))
#' plot(mod.disc1, features = c(1, 3))
#' plot(mod.disc2, xlim = c(4, 8), features = c(1, 3))
#'
#'
#' k <- sample(nrow(iris), 10)
#' newdata1 <- iris.disc[k, 1:4]
#' newdata2 <- iris[k, 1:4]
#' y <- iris[k, 5]
#'
#' upd1 <- update(mod.disc1, newdata1, y)
#' upd2 <- update(mod.disc2, newdata2, y)
#' all.equal(upd1, upd2)
#'

update.nb2 <- function(object, newdata, y, ...) {

  laplace <- 0

  if (missing(object) || is.null(object)) {
    return(nb2(newdata, y, ...))
  }

  newdata <- as.data.frame(newdata)
  attribs <- match(names(object$tables), names(newdata))
  # isnumeric <- sapply(newdata, is.numeric)
  # newdatamat <- data.matrix(newdata)

  newobj <- object

  ## new levels
  ylevs <- if (is.factor(y)) levels(y) else {
    yl <- unique(y)
    ind <- sort.list(yl)
    yl <- as.character(yl)
    unique(yl[ind])
  }
  levelstest <- ylevs %in% object$levels

  if (!all(levelstest)) {
    newlevels <- ylevs[!levelstest]
    zzztest <- ("ZZZ_TEMP" == object$levels)
    if (any(zzztest)) {
      alllevels <- c(object$levels[-which(zzztest)], newlevels)
      newobj$apriori <- newobj$apriori[-which(zzztest)]
      newobj$apriori.list <- sapply(object$apriori.list, function(x) x[-which(zzztest)],
                                    simplify = FALSE)
      newobj$tables <- lapply(newobj$tables, function(x) {
        x[-which(zzztest), , drop = FALSE]
      })
      # newobj$levels <- newobj$levels[-which(zzztest)]
    } else {
      alllevels <- c(object$levels, newlevels)
    }
    neworder <- order(alllevels)
    alllevels <- alllevels[neworder]
    newobj$levels <- alllevels
    adds <- rep(0, length(newlevels))
    names(adds) <- newlevels
    temp <- as.table(c(newobj$apriori, adds))
    newobj$apriori <- temp[neworder]
    newobj$apriori.list <- sapply(newobj$apriori.list, function(x) {
      temp <- as.table(c(x, adds))
      temp[neworder]
    }, simplify = FALSE)
    newobj$tables <- lapply(newobj$tables, function(x) {
      tmpmat <- matrix(0, nrow = length(newlevels), ncol = ncol(x),
                      dimnames = list(newlevels))
      tmptab <- rbind(x, tmpmat)
      tmptab[neworder, ]
    })
    newobj$yc <- lapply(newobj$yc, function(x) {
      tmplist <- rep(list(NULL), length(newlevels))
      c(x, tmplist)[neworder]
    })
    y <- factor(y, levels = alllevels)
  } else {
    y <- factor(y, levels = object$levels)
  }


  # update apriori
  ytab <- table(y)
  newobj$apriori <- newobj$apriori + ytab
  names(dimnames(newobj$apriori)) <- names(dimnames(object$apriori))
  oldapriorilist <- newobj$apriori.list
  newobj$apriori.list <- newapriorilist <- sapply(newobj$apriori.list,
    function(x) x + ytab, simplify = FALSE)


  # discretization if necessary
  if (!is.null(object$disc)) {
    if (object$disc$discretize == "fixed") {
      discNames <- names(object$disc$params$breaks)
      out <- sapply(discNames, simplify = FALSE, function(dN) {
        return(cut(newdata[, dN], breaks = object$disc$params$breaks[[dN]]))
      })
      newdata[discNames] <- out
      # params <- list(breaks = lapply(x[discNames], attr, "breaks"))
    }
  }


  # update tables
  newdatasplit <- split(newdata, y)

  res <- lapply(seq_len(ncol(newdata)), function(i) {
    out <- list(tables = matrix(0, nrow = length(newobj$levels), ncol = 2),
                yc = list())

    ndi <- newdata[, i]
    if (is.numeric(ndi)) {
      ### update means
      newmeans <- unlist(sapply(newdatasplit, function(x) mean(x[, i]), simplify = FALSE))
      newmeans[is.nan(newmeans)] <- 0
      oldmeans <- newobj$tables[[i]][, 1]
      oldmeans[is.na(oldmeans)] <- 0
      updmean <- oldmeans * (oldapriorilist[[i]] / newapriorilist[[i]]) +
        newmeans * (ytab / newapriorilist[[i]])
      updmean[is.nan(updmean)] <- NA
      out$tables[, 1] <- updmean

      ### update standard deviations
      updvar <- lapply(seq_len(length(newdatasplit)), function(j) {
        x = newdatasplit[[j]]
        if (nrow(x) == 0) newobj$yc[[i]][[j]] else
          updateVarYC(x[, i], newobj$yc[[i]][[j]])
      })
      out$tables[, 2] <- sqrt(sapply(updvar, function(x)
        if (is.null(x)) NA else x[["var"]]))
      out$yc <- updvar

    } else {
      oldtab <- newobj$tables[[i]]
      oldtab.na <- is.na(oldtab)
      if (any(oldtab.na))
        oldtab[oldtab.na] <- 0
      updtab <- oldtab * as.numeric(oldapriorilist[[i]] / newapriorilist[[i]]) +
        table(y, ndi) / as.numeric(newapriorilist[[i]])
      updtab.na <- is.na(updtab)
      if (any(updtab.na))
        updtab[updtab.na] <- 0
      out$tables <- (updtab + laplace)/(rowSums(updtab) + laplace * nlevels(ndi))
      outtab.na <- is.na(out$tables)
      if (any(outtab.na))
        out$tables[outtab.na] <- 0
    }

    return(out)
  })

  for (i in 1:length(newobj$tables)) {
    # ncells <- prod(dim(newobj$tables[[i]])) # 1:(ncells * length(newobj$levels))
    newobj$tables[[i]][] <- res[[i]]$tables
    names(dimnames(newobj$tables[[i]])) <- names(dimnames(object$tables[[i]]))
    if (!is.null(newobj$apriori.list)) {
      names(dimnames(newobj$apriori.list[[i]])) <- names(dimnames(newobj$apriori))
    }
    newobj$yc[[i]] <- res[[i]]$yc
  }

  return(structure(newobj, class = c("nb2", "naiveBayes")))

}



#' Prediction for Extended Naive Bayes (nb2)
#'
#' Modified and extended copy of \link[e1071]{predict.naiveBayes} from
#' \code{e1071}.
#'
#' @param object
#' An object of class "nb2".
#'
#' @param newdata
#' A dataframe with new predictors (with possibly fewer columns than
#' the training data). Note that the column names of newdata are matched
#' against the training data ones.
#'
#' @param type
#' If "raw", the conditional a-posterior probabilities for each class are
#' returned, and the class with maximal probability else.
#'
#' @param threshold
#' Value replacing cells with probabilities within eps range.
#'
#' @param eps
#' double for specifying an epsilon-range to apply laplace smoothing
#' (to replace zero or close-zero probabilities by theshold.)
#'
#' @param ...
#' Currently ignored.
#'
#' @author
#' \itemize{
#'   \item David Meyer \email{David.Meyer@@R-project.org}: original \link[e1071]{predict.naiveBayes}
#'   \item Philipp Aschersleben \email{aschersleben@@statistik.tu-dortmund.de}:
#'   modifications and extensions for the NBCD package
#' }
#'
#' @export
#'
#' @seealso \code{\link{nb2}}, \code{\link{update.nb2}}, \code{\link{plot.nb2}},
#' \code{\link{print.nb2}}
#'
#' @examples
#' mod <- nb2(iris[, 1:4], iris[, 5])
#' predict(mod, iris)
#'
#' discParam <- list(Sepal.L = c(5, 6, 7), Sepal.W = c(2.5, 3.5))
#' mod2 <- nb2(iris[, 1:4], iris[, 5], discretize = "fixed", discParam = discParam)
#' predict(mod2, iris)
#'
#' c(mean(predict(mod, iris) == iris$Species),
#'   mean(predict(mod2, iris) == iris$Species))
#'

predict.nb2 <- function(object, newdata, type = c("class", "raw"),
                        threshold = 0.001, eps = 0, ...) {
  type <- match.arg(type)
  newdata <- as.data.frame(newdata)
  attribs <- match(names(object$tables), names(newdata))

  # discretization if necessary
  disc <- object$disc
  if (!is.null(disc)) {
    if (disc$discretize == "fixed") {
      discNames <- names(disc$params$breaks)
      nd.disc <- discNames %in% names(newdata)
      if (any(nd.disc)) {
        out <- sapply(discNames[nd.disc], simplify = FALSE, function(dN) {
          return(cut(newdata[, dN], breaks = disc$params$breaks[[dN]]))
        })
        newdata[discNames[nd.disc]] <- out
        # params <- list(breaks = lapply(x[discNames], attr, "breaks"))
      }
    }
  }

  isnumeric <- sapply(newdata, is.numeric)
  newdata <- data.matrix(newdata)

  L <- sapply(1:nrow(newdata), function(i) {
    ndata <- newdata[i, ]
    apl <- is.null(object$apriori.list)
    if (apl) {
      la <- log(probnorm(object$apriori))
    } else la <- 0
    L <- la + apply(log(sapply(seq_along(attribs),
           function(v) {
             if (apl) {
               apr <- 1
             } else apr <- probnorm(object$apriori.list[[v]])^(1 / length(attribs))
             nd <- ndata[attribs[v]]
             if (is.na(nd)) rep(1, length(object$apriori)) else {
               prob <- if (isnumeric[attribs[v]]) {
                 msd <- object$tables[[v]]
                 msd[, 2][msd[, 2] <= eps] <- threshold
                 dnorm(nd, msd[, 1], msd[, 2])
               } else object$tables[[v]][, nd]
               prob[prob <= eps] <- threshold
               apr * prob
             }
           })), 1, sum)
    if (type == "class")
      L
    else {
      sapply(L, function(lp) {
        1/sum(exp(L - lp))
      })
    }
  })
  if (type == "class")
    factor(object$levels[apply(L, 2, which.max)], levels = object$levels)
  else t(L)
}



#' Plot Naive Bayes Models
#'
#' @param x \code{[nb2]}\cr
#' A Naive Bayes Model
#'
#' @param xlim,ylim \code{[numeric(2)]}\cr
#'
#' @param features \code{[numeric(2) | character(2)]}\cr
#' If missing, the first two features are plotted.
#'
#' @param sdevs \code{[numeric(1)]}\cr
#' If missing(xlim) & missing(ylim): How many standard deviations from the mean
#' in x and y direction.
#'
#' @param gridsize \code{[numeric(1)]}\cr
#' Size of grid to be predicted and plotted.
#'
#' @param data \code{[list | data.frame]}\cr
#' List with elements "x" and "class" or data.frame (then see argument
#' class.name).
#'
#' @param class.name \code{[character]}\cr
#' If data is data.frame, then give the class' column name here.
#' Default: \code{"class"}.
#'
#' @param gg \code{[logical]}\cr
#'
#' @param ...
#' Currently ignored.
#'
#' @export
#'
#' @seealso \code{\link{nb2}}, \code{\link{update.nb2}},
#' \code{\link{print.nb2}}, \code{\link{predict.nb2}}
#'
#' @examples
#' mod <- nb2(iris[, 1:4], iris[, 5])
#' plot(mod, features = 2:3, xlim = c(1.5, 4.5), ylim = c(0.5, 7.5),
#'      gridsize = 75, gg = FALSE)
#' points(iris[, 2:3], col = iris[, 5], pch = 4, lwd = 3)
#'
#' require(ggplot2)
#' mod1 <- nb2(iris[, 1:4], iris[, 5])
#' preds <- predict(mod1, iris[, c(1, 3)])
#' iris2 <- cbind(iris, err = (preds != iris$Species))
#' plot(mod1, features = c(1, 3), xlim = c(4, 8), ylim = c(0, 8)) +
#'   geom_point(data = subset(iris2, iris2$err),
#'              aes(x = Sepal.Length, y = Petal.Length, shape = Species),
#'              size = 3.5, col = "white", show.legend = FALSE) +
#'   geom_point(aes(x = Sepal.Length, y = Petal.Length, shape = Species),
#'              data = iris2, size = 2) +
#'   guides(colour = FALSE, shape = FALSE)
#'
#' discParam <- list(Sepal.L = c(5, 6, 7), Sepal.W = c(2.5, 3.5))
#' mod2 <- nb2(iris[, 1:4], iris[, 5], discretize = "fixed",
#'             discParams = discParam)
#' preds <- predict(mod2, iris[, c(1, 3)])
#' iris2 <- cbind(iris, err = (preds != iris$Species))
#' plot(mod2, features = c(1, 3), xlim = c(4, 8), ylim = c(0, 8)) +
#'   geom_point(aes(x = Sepal.Length, y = Petal.Length, shape = Species),
#'              data = subset(iris2, iris2$err), size = 3.5, col = "white",
#'              show.legend = FALSE) +
#'   geom_point(aes(x = Sepal.Length, y = Petal.Length, shape = Species),
#'   data = iris2, size = 2) +
#'   guides(colour = FALSE, shape = FALSE)
#'

plot.nb2 <- function(x, xlim, ylim, features, sdevs = 2, gridsize, data,
                     class.name = "class", gg = TRUE, ...) {

  std.gridsize <- 25

  tabs <- x$tables
  if (missing(features)) {
    feats <- 1:2
  } else if (!is.numeric(features)) {
    feats <- which(names(x$tables) == features[1])
    feats[2] <- which(names(x$tables) == features[2])
  } else feats <- features
  feat.names <- names(tabs)[feats]

  if (!is.null(x$disc)) {
    disc.names <- names(x$disc$params$breaks)
    disc.breaks <- x$disc$params$breaks
    disc <- TRUE
  } else disc <- FALSE

  if (missing(xlim)) {
    if (disc && feat.names[1] %in% disc.names) {
      xlim <- range(disc.breaks[[(feat.names[1])]][is.finite(disc.breaks[[(feat.names[1])]])])
    } else if (length(x$yc[[(feat.names[1])]]) == 0) {
      xlim <- NULL
    } else if (!missing(data)) {
      if ((checkmate::testList(data) && feat.names[1] %in% names(data$x)) ||
         (checkmate::testDataFrame(data) && feat.names[1] %in% names(data))) {
        if (is.data.frame(data)) {
          xlim <- range(data[, feat.names[1]], na.rm = TRUE)
        } else {
          xlim <- range(data$x[, feat.names[1]], na.rm = TRUE)
        }
      }
    } else {
      tab1 <- tabs[[(feats[1])]]
      xlim <- range(tab1[, 1] + t(sdevs * c(-1, 1) %*% t(tab1[, 2])), na.rm = TRUE)
    }
  }
  xlim <- suppressWarnings(sort(xlim))
  if (missing(ylim)) {
    if (disc && feat.names[2] %in% disc.names) {
      ylim <- range(disc.breaks[[(feat.names[2])]][is.finite(disc.breaks[[(feat.names[2])]])])
    } else if (length(x$yc[[(feat.names[2])]]) == 0) {
      ylim <- NULL
    } else if (!missing(data)) {
      if ((checkmate::testList(data) && feat.names[2] %in% names(data$x)) ||
          (checkmate::testDataFrame(data) && feat.names[2] %in% names(data))) {
        if (is.data.frame(data)) {
          ylim <- range(data[, feat.names[2]], na.rm = TRUE)
        } else {
          ylim <- range(data$x[, feat.names[2]], na.rm = TRUE)
        }
      }
    } else {
      tab2 <- tabs[[(feats[2])]]
      ylim <- range(tab2[, 1] + t(sdevs * c(-1, 1) %*% t(tab2[, 2])), na.rm = TRUE)
    }
  }
  ylim <- suppressWarnings(sort(ylim))

  if (disc && feat.names[1] %in% disc.names) {
    db1fin <- disc.breaks[[1]][is.finite(disc.breaks[[1]])]
    db1fin <- c(xlim[1], db1fin[db1fin > xlim[1] & db1fin < xlim[2]], xlim[2])
    xvals <- diff(db1fin) / 2 + (db1fin)[-length(db1fin)]
    gridsize.x <- length(xvals)
    width.x <- diff(db1fin)
  } else if (length(x$yc[[(feat.names[1])]]) == 0) {
    xvals <- dimnames(x$tables[[(feat.names[1])]])[[(feat.names[1])]]
    gridsize.x <- length(xvals)
    width.x <- NULL
    # stop("Not yet implemented for factorial features.")
  } else {
    gridsize.x <- if (missing(gridsize)) std.gridsize else gridsize
    xl <- seq(xlim[1], xlim[2], length.out = gridsize.x + 1)
    xvals <- diff(xl) / 2 + xl[-length(xl)]
    width.x <- NULL
  }
  if (disc && feat.names[2] %in% disc.names) {
    db2fin <- disc.breaks[[2]][is.finite(disc.breaks[[2]])]
    db2fin <- c(min(ylim), db2fin[db2fin > min(ylim) & db2fin < max(ylim)], max(ylim))
    yvals <- diff(db2fin) / 2 + (db2fin)[-length(db2fin)]
    gridsize.y <- length(yvals)
    width.y <- diff(db2fin)
  } else if (length(x$yc[[(feat.names[2])]]) == 0) {
    yvals <- dimnames(x$tables[[(feat.names[2])]])[[(feat.names[2])]]
    gridsize.y <- length(yvals)
    width.y <- NULL
    # stop("Not yet implemented for factorial features.")
  } else {
    gridsize.y <- if (missing(gridsize)) std.gridsize else gridsize
    yl <- seq(ylim[1], ylim[2], length.out = gridsize.y + 1)
    yvals <- diff(yl) / 2 + yl[-length(yl)]
    width.y <- NULL
  }
  ndata <- expand.grid(xvals, yvals, KEEP.OUT.ATTRS = FALSE)
  names(ndata) <- feat.names

  if ("ZZZ_TEMP" %in% x$levels) {
    preds <- matrix(1, nrow(ndata), 1)
    imgvals <- matrix(1, gridsize.x, gridsize.y)
  } else {
    preds <- predict(x, newdata = ndata, type = "raw")
    imgvals <- matrix(apply(preds, 1, max, na.rm = TRUE), gridsize.x, gridsize.y)
  }

  if (gg) {
    grid <- ndata
    target = class.name
    grid[, target] <- factor(apply(preds, 1, which.max), levels = 1:length(x$levels),
                             labels = x$levels)
                      # predict(x, newdata = grid)
    grid$.prob.pred.class = apply(preds, 1, max)
    if (!is.null(width.x) & !is.null(width.y)) {
      aes1 <- aes_string(fill = target, width = "width.x", height = "width.y", alpha = ".prob.pred.class")
      grid <- cbind(grid, width.x = rep(width.x, gridsize.y), width.y = rep(width.y, each = gridsize.x))
    } else if (!is.null(width.x)) {
      aes1 <- aes_string(fill = target, width = "width.x", alpha = ".prob.pred.class")
      grid <- cbind(grid, width.x = rep(width.x, gridsize.y))
    } else if (!is.null(width.y)) {
      aes1 <- aes_string(fill = target, height = "width.y", alpha = ".prob.pred.class")
      grid <- cbind(grid, width.y = rep(width.y, each = gridsize.x))
    } else {
      aes1 <- aes_string(fill = target, alpha = ".prob.pred.class")
    }

    # g <- guide_legend("class")
    p <- ggplot(grid, aes_string(x = feat.names[1L], y = feat.names[2L]))
    p <- p + geom_tile(aes1, grid, show.legend = TRUE,
                       width = if (is.null(xlim)) rep(0.75, length(xvals)) else NULL,
                       height = if (is.null(ylim)) rep(0.75, length(yvals)) else NULL)
    p <- p + scale_alpha(limits = c(1 / length(x$levels), 1))
    p <- p + guides(alpha = FALSE)
    # p <- p + guides(alpha = FALSE, fill = g)
    if (!is.null(xlim)) p <- p + xlim(xlim)
    if (!is.null(ylim)) p <- p + ylim(ylim)
    if (!missing(data)) {
      if (checkmate::testList(data) &
          checkmate::testSubset(c("x", "class"), names(data))) {
        data <- data.frame(data$x, class = data$class)
      }
      if (checkmate::testDataFrame(data) &
          checkmate::testSubset(class.name, colnames(data))) {
        if (!is.factor(data[, class.name]))
          data[, class.name] <- factor(data[, class.name])
        preds2 <- predict(x, data)
        data$err <- (preds2 != data[, class.name])
        p <- p + geom_point(data = subset(data, !data$err),
                     aes_string(x = feat.names[1L], y = feat.names[2L], shape = class.name),
                     size = 2.5, col = "black", show.legend = TRUE) +
          geom_point(data = subset(data, data$err),
                            aes_string(x = feat.names[1L], y = feat.names[2L], shape = class.name),
                            size = 4, col = "white", show.legend = FALSE) +
          geom_point(data = subset(data, data$err),
                     aes_string(x = feat.names[1L], y = feat.names[2L], shape = class.name),
                     size = 2.5, col = "black", show.legend = FALSE) +
          scale_shape_discrete(drop = FALSE) +
          # guides(colour = FALSE, shape = guide_legend(class.name))
          guides(colour = FALSE)
      } else {
        warning("No points included. Did you set argument class.name?")
      }
    }
    return(p)
  } else {
    image(xvals, yvals, imgvals, col = grey.colors(std.gridsize))
    points(ndata, col = rgb(t(col2rgb(apply(preds, 1, which.max), FALSE)),
                            alpha = 0.25 * 255, maxColorValue = 255), pch = 15)
    return(invisible(NULL))
  }

}



#' Printing nb2 objects.
#'
#' Modified copy of \link[e1071]{print.naiveBayes} from \code{e1071}.
#'
#' @param x \code{[nb2]}\cr
#'   nb2 object to print.
#'
#' @param ...
#'   Currently ignored.
#'
#' @param digits \code{[numeric]}\cr
#'
#' @param print.apl \code{[logical]}\cr
#'   Print apriori.list? Else: apriori.
#'
#' @return
#'   NULL
#'
#' @author
#' \itemize{
#'   \item David Meyer \email{David.Meyer@@R-project.org}: original \link[e1071]{print.naiveBayes}
#'   \item Philipp Aschersleben \email{aschersleben@@statistik.tu-dortmund.de}:
#'   modifications and extensions for the NBCD package
#' }
#'
#' @export
#'
#' @seealso \code{\link{nb2}}, \code{\link{update.nb2}}, \code{\link{plot.nb2}},
#' \code{\link{predict.nb2}}
#'

print.nb2 <- function(x, ..., digits = getOption("digits"), print.apl = FALSE) {
  cat("\nNaive Bayes Classifier for Discrete Predictors\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nA-priori probabilities:\n")
  if (is.null(x$apriori.list) | !print.apl) {
    xapr <- x$apriori
    zzz <- which(names(xapr) == "ZZZ_TEMP")
    if (length(zzz) > 0) xapr <- xapr[-zzz, drop = FALSE]
    xapr[] <- paste(format(xapr / sum(xapr), digits = digits),
                         paste0("(", xapr, "L)"))
    print(xapr)
  } else {
    invisible(mapply(function(y, z) {
      zzz <- which(names(z) == "ZZZ_TEMP")
      if (length(zzz) > 0) z <- z[-zzz, drop = FALSE]
      z[] <- paste(format(z / sum(z), digits = digits), paste0("(", z, "L)"))
      names(dimnames(z)) <- paste(names(dimnames(z)), paste0("(", y, ")"))
      print(z); cat("\n")
    }, names(x$apriori.list), x$apriori.list))
  }

  cat("\nConditional probabilities:\n")
  for (i in x$tables) {
    if (exists("zzz") && length(zzz) > 0) print(i[-zzz, , drop = FALSE]) else print(i)
    cat("\n")
  }
}


