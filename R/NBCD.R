#' Make or update NBCD model.
#'
#' @param new.obs \code{[list]}\cr
#'   Has elements "x" \code{[data.frame]} and "class" \code{[factor]}
#'   (and "time", if missing: +1)
#'
#' @param model \code{[list]}\cr
#'   Model from makeNBCDmodel or empty.
#'
#' @param max.waiting.time \code{[numeric | list]}\cr
#'   Can be \code{numeric(1)}, then all features with same waiting time. In the
#'   other case, you have to give an individual waiting time for all features
#'   as an (un)named vector or list.
#'
#' @param init.obs \code{[numeric(1)]}\cr
#'
#' @param verbose \code{[logical(1)]}\cr
#'
#' @param ... \code{[any]}\cr
#'   Arguments passed to update.nb2(...).
#'
#' @param waiting.time \code{[character(1)]}\cr
#'   Default is "fixed". If "auto", then calculate waiting time between two
#'   models (for factorial variables) with getWaitingTime2(). If "fixed", use
#'   max.waiting.time.
#'
#' @param min.waiting.time \code{[numeric(1)]}\cr
#'   Minimal waiting time, if waiting.time = "auto".
#'
#' @param wait.all.classes \code{[logical(1)]}\cr
#'   If waiting.time = "auto", should algo wait for all classes to be appeared?
#'
#' @param k.waiting.time \code{[numeric(1)]}
#'   If \code{waiting.time = "auto"}, the calculated waiting time willl be
#'   multiplied by this argument's value.
#'
#' @return \code{[NBCD]}\cr
#' NBCD model
#'
#' @export
#'

makeNBCDmodel <- function(new.obs, model, max.waiting.time, init.obs,
                          verbose = FALSE, ...,
                          waiting.time = c("fixed", "auto"),
                          min.waiting.time = 3,
                          wait.all.classes = TRUE,
                          k.waiting.time = 1) {

  waiting.time <- match.arg(waiting.time)

  asscoll <- checkmate::makeAssertCollection()
  checkmate::assertNumber(min.waiting.time, lower = 0, add = asscoll)
  checkmate::assertFlag(wait.all.classes, add = asscoll)
  checkmate::assertNumber(k.waiting.time, lower = 0, finite = TRUE, add = asscoll)
  checkmate::assertList(new.obs, add = asscoll)
  checkmate::assertSubset(c("x", "class"), names(new.obs), add = asscoll)
  checkmate::assertDataFrame(new.obs$x, add = asscoll)
  checkmate::reportAssertions(asscoll)
  checkmate::assert(checkNumeric(max.waiting.time, lower = 0, finite = TRUE,
                                 len = ncol(new.obs$x), any.missing = FALSE),
                    checkNumber(max.waiting.time, finite = TRUE),
                    checkList(max.waiting.time, any.missing = FALSE, len = ncol(new.obs$x)))
  checkmate::assert(checkNumeric(init.obs, lower = 0, finite = TRUE,
                                 len = ncol(new.obs$x), any.missing = FALSE),
                    checkNumber(init.obs, finite = TRUE),
                    checkList(init.obs, any.missing = FALSE, len = ncol(new.obs$x)))

  if (nrow(new.obs$x) > 1) {
    if (verbose)
      cat("split new.obs")
    for (i in seq_len(nrow(new.obs$x))) {
      nobi <- list(x = new.obs$x[i, ], class = new.obs$class[i])
      if (!is.null(new.obs$time))
        nobi <- c(nobi, time = new.obs$time[i])
      model <- makeNBCDmodel(new.obs = nobi, model = model,
                             max.waiting.time = max.waiting.time,
                             init.obs = init.obs, verbose = verbose, ...,
                             waiting.time = waiting.time,
                             min.waiting.time = min.waiting.time,
                             wait.all.classes = wait.all.classes,
                             k.waiting.time = k.waiting.time)
    }
    return(model)
  }

  if (verbose) {
    if (inherits(try(is.environment(.NBCD), silent = TRUE), "try-error"))
      .NBCD <- new.env()
    if (inherits(try(get("cat.flag", envir = .NBCD), silent = TRUE), "try-error"))
      assign("cat.flag", TRUE, envir = .NBCD)
    if (get("cat.flag", pos = .NBCD))
      cat("makeNBCDmodel: ")
  }

  if (missing(model) || length(model) == 0) {
    model <- initModel()
  }
  curr.model <- model$current
  old.model <- model$old

  if (checkmate::testNumber(max.waiting.time)) {
    mwt.list <- vector("list", ncol(new.obs$x))
    names(mwt.list) <- names(new.obs$x)
    max.waiting.time <- lapply(mwt.list, function(x) max.waiting.time)
  } else {
    # if (checkmate::testNumeric(max.waiting.time)) {
    max.waiting.time <- as.list(max.waiting.time)
    if (is.null(names(max.waiting.time))) {
      names(max.waiting.time) <- names(new.obs$x)
    } else {
      names(max.waiting.time) <- match.arg(names(max.waiting.time), names(new.obs$x), several.ok = TRUE)
    }
  # } else {
  #   names(max.waiting.time) <- match.arg(names(max.waiting.time), names(new.obs$x), several.ok = TRUE)
  }

  if (checkmate::testNumber(init.obs)) {
    init.list <- vector("list", ncol(new.obs$x))
    names(init.list) <- names(new.obs$x)
    init.obs <- lapply(init.list, function(x) init.obs)
  } else {
    # if (checkmate::testNumeric(init.obs)) {
    init.obs <- as.list(init.obs)
    if (is.null(names(init.obs))) {
      names(init.obs) <- names(new.obs$x)
    } else {
      names(init.obs) <- match.arg(names(init.obs), names(new.obs$x), several.ok = TRUE)
    }
    # } else {
    #   names(init.obs) <- match.arg(names(init.obs), names(new.obs$x), several.ok = TRUE)
  }

  # args <- as.list(match.call())
  # args <- mget(names(formals()), sys.frame(sys.nframe()))
  # args[[1]] <- NULL
  tn <- names(formals())
  tn <- tn[-which(tn == "...")]
  args <- c(mget(tn, sys.frame(sys.nframe())), list(...))
  args$new.obs <- args$model <- NULL
  # args <- lapply(args, function(a) if (is.symbol(a)) eval(a) else a)
  model$args <- args

  # Fuege neue Beobachtungen dem Modell hinzu (wird angelegt, falls noch nicht vorhanden)
  if (is.null(curr.model$general$nb2)) {
    curr.model$general$nb2 <- nb.mod <- nb2(x = new.obs$x, y = new.obs$class, ...)
  } else {
    curr.model$general$nb2 <- nb.mod <- update(object = curr.model$general$nb2,
                                               newdata = new.obs$x, y = new.obs$class, ...)
  }

  # Falls Listeneintraege fuer die Variablen noch nicht in curr.model enthalten
  missing.vars <- names(nb.mod$tables)[!(names(nb.mod$tables) %in% names(curr.model))]
  for (i in missing.vars) {
    curr.model[[i]] <- vector("list", 5)
    names(curr.model[[i]]) <- c("nobs", "time", "wait", "type", "init")
    curr.model[[i]]$type <- if (is.numeric(new.obs$x[[i]]) & !(i %in% names(model$args$discParams)))
      "numeric" else "factor"
    curr.model[[i]]$init <- TRUE
  }

  curr.model <- lapply(curr.model, function(x) {
    if (x$nobs == 0 || is.null(x$nobs))
      x$nobs <- 1 else x$nobs <- x$nobs + 1
    a <- as.character(new.obs$class)
    xtime <- x$time[[a]]
    if (is.null(new.obs$time)) {
      if (xtime$time1 == 0 || is.null(xtime$time1)) {
        xtime$time1 <- xtime$time2 <- 1
      } else {
        xtime$time2 <- xtime$time2 + 1
      }
    } else {
      xtime$time2 <- new.obs$time
      if (is.null(xtime$time1) || is.na(xtime$time1))
        xtime$time1 <- new.obs$time
    }
    x$time[[as.character(new.obs$class)]] <- xtime
    return(x)
  })
  curr.model$general$time.last <- max(unlist(curr.model$general$time), na.rm = TRUE)

  if (verbose) {
    if (get("cat.flag", pos = .NBCD)) {
      cat("Obs. in model:", "")
      assign("cat.flag", FALSE, envir = .NBCD)
    }
    cat(paste0(curr.model$general$nobs, ", "))
  }

  model$current <- curr.model

  # init.flag <- old.flag <- TRUE
  # Fuer alle Variablen:
  for (i in names(nb.mod$tables)) {

    # Pruefe, ob init.time / waiting.time erreicht
    if ((curr.model[[i]]$init & curr.model[[i]]$nobs >= init.obs[[i]]) |
        ((curr.model[[i]]$nobs >= ifelse(is.null(old.model[[i]]), Inf, old.model[[i]]$wait)) &
         if (curr.model[[i]]$type == "numeric") TRUE else
           (!wait.all.classes | isTRUE(all(sapply(
             curr.model$general$nb2$tables, function(x) all(rowSums(x) > 0))
           ))))) {

      if (curr.model[[i]]$init & verbose)
        cat("Init", i, "finished. ")
      if (curr.model[[i]]$init) {
        curr.model[[i]]$init <- FALSE
        miss.class <- (curr.model$general$nb2$apriori == 0)
        if (any(miss.class)) {
          mc <- which(miss.class)
          for (j in names(nb.mod$tables)) {
            rownames(curr.model$general$nb2$tables[[j]])[mc] <- "ZZZ_TEMP"
            names(curr.model$general$nb2$apriori.list[[j]])[mc] <- "ZZZ_TEMP"
          }
          names(curr.model$general$nb2$apriori)[mc] <- "ZZZ_TEMP"
          curr.model$general$nb2$levels[mc] <- "ZZZ_TEMP"
        }
      }

      # if (is.null(old.model[[i]]))
      #   old.flag <- FALSE

      old.model[[i]] <- curr.model[[i]]

      # Modell zur Vorhersage (PredictionModel) fuer die naechsten Beobachtungen erstellen
      # if (old.flag)
        old.model[[i]]$pred.mod <- setPredictionModel(old.model = model$old,
                                                      new.model = curr.model,
                                                      var.name = i)

      # Setze neue waiting.time
      if (curr.model[[i]]$type == "numeric") {
        curr.model[[i]]$wait <- max.waiting.time[[i]]
      } else {
        if (waiting.time == "fixed") {
          curr.model[[i]]$wait <- max.waiting.time[[i]]
        } else {
          wt <- max(min.waiting.time, getWaitingTime2(curr.model$general$nb2$tables[[i]]))
          if (!missing(k.waiting.time))
            wt <- wt * k.waiting.time
          curr.model[[i]]$wait <- min(wt, max.waiting.time[[i]])
        }
      }

      old.model[[i]]$nb2 <- curr.model$general$nb2
      old.model[[i]]$wait <- curr.model[[i]]$wait
      old.model[[i]]$init <- curr.model[[i]]$init

      model <- resetModel(model = model, var.name = i)
      model$old <- old.model
      if (verbose) {
        cat("reset model for", i, "")
        assign("cat.flag", TRUE, envir = .NBCD)
      }
    }
  }

  if (model$current$general$init & all(!sapply(model$current, "[[", "init")[-1]))
    model$current$general$init <- FALSE

  if (verbose && get("cat.flag", pos = .NBCD))
    cat("\n")

  return(structure(model, class = "NBCD"))
}



#' Initialize NBCD model.
#'
#' @return \code{[NBCD]}\cr
#' Empty model.
#'

initModel <- function() {
  list(current = list(general = list(nb2 = NULL, time = list(),
                                     nobs = 0, init = TRUE)),
       old = NULL,
       args = NULL)
}



#' Reset Model when waiting.time reached
#'
#' @param model \code{[NBCD]}\cr
#'   "model$current"
#'
#' @param var.name \code{[character]}\cr
#'
#'
#' @return \code{[NBCD]}\cr
#'   NBCD model with resetted "current" part.
#'

resetModel <- function(model, var.name) {
  model$current$general$nb2$apriori.list[[var.name]][] <- 0
  model$current$general$nb2$tables[[var.name]][] <- NA
  model$current$general$nb2$yc[[var.name]] <- lapply(model$current$general$nb2$yc[[var.name]], function(x) NULL)
  model$current$general$nobs <- 0
  model$current$general$init <- FALSE
  model$current$general$time <- lapply(model$current$general$time,
                                           lapply, function(x) NA)
  model$current[[var.name]]$nobs <- 0
  model$current[[var.name]]$time <- lapply(model$current[[var.name]]$time,
                                           lapply, function(x) NA)
  return(model)
}



#' Create a prediction model from parts of an NBCD model.
#'
#' @param old.model \code{[list]}\cr
#'   "old" part of NBCD model.
#'
#' @param new.model \code{[list]}\cr
#'   "current" part of NBCD model.
#'
#' @param var.name \code{[character]}\cr
#'   Name of variable for which to create prediction model.
#'
#' @param n.models \code{[numeric(1)]}\cr
#'   Number of past model values for the linear prediction.
#'
#' @return \code{[list]}\cr
#'   To be used as "pred.mod" part of "old" part of NBCD model.
#'   For functions makeNBCDmodel(...) and predict.NBCD(...).
#'

setPredictionModel <- function(old.model, new.model, var.name, n.models = Inf) {

  checkmate::assertNumber(n.models, lower = 2)

  old.nb2.tabs <- old.model[[var.name]]$nb2$tables[[var.name]]
  new.nb2.tabs <- new.model$general$nb2$tables[[var.name]]
  if (is.null(old.nb2.tabs)) {
    old.nb2.tabs <- new.nb2.tabs
    old.model <- new.model
    old <- FALSE
  } else old <- TRUE
  class.names <- rownames(new.nb2.tabs)
  if ("ZZZ_TEMP" %in% class.names)
    class.names <- class.names[-which(class.names == "ZZZ_TEMP")]
  class.names.old <- rownames(old.nb2.tabs)
  if ("ZZZ_TEMP" %in% class.names.old) {
    if (!("ZZZ_TEMP" %in% class.names)) {
      newname <- class.names[!(class.names %in% class.names.old)]
      for (i in newname) old.nb2.tabs <- rbind(old.nb2.tabs, old.nb2.tabs["ZZZ_TEMP", ])
      rownames(old.nb2.tabs) <- c(class.names.old, newname)
    }
  }

  out <- list()

  if (new.model[[var.name]]$type == "numeric") {
    for (i in class.names) {
      if (old) {
        old.time <- suppressWarnings(mean(c(old.model[[var.name]]$time[[i]]$time1,
                                            old.model[[var.name]]$time[[i]]$time2)))
      } else {
        old.time <- old.nb2.tabs <- NULL
      }
      new.time <- mean(c(new.model[[var.name]]$time[[i]]$time1,
                         new.model[[var.name]]$time[[i]]$time2))
      old.lm.data <- old.model[[var.name]]$pred.mod[[i]]$data
      lm.data <- data.frame(mean = c(old.nb2.tabs[i, 1], new.nb2.tabs[i, 1]),
                            sd = c(old.nb2.tabs[i, 2], new.nb2.tabs[i, 2]),
                            time = c(old.time, new.time))
      if (!is.null(old.lm.data))
        lm.data <- rbind(tail(head(old.lm.data, -1), n.models), lm.data)
      out[[i]]$mean <- lm(mean ~ time + 1, data = lm.data)
      out[[i]]$sd <- lm(sd ~ time + 1, data = lm.data)
      out[[i]]$data <- lm.data
    }
  } else {
    out <- lapply(vector("list", length(class.names)),
                  function(x) vector("list", ncol(new.nb2.tabs)))
    names(out) <- class.names
    for (i in class.names) {
      for (j in 1:ncol(new.nb2.tabs)) {
        if (old) {
          old.time <- suppressWarnings(mean(c(old.model[[var.name]]$time[[i]]$time1,
                                              old.model[[var.name]]$time[[i]]$time2)))
        } else {
          old.time <- old.nb2.tabs <- NULL
        }
        new.time <- mean(c(new.model[[var.name]]$time[[i]]$time1,
                           new.model[[var.name]]$time[[i]]$time2))
        old.lm.data <- old.model[[var.name]]$pred.mod[[i]]$data
        lm.data <- data.frame(y = c(old.nb2.tabs[i, j], new.nb2.tabs[i, j]),
                              time = c(old.time, new.time))
        if (!is.null(old.lm.data))
          lm.data <- rbind(tail(head(old.lm.data, -1), n.models), lm.data)
        if (all(is.na(lm.data$time))) {
          out[[i]][[j]] <- old.model[[var.name]]$pred.mod[[i]][[j]]
        } else
          out[[i]][[j]] <- lm(y ~ time + 1, lm.data)
        out[[i]][[j]]$data <- lm.data
      }
    }
  }
  return(out)
}



#' Get a prediction model from NBCD model.
#'
#' @param model \code{[NBCD]}\cr
#'   NBCD model.
#'
#' @param pred.time \code{[numeric(1)]}\cr
#'   Time point of prediction.
#'
#' @param use.lm \code{[logical(1)]}\cr
#'   Use lm models to forecast the movements of mean and sd? (Default TRUE)
#'   If FALSE, use mean and sd of old.model.
#'
#' @param n.models \code{[numeric(1)]}\cr
#'   Number of past model values for the linear prediction.
#'
#' @export
#'
#' @return \code{[nb2]}\cr
#'   nb2 model
#'

getPredictionModel <- function(model, pred.time, use.lm = TRUE, n.models = Inf) {

  if (is.logical(use.lm)) {
    use.lm = if (isTRUE(use.lm)) "mean" else "none"
  } else use.lm = match.arg(use.lm, c("none", "mean", "both"))

  if (missing(pred.time)) {
    pred.time <- model$current$general$time.last
    message("Missing pred.time argument. ",
            "Predicition for \"last.time\" = ", pred.time)
  }

  asscoll <- checkmate::makeAssertCollection()
  checkmate::assertNumber(pred.time, na.ok = FALSE, lower = 0, finite = TRUE, add = asscoll)
  checkmate::assertNumber(n.models, lower = 2, add = asscoll)
  checkmate::reportAssertions(asscoll)

  pred.model <- model$old
  var.names <- names(pred.model)
  pred.mod.list <- sapply(var.names, function(x) pred.model[[x]]$pred.mod,
                          simplify = FALSE)
  ndata <- data.frame(time = pred.time)

  out.model <- pred.model[[1]]$nb2

  for (i in var.names) {                           # ueber die features
    class.names <- rownames(out.model$tables[[i]])
    if ("ZZZ_TEMP" %in% class.names)
      class.names <- class.names[-which(class.names == "ZZZ_TEMP")]
    for (j in class.names) {   # ueber die target classes
      tmp.pred.ij <- pred.mod.list[[i]][[j]]

      if (pred.model[[i]]$type == "numeric") {
        # bei numerischer Variable
        # tmp.pred.ij ist list(mean, sd) mit lm-Modellen
        if (use.lm == "both") {
          if (is.finite(n.models)) {
            lm.data <- tail(pred.model[[i]]$pred.mod[[j]]$data, n.models)
            mean.mod <- lm(mean ~ time + 1, data = lm.data)
            sd.mod <- lm(sd ~ time + 1, data = lm.data)
          } else {
            mean.mod <- tmp.pred.ij$mean
            sd.mod <- tmp.pred.ij$sd
          }
          mean.pred <- predict(mean.mod, newdata = ndata)
          sd.pred <- predict(sd.mod, newdata = ndata)
          ms.pred <- c(mean.pred, sd.pred)
        } else if (use.lm == "mean") {
          if (is.finite(n.models)) {
            lm.data <- tail(pred.model[[i]]$pred.mod[[j]]$data, n.models)
            mean.mod <- lm(mean ~ time + 1, data = lm.data)
          } else {
            mean.mod <- tmp.pred.ij$mean
          }
          mean.pred <- predict(mean.mod, newdata = ndata)
          sd.pred <- pred.model[[i]]$nb2$tables[[i]][j, 2]
          ms.pred <- c(mean.pred, sd.pred)
        } else {
          ms.pred <- pred.model[[i]]$nb2$tables[[i]][j, ]
        }
        out.model$tables[[i]][j, ] <- ms.pred

      } else {
        # bei kategorieller Variable
        # tmp.pred.ij ist vector("list", number.of.feature.classes) mit lm-Modellen
        if (use.lm != "none") {
          for (k in seq_along(tmp.pred.ij)) {
            if (is.finite(n.models)) {
              lm.data <- tail(pred.model[[i]]$pred.mod[[j]][[k]]$data, n.models)
              lm.mod <- lm(y ~ time + 1, data = lm.data)
            } else {
              lm.mod <- tmp.pred.ij[[k]]
            }
            out.model$tables[[i]][j, k] <- predict(lm.mod, newdata = ndata)
          }
        } else {
          out.model$tables[[i]][j, ] <- pred.model[[i]]$nb2$tables[[i]][j, ]
        }
      }
    }
    out.model$apriori.list[[i]] <- pred.model[[i]]$nb2$apriori.list[[i]]
  }
  return(out.model) # forward to predict.nb2
}



#' Predict class values for new data with NBCD model.
#'
#' @param object \code{[NBCD]}\cr
#'   NBCD model.
#'
#' @param newdata \code{[data.frame]}\cr
#'   Containing variables.
#'
#' @param time \code{[numeric(1)]}\cr
#'   Time for which to predict.
#'
#' @param use.lm \code{[logical(1)]}\cr
#'   Use lm models to forecast the movements of mean and sd? (Default TRUE)
#'   If FALSE, use mean and sd of old.model.
#'
#' @param ...
#'   Arguments passed to predict.nb2(...).
#'
#' @param n.models \code{[numeric(1)]}\cr
#'   Number of past model values for the linear prediction.
#'
#' @return
#'   Results from predict.nb2(...).
#'
#' @export
#'

predict.NBCD <- function(object, newdata, time, use.lm, n.models = Inf, ...) {
  checkmate::assertNumber(time)
  checkmate::assertNumber(n.models)
  predict(
    suppressWarnings(getPredictionModel(object, pred.time = time, use.lm = use.lm,
                                        n.models = n.models)),
    newdata = newdata, ...)
}



#' Plot predicted class values of NBCD model.
#'
#' @param x \code{[NBCD]}\cr
#'   NBCD model.
#'
#' @param time \code{[numeric(1)]}\cr
#'   Time for which to predict (passed to predict.NBCD(...)).
#'
#' @param use.lm \code{[logical(1)]}\cr
#'   Use lm models to forecast the movements of mean and sd? (Default TRUE)
#'   If FALSE, use mean and sd of old.model.
#'   (passed to predict.NBCD(...))
#'
#' @param ... \code{[any]}\cr
#'   Arguments passed to plot.nb2(...).
#'
#' @param n.models \code{[numeric(1)]}\cr
#'   Number of past model values for the linear prediction.
#'
#' @param .verb \code{[logical(1)]}\cr
#'   Show messages? Default TRUE.
#'
#' @export
#'

plot.NBCD <- function(x, time, use.lm = FALSE, ..., n.models = Inf,
                      .verb = TRUE) {
  if (missing(time)) {
    time <- x$current$general$time.last
    if (.verb) message("Missing time argument. ",
                       "Predicition for \"last.time\" = ", time)
  }
  p <- suppressWarnings(getPredictionModel(x, pred.time = time, use.lm = use.lm,
                                           n.models = n.models))
  plot(p, ...)
}



#' Print information about NBCD model.
#'
#' @param x \code{[NBCD]}\cr
#'   NBCD model.
#'
#' @param size \code{[character]}\cr
#'   small or big
#'
#' @param ... \code{[any]}\cr
#'   Currently ignored.
#'
#' @param use.lm \code{[logical]}\cr
#'   use lm for prediction?
#'
#' @param len \code{[numeric(1)]}\cr
#'   max. number of printed names/numbers/...
#'
#' @export
#'

print.NBCD <- function(x, size = c("small", "big"), ..., use.lm = FALSE, time,
                       len = 3) {

  size <- match.arg(size)
  # general:
  #   names of variables + "type" (general)
  #   arguments passed to makeNBCDmodel (--> stored in model$args)
  #   ...
  #
  # current:
  #   general
  #     $init
  #     $nobs in current model (apriori.list) + time ($args$waiting.time == "auto"?)
  #
  # old:
  #   $nobs in pred.model
  #   $wait
  #

  var.names <- names(x$current)[-1]
  var.types <- sapply(var.names, function(y) x$current[[y]]$type)
  classes <- x$current$general$nb2$levels
  args <- x$args

  initflag <- x$current$general$init
  init <- if (initflag) "# INIT #" else ""
  all.nobs <- sum(x$current$general$nb2$apriori)
  curr.nobs <- sapply(x$current, "[[", "nobs")[-1]
  if (missing(time)) time <- x$current$general$time.last

  old.nobs <- sapply(x$old, "[[", "nobs")
  waits <- sapply(x$old, "[[", "wait")

  if (size == "big") {

    cat("\n# # # NBCD Model # # # \n\n")

    cat("Variables: \n")
    for (i in var.names) {
      cat("  >", i, paste0("(", var.types[i]))
      if (var.types[i] == "factor") {
        cat("; levels: ")
        cat(pastehead(colnames(x$current$general$nb2$tables[[i]]), len = len))
      }
      cat(") \n")
    }

    cat("\nClasses: \n  >", pastehead(classes, len = len), "\n\n")

    cat("Number of Observations: \n")
    cat("  > total:", all.nobs, init, "\n")
    cat("  > new model:", pastehead(paste0(curr.nobs, " (", names(curr.nobs), ")"), len = len), "\n")
    if (!is.null(x$old)) {
      cat("  > pred. model:", pastehead(paste0(old.nobs, " (", names(old.nobs), ")"), len = len), "\n\n")

      cat("Waiting Time until Model Update: \n")
      cat("  >", pastehead(paste0(ceiling(waits), " (", names(waits), ")"), len = len), "\n\n")

      cat(paste0("Prediction Model for Current Time (", time,") ", ifelse(use.lm, "w/", "w/o"), " \"lm\" Usage:"), "\n")
      predmod <- capture.output(print(getPredictionModel(x, time, use.lm = use.lm), print.apl = TRUE))
      cat(paste("  >", head(predmod[-(1:6)], -1)), sep = "\n")
    }

    cat("\nArguments Passed to Function: \n")
    cat("  > max.waiting.time:", pastehead(x$args$max.waiting.time, len = len), "\n")
    cat("  > init.obs:", pastehead(x$args$init.obs, len = len), "\n")
    cat("  > waiting.time:", x$args$waiting.time, "\n")
    cat("  > min.waiting.time:", x$args$min.waiting.time, "\n")
    cat("  > wait.all.classes:", x$args$wait.all.classes, "\n")

  } else {

    cat("NBCD Model:", length(var.names), "Variables",
        paste0("(", pastehead(var.names, len = len), ")"),
        "with", length(classes), "Classes",
        paste0("(", pastehead(classes, len = len), ")"), "\n")

    if (initflag) {
      cat(" ", init," ")
    } else {
      cat("           ", all.nobs, "Observations,", pastehead(curr.nobs, len = len),
          "in Current and", pastehead(old.nobs, len = len), "in Pred. Model\n")
      cat("            Waiting Times:", pastehead(ceiling(waits), len = len))
    }
    cat("\n")

  }

  return(invisible(x))
}
