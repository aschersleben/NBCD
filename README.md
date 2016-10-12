# NBCD: Naive Bayes for Online Classification with Concept Drift

This package provides an online classification method based on Naive Bayes, that
is able to handle concept drift. Furthermore, it comes with extended Naive Bayes
functions, that can be printed, plotted, predicted and updated (see below). The
same holds for the NBCD method.

[![Travis-CI Build Status](https://travis-ci.org/aschersleben/NBCD.svg?branch=master)](https://travis-ci.org/aschersleben/NBCD)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/aschersleben/NBCD?branch=master&svg=true)](https://ci.appveyor.com/project/aschersleben/NBCD)



## Installation:
```r
devtools::install_github("aschersleben/NBCD")
library("NBCD")
```


## Simple example:
We use the well-known `iris` dataset and add a "concept drift":
```r
set.seed(1234)
iris2 <- iris[sample(150), ]
iris2$Sepal.Width <- iris2$Sepal.Width + seq(1, 30, len = 150)
model <- NULL
for (i in 1:120) {
  model <- makeNBCDmodel(list(x = iris2[i, 1:4], class = iris2[i, 5], time = i), model = model,
                         discretize = "fixed", discParams = list(Sepal.Length = 4:8),
                         init.obs = 20, max.waiting.time = 20, waiting.time = "auto")
}
print(model)
```

For plotting, the NBCD package uses `ggplot2`.
```r
plot(model, ylim = c(25, 35))
plot(model, ylim = c(25, 35), use.lm = TRUE, time = 150)
```

You can directly add data to the plots, predictions are included automatically:
```r
plot(model, ylim = c(20, 40), use.lm = FALSE,
     data = iris2[140:150, ], class.name = "Species")
plot(model, ylim = c(20, 40), use.lm = TRUE, time = 150,
     data = iris2[140:150, ], class.name = "Species")
```


## About the NBCD method

[...]


## Extended Naive Bayes

This package also includes an extended version of `e1071::naiveBayes` called
`nb2`. It can be updated with new observations and includes an automated
discretization.

At the first look, there is no difference to the `e1071` function:
```r
mod <- nb2(iris[, 1:4], iris[, 5])
print(mod)
```

But you can not only print but also plot the model:
```r
plot(mod)
plot(mod, data = iris, class.name = "Species")
```

Easy discretization:
```r
discParam <- list(Sepal.L = 4:8, Sepal.W = 1:5)
mod2 <- nb2(iris[, 1:4], iris[, 5], discretize = "fixed", discParams = discParam)
print(mod2)
plot(mod2, data = iris, class.name = "Species")
```

Easy updates (= adding new observations to the model without re-computing):
```r
mod.upd <- update(mod, newdata = iris[1:50, 1:4], y = iris$Species[1:50])
print(mod.upd)
```


## Concept Drift

Read about concept drift in Webb et al. (2016, [DOI:10.1007/s10618-015-0448-4](http://dx.doi.org/DOI:10.1007/s10618-015-0448-4)).
