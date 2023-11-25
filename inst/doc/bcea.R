## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, warning=FALSE, message=FALSE--------------------------------------
library(BCEA)

## -----------------------------------------------------------------------------
data(Smoking)

## -----------------------------------------------------------------------------
treats <- c("No intervention", "Self-help", "Individual counselling", "Group counselling")

## -----------------------------------------------------------------------------
bcea_smoke <- bcea(eff, cost, ref = 4, interventions = treats, Kmax = 500)

## ----fig.width=10, fig.height=10----------------------------------------------
library(ggplot2)
library(purrr)

plot(bcea_smoke)

## -----------------------------------------------------------------------------
ceplane.plot(bcea_smoke, comparison = 2, wtp = 250)

eib.plot(bcea_smoke)

contour(bcea_smoke)

ceac.plot(bcea_smoke)

ib.plot(bcea_smoke)

## -----------------------------------------------------------------------------
plot(bcea_smoke,
     graph = "ggplot2",
     wtp = 250,
     line = list(color = "red", size = 1),
     point = list(color = c("plum", "tomato", "springgreen"), shape = 3:5, size = 2),
     icer = list(color = c("red", "orange", "black"), size = 5))

