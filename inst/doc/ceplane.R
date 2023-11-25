## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>",
fig.width = 6
)

## ----setup, results='hide', message=FALSE, warning=FALSE, echo=FALSE----------
library(BCEA)
library(dplyr)
library(reshape2)
library(ggplot2)
library(purrr)

## -----------------------------------------------------------------------------
data("Vaccine")

he <- bcea(eff, cost)

## -----------------------------------------------------------------------------
ceplane.plot(he, graph = "base")
ceplane.plot(he, graph = "ggplot2")
# ceac.plot(he, graph = "plotly")

## -----------------------------------------------------------------------------
ceplane.plot(he,
             graph = "ggplot2",
             title = "my title",
             line = list(color = "green", size = 3),
             point = list(color = "blue", shape = 10, size = 5),
             icer = list(color = "orange", size = 5),
             area = list(fill = "grey"),
             theme = theme_linedraw())

## -----------------------------------------------------------------------------
ceplane.plot(he,
             graph = "ggplot2",
             point = list(size = NA),
             icer = list(size = 5))

## -----------------------------------------------------------------------------
data("Smoking")

he <- bcea(eff, cost, ref = 4)
# str(he)

## -----------------------------------------------------------------------------
ceplane.plot(he)
ceplane.plot(he, graph = "ggplot2")

## -----------------------------------------------------------------------------
ceplane.plot(he,
             graph = "ggplot2",
             title = "my title",
             line = list(color = "red", size = 1),
             point = list(color = c("plum", "tomato", "springgreen"), shape = 3:5, size = 2),
             icer = list(color = c("red", "orange", "black"), size = 5))

## -----------------------------------------------------------------------------
ceplane.plot(he, pos = FALSE) # bottom right
ceplane.plot(he, pos = c(0, 0))
ceplane.plot(he, pos = c(0, 1))
ceplane.plot(he, pos = c(1, 0))
ceplane.plot(he, pos = c(1, 1))

