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
ceplane.plot(he, graph = "ggplot2", text = list(size = 20))
ceplane.plot(he, graph = "ggplot2", text = list(size = rel(2)))  # relative scaling, double size

# equivalent but more flexible and direct
ceplane.plot(he, graph = "ggplot2") +
  theme(axis.text = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

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

## -----------------------------------------------------------------------------
ceplane.plot(he, graph = "ggplot2")

## -----------------------------------------------------------------------------
ceplane.plot(he, graph = "ggplot2", wtp = 10000)

## -----------------------------------------------------------------------------
ceplane.plot(he, graph = "ggplot2", wtp = list(value = 10000))
ceplane.plot(he, graph = "ggplot2", wtp = list(value = 10000, colour = "blue"))
ceplane.plot(he, graph = "ggplot2", wtp = list(colour = "blue"))
ceplane.plot(he, graph = "ggplot2", wtp = list(size = 5))

## -----------------------------------------------------------------------------
ceplane.plot(he, graph = "ggplot2", wtp = list(size = 0))

## -----------------------------------------------------------------------------
data("Vaccine")
he <- bcea(eff, cost, ref=2)

ceplane.plot(he, graph = "ggplot2")
ceplane.plot(he, graph = "ggplot2", wtp = 1000, xlim = c(-0.0005, 0.0015))
ceplane.plot(he, graph = "ggplot2", wtp = 10000, xlim = c(-0.0005, 0.0015))
ceplane.plot(he, graph = "ggplot2", wtp = 25000, xlim = c(-0.0005, 0.0015))
ceplane.plot(he, graph = "ggplot2", wtp = 50000, xlim = c(-0.0005, 0.0015))

## -----------------------------------------------------------------------------
ceplane.plot(he, graph = "ggplot2", wtp = 1000, xlim = c(-0.0005, 0.0015), label.pos = TRUE)
ceplane.plot(he, graph = "ggplot2", wtp = 1000, xlim = c(-0.0005, 0.0015), label.pos = FALSE)

## -----------------------------------------------------------------------------
ceplane.plot(he, graph = "ggplot2", wtp = list(y = 8))

## -----------------------------------------------------------------------------
ceplane.plot(he)  # default

## -----------------------------------------------------------------------------
ceplane.plot(he, wtp = 1000, xlim = c(-0.0005, 0.0015))
ceplane.plot(he, wtp = 10000, xlim = c(-0.0005, 0.0015))
ceplane.plot(he, wtp = 25000, xlim = c(-0.0005, 0.0015))
ceplane.plot(he, wtp = 50000, xlim = c(-0.0005, 0.0015))

## -----------------------------------------------------------------------------
ceplane.plot(he, wtp = 1000, xlim = c(-0.0005, 0.0015), label.pos = TRUE)
ceplane.plot(he, wtp = 1000, xlim = c(-0.0005, 0.0015), label.pos = FALSE)

## -----------------------------------------------------------------------------
##TODO: not yet implemented ggplot syntax for base R

# ceplane.plot(he, wtp = list(value = 10000))
# ceplane.plot(he, wtp = list(value = 10000, colour = "blue"))
# ceplane.plot(he, wtp = list(colour = "blue"))
# ceplane.plot(he, wtp = list(y = 8))
# ceplane.plot(he, wtp = list(size = 5))
# 
# # to hide text
# ceplane.plot(he, wtp = list(size = 0))

