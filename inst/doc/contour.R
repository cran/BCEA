## ---- include = FALSE---------------------------------------------------------
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
he <- bcea(eff, cost, ref = 2)

## -----------------------------------------------------------------------------
contour(he, graph = "base")
contour(he, graph = "ggplot2")
# ceac.plot(he, graph = "plotly")

## -----------------------------------------------------------------------------
contour(he, levels = c(0.2, 0.8))
contour(he, graph = "ggplot2", contour = list(breaks = c(0.2, 0.8)))

## -----------------------------------------------------------------------------
contour(he,
        graph = "ggplot2",
        title = "my title",
        point = list(color = "blue", shape = 2, size = 5),
        contour = list(size = 2))

## -----------------------------------------------------------------------------
contour(he,
        graph = "base",
        title = "my title",
        point = list(color = "blue", shape = 2, size = 2),
        contour = list(size = 2))

## -----------------------------------------------------------------------------
contour2(he, graph = "base")
contour2(he, graph = "ggplot2")
# ceac.plot(he, graph = "plotly")

## -----------------------------------------------------------------------------
contour2(he,
         graph = "ggplot2",
         title = "my title",
         point = list(color = "blue", shape = 10, size = 5),
         contour = list(size = 2))

## -----------------------------------------------------------------------------
contour2(he,
         graph = "base",
         title = "my title",
         point = list(color = "blue", shape = 2, size = 3),
         contour = list(size = 4))

## -----------------------------------------------------------------------------
data("Smoking")
he <- bcea(eff, cost, ref = 4)
# str(he)

## -----------------------------------------------------------------------------
contour(he)
contour(he, graph = "ggplot2")

## -----------------------------------------------------------------------------
contour(he, scale = 0.9)
contour(he, graph = "ggplot2", scale = 0.9)  ##TODO: what is the equivalent ggplot2 argument?

## -----------------------------------------------------------------------------
contour(he, nlevels = 10)
contour(he, graph = "ggplot2", contour = list(bins = 10))

## -----------------------------------------------------------------------------
contour(he, levels = c(0.2, 0.8))
contour(he, graph = "ggplot2", contour = list(breaks = c(0.2, 0.8)))

## -----------------------------------------------------------------------------
contour(he,
        graph = "ggplot2",
        title = "my title",
        line = list(color = "red", size = 1),
        point = list(color = c("plum", "tomato", "springgreen"), shape = 3:5, size = 2),
        icer = list(color = c("red", "orange", "black"), size = 5),
        contour = list(size = 2))

## -----------------------------------------------------------------------------
contour(he,
        graph = "base",
        title = "my title",
        line = list(color = "red", size = 1),
        point = list(color = c("plum", "tomato", "springgreen"), shape = 3:5, size = 2),
        icer = list(color = c("red", "orange", "black"), size = 5),
        contour = list(size = 4))

## -----------------------------------------------------------------------------
contour2(he, wtp = 250)
contour2(he, wtp = 250, graph = "ggplot2")

## -----------------------------------------------------------------------------
contour2(he, wtp = 250,
         graph = "ggplot2",
         title = "my title",
         line = list(color = "red", size = 1),
         point = list(color = c("plum", "tomato", "springgreen"), shape = 3:5, size = 2),
         icer = list(color = c("red", "orange", "black"), size = 5),
         contour = list(size = 2))

## -----------------------------------------------------------------------------
contour2(he, wtp = 250,
         graph = "base",
         title = "my title",
         line = list(color = "red", size = 1),
         point = list(color = c("plum", "tomato", "springgreen"), shape = 3:5, size = 2),
         icer = list(color = c("red", "orange", "black"), size = 5),
         contour = list(size = 4))

## -----------------------------------------------------------------------------
contour(he, pos = FALSE)    # bottom right
contour(he, pos = c(0, 0))
contour(he, pos = c(0, 1))
contour(he, pos = c(1, 0))
contour(he, pos = c(1, 1))

