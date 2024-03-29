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
# str(he)

ceac.plot(he)

## -----------------------------------------------------------------------------
ceac.plot(he, graph = "base")
ceac.plot(he, graph = "ggplot2")
# ceac.plot(he, graph = "plotly")

## -----------------------------------------------------------------------------
ceac.plot(he,
          graph = "ggplot2",
          title = "my title",
          line = list(color = "green"),
          theme = theme_dark())

## -----------------------------------------------------------------------------
data("Smoking")

he <- bcea(eff, cost, ref = 4)
# str(he)

## -----------------------------------------------------------------------------
ceac.plot(he)

ceac.plot(he,
          graph = "base",
          title = "my title",
          line = list(color = "green"))

## -----------------------------------------------------------------------------
ceac.plot(he,
          graph = "ggplot2",
          title = "my title",
          line = list(color = "green"))

## -----------------------------------------------------------------------------
ceac.plot(he, pos = FALSE) # bottom right
ceac.plot(he, pos = c(0, 0))
ceac.plot(he, pos = c(0, 1))
ceac.plot(he, pos = c(1, 0))
ceac.plot(he, pos = c(1, 1))

## -----------------------------------------------------------------------------
ceac.plot(he, graph = "ggplot2", pos = c(0, 0))
ceac.plot(he, graph = "ggplot2", pos = c(0, 1))
ceac.plot(he, graph = "ggplot2", pos = c(1, 0))
ceac.plot(he, graph = "ggplot2", pos = c(1, 1))

## -----------------------------------------------------------------------------
mypalette <- RColorBrewer::brewer.pal(3, "Accent")

ceac.plot(he,
          graph = "base",
          title = "my title",
          line = list(color = mypalette),
          pos = FALSE)

ceac.plot(he,
          graph = "ggplot2",
          title = "my title",
          line = list(color = mypalette),
          pos = FALSE)

## -----------------------------------------------------------------------------
he <- multi.ce(he)

## -----------------------------------------------------------------------------
ceac.plot(he, graph = "base")

ceac.plot(he,
          graph = "base",
          title = "my title",
          line = list(color = "green"),
          pos = FALSE)

mypalette <- RColorBrewer::brewer.pal(4, "Dark2")

ceac.plot(he,
          graph = "base",
          title = "my title",
          line = list(color = mypalette),
          pos = c(0,1))

## -----------------------------------------------------------------------------
ceac.plot(he,
          graph = "ggplot2",
          title = "my title",
          line = list(color = mypalette),
          pos = c(0,1))

## -----------------------------------------------------------------------------
ceac.plot(he,
          graph = "ggplot2",
          title = "my title",
          line = list(size = 2))

## -----------------------------------------------------------------------------
ceac.plot(he,
          graph = "ggplot2",
          title = "my title",
          line = list(size = c(1,2,3)))

## ----echo=FALSE---------------------------------------------------------------
# create output docs
# rmarkdown::render(input = "vignettes/ceac.Rmd", output_format = "pdf_document", output_dir = "vignettes")
# rmarkdown::render(input = "vignettes/ceac.Rmd", output_format = "html_document", output_dir = "vignettes")

