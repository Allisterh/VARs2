# VARs

This repository is based on the  **macrometrics** toolbox by [Gabriel Züllig](https://gabrielzuellig.ch/macrometrics/).

## Changes:

- The folder `_tbx` is renamed to `toolbox`.

- The use of an R Project supersedes the need to set the working directoy with `setwd()`.

- The scripts have been formatted with the [styler](https://styler.r-lib.org/) package.

- Data is read into R with functions from the `tidyverse`, hence `base::read.csv()` is replaced by `readr::read_csv()`.
  Note that the latter function returns a tibble object instead of a data.frame, hence a single column-vector has to be transposed into a row-fector with either the transpose t() function or with dplyr::pull().