
<!-- README.md is generated from README.Rmd. Please edit that file -->

# doSurvival

<!-- badges: start -->
<!-- badges: end -->

The goal of doSurvival is to do radiomics-based survival analysis

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("yzhizai/doSurvival")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(doSurvival)
## basic example code
radFile <- file.choose()
cliFile <- file.choose()
dt.radiomics <- read_csv(radFile) %>% select(-1)
dt.clinics <- read_csv(cliFile) %>% select(-1)

mod_rad <- survRadiomics$new(time = 'PFS', 
                             event = 'Progress')
mod_rad <- mod_rad$run(dt.radiomics = dt.radiomics, 
                       dt.clinics = dt.clinics)
mod_rad$figure(outName = 'output.pptx')



mod_nomo <- survNomogram$new(time = 'PFS', 
                             event = 'Progress')
mod_nomo <- mod_nomo$run(dt = dt.clinics, restep = T)
mod_nomo$figure(dt = dt.clinics, outName = 'nomogram.pptx')

out_nomo <- mod_nomo$predict(dt = dt.clinics)
out_nomo$kmplot('km.pptx')
```
