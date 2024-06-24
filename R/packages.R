# load packages

# # NOTE: patchwork is bugging out with recent ggplot, use 3.4.4:
# remotes::install_version("ggplot2",
#                          version = "3.4.4",
#                          repos = "http://cran.us.r-project.org")
# # and older tidyterra bc dependency
# remotes::install_version("tidyterra",
#                          version = "0.4.0",
#                          dependencies = FALSE,
#                          repos = "http://cran.us.r-project.org")
# # install specific version of greta.dynamics:
# remotes::install_github("greta-dev/greta.dynamics@greta_2")

library(tidyverse)
library(readxl)
library(greta)
library(lme4)
library(terra)
library(tidyterra)
library(greta.dynamics)
library(tidygeocoder)
library(future)
library(future.apply)
library(DHARMa)
library(patchwork)
