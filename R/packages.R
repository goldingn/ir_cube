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
# # greta and greta.dynamics versions. The greta_2 branch this used to pin no
# # longer exists upstream; it was merged into main. Sampling a model whose
# # likelihood depends on iterate_dynamic_function() fails under greta 0.6.0
# # with a tf.while_loop shape error, so pin the last working pair:
# remotes::install_version("tensorflow", version = "2.16.0")
# remotes::install_github("greta-dev/greta@v0.5.0")
# remotes::install_github("greta-dev/greta.dynamics@v0.2.2")

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
library(future.callr)
library(DHARMa)
library(Hmisc)
library(patchwork)
library(extraDistr)
library(ggtext)
library(geodata)
library(sf)
