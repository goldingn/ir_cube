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
# # longer exists upstream; it was merged into main, of which 0.2.3 is the
# # release.
# #
# # WARNING: sampling this model currently fails on every version combination
# # tried here, with the likelihood on iterate_dynamic_function() output:
# #   greta 0.6.0 + greta.dynamics 0.2.2 or 0.2.3
# #     ValueError: Input tensor enters the loop with shape (1, ...) but has
# #     shape (None, ...) after one iteration
# #   greta 0.5.0 + greta.dynamics 0.2.2 (tensorflow R 2.16.0)
# #     ValueError: Shape must be rank 3 but is rank 2 for node while/Tile
# # This is not specific to the validation code: the model definition on master
# # fails identically. A working combination has not been identified; whichever
# # versions were installed before the R 4.6.1 upgrade are the ones to restore.
# remotes::install_github("greta-dev/greta.dynamics")

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
