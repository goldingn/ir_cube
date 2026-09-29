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
# # greta. The model needs greta >= 43f9c52 (0.6.0.9000), where a
# # subassignment into a greta array fills it column-major, as R does
# # (greta-dev/greta#844, fixed in #847); the check below tests that behaviour.
# # greta.dynamics is no longer needed: the selection recursion is a closed-form
# # op (R/dynamical_model.R). Tested with greta at 282944f, TensorFlow 2.21.0
# # and TensorFlow Probability 0.25.0 under python 3.12, installed in a separate
# # library and conda environment (see doc/cv_run_plan.md, section 1):
# remotes::install_github("greta-dev/greta@282944f",
#                         lib = "~/R/greta06-lib")
# # conda create -n greta06-env python=3.12
# # <env>/bin/python -m pip install "tensorflow==2.21.*" \
# #   "tensorflow_probability[tf]==0.25.*"
# # and run R with
# # R_LIBS=~/R/greta06-lib \
# #   RETICULATE_PYTHON=~/.local/share/r-miniconda/envs/greta06-env/bin/python

library(tidyverse)
library(readxl)
library(greta)

# greta must fill subassignments column-major. This initialises python, so it
# comes before terra and sf are attached (doc/cv_run_plan.md, section 2)
local({
  g <- greta::zeros(4, 2)
  g[c(1, 3), ] <- matrix(1:4, 2)
  stopifnot(
    "greta fills subassignments row-major; install greta >= 43f9c52 (greta-dev/greta#847)" =
      identical(unname(greta::calculate(g)[[1]]), rbind(c(1, 3), 0, c(2, 4), 0))
  )
})

library(lme4)
library(terra)
library(tidyterra)
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
