# Second-stage geostatistical correction to the dynamical model (#21).
#
# The dynamical model captures broad trends, but its residuals carry
# fine-scale spatiotemporal structure. Rather than add a space-time random
# effect inside the dynamical model (infeasible, and confounded with the process
# covariates), this fits a Gaussian model to the residuals on the logit scale,
# conditioning on the dynamical model without ever updating it: a cut posterior.
#
#   z_i = m_i + omega(s_i) + xi(s_i, t_i) + u_j(i) + e_i,   e_i ~ N(0, v_i)
#
# z and v are the empirical logit and its variance, inflated for beta-binomial
# noise by the replicate-based rho (stage A of the issue). m is the dynamical
# model's logit prediction. omega corrects the initial conditions, xi
# accumulates an AR(1) field of annual selection anomalies eta, so that
# forecasts of the correction plateau, and u is the pixel-year deviation from
# the smooth fields.
#
# Optional terms (doc/two_stage_plan.md, "Error structure"), off by default so
# that the defaults reproduce the model above exactly:
#
#   + p_c(i)  static pixel effect, iid N(0, sigma_p^2) per pixel (raster cell):
#             a persistent local deviation that omega otherwise absorbs as a
#             node spike. Part of the prediction target, like u.
#   + s_k(i)  survey effect, iid N(0, sigma_s^2) per survey (citation x country
#             x year, see survey_id()): shared batch / measurement error between
#             the assays of one survey. Not part of the prediction target: it is
#             left out of maps, and added (as a fresh draw per held-out survey)
#             only to predict held-out assays, like the assay noise e.
#
# Hyperparameters are estimated by penalised maximum marginal
# likelihood (TMB, template in tmb/two_stage_correction.cpp); given them, the
# latent posterior is exactly Gaussian with precision equal to the random
# effects Hessian, which is what prediction samples from.
#
# Fitted separately per insecticide. Functions only: source from the repo root.

suppressMessages({
  library(TMB)
  library(Matrix)
  library(fmesher)
  library(sf)
  library(dplyr)
})

correction_template <- "tmb/two_stage_correction.cpp"

# compile the template if its shared object is missing or stale, and load it
load_correction_template <- function(path = correction_template) {
  dll_path <- TMB::dynlib(sub("\\.cpp$", "", path))
  if (!file.exists(dll_path) || file.mtime(dll_path) < file.mtime(path)) {
    TMB::compile(path, flags = "-O2", framework = "TMBad")
  }
  dll_name <- basename(sub("\\.cpp$", "", path))
  if (!dll_name %in% names(getLoadedDLLs())) {
    dyn.load(dll_path)
  }
  dll_name
}


# stage A response ---------------------------------------------------------

# empirical logit of bioassay mortality, and its approximate sampling variance
# inflated by the beta-binomial design effect 1 + (n - 1) rho. rho is the
# replicate-based estimate (R/estimate_bioassay_rho.R), fixed rather than
# estimated, so v carries all of the within-pixel-year assay noise and u can be
# indexed by pixel-year
empirical_logit <- function(died, n, rho) {
  z <- log((died + 0.5) / (n - died + 0.5))
  v <- (1 / (died + 0.5) + 1 / (n - died + 0.5)) * (1 + (n - 1) * rho)
  list(z = z, v = v)
}


# coordinates and mesh ---------------------------------------------------------

# Albers equal-area conic for Africa, in km. Equal-area so that the Matern
# range means the same distance everywhere on the continent, and in km so that
# the range priors can be stated in km
africa_equal_area_crs <- paste(
  "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25",
  "+x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"
)

project_km <- function(lon, lat, crs = africa_equal_area_crs) {
  xy <- sf::sf_project(
    from = "+proj=longlat +datum=WGS84 +no_defs",
    to = crs,
    pts = cbind(lon, lat)
  )
  colnames(xy) <- c("x_km", "y_km")
  xy
}

# projected coordinates of a data frame of observations or prediction points.
# Uses x_km and y_km if present (e.g. simulated data), otherwise projects lon
# and lat
coords_km <- function(df) {
  if (all(c("x_km", "y_km") %in% names(df))) {
    cbind(x_km = df$x_km, y_km = df$y_km)
  } else {
    project_km(df$lon, df$lat)
  }
}

# fmesher mesh, fine near the data and coarse elsewhere.
#
# The issue asks for edges of at most about range / 4 near the data. A uniform
# fine mesh over the data's extent cannot meet the node budget: the bioassay
# sites span most of sub-Saharan Africa, so ~2500 nodes spread evenly over a
# non-convex hull of the sites gives ~150 km edges everywhere (checked on the
# Deltamethrin and Fenitrothion data). Instead the mesh has a node at every
# site, with sites closer than `cutoff` merged, and triangles of at most
# max_edge_inner km between them. Resolution then follows data density: where
# sites are dense the edges are the (>= cutoff) site spacing, where they are
# sparse they are up to max_edge_inner, which matters little because there is
# no information there to resolve. The inner boundary is a non-convex hull with
# a buffer, and an outer extension with coarse triangles keeps the SPDE's
# boundary effects away from the data.
#
# If the mesh exceeds max_nodes, the cutoff is increased (coarsening the
# densest clusters first) until it doesn't, and the final value is reported.
# With the defaults, Deltamethrin (the most sampled insecticide, ~4000 distinct
# sites) gives ~2400 nodes with 90% of sites within ~20 km of a node
build_correction_mesh <- function(coords_km,
                                  max_edge_inner = 250,
                                  max_edge_outer = 800,
                                  cutoff = 30,
                                  inner_buffer = 200,
                                  outer_buffer = 1500,
                                  hull_resolution = c(60, 60),
                                  max_nodes = 2500,
                                  cutoff_growth = 1.2,
                                  verbose = TRUE) {

  coords_km <- unique(as.matrix(coords_km))
  # a smooth hull (coarse resolution) avoids many short boundary segments,
  # each of which would add nodes
  inner <- suppressWarnings(
    fmesher::fm_nonconvex_hull(coords_km,
                               convex = inner_buffer,
                               resolution = hull_resolution,
                               format = "fm")
  )

  repeat {
    mesh <- fmesher::fm_mesh_2d(
      loc = coords_km,
      boundary = list(inner),
      max.edge = c(max_edge_inner, max_edge_outer),
      cutoff = cutoff,
      offset = c(-0.01, outer_buffer)
    )
    if (mesh$n <= max_nodes) break
    cutoff <- cutoff * cutoff_growth
  }

  if (verbose) {
    message(sprintf("correction mesh: %i nodes, cutoff %.0f km",
                    mesh$n, cutoff))
  }
  attr(mesh, "cutoff") <- cutoff
  mesh
}

# FEM matrices in the form R_inla::spde_t expects: c0 is the lumped (diagonal)
# mass matrix, as in INLA, so that Q is sparse
correction_fem <- function(mesh) {
  fem <- fmesher::fm_fem(mesh, order = 2)
  as_dgc <- function(x) as(as(as(x, "dMatrix"), "generalMatrix"),
                           "CsparseMatrix")
  list(M0 = as_dgc(fem$c0),
       M1 = as_dgc(fem$g1),
       M2 = as_dgc(fem$g2))
}

# Matern precision at the nodes with marginal SD sigma (alpha = 2, d = 2),
# matching matern_precision() in the template
matern_precision_r <- function(fem, kappa, sigma) {
  tau_spde <- 1 / (sigma * kappa * sqrt(4 * pi))
  tau_spde ^ 2 * (kappa ^ 4 * fem$M0 + 2 * kappa ^ 2 * fem$M1 + fem$M2)
}


# priors ------------------------------------------------------------------------

# PC priors: P(range < 50 km) = 0.05 and P(sigma > 1) = 0.05 for both Matern
# fields, P(tau > 1) = 0.05, and persistence 1 / (1 - phi) ~ lognormal(log 5,
# 0.62), i.e. a median of five years with a 95% interval of roughly 1.5-17
# P(sigma_p > 1) = 0.05 and P(sigma_s > 1) = 0.05 for the optional pixel and
# survey effects
correction_priors <- function(range0 = 50,
                              alpha_range = 0.05,
                              sigma0 = 1,
                              alpha_sigma = 0.05,
                              tau0 = 1,
                              alpha_tau = 0.05,
                              persistence_median = 5,
                              persistence_sdlog = 0.62,
                              sigma_p0 = 1,
                              alpha_sigma_p = 0.05,
                              sigma_s0 = 1,
                              alpha_sigma_s = 0.05) {
  list(
    pc_omega = c(range0, alpha_range, sigma0, alpha_sigma),
    pc_eta = c(range0, alpha_range, sigma0, alpha_sigma),
    pc_tau = c(tau0, alpha_tau),
    persistence_prior = c(log(persistence_median), persistence_sdlog),
    pc_sigma_p = c(sigma_p0, alpha_sigma_p),
    pc_sigma_s = c(sigma_s0, alpha_sigma_s)
  )
}


# surveys -----------------------------------------------------------------------

# Survey identifier for the survey effect s: citation x country x year. The
# citation alone does not identify a study: aggregated sources ("Ministry of
# Health", "PMI 2016", "VectorBase", personal communications) span countries,
# so the country splits them (110 of 1434 citation-years span more than one
# country, with 9,902 assays). Assays with no citation (12 in the cleaned data,
# all Madagascar 2023) fall back to spatial clusters within a year: pixels
# within cluster_km of each other (single linkage) in the same country-year.
survey_id <- function(citation, country, year, coords_km, cluster_km = 25) {
  id <- paste(citation, country, year, sep = " | ")
  missing <- is.na(citation) | citation == ""
  if (any(missing)) {
    group <- paste(country, year)[missing]
    xy <- coords_km[missing, , drop = FALSE]
    cluster <- integer(sum(missing))
    for (g in unique(group)) {
      rows <- which(group == g)
      cluster[rows] <- if (length(rows) == 1) 1L else
        stats::cutree(stats::hclust(stats::dist(xy[rows, , drop = FALSE]),
                                    method = "single"), h = cluster_km)
    }
    id[missing] <- paste("no citation", group, "cluster", cluster, sep = " | ")
  }
  id
}


# fitting -----------------------------------------------------------------------

# sparse projection from points to mesh nodes
mesh_basis <- function(mesh, coords) {
  A <- fmesher::fm_basis(mesh, loc = coords)
  as(as(A, "generalMatrix"), "CsparseMatrix")
}

# observation -> latent design matrices. xi(., t0) = 0 and the columns of x are
# years t0 + 1, ..., T, so an observation in year t touches column t - t0 and
# observations in year t0 touch none
correction_design <- function(mesh, mesh_xi, coords, year, t0, T) {
  A_omega <- mesh_basis(mesh, coords)
  A_xi_space <- mesh_basis(mesh_xi, coords)
  n_nodes <- mesh_xi$n
  n_years <- T - t0
  year_col <- year - t0
  trip <- summary(A_xi_space)
  keep <- year_col[trip$i] > 0 & year_col[trip$i] <= n_years
  A_xi <- Matrix::sparseMatrix(
    i = trip$i[keep],
    j = trip$j[keep] + (year_col[trip$i[keep]] - 1) * n_nodes,
    x = trip$x[keep],
    dims = c(length(year), n_nodes * n_years)
  )
  list(A_omega = A_omega, A_xi = A_xi)
}

# assemble the TMB object. Kept separate from fit_correction() so that the
# checks can rebuild it with a perturbed m and fixed hyperparameters
correction_adfun <- function(data, parameters, variant, fix_hyper = FALSE,
                             silent = TRUE) {
  dll <- load_correction_template()
  hyper_names <- c("log_sigma_omega", "log_kappa_omega", "log_sigma_eta",
                   "log_kappa_eta", "logit_phi", "log_tau")
  map <- list()
  random <- c("w_omega", "u")
  # the optional pixel and survey effects: mapped to zero (and dropped from the
  # parameter vector, so the model is exactly the one without them) when off
  for (term in c("p", "s")) {
    if (data[[paste0("include_", term)]] == 1) {
      random <- c(random, term)
      hyper_names <- c(hyper_names, paste0("log_sigma_", term))
    } else {
      map[[term]] <- factor(rep(NA, length(parameters[[term]])))
      map[[paste0("log_sigma_", term)]] <- factor(NA)
    }
  }
  if (variant == "omega_u") {
    # xi and its hyperparameters drop out of the model entirely
    map$x <- factor(rep(NA, length(parameters$x)))
    map$log_sigma_eta <- factor(NA)
    map$log_kappa_eta <- factor(NA)
    map$logit_phi <- factor(NA)
  } else {
    random <- c(random, "x")
  }
  if (fix_hyper) {
    for (name in hyper_names) map[[name]] <- factor(NA)
  }
  TMB::MakeADFun(data = data,
                 parameters = parameters,
                 map = map,
                 random = random,
                 DLL = dll,
                 silent = silent)
}

# Fit the correction model for one insecticide.
#
# train needs lon and lat (or x_km and y_km), year, cell, m (the dynamical
# model's posterior mean logit prediction), and either z and v, or died,
# mosquito_number and rho. t0 is the dynamical model's start year and T the
# last year for which xi is represented (default: the last data year); xi is
# forecast beyond T at prediction time.
#
# mesh is used for omega and mesh_xi for xi; either is built from the training
# locations if not supplied. mesh_xi defaults to a coarser mesh (at most 600
# nodes) because the cost of the fit is dominated by the Cholesky factor of the
# space-time block, whose fill-in grows faster than linearly in
# nodes x years. On the simulation check (20 years), xi on the 1441-node omega
# mesh gave a factor with 29M non-zeros, a 6 minute fit and 3.3 GB peak memory,
# against 2 minutes and 2.3 GB on a 587-node xi mesh; the real data (up to 2500
# nodes, ~30 years) would be several times worse again. Pass mesh_xi = mesh to
# use one mesh for both.
#
# Speed depends heavily on the BLAS used by CHOLMOD's supernodal factorisation:
# with R's reference BLAS the same fits are ~10x slower than with OpenBLAS.
fit_correction <- function(train,
                           variant = c("omega_xi_u", "omega_u"),
                           t0,
                           T = max(train$year),
                           mesh = NULL,
                           mesh_args = list(),
                           mesh_xi = NULL,
                           mesh_xi_args = list(max_nodes = 600),
                           priors = correction_priors(),
                           start = list(),
                           control = list(eval.max = 1000, iter.max = 500),
                           sdreport = FALSE,
                           silent = TRUE,
                           pixel_effect = FALSE,
                           survey_effect = FALSE) {

  variant <- match.arg(variant)
  if (survey_effect && !"survey" %in% names(train)) {
    stop("survey_effect = TRUE needs a survey column in train (see survey_id())")
  }
  timing <- list(start = Sys.time())

  if (any(train$year < t0) || any(train$year > T)) {
    stop("training years must lie in [t0, T]")
  }
  if (!all(c("z", "v") %in% names(train))) {
    stage_a <- empirical_logit(train$died, train$mosquito_number, train$rho)
    train$z <- stage_a$z
    train$v <- stage_a$v
  }

  coords <- coords_km(train)
  if (is.null(mesh)) {
    mesh <- do.call(build_correction_mesh, c(list(coords), mesh_args))
  }
  if (is.null(mesh_xi)) {
    mesh_xi <- if (variant == "omega_xi_u") {
      do.call(build_correction_mesh, c(list(coords), mesh_xi_args))
    } else {
      mesh
    }
  }
  fem <- correction_fem(mesh)
  fem_xi <- correction_fem(mesh_xi)
  n_nodes <- mesh$n
  # keep at least one column so the parameter exists for the omega_u variant
  n_years <- max(T - t0, 1)

  design <- correction_design(mesh, mesh_xi, coords, train$year, t0, T)

  # pixel-years with data, each with its own u
  pixel_years <- train %>%
    distinct(cell, year) %>%
    mutate(u_index = row_number())
  u_index <- pixel_years$u_index[match(paste(train$cell, train$year),
                                       paste(pixel_years$cell,
                                             pixel_years$year))]

  # static pixels and surveys with data, each with its own p and s (a single
  # unused level when the term is off)
  pixels <- tibble(cell = unique(train$cell)) %>% mutate(p_index = row_number())
  p_index <- match(train$cell, pixels$cell)
  surveys <- if (survey_effect) {
    tibble(survey = unique(train$survey)) %>% mutate(s_index = row_number())
  } else {
    tibble(survey = NA_character_, s_index = 1L)
  }
  s_index <- if (survey_effect) match(train$survey, surveys$survey) else
    rep(1L, nrow(train))

  data <- list(
    z = train$z,
    v = train$v,
    m = train$m,
    A_omega = design$A_omega,
    A_xi = design$A_xi,
    u_index = as.integer(u_index - 1),
    spde = fem,
    spde_xi = fem_xi,
    include_xi = as.integer(variant == "omega_xi_u"),
    pc_omega = priors$pc_omega,
    pc_eta = priors$pc_eta,
    pc_tau = priors$pc_tau,
    persistence_prior = priors$persistence_prior,
    include_p = as.integer(pixel_effect),
    p_index = as.integer(p_index - 1),
    pc_sigma_p = if (is.null(priors$pc_sigma_p)) c(1, 0.05) else
      priors$pc_sigma_p,
    include_s = as.integer(survey_effect),
    s_index = as.integer(s_index - 1),
    pc_sigma_s = if (is.null(priors$pc_sigma_s)) c(1, 0.05) else
      priors$pc_sigma_s
  )

  # start the ranges at a few hundred km, well inside the prior, and the
  # SDs small: the correction should be modest if the dynamical model is good
  defaults <- list(
    log_sigma_omega = log(0.5),
    log_kappa_omega = log(sqrt(8) / 300),
    log_sigma_eta = log(0.2),
    log_kappa_eta = log(sqrt(8) / 300),
    logit_phi = qlogis(0.8),
    log_tau = log(0.3),
    log_sigma_p = log(0.3),
    log_sigma_s = log(0.3)
  )
  start <- modifyList(defaults, start)
  parameters <- c(
    list(w_omega = rep(0, n_nodes),
         x = matrix(0, mesh_xi$n, n_years),
         u = rep(0, nrow(pixel_years))),
    start[c("log_sigma_omega", "log_kappa_omega", "log_sigma_eta",
            "log_kappa_eta", "logit_phi", "log_tau")],
    list(p = rep(0, if (pixel_effect) nrow(pixels) else 1),
         s = rep(0, if (survey_effect) nrow(surveys) else 1)),
    start[c("log_sigma_p", "log_sigma_s")]
  )

  obj <- correction_adfun(data, parameters, variant, silent = silent)
  timing$setup <- Sys.time()

  opt <- nlminb(obj$par, obj$fn, obj$gr, control = control)
  if (opt$convergence != 0) {
    warning("nlminb did not converge: ", opt$message)
  }
  timing$optimise <- Sys.time()

  # latent mode and its Hessian at the optimum. Evaluating fn at opt$par
  # leaves the inner problem solved at exactly these hyperparameters
  obj$fn(opt$par)
  par_full <- obj$env$last.par
  random <- obj$env$random
  mode <- par_full[random]
  block <- names(par_full)[random]
  blocks <- split(seq_along(random), factor(block, levels = unique(block)))

  H <- obj$env$spHess(par_full, random = TRUE)
  H <- Matrix::forceSymmetric(as(H, "CsparseMatrix"), uplo = "L")
  H_chol <- Matrix::Cholesky(H, perm = TRUE, LDL = FALSE, super = TRUE)
  timing$hessian <- Sys.time()

  # full observation -> latent design in the latent (Hessian) order, used for
  # the cut-posterior shift of the mode
  A_u <- Matrix::sparseMatrix(i = seq_len(nrow(train)), j = u_index, x = 1,
                              dims = c(nrow(train), nrow(pixel_years)))
  A_blocks <- list(w_omega = design$A_omega, x = design$A_xi, u = A_u)
  if (pixel_effect) {
    A_blocks$p <- Matrix::sparseMatrix(i = seq_len(nrow(train)), j = p_index,
                                       x = 1, dims = c(nrow(train),
                                                       nrow(pixels)))
  }
  if (survey_effect) {
    A_blocks$s <- Matrix::sparseMatrix(i = seq_len(nrow(train)), j = s_index,
                                       x = 1, dims = c(nrow(train),
                                                       nrow(surveys)))
  }
  A_latent <- do.call(cbind, A_blocks[names(blocks)])

  report <- obj$report(par_full)
  hyper <- list(
    sigma_omega = report$sigma_omega,
    range_omega = report$range_omega,
    kappa_omega = sqrt(8) / report$range_omega,
    tau = report$tau
  )
  if (pixel_effect) hyper$sigma_p <- report$sigma_p
  if (survey_effect) hyper$sigma_s <- report$sigma_s
  if (variant == "omega_xi_u") {
    hyper <- c(hyper, list(
      sigma_eta = report$sigma_eta,
      range_eta = report$range_eta,
      kappa_eta = sqrt(8) / report$range_eta,
      phi = report$phi,
      persistence = report$persistence
    ))
  }

  sd_report <- if (sdreport) TMB::sdreport(obj) else NULL

  timing$end <- Sys.time()
  timings <- c(
    setup = difftime(timing$setup, timing$start, units = "secs"),
    optimise = difftime(timing$optimise, timing$setup, units = "secs"),
    hessian = difftime(timing$hessian, timing$optimise, units = "secs"),
    total = difftime(timing$end, timing$start, units = "secs")
  )

  structure(
    list(
      variant = variant,
      t0 = t0,
      T = T,
      n_years = n_years,
      mesh = mesh,
      fem = fem,
      mesh_xi = mesh_xi,
      fem_xi = fem_xi,
      obj = obj,
      opt = opt,
      par_full = par_full,
      par_list = obj$env$parList(par = par_full),
      tmb_data = data,
      hyper = hyper,
      sd_report = sd_report,
      mode = mode,
      blocks = blocks,
      H = H,
      H_chol = H_chol,
      A_latent = A_latent,
      precision_obs = 1 / train$v,
      m_ref = train$m,
      pixel_years = pixel_years,
      pixel_effect = pixel_effect,
      survey_effect = survey_effect,
      pixels = pixels,
      surveys = surveys,
      n_obs = nrow(train),
      timings = as.numeric(timings, units = "secs") |> setNames(names(timings))
    ),
    class = "correction_fit"
  )
}


# prediction --------------------------------------------------------------------

# Cut-posterior shift of the latent mode when the offset changes from m_ref to
# m_new at the training observations. Given the hyperparameters the mode solves
# H theta = A' D (z - m), so it moves by H^-1 A' D (m_ref - m_new). m_new may
# be a matrix with one column per dynamical draw; the single Cholesky factor of
# H is reused for all of them
correction_mode_shift <- function(fit, m_new) {
  m_new <- as.matrix(m_new)
  rhs <- Matrix::crossprod(fit$A_latent,
                           fit$precision_obs * (fit$m_ref - m_new))
  as.matrix(Matrix::solve(fit$H_chol, rhs, system = "A"))
}

# n draws from N(0, H^-1) using the permuted Cholesky factor P H P' = L L':
# P' L^-T e has covariance P' L^-T L^-1 P = H^-1
sample_latent_deviation <- function(H_chol, n_latent, n) {
  e <- matrix(rnorm(n_latent * n), n_latent, n)
  deviation <- Matrix::solve(H_chol, e, system = "Lt")
  as.matrix(Matrix::solve(H_chol, deviation, system = "Pt"))
}

# Joint predictive draws of lambda = m + omega + xi + u (logit scale, excluding
# the assay noise e) at the rows of new.
#
# new needs lon and lat (or x_km and y_km), year, cell and m (the dynamical
# model's posterior mean at those pixel-years). m_draws_train (K x n_obs) and
# m_draws_new (K x n_new) are the matching dynamical posterior draws; with them,
# draw d uses dynamical draw ((d - 1) mod K) + 1, shifts the latent mode by the
# cut-posterior formula, and adds that draw's m at the new points. Without them
# m is treated as known (no mechanistic uncertainty).
#
# Latent draws are joint samples with precision H, so omega, xi and u keep
# their posterior correlations. xi beyond T is forecast by running the AR(1)
# for eta forward from each draw's eta_T = x_T - x_{T-1} with fresh Matern
# innovations, and accumulating. u is taken from the joint draw at pixel-years
# with data, and drawn from N(0, tau^2) (once per distinct pixel-year) elsewhere.
# The static pixel effect p (if fitted) is handled the same way per pixel: the
# joint draw at pixels with data, N(0, sigma_p^2) once per distinct new pixel
# elsewhere.
#
# The survey effect s (if fitted) is batch error, not part of the target, so
# `survey` says what to do with it:
#   "none"       leave it out: the prediction of the population fraction, what
#                a map shows (the default);
#   "fresh"      add a fresh N(0, sigma_s^2) draw per distinct new$survey,
#                shared by the rows of that survey, whether or not the survey
#                has training data: the predictive distribution of a new assay
#                from an unknown survey. This is what the held-out scoring uses:
#                the scoring code turns logit draws into draws of p and scores
#                them with a beta-binomial at the external per-type rho, so the
#                survey draw has to be inside the draws to widen the predictive
#                distribution the way the assay noise does;
#   "posterior"  as "fresh", but surveys with training data take their joint
#                posterior draw (a sensitivity check: it uses the survey's
#                other assays, which a map user does not have).
# Several values give a named list of matrices from the same latent draws, so
# they differ only in the survey term.
#
# Returns an n_draws x nrow(new) matrix (or a list of them). Draws are generated
# in batches of batch_size to bound memory
predict_correction <- function(fit,
                               new,
                               m_draws_train = NULL,
                               m_draws_new = NULL,
                               n_draws = 1000,
                               batch_size = 100,
                               survey = "none") {

  survey <- unique(survey)
  stopifnot(length(survey) >= 1,
            all(survey %in% c("none", "fresh", "posterior")))
  include_s <- isTRUE(fit$survey_effect)
  include_p <- isTRUE(fit$pixel_effect)
  if (include_s && any(survey != "none") && !"survey" %in% names(new)) {
    stop("survey draws need a survey column in new")
  }

  if (xor(is.null(m_draws_train), is.null(m_draws_new))) {
    stop("supply both m_draws_train and m_draws_new, or neither")
  }
  use_m_draws <- !is.null(m_draws_train)
  if (use_m_draws) {
    stopifnot(ncol(m_draws_train) == fit$n_obs,
              ncol(m_draws_new) == nrow(new),
              nrow(m_draws_train) == nrow(m_draws_new))
    n_m_draws <- nrow(m_draws_train)
  }

  n_new <- nrow(new)
  n_nodes_xi <- fit$mesh_xi$n
  n_latent <- length(fit$mode)
  hyper <- fit$hyper
  include_xi <- fit$variant == "omega_xi_u"

  coords <- coords_km(new)
  A_new <- mesh_basis(fit$mesh, coords)
  outside <- Matrix::rowSums(A_new) < 0.5
  if (any(outside)) {
    warning(sum(outside), " prediction points lie outside the mesh; ",
            "their spatial fields are set to zero")
  }

  # xi design: years within (t0, T] are read off x, later years are forecast,
  # years at or before t0 have xi = 0
  in_range <- new$year > fit$t0 & new$year <= fit$T
  forecast <- new$year > fit$T
  if (include_xi) {
    A_xi_new <- correction_design(fit$mesh, fit$mesh_xi, coords,
                                  ifelse(in_range, new$year, fit$t0),
                                  fit$t0, fit$T)$A_xi
    A_xi_space_new <- mesh_basis(fit$mesh_xi, coords)
    horizon <- new$year[forecast] - fit$T
    max_horizon <- max(c(0, horizon))
    if (max_horizon > 0) {
      Q_eta <- matern_precision_r(fit$fem_xi, hyper$kappa_eta,
                                  hyper$sigma_eta)
      Q_eta_chol <- Matrix::Cholesky(Matrix::forceSymmetric(Q_eta),
                                     perm = TRUE, LDL = FALSE,
                                     super = TRUE)
    }
  }

  # u: pixel-years with data take the latent draw, others a fresh iid draw
  new_key <- paste(new$cell, new$year)
  u_train <- fit$pixel_years$u_index[match(new_key,
                                           paste(fit$pixel_years$cell,
                                                 fit$pixel_years$year))]
  has_u <- !is.na(u_train)
  unseen_keys <- unique(new_key[!has_u])
  unseen_index <- match(new_key[!has_u], unseen_keys)

  idx_w <- fit$blocks$w_omega
  idx_x <- fit$blocks$x
  idx_u <- fit$blocks$u

  # p: pixels with data take the latent draw, others a fresh iid draw
  if (include_p) {
    p_train <- fit$pixels$p_index[match(new$cell, fit$pixels$cell)]
    has_p <- !is.na(p_train)
    unseen_cells <- unique(new$cell[!has_p])
    unseen_cell_index <- match(new$cell[!has_p], unseen_cells)
    idx_p <- fit$blocks$p
  }
  # s: one fresh draw per distinct new survey; for "posterior", surveys with
  # training data take the latent draw
  if (include_s && any(survey != "none")) {
    new_surveys <- unique(new$survey)
    new_survey_index <- match(new$survey, new_surveys)
    s_train <- fit$surveys$s_index[match(new_surveys, fit$surveys$survey)]
    idx_s <- fit$blocks$s
  }

  draws <- lapply(setNames(survey, survey), function(x)
    matrix(NA_real_, n_draws, n_new))
  batches <- split(seq_len(n_draws), ceiling(seq_len(n_draws) / batch_size))

  for (batch in batches) {
    nb <- length(batch)

    # latent draws: mode (shifted per dynamical draw) plus a joint deviation
    theta <- sample_latent_deviation(fit$H_chol, n_latent, nb) + fit$mode
    if (use_m_draws) {
      k <- ((batch - 1) %% n_m_draws) + 1
      theta <- theta + correction_mode_shift(fit, t(m_draws_train[k, ,
                                                                   drop = FALSE]))
      m_new <- t(m_draws_new[k, , drop = FALSE])
    } else {
      m_new <- matrix(new$m, n_new, nb)
    }

    lambda <- m_new + as.matrix(A_new %*% theta[idx_w, , drop = FALSE])

    if (include_xi) {
      x_draw <- theta[idx_x, , drop = FALSE]
      lambda[in_range, ] <- lambda[in_range, ] +
        as.matrix(A_xi_new[in_range, , drop = FALSE] %*% x_draw)

      if (max_horizon > 0) {
        # node values of x_T and x_{T-1} (x_t0 = 0) in each draw
        last <- (fit$n_years - 1) * n_nodes_xi + seq_len(n_nodes_xi)
        xi_nodes <- x_draw[last, , drop = FALSE]
        eta_nodes <- if (fit$n_years > 1) {
          xi_nodes - x_draw[last - n_nodes_xi, , drop = FALSE]
        } else {
          xi_nodes
        }
        for (h in seq_len(max_horizon)) {
          innovation <- sample_latent_deviation(Q_eta_chol, n_nodes_xi, nb)
          eta_nodes <- hyper$phi * eta_nodes +
            sqrt(1 - hyper$phi ^ 2) * innovation
          xi_nodes <- xi_nodes + eta_nodes
          rows <- which(forecast)[horizon == h]
          if (length(rows) > 0) {
            lambda[rows, ] <- lambda[rows, ] +
              as.matrix(A_xi_space_new[rows, , drop = FALSE] %*% xi_nodes)
          }
        }
      }
    }

    if (any(has_u)) {
      lambda[has_u, ] <- lambda[has_u, ] +
        theta[idx_u[u_train[has_u]], , drop = FALSE]
    }
    if (any(!has_u)) {
      u_fresh <- matrix(rnorm(length(unseen_keys) * nb, 0, hyper$tau),
                        length(unseen_keys), nb)
      lambda[!has_u, ] <- lambda[!has_u, ] +
        u_fresh[unseen_index, , drop = FALSE]
    }

    if (include_p) {
      if (any(has_p)) {
        lambda[has_p, ] <- lambda[has_p, ] +
          theta[idx_p[p_train[has_p]], , drop = FALSE]
      }
      if (any(!has_p)) {
        p_fresh <- matrix(rnorm(length(unseen_cells) * nb, 0, hyper$sigma_p),
                          length(unseen_cells), nb)
        lambda[!has_p, ] <- lambda[!has_p, ] +
          p_fresh[unseen_cell_index, , drop = FALSE]
      }
    }

    for (mode in survey) {
      lambda_mode <- lambda
      if (include_s && mode != "none") {
        s_draw <- matrix(rnorm(length(new_surveys) * nb, 0, hyper$sigma_s),
                         length(new_surveys), nb)
        if (mode == "posterior" && any(!is.na(s_train))) {
          seen <- !is.na(s_train)
          s_draw[seen, ] <- theta[idx_s[s_train[seen]], , drop = FALSE]
        }
        lambda_mode <- lambda_mode + s_draw[new_survey_index, , drop = FALSE]
      }
      draws[[mode]][batch, ] <- t(lambda_mode)
    }
  }

  if (length(survey) == 1) draws[[1]] else draws
}

# The named mesh configurations of the mesh-resolution experiment
# (doc/two_stage_plan.md, "Mesh resolution results"), as arguments to
# build_correction_mesh() for omega and xi; xi = "omega" puts xi on the omega
# mesh. Mirrors mesh_configs in R/run_two_stage_folds.R. Cross-validation chose
# omega5000_xi2500: omega on a finer mesh (15 km cutoff, inner edges of at most
# 150 km, at most 5000 nodes) and xi on the base omega mesh
correction_mesh_configs <- function() {
  list(
    base = list(omega = list(), xi = list(max_nodes = 600)),
    xi1200 = list(omega = list(), xi = list(max_nodes = 1200)),
    xifull = list(omega = list(), xi = "omega"),
    omega5000 = list(omega = list(cutoff = 15, max_edge_inner = 150,
                                  max_nodes = 5000),
                     xi = list(max_nodes = 1200)),
    omega5000_xi2500 = list(omega = list(cutoff = 15, max_edge_inner = 150,
                                         max_nodes = 5000),
                            xi = list(max_nodes = 2500))
  )
}

# the omega and xi meshes of a named configuration, for coordinates in km
build_correction_meshes <- function(coords_km, config = "omega5000_xi2500",
                                    verbose = FALSE) {
  configs <- correction_mesh_configs()
  if (!config %in% names(configs)) {
    stop("unknown mesh configuration: ", config, "; one of ",
         paste(names(configs), collapse = ", "))
  }
  spec <- configs[[config]]
  mesh <- do.call(build_correction_mesh,
                  c(list(coords_km), spec$omega, verbose = verbose))
  mesh_xi <- if (identical(spec$xi, "omega")) mesh else
    do.call(build_correction_mesh, c(list(coords_km), spec$xi,
                                     verbose = verbose))
  list(omega = mesh, xi = mesh_xi)
}
