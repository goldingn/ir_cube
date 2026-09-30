# The dynamical model, defined once.
#
# build_dynamical_model() builds the greta model from the data, the covariates
# and a training set. fit_model.R calls it with the full data, fit_fold()
# (R/fit_validation_fold.R) with a cross-validation training fold; the only
# difference between the two is which records enter the likelihood.
#
# The priors and the deterministic transforms from the sampled parameters to the
# selection coefficients and the initial states are in dynamical_variables() and
# dynamical_terms(). dynamical_terms() is written to work both on greta arrays
# and on one posterior draw as plain R arrays, so that the plain-R predictions
# in R/dynamical_predictions.R and R/two_stage_map_functions.R use the same
# definitions. A new model term goes in those two functions (and, if it changes
# the recursion, in the state computation), switched on through
# dynamical_model_options().
#
# Source from the repo root, after R/packages.R and R/functions.R.

source("R/greta_setup.R")
source("R/model_covariates.R")
source("R/windowed_hmc.R")


# model options ------------------------------------------------------------

# Switches for model terms. The defaults are the model for the refit: rho per
# type, the mortality floor and both initial-state covariates on; in the
# selection design, population by the encounter transform of density
# (d_half = 50 per km2, not yet confirmed) times g_dom, and the crops times
# g_ag, both linear, 0 in 1995 and 1 in 2025, nets as MITN run 06 use times
# the pyrethroid-only share (w = 0.25), and no hinges (selection_design());
# and reversion estimated, one rate per class (#24). The fits before these terms were added had rho = "class", mortality_floor = FALSE,
# init_covariates = NULL (see fold_options()). build_dynamical_model() refuses
# settings that are not implemented:
#   rho               "class": one overdispersion per insecticide class;
#                     "type": one per insecticide type, types nested in
#                     class (#20)
#   mortality_floor   TRUE for an estimated floor on bioassay mortality, the
#                     mortality of a fully resistant population (#14)
#   init_covariates   names of static covariates of the initial state, from
#                     init_covariate_names(selection_columns)
#                     (R/model_covariates.R), or NULL
#                     for none (#19)
#   selection_columns how the selection design matrix is built:
#                     selection_design() (R/model_covariates.R), the
#                     population transform, any hinge columns and the trend
#                     on the population and crop columns (#23). The
#                     model takes the matrix as built; build_dynamical_model()
#                     checks its columns against this
#   reversion         reversion to susceptibility (#24): FALSE for none,
#                     "estimated" for one rate per class, or fixed values of
#                     kappa (<= 0, one, or one per class) for sensitivity
#                     analysis. See reversion_kappa().
#   init_centred      how the hierarchy on the initial state is
#                     parameterised (#25): "country", the countries' logit
#                     relative initial states sampled directly (centred), the
#                     regions' deviations non-centred; "none", N(0, 1)
#                     deviations scaled by the sds, as in the fits before.
#                     The same model; see dynamical_variables()
dynamical_model_options <- function(rho = c("type", "class"),
                                    mortality_floor = TRUE,
                                    init_covariates =
                                      init_covariate_names(selection_columns),
                                    selection_columns = selection_design(),
                                    reversion = "estimated",
                                    init_centred = c("country", "none")) {
  list(rho = match.arg(rho),
       mortality_floor = mortality_floor,
       init_covariates = init_covariates,
       selection_columns = selection_columns,
       reversion = reversion,
       init_centred = match.arg(init_centred))
}

# Options saved before init_centred was added were non-centred.
complete_dynamical_model_options <- function(options) {
  if (is.null(options$init_centred)) {
    options$init_centred <- "none"
  }
  options
}

check_dynamical_model_options <- function(options) {
  defaults <- dynamical_model_options()
  # init_covariate_centre is set by build_dynamical_model()
  stopifnot(setequal(setdiff(names(options), "init_covariate_centre"),
                     names(defaults)),
            is.null(options$init_covariate_centre) ||
              identical(colnames(options$init_covariate_centre),
                        options$init_covariates))
  options$init_covariate_centre <- NULL
  implemented <- list(rho = list("class", "type"),
                      mortality_floor = list(FALSE, TRUE),
                      init_covariates = list(NULL),
                      reversion = list(FALSE),
                      init_centred = list("none", "country"))
  reversion <- options$reversion
  if (identical(reversion, "estimated") ||
      (is.numeric(reversion) && length(reversion) > 0 &&
         all(is.finite(reversion)) && all(reversion <= 0))) {
    implemented$reversion <- list(FALSE, reversion)
  }
  # any design selection_design() builds, including one saved before the
  # trend was added (complete_selection_design())
  design <- options$selection_columns
  if (identical(complete_selection_design(design),
                do.call(selection_design, complete_selection_design(design)))) {
    implemented$selection_columns <- list(design)
  }
  if (!is.null(options$init_covariates)) {
    stopifnot(is.character(options$init_covariates),
              !anyDuplicated(options$init_covariates),
              all(options$init_covariates %in%
                    init_covariate_names(options$selection_columns)))
    implemented$init_covariates <- list(NULL, options$init_covariates)
  }
  for (name in names(defaults)) {
    allowed <- if (name %in% names(implemented)) implemented[[name]] else
      list(defaults[[name]])
    if (!any(vapply(allowed, identical, logical(1), options[[name]]))) {
      stop("dynamical model option '", name, "' is not implemented yet")
    }
  }
  invisible(options)
}


# fixed quantities ---------------------------------------------------------

# Prior and minimum values for the initial fraction susceptible, per type. The
# initial state is modelled on the logit of its relative position between
# init_frac_min and 1. More flexibility for DDT, less for the others.
init_frac_constants <- function(types) {
  prior <- ifelse(types == "DDT", 0.9, 0.95)
  min <- ifelse(types == "DDT", 0.75, 0.9)
  list(prior = prior,
       min = min,
       # mean logit proportion of the distance from the minimum to 1
       relative_prior = (prior - min) / (1 - min))
}

# Lookups from the full data (never a fold, so that every fold and the full fit
# index the same countries and cells): the region of each country, and the
# country whose initial state each cell takes, i.e. the country of the cell's
# first record.
dynamical_lookups <- function(df) {
  country_region_index <- df %>%
    group_by(country_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(country_id, region_id) %>%
    arrange(country_id) %>%
    pull(region_id)

  cell_country_lookup <- df %>%
    group_by(cell_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(cell_id, country_id) %>%
    arrange(cell_id) %>%
    pull(country_id)

  list(country_region_index = country_region_index,
       cell_country_lookup = cell_country_lookup)
}


# parameters -----------------------------------------------------------------

# The model's variables (the free parameters), with their priors, as a named
# list of greta arrays. Every element is passed to model(), so every one has
# named columns in the draws.
dynamical_variables <- function(n_covs, n_classes, n_types, n_regions,
                                n_countries, types,
                                options = dynamical_model_options(),
                                country_region_index = NULL) {

  init <- init_frac_constants(types)

  variables <- list(
    # initial fractions susceptible: a prior logit-mean per type, and IID
    # deviations by region and by country within region
    init_region_sd = normal(0, 1, truncation = c(0, Inf), dim = n_types),
    init_country_sd = normal(0, 1, truncation = c(0, Inf), dim = n_types),
    init_region_raw = normal(0, 1, dim = c(n_regions, n_types)),
    init_country_raw = normal(0, 1, dim = c(n_countries, n_types)),
    # hierarchical regression coefficients: overall -> class -> type
    beta_overall = normal(0, 1, dim = n_covs),
    beta_class_raw = normal(0, 1, dim = c(n_covs, n_classes)),
    beta_type_raw = normal(0, 1, dim = c(n_covs, n_types)),
    sigma_overall = normal(0, 1, dim = n_covs, truncation = c(0, Inf)),
    sigma_class = normal(0, 1, dim = n_covs, truncation = c(0, Inf)),
    logit_init_mean = normal(qlogis(init$relative_prior), 1, dim = n_types)
  )

  rho <- switch(
    options$rho,
    # observation overdispersion per class. With rho = 0.067 the 95% interval
    # of a betabinomial at p = 0.5 is 0.5 wide, a range to treat as unlikely:
    # the half-normal sd of 0.025 puts P(rho < 0.067) at about 0.99 (see the
    # history of fit_model.R for the calculation)
    class = list(
      rho_classes = normal(0, 0.025, truncation = c(0, 1), dim = n_classes)
    ),
    # Per type, nested in class, on the logit scale and non-centred, as the
    # replicate-assay estimate in R/fig_illustrate_bioassay_variability.R:
    #   logit rho_type = rho_mu + rho_sigma_class z_class + rho_sigma_type z_type
    # with the same priors. The prior centre is that of the replicate-assay
    # rho (0.15); rho here also absorbs misfit of the model, and the class-level
    # rho above reached posterior means of 0.24-0.38 despite its half-normal
    # prior, so no stronger prior is put on it.
    type = list(
      rho_mu = normal(qlogis(0.15), 1),
      rho_sigma_class = normal(0, 0.5, truncation = c(0, Inf)),
      rho_sigma_type = normal(0, 0.5, truncation = c(0, Inf)),
      rho_class_raw = normal(0, 1, dim = n_classes),
      rho_type_raw = normal(0, 1, dim = n_types)
    )
  )

  # Floor on bioassay mortality. Predicted mortality is the fraction susceptible
  # q_t (a susceptible mosquito dies at the discriminating dose, a resistant
  # one survives), so a floor f makes it f + (1 - f) q_t: f is the mortality
  # of a fully resistant population, from resistance mechanisms that give
  # finite protection at the discriminating dose, and from deaths by handling
  # rather than insecticide. Beta(1, 9): the density is highest at 0 (no floor,
  # the model without it), with mean 0.1, P(f < 0.2) = 0.87 and
  # P(f < 0.3) = 0.96. WHO tests with control mortality above 20% are
  # discarded and those at 5-20% Abbott-corrected, so handling mortality in the
  # data should be below 0.2, and the discriminating doses are set to kill
  # susceptible mosquitoes with a margin, not resistant ones.
  floor <- if (isTRUE(options$mortality_floor)) {
    list(mortality_floor = beta(1, 9))
  }

  # Coefficients of the standardised initial-state covariates on the logit
  # relative initial state, per type. Independent N(0, 1) rather than
  # hierarchical like the selection effects: the rest of the initial state
  # (logit_init_mean, the region and country sds) is estimated per type with
  # no pooling by class either, and a coefficient of 1 moves the initial state
  # by less than the country deviations do (their sds were 1-3 in the last
  # fit).
  n_init_covs <- length(options$init_covariates)
  init_covariates <- if (n_init_covs > 0) {
    list(init_coef = normal(0, 1, dim = c(n_init_covs, n_types)))
  }

  # Reversion to susceptibility (#24), estimated: a per-year rate per class,
  # constrained to move towards susceptibility; see reversion_kappa() for the
  # sign convention. Half-normal with sd 0.3 on the logit scale per year: at
  # 0.1 the odds of resistance halve in 7 years without selection, and at the
  # 97.5% prior quantile (0.67) in one year. The simulation of #24 found the
  # rate recovered, with the data dominating a prior of sd 0.1, and trading
  # off mildly with the selection effects and the mortality floor.
  reversion <- if (identical(options$reversion, "estimated")) {
    list(reversion_rate = normal(0, 0.3, truncation = c(0, Inf),
                                 dim = n_classes))
  }

  # Centred countries (#25): init_country_level, each country's logit
  # relative initial state, with prior N(its region's, init_country_sd),
  # replaces init_country_raw; the region's is logit_init_mean plus its
  # non-centred deviation. The same model as the non-centred one, which puts
  # N(0, 1) on init_country_raw, the country's level less its region's and
  # divided by init_country_sd. With initial-state covariates, the level is
  # the country's at their mean over its modelled cells: the covariates are
  # standardised over the whole mask, and the data cells lie mostly above its
  # mean, so a level at 0 trades off against their coefficients. Most countries have data for most types, which
  # pins their initial states: with non-centred deviations those of every
  # country in a region then move together against their region's and
  # logit_init_mean, a ridge HMC mixes along slowly. Centring the countries'
  # deviations from logit_init_mean alone leaves the ridge against
  # logit_init_mean; centring the regions too puts them in a funnel with
  # init_region_sd, which 5 regions barely identify (a chain stuck for 2,000
  # iterations at init_region_sd 0.06).
  if (identical(options$init_centred, "country")) {
    stopifnot(length(country_region_index) == n_countries)
    region_level <- sweep(sweep(variables$init_region_raw, 2,
                                variables$init_region_sd, FUN = "*"),
                          2, variables$logit_init_mean, FUN = "+")
    country_mean <- region_level[country_region_index, ]
    # the level is at the country's mean initial-state covariates, not at 0
    # (see init_covariate_shift())
    shift <- init_covariate_shift(init_covariates$init_coef, options)
    if (!is.null(shift)) {
      country_mean <- country_mean + shift
    }
    country_sd <- sweep(zeros(n_countries, n_types), 2,
                        variables$init_country_sd, FUN = "+")
    variables$init_country_raw <- NULL
    variables$init_country_level <- normal(country_mean, country_sd)
  }

  c(variables, rho, floor, init_covariates, reversion)
}

# The non-centred deviations of the countries' initial states
# (init_country_raw) from the centred levels (#25), for draws x dim arrays of
# each variable (as variable_draws() gives); `v` is returned with the levels
# replaced. The prediction code draws the deviations of countries without data
# on the non-centred scale.
init_noncentred_draws <- function(v, country_region_index,
                                  options = NULL) {
  if (is.null(v$init_country_level)) {
    return(v)
  }
  centre <- options$init_covariate_centre
  n_draws <- dim(v$init_region_sd)[1]
  n_types <- length(v$init_region_sd) / n_draws
  mean <- array(v$logit_init_mean, c(n_draws, n_types))
  region_sd <- array(v$init_region_sd, c(n_draws, n_types))
  country_sd <- array(v$init_country_sd, c(n_draws, n_types))
  v$init_country_raw <- v$init_country_level
  for (k in seq_len(n_types)) {
    region_level <- v$init_region_raw[, country_region_index, k,
                                      drop = FALSE][, , 1] * region_sd[, k] +
      mean[, k]
    if (!is.null(centre)) {
      region_level <- region_level +
        matrix(v$init_coef[, , k], n_draws) %*% t(centre)
    }
    v$init_country_raw[, , k] <- (v$init_country_level[, , k] - region_level) /
      country_sd[, k]
  }
  v$init_country_level <- NULL
  v
}

# Reversion to susceptibility (#24). A fitness cost c of resistance, paid
# whatever the selection, makes the relative fitness of the resistant phenotype
# w = (1 - c)(1 + x' exp(beta)), so on the logit scale
#   logit q_t = logit q_0 - sum_{s <= t} log(1 + x_s' exp(beta)) - t kappa,
#   kappa = log(1 - c) <= 0,
# i.e. without selection, logit resistance changes by kappa per year and logit
# susceptibility by -kappa. The selection terms stay strictly positive. kappa
# is per class (kdr gives resistance to both DDT and the pyrethroids), expanded
# here to types: n_types, -reversion_rate when estimated, the fixed values
# when options$reversion is numeric, and NULL for none. For greta arrays or one
# draw in plain R.
reversion_kappa <- function(v, classes_index, options) {
  reversion <- options$reversion
  if (isFALSE(reversion)) {
    return(NULL)
  }
  if (identical(reversion, "estimated")) {
    return(-v$reversion_rate[classes_index])
  }
  n_classes <- max(classes_index)
  stopifnot(length(reversion) %in% c(1, n_classes))
  rep_len(reversion, n_classes)[classes_index]
}

# The effect of the initial-state covariates at each country's mean over its
# modelled cells, options$init_covariate_centre (countries x covariates), as
# countries x types: centre %*% init_coef, or NULL without covariates or a
# centre (#25). For greta arrays or plain R.
init_covariate_shift <- function(init_coef, options) {
  centre <- options$init_covariate_centre
  if (is.null(init_coef) || is.null(centre)) {
    return(NULL)
  }
  if (!inherits(init_coef, "greta_array")) {
    init_coef <- matrix(init_coef, ncol(centre))
  }
  centre %*% init_coef
}

# Deterministic transforms from the variables `v` to the quantities the dynamics
# need:
#   beta_type           n_covs x n_types, log selection effect of each covariate
#   logit_init_country  n_countries x n_types, logit of the initial fraction
#                       susceptible q_0 in each country
#   rho_types           n_types, the observation overdispersion of each type
#                       (with rho = "class", its class's)
#   mortality_floor     the floor on bioassay mortality, or NULL for none
#   logit_init_relative n_countries x n_types, the logit relative initial state
#                       (above init_frac_min) of each country, to which the
#                       initial-state covariates are added at each cell
#   init_coef           n_init_covs x n_types, their coefficients, or NULL
#   kappa_type          n_types, the per-year change in logit resistance from
#                       reversion (<= 0), or NULL for none (reversion_kappa())
# and some intermediate quantities. logit_init_country is the initial state
# without covariates, i.e. at a cell whose covariates are all 0 (the mean).
#
# `v` is a named list of either greta arrays or plain R arrays for a single
# posterior draw (dimensions as in dynamical_variables(), vectors as vectors or
# one-column matrices). `country_region_index` maps the rows of
# v$init_country_raw to the rows of v$init_region_raw; for prediction at
# countries without data, pass raw deviations and an index for all countries.
dynamical_terms <- function(v, classes_index, country_region_index, types,
                            options = dynamical_model_options()) {

  is_greta <- inherits(v$beta_overall, "greta_array")
  if (is_greta) {
    inv_logit <- greta::ilogit
  } else {
    v <- lapply(v, function(x) {
      if (length(dim(x)) <= 1 || (is.matrix(x) && ncol(x) == 1)) c(x) else x
    })
    inv_logit <- stats::plogis
  }

  # selection effects: doubly hierarchical
  beta_class_sigma <- sweep(v$beta_class_raw, 1, v$sigma_overall, FUN = "*")
  beta_class <- sweep(beta_class_sigma, 1, v$beta_overall, FUN = "+")
  beta_type_sigma <- sweep(v$beta_type_raw, 1, v$sigma_class, FUN = "*")
  beta_type <- beta_class[, classes_index] + beta_type_sigma

  # initial state: logit relative position above init_frac_min, the prior mean
  # plus region and country deviations
  init_region_effect <- sweep(v$init_region_raw, 2, v$init_region_sd,
                              FUN = "*")
  if (!is.null(v$init_country_level)) {
    # centred countries (#25), with levels at the mean covariates
    shift <- init_covariate_shift(v$init_coef, options)
    logit_init_relative <- if (is.null(shift)) v$init_country_level else
      v$init_country_level - shift
    init_country_effect <- logit_init_relative -
      sweep(init_region_effect, 2, v$logit_init_mean,
            FUN = "+")[country_region_index, ]
  } else {
    init_country_effect <- sweep(v$init_country_raw, 2, v$init_country_sd,
                                 FUN = "*")
    init_country_overall_effect <- init_country_effect +
      init_region_effect[country_region_index, ]
    logit_init_relative <- sweep(init_country_overall_effect, 2,
                                 v$logit_init_mean, FUN = "+")
  }

  init_min <- init_frac_constants(types)$min
  logit_init_country <- logit_init_from_relative(
    logit_init_relative,
    matrix(init_min, nrow(logit_init_relative), length(types), byrow = TRUE))

  # observation overdispersion per type
  rho_types <- switch(
    options$rho,
    class = v$rho_classes[classes_index],
    type = {
      logit_rho_class <- v$rho_mu + v$rho_sigma_class * v$rho_class_raw
      inv_logit(logit_rho_class[classes_index] +
                  v$rho_sigma_type * v$rho_type_raw)
    }
  )

  list(beta_type = beta_type,
       logit_init_country = logit_init_country,
       rho_types = rho_types,
       mortality_floor = v$mortality_floor,
       logit_init_relative = logit_init_relative,
       init_coef = v$init_coef,
       kappa_type = reversion_kappa(v, classes_index, options),
       # intermediate quantities, which the figure scripts read
       beta_class = beta_class,
       init_region_effect = init_region_effect,
       init_country_effect = init_country_effect)
}

# The logit initial state from its logit relative position l above the minimum
# init_frac_min (`min`, conformable with l): q_0 = min + range * ilogit(l). Its
# logit is computed without forming q_0, which rounds to 1 in double precision
# when l is large: with 1 - q_0 = range * ilogit(-l),
#   logit q_0 = log(min + range * ilogit(l)) - log(range) + softplus(l)
# For greta arrays or plain R.
logit_init_from_relative <- function(l, min) {
  if (inherits(l, "greta_array")) {
    inv_logit <- greta::ilogit
    softplus <- greta::log1pe
  } else {
    inv_logit <- stats::plogis
    softplus <- function(x) -stats::plogis(-x, log.p = TRUE)
  }
  range <- 1 - min
  log(min + range * inv_logit(l)) - log(range) + softplus(l)
}


# Initial values for the model's `variables` from a cached set
# (dynamical_inits_file, posterior means from an earlier fit): the cached values
# of the variables the model has (where the dimensions match), for
# rho = "type" a start for rho_mu from the cached class-level rho, and starts
# for the floor, the initial-state coefficients and the reversion rate. Other
# variables start where greta puts them. `columns` are the columns of the
# selection design matrix, to match the selection coefficients to the cached
# fit's (cached_columns: the "columns" attribute of the cached set, or else
# the design without the trend, of the fits before it).
dynamical_inits <- function(cached, variables, columns = NULL,
                            cached_columns = attr(cached, "columns") %||%
                              selection_column_names(
                                selection_design_untrended()),
                            country_region_index = NULL,
                            init_covariate_centre = NULL) {
  cached_columns <- force(cached_columns)
  cached <- unclass(cached)
  # the centred country levels (#25) from cached non-centred deviations,
  # which need the cached logit_init_mean
  if ("init_country_level" %in% names(variables) &&
      is.null(cached$init_country_level)) {
    stopifnot(!is.null(cached$init_country_raw),
              !is.null(cached$logit_init_mean),
              !is.null(country_region_index))
    region_level <- sweep(
      sweep(as.matrix(cached$init_region_raw), 2, c(cached$init_region_sd),
            FUN = "*"),
      2, c(cached$logit_init_mean), FUN = "+")
    cached$init_country_level <-
      sweep(as.matrix(cached$init_country_raw), 2,
            c(cached$init_country_sd), FUN = "*") +
      region_level[country_region_index, , drop = FALSE]
    if (!is.null(init_covariate_centre) && !is.null(cached$init_coef)) {
      cached$init_country_level <- cached$init_country_level +
        init_covariate_centre %*% as.matrix(cached$init_coef)
    }
  }
  # the selection coefficients, by covariate (row): the cached values where
  # the column is in the cached fit's design, and weak selection (a log effect
  # of -4, no deviations) for new columns. The cached pop coefficients put on
  # log population, which is near 0.8 where raw population is near 0, drove p
  # to 0 and stalled a short run at its initial values
  if (!is.null(columns)) {
    selection_starts <- c(beta_overall = -4, beta_class_raw = 0,
                          beta_type_raw = 0, sigma_overall = 0.5,
                          sigma_class = 0.5)
    for (name in intersect(names(selection_starts), names(cached))) {
      old <- as.matrix(cached[[name]])
      new <- matrix(selection_starts[[name]], length(columns), ncol(old))
      shared <- match(columns, cached_columns)
      new[!is.na(shared), ] <- old[shared[!is.na(shared)], ]
      cached[[name]] <- new
    }
  }
  out <- cached[intersect(names(cached), names(variables))]
  # only where the dimensions match (a different selection design changes the
  # number of covariates)
  matches <- vapply(names(out), function(name) {
    identical(as.integer(dim(out[[name]])),
              as.integer(dim(variables[[name]])))
  }, logical(1))
  out <- out[matches]
  if ("rho_mu" %in% names(variables) && !is.null(cached$rho_classes)) {
    out$rho_mu <- mean(qlogis(c(cached$rho_classes)))
  }
  # the other new terms start near the model without them. Left to greta, the
  # reversion rate starts around 1 per year, which drives p to 1 at most
  # assays and stalled a short run at its initial values
  starts <- c(mortality_floor = 0.02, init_coef = 0, reversion_rate = 0.01)
  for (name in intersect(names(starts), setdiff(names(variables),
                                                names(out)))) {
    out[[name]] <- array(starts[[name]], dim(variables[[name]]))
  }
  do.call(greta::initials, out)
}


# sampling -------------------------------------------------------------------

# The sampler settings for the dynamical model, used by fit_fold()
# (R/fit_validation_fold.R) and fit_model.R. The arguments override single
# settings, e.g. for a smoke test. The defaults, and the evidence for them, are
# in doc/cv_run_plan.md (section 3, sampling settings): windowed_hmc() with
# target acceptance 0.65, the number of leapfrog steps
# redrawn every 10 iterations, 4 chains, 2,000 warmup and 5,000 samples.
#   sampler        "hmc", greta's hmc(), or "windowed", windowed_hmc()
#                  (R/windowed_hmc.R), which adapts the mass matrix in windows
#   Lmin, Lmax     range of the number of leapfrog steps, drawn afresh for each
#                  burst of iterations
#   accept_target  target acceptance of the step-size adaptation ("windowed"
#                  only; hmc() fixes it at 0.5)
#   pb_update      iterations per burst while sampling, so how often the
#                  number of leapfrog steps is redrawn. With it fixed for a
#                  burst, a parameter whose trajectory returns near its start
#                  hardly moves for the whole burst
#   one_by_one     TRUE to redraw it every iteration (a burst of one)
dynamical_mcmc_settings <- function(n_chains = 4,
                                    warmup = 2000,
                                    n_samples = 5000,
                                    sampler = c("windowed", "hmc"),
                                    Lmin = 15,
                                    Lmax = 30,
                                    accept_target = 0.65,
                                    pb_update = 10,
                                    one_by_one = FALSE) {
  list(n_chains = n_chains,
       warmup = warmup,
       n_samples = n_samples,
       sampler = match.arg(sampler),
       Lmin = Lmin,
       Lmax = Lmax,
       accept_target = accept_target,
       pb_update = pb_update,
       one_by_one = one_by_one)
}

# Sample the model `m` with `settings` (dynamical_mcmc_settings()), with all
# chains starting from `inits_one` (dynamical_inits()). mcmc() matches the
# initial values to the variables by name, so `variables` (the model's
# variables, as from build_dynamical_model()) are put in the calling frame.
run_dynamical_mcmc <- function(m, variables, inits_one,
                               settings = dynamical_mcmc_settings()) {
  list2env(variables, environment())
  inits <- replicate(settings$n_chains, inits_one, simplify = FALSE)
  sampler <- switch(
    settings$sampler %||% "hmc",
    hmc = {
      # greta's hmc() fixes the target acceptance at 0.5
      hmc(Lmin = settings$Lmin, Lmax = settings$Lmax)
    },
    windowed = windowed_hmc(Lmin = settings$Lmin, Lmax = settings$Lmax,
                            accept_target = settings$accept_target)
  )
  mcmc(m,
       chains = settings$n_chains,
       initial_values = inits,
       warmup = settings$warmup,
       sampler = sampler,
       n_samples = settings$n_samples,
       pb_update = settings$pb_update %||% 50,
       one_by_one = settings$one_by_one %||% FALSE)
}


# The cached initial values for the fits: posterior means of every variable
# of the default model, from a 4-chain fit to the interpolation fold
# (September 2026, before reversion was added; reversion_rate starts at
# 0.01), with the columns of its selection design as attribute "columns".
# Unlike temporary/inits.RDS, from the fits before the refit, it has
# logit_init_mean, which the centred hierarchy needs. Not in git: remake it
# from a fit's draws as in fit_model.R.
dynamical_inits_file <- "temporary/inits_refit.RDS"


# Bioassay mortality from the fraction susceptible q, with the floor f (a
# scalar, or one per draw conformable with q; NULL for none): f + (1 - f) q.
# For greta arrays or plain R.
floored_mortality <- function(q, floor) {
  if (is.null(floor)) {
    return(q)
  }
  floor + (1 - floor) * q
}

# The same on the logit scale: the logit of f + (1 - f) q from l = logit q,
# computed without forming q, which rounds to 1 when l is large. With
# 1 - p = (1 - f) ilogit(-l),
#   logit p = log(f + (1 - f) ilogit(l)) - log(1 - f) + softplus(l)
# which is l when f = 0. Plain R.
floored_logit_mortality <- function(l, floor) {
  if (is.null(floor)) {
    return(l)
  }
  log(floor + (1 - floor) * plogis(l)) - log1p(-floor) -
    plogis(-l, log.p = TRUE)
}


# the model ------------------------------------------------------------------

# Build the greta model with the likelihood over `train_df`.
#
#   train_df          the records in the likelihood (the full data, or a fold)
#   df                the full data, for the dimensions and lookups
#   x_cell_years      covariates, one row per (cell_id, year_id), cell-major
#   cell_years_index  data frame of the cell_id and year_id of each row
#   classes_index     class of each type
#   types             type names, in type_id order
#   x_cells_init      initial-state covariates, one row per cell_id, with
#                     named columns (init_covariate_matrix()); needed when
#                     options$init_covariates is set
#
# Returns a list with the model, its variables (as passed to model()), the
# derived terms, and two functions for predictions after sampling:
# mortality(rows) returns the greta array of predicted bioassay mortality at the
# (cell_id, type_id, year_id) of `rows`, and all_states() the predicted
# mortality at every cell, type and year. Mortality is the state (the fraction
# susceptible) with the mortality floor applied, if there is one.
build_dynamical_model <- function(train_df,
                                  df,
                                  x_cell_years,
                                  cell_years_index,
                                  classes_index,
                                  types,
                                  options = dynamical_model_options(),
                                  x_cells_init = NULL) {

  options <- complete_dynamical_model_options(options)
  check_dynamical_model_options(options)
  check_greta_fill()
  if (!is.null(colnames(x_cell_years))) {
    stopifnot(identical(colnames(x_cell_years),
                        selection_column_names(options$selection_columns)))
  }

  n_covs <- ncol(x_cell_years)
  n_unique_cells <- max(df$cell_id)
  n_times <- max(cell_years_index$year_id)
  n_classes <- max(df$class_id)
  n_types <- length(types)
  n_regions <- max(df$region_id)
  n_countries <- max(df$country_id)
  stopifnot(
    length(classes_index) == n_types,
    # the layout the state computation relies on
    identical(as.integer(cell_years_index$cell_id),
              rep(seq_len(n_unique_cells), each = n_times)),
    identical(as.integer(cell_years_index$year_id),
              rep(seq_len(n_times), n_unique_cells))
  )

  lookups <- dynamical_lookups(df)
  x_init <- select_init_covariates(x_cells_init, options, n_unique_cells)
  # the centred country levels are at each country's mean initial-state
  # covariates over its modelled cells (all of them, whatever the fold; the
  # overall mean for a country with none), recorded in the options for the
  # plain-R predictions
  options$init_covariate_centre <- NULL
  if (identical(options$init_centred, "country") && !is.null(x_init)) {
    x_cells <- x_init[seq_len(n_unique_cells), , drop = FALSE]
    centre <- t(vapply(seq_len(n_countries), function(country) {
      rows <- lookups$cell_country_lookup == country
      if (any(rows)) colMeans(x_cells[rows, , drop = FALSE]) else
        colMeans(x_cells)
    }, numeric(ncol(x_cells))))
    colnames(centre) <- colnames(x_cells)
    options$init_covariate_centre <- centre
  }

  variables <- dynamical_variables(n_covs = n_covs,
                                   n_classes = n_classes,
                                   n_types = n_types,
                                   n_regions = n_regions,
                                   n_countries = n_countries,
                                   types = types,
                                   options = options,
                                   country_region_index =
                                     lookups$country_region_index)
  terms <- dynamical_terms(variables,
                           classes_index = classes_index,
                           country_region_index = lookups$country_region_index,
                           types = types,
                           options = options)

  # predicted mortality (the fraction susceptible) at the (cell_id, type_id,
  # year_id) of `rows`, computing the states only for the cell-type pairs there
  mortality <- function(rows) {
    pairs <- distinct(tibble(cell_id = as.integer(rows$cell_id),
                             type_id = as.integer(rows$type_id)))
    states <- closed_form_states(terms, x_cell_years, pairs$cell_id,
                                 pairs$type_id, lookups$cell_country_lookup,
                                 n_times, types, x_init)
    pair_index <- match(paste(rows$cell_id, rows$type_id),
                        paste(pairs$cell_id, pairs$type_id))
    floored_mortality(states[cbind(pair_index, rows$year_id)],
                      terms$mortality_floor)
  }

  # the predicted mortality at every cell, type and year, as n_unique_cells x
  # n_types x n_times, for predictions after sampling. (Creating this before model() would
  # put it in the model's graph.)
  all_states <- function() {
    pairs <- expand.grid(cell_id = seq_len(n_unique_cells),
                         type_id = seq_len(n_types))
    states <- closed_form_states(terms, x_cell_years, pairs$cell_id,
                                 pairs$type_id, lookups$cell_country_lookup,
                                 n_times, types, x_init)
    dim(states) <- c(n_unique_cells, n_types, n_times)
    floored_mortality(states, terms$mortality_floor)
  }

  # likelihood
  population_mortality_vec <- mortality(train_df)
  rho <- terms$rho_types[train_df$type_id]
  distribution(train_df$died) <- betabinomial_p_rho(
    N = train_df$mosquito_number,
    p = population_mortality_vec,
    rho = rho)

  # model() takes the target names from its call, so it is called with the
  # variables' names as symbols
  target_env <- list2env(variables, parent = environment())
  m <- do.call("model", lapply(names(variables), as.name), envir = target_env)

  list(model = m,
       variables = variables,
       terms = terms,
       population_mortality_vec = population_mortality_vec,
       mortality = mortality,
       all_states = all_states,
       lookups = lookups,
       options = options)
}


# The columns of the initial-state covariates the options ask for, as a
# cells x covariates matrix, or NULL for none.
select_init_covariates <- function(x_cells_init, options, n_cells = NULL) {
  if (is.null(options$init_covariates)) {
    return(NULL)
  }
  if (is.null(x_cells_init)) {
    stop("the initial-state covariates ", toString(options$init_covariates),
         " are needed (see init_covariate_matrix())")
  }
  stopifnot(all(options$init_covariates %in% colnames(x_cells_init)),
            is.null(n_cells) || nrow(x_cells_init) >= n_cells)
  x <- x_cells_init[, options$init_covariates, drop = FALSE]
  stopifnot(!anyNA(x))
  x
}

# The logit relative initial state (above init_frac_min) of rows with countries
# `country` and types `type`, and covariates x_init (rows x covariates, NULL
# for none), from dynamical_terms(): the country's value plus the covariate
# effects of the type. For greta arrays (a column vector) or one draw in plain
# R.
logit_init_relative_rows <- function(terms, country, type, x_init = NULL) {
  n_countries <- nrow(terms$logit_init_relative)
  l <- terms$logit_init_relative[(type - 1) * n_countries + country]
  if (!is.null(x_init)) {
    stopifnot(nrow(x_init) == length(country))
    coef_rows <- t(terms$init_coef)[type, , drop = FALSE]
    # a sum over columns rather than rowSums(), which Matrix (attached by
    # lme4) masks with a version that does not dispatch to greta
    for (j in seq_len(ncol(x_init))) {
      l <- l + x_init[, j] * coef_rows[, j]
    }
  }
  l
}


# the selection recursion ---------------------------------------------------

# Haploid selection: with q the fraction susceptible and w the relative fitness
# of the resistant phenotype,
#   q_t = q_{t-1} / (q_{t-1} + (1 - q_{t-1}) w_t),
# which divides the odds of susceptibility by w_t each year, so it is exactly
# additive on the logit scale:
#   logit q_t = logit q_0 - sum_{s <= t} log w_s,  w_s = 1 + x_s' exp(beta).
# The state recorded for year t has had the fitness of years 1..t applied (as
# greta.dynamics recorded it). This is computed by one greta op (#25).
#
# log w is computed from the linear predictor as
#   m + log(exp(-m) + x' exp(beta - m)),  m = max(0, max_k beta_k)
# per type, with the gradient stopped through m. That is exact for any m, and
# cannot overflow in float64. The covariates are all non-negative.

# The TensorFlow side. Arguments are tensors with a leading batch dimension B
# (greta's), then the constants:
#   beta_type            (B, n_covs, n_types)
#   logit_init           (B, J, 1) logit q_0 of each pair
#   kappa_type           (B, n_types, 1) reversion kappa per type (optional)
#   x_pairs              (J, n_times, n_covs) covariates of each pair's cell
#   pair_type            (J) 0-based type of each pair
# Returns (B, J, n_times): the fraction susceptible for each pair and year.
tf_closed_form_states <- function(beta_type, logit_init, kappa_type = NULL,
                                  x_pairs, pair_type) {
  tf <- tensorflow::tf
  dtype <- beta_type$dtype
  # reshaped to a vector, since reticulate passes a length-one R vector as a
  # scalar
  pair_type <- tf$reshape(tf$constant(pair_type, dtype = tf$int32), list(-1L))
  n_times <- dim(x_pairs)[2]
  x_pairs <- tf$constant(x_pairs, dtype = dtype)

  # log fitness, (B, J, n_times)
  m <- tf$stop_gradient(tf$maximum(tf$reduce_max(beta_type, axis = 1L,
                                                  keepdims = TRUE),
                                   tf$constant(0, dtype = dtype)))
  effect <- tf$gather(tf$exp(beta_type - m), pair_type, axis = 2L)
  m <- tf$transpose(tf$gather(m, pair_type, axis = 2L), c(0L, 2L, 1L))
  selection <- tf$einsum("jtp,bpj->bjt", x_pairs, effect)
  log_w <- m + tf$math$log(tf$exp(-m) + selection)

  logit_q <- logit_init - tf$cumsum(log_w, axis = 2L)

  # reversion: - t kappa in year t (see reversion_kappa())
  if (!is.null(kappa_type)) {
    years <- tf$reshape(tf$range(1, n_times + 1, dtype = dtype),
                        c(1L, 1L, -1L))
    logit_q <- logit_q - tf$gather(kappa_type, pair_type, axis = 1L) * years
  }

  tf$sigmoid(logit_q)
}

# The greta side: the fraction susceptible for the cell-type pairs
# (pair_cell, pair_type), as a J x n_times greta array. `terms` is the output of
# dynamical_terms(); `x_cell_years` has one row per (cell, year), cell-major;
# each cell takes the initial state of country cell_country_lookup[cell], plus
# the effects of its initial-state covariates, row `cell` of x_init (NULL for
# none).
closed_form_states <- function(terms, x_cell_years, pair_cell, pair_type,
                               cell_country_lookup, n_times, types,
                               x_init = NULL) {
  n_covs <- ncol(x_cell_years)
  stopifnot(nrow(x_cell_years) %% n_times == 0,
            length(pair_cell) == length(pair_type))

  # (cells, years, covariates) and each pair's slice of it
  x_cells <- aperm(array(x_cell_years,
                         c(n_times, nrow(x_cell_years) / n_times, n_covs)),
                   c(2, 1, 3))
  x_pairs <- x_cells[pair_cell, , , drop = FALSE]
  pair_country <- cell_country_lookup[pair_cell]
  stopifnot(!anyNA(pair_country))

  # the initial state of each pair
  l <- logit_init_relative_rows(
    terms, pair_country, pair_type,
    if (!is.null(x_init)) x_init[pair_cell, , drop = FALSE])
  logit_init <- logit_init_from_relative(
    l, init_frac_constants(types)$min[pair_type])

  # the TensorFlow function is found in this small environment, which is saved
  # with the node, so a reloaded draws object can still calculate() through it
  op_env <- new.env(parent = globalenv())
  op_env$tf_closed_form_states <- tf_closed_form_states

  # the reversion kappa of each type, if any, is a greta array when estimated
  # and data when fixed
  kappa <- if (!is.null(terms$kappa_type)) {
    list(if (inherits(terms$kappa_type, "greta_array")) terms$kappa_type else
      as_data(terms$kappa_type))
  }

  do.call(greta:::op, c(
    list("closed_form_states",
         terms$beta_type,
         logit_init),
    kappa,
    list(operation_args = list(
           x_pairs = x_pairs,
           pair_type = as.integer(pair_type - 1)),
         tf_operation = "tf_closed_form_states",
         tf_function_env = op_env,
         dim = c(length(pair_cell), n_times))))
}
