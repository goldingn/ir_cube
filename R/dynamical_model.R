# The dynamical model, defined once.
#
# build_dynamical_model() builds the greta model; fit_model.R calls it with the
# full data, fit_fold() (R/fit_validation_fold.R) with a cross-validation
# training fold. The priors are in dynamical_variables(), and the transforms
# from the parameters to the selection effects and initial states in
# dynamical_terms(), which works both on greta arrays and on one posterior draw
# in plain R, so the plain-R predictions (R/dynamical_predictions.R) use the
# same definitions. New terms are switched on through dynamical_model_options().
#
# The strength of selection and the mortality floor vary in space by two
# latent smooths, for the pyrethroids and DDT (#47; R/latent_smooth.R):
#   logit q_t(x) = logit q_0(x) - exp(u_s(x)) sum_{s <= t} log w_s - t kappa
#   p_t(x) = f(x) + (1 - f(x)) q_t(x),  f(x) = 1 - exp(-lambda(x))
# with u_s(x) the log multiplier of selection and log lambda(x) a log hazard
# per insecticide type plus the smooth u_f(x). The other classes have neither
# smooth: selection by the covariates alone and a flat floor.
#
# Source from the repo root, after R/packages.R and R/functions.R.

source("R/greta_setup.R")
source("R/model_covariates.R")
source("R/latent_smooth.R")
source("R/windowed_hmc.R")


# model options ------------------------------------------------------------

# Switches for model terms. The defaults (population d_half 270 in
# selection_design(); #37) are not those of the fits before #37 (d_half 50).
# A fit's own options are saved with it (model_options), with the basis of
# the latent smooths (options$smooth, smooth_basis()), and the scripts that
# use a fit take them from there.
#   init_covariates   names of static covariates of the initial state, from
#                     init_covariate_names(selection_columns), or NULL for
#                     none (#19)
#   selection_columns how the selection design matrix is built
#                     (selection_design(), R/model_covariates.R; #23);
#                     build_dynamical_model() checks the matrix's columns
#                     against it
#   reversion         reversion to susceptibility (#24): "estimated" for one
#                     rate per class, or FALSE for none
dynamical_model_options <- function(init_covariates =
                                      init_covariate_names(selection_columns),
                                    selection_columns = selection_design(),
                                    reversion = "estimated") {
  list(init_covariates = init_covariates,
       selection_columns = selection_columns,
       reversion = reversion)
}

check_dynamical_model_options <- function(options) {
  reversion <- options$reversion
  init_covariates <- options$init_covariates
  stopifnot(
    # init_covariate_centre and smooth are set by build_dynamical_model()
    setequal(setdiff(names(options), c("init_covariate_centre", "smooth")),
             names(dynamical_model_options())),
    is.null(options$init_covariate_centre) ||
      identical(colnames(options$init_covariate_centre), init_covariates),
    isFALSE(reversion) || identical(reversion, "estimated"),
    is.null(init_covariates) ||
      (is.character(init_covariates) && !anyDuplicated(init_covariates) &&
         all(init_covariates %in%
               init_covariate_names(options$selection_columns))))
  # errors on a design selection_design() does not build
  complete_selection_design(options$selection_columns)
  invisible(options)
}


# fixed quantities ---------------------------------------------------------

# The number of selection effects (rows of beta_class and beta_type, the first
# columns of the selection design: net use, IRS and population;
# selection_column_names()) sampled centred at both levels (#48). A
# hierarchical effect b ~ N(m, s) can be sampled non-centred, as z with
# b = m + s z, or centred, as b itself; the posterior is the same. Where the
# data pin b down, the non-centred form is a curved ridge in (s, z) that HMC
# follows only with small steps; where they say little, the centred form is a
# funnel. The data pin these effects down; the crop effects stay non-centred
n_centred_rows <- 3

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

  # the names of the types, classes, regions and countries, in id order
  in_id_order <- function(name, id) df[[name]][match(seq_len(max(df[[id]])),
                                                     df[[id]])]
  levels <- list(types = in_id_order("insecticide_type", "type_id"),
                 classes = in_id_order("insecticide_class", "class_id"),
                 regions = in_id_order("region", "region_id"),
                 countries = in_id_order("country_name", "country_id"))

  list(country_region_index = country_region_index,
       cell_country_lookup = cell_country_lookup,
       # the mask cell of each cell_id
       cells = in_id_order("cell", "cell_id"),
       levels = levels)
}


# parameters -----------------------------------------------------------------

# The model's variables (the free parameters), with their priors, as a named
# list of greta arrays. Every element is passed to model(), so every one has
# named columns in the draws.
dynamical_variables <- function(n_covs, n_classes, n_types, n_regions,
                                n_countries, types, classes_index,
                                options = dynamical_model_options(),
                                country_region_index = NULL) {

  init <- init_frac_constants(types)
  noncentred <- seq(n_centred_rows + 1, n_covs)

  variables <- list(
    # initial fractions susceptible: a prior logit-mean per type, and IID
    # deviations by region and by country within region
    init_region_sd = normal(0, 1, truncation = c(0, Inf), dim = n_types),
    init_country_sd = normal(0, 1, truncation = c(0, Inf), dim = n_types),
    init_region_raw = normal(0, 1, dim = c(n_regions, n_types)),
    # hierarchical regression coefficients: overall -> class -> type, with the
    # standard normal deviations of the non-centred rows at each level
    beta_overall = normal(0, 1, dim = n_covs),
    beta_class_raw = normal(0, 1, dim = c(length(noncentred), n_classes)),
    beta_type_raw = normal(0, 1, dim = c(length(noncentred), n_types)),
    sigma_overall = normal(0, 1, dim = n_covs, truncation = c(0, Inf)),
    sigma_class = normal(0, 1, dim = n_covs, truncation = c(0, Inf)),
    logit_init_mean = normal(qlogis(init$relative_prior), 1, dim = n_types)
  )

  # The centred rows of the selection effects (n_centred_rows): the class and
  # type effects themselves, with the priors the non-centred rows imply,
  #   beta_class[c, ] ~ N(beta_overall[c], sigma_overall[c])
  #   beta_type[c, ] ~ N(beta_class[c, class of type], sigma_class[c])
  centred <- seq_len(n_centred_rows)
  class_mean <- sweep(zeros(n_centred_rows, n_classes), 1,
                      variables$beta_overall[centred], FUN = "+")
  class_sd <- sweep(zeros(n_centred_rows, n_classes), 1,
                    variables$sigma_overall[centred], FUN = "+")
  variables$beta_class_centred <- normal(class_mean, class_sd)
  type_sd <- sweep(zeros(n_centred_rows, n_types), 1,
                   variables$sigma_class[centred], FUN = "+")
  variables$beta_type_centred <- normal(
    variables$beta_class_centred[, classes_index], type_sd)

  # Observation overdispersion per type, nested in class, on the logit scale,
  # as the replicate-assay estimate in R/fig_illustrate_bioassay_variability.R
  # (#20):
  #   logit rho_class = rho_mu + rho_sigma_class z_class
  #   logit rho_type ~ N(logit rho_class, rho_sigma_type)
  # with the same priors, the type level centred (as the data pin it down)
  # and the class level not. The prior centre is that of the replicate-assay
  # rho (0.15); rho here also absorbs misfit of the model, and the class-level
  # rho of the fits before the refit reached posterior means of 0.24-0.38
  # despite a half-normal prior with sd 0.025, so no stronger prior is put on
  # it.
  rho <- list(
    rho_mu = normal(qlogis(0.15), 1),
    rho_sigma_class = normal(0, 0.5, truncation = c(0, Inf)),
    rho_sigma_type = normal(0, 0.5, truncation = c(0, Inf)),
    rho_class_raw = normal(0, 1, dim = n_classes)
  )
  logit_rho_class <- rho$rho_mu + rho$rho_sigma_class * rho$rho_class_raw
  rho$logit_rho_type <- normal(logit_rho_class[classes_index],
                               rho$rho_sigma_type)

  # Floor on bioassay mortality (#14, #47): predicted mortality
  # f + (1 - f) q_t, f the mortality of a fully resistant population
  # (mechanisms with finite protection at the discriminating dose, handling
  # deaths), one floor per insecticide type, 1 - exp(-lambda(x)), with the
  # floor smooth for the pyrethroids and DDT (floor_log_hazard(),
  # R/latent_smooth.R). floor_flat is each class's mean hazard, its types'
  # deviations about it with sd floor_sd_type. The half-normal prior of
  # floor_flat keeps the floors near 0 unless the data demand them, so that
  # resistance rises by selection; looser ones let the floors rise to 0.4-0.6
  # and the chains split into modes.
  floor <- list(
    floor_flat = normal(0, floor_hazard_scale, dim = n_classes,
                        truncation = c(0, Inf)),
    floor_sd_type = normal(0, floor_type_sd_scale, truncation = c(0, Inf)),
    floor_type_raw = normal(0, 1, dim = n_types))

  # The latent smooths of selection and of the floor (#47;
  # R/latent_smooth.R): standard normal raw weights of the basis functions
  # (non-centred), and each smooth's marginal sd, half-normal. The range is
  # fixed (smooth_range)
  n_basis <- nrow(options$smooth$indices)
  smooths <- list(
    smooth_raw_selection = normal(0, 1, dim = n_basis),
    smooth_sd_selection = normal(0, smooth_sd_scale, truncation = c(0, Inf)),
    smooth_raw_floor = normal(0, 1, dim = n_basis),
    smooth_sd_floor = normal(0, smooth_sd_scale, truncation = c(0, Inf)))

  # Coefficients of the standardised initial-state covariates on the logit
  # relative initial state, per type. Independent N(0, 1) rather than
  # hierarchical like the selection effects: the rest of the initial state
  # (logit_init_mean, the region and country sds) is estimated per type with
  # no pooling by class either, and a coefficient of 1 moves the initial state
  # by less than the country deviations do (their sds were 1-3 in the last
  # fit).
  n_init_covs <- length(options$init_covariates)
  init_covariates <- if (n_init_covs > 0) {
    # constrained to be <= 0, so each covariate can only make the initial state
    # more resistant, as the selection effects can only make selection faster.
    # The bioassays sit at the populated end of these covariates, so an
    # unconstrained slope extrapolates a correlative effect without a mechanism
    # across most of the map (positive for population in the first refit, and
    # traded against population-driven selection)
    list(init_coef = normal(0, 1, dim = c(n_init_covs, n_types),
                            truncation = c(-Inf, 0)))
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
  # relative initial state, with prior N(its region's, init_country_sd); the
  # region's is logit_init_mean plus its non-centred deviation. The same model
  # as the non-centred one of the fits before the refit. With
  # initial-state covariates, the level is the country's at their mean over
  # its modelled cells: the covariates are standardised over the whole mask
  # and the data cells lie mostly above its mean, so a level at 0 trades off
  # against their coefficients. Most countries have data for most types, which
  # pins their initial states: non-centred, those of every country in a region
  # then move together against their region's and logit_init_mean, a ridge
  # HMC mixes along slowly. Centring the regions too puts them in a funnel
  # with init_region_sd, which 5 regions barely identify.
  stopifnot(length(country_region_index) == n_countries)
  region_level <- sweep(sweep(variables$init_region_raw, 2,
                              variables$init_region_sd, FUN = "*"),
                        2, variables$logit_init_mean, FUN = "+")
  country_mean <- region_level[country_region_index, ]
  shift <- init_covariate_shift(init_covariates$init_coef, options)
  if (!is.null(shift)) {
    country_mean <- country_mean + shift
  }
  country_sd <- sweep(zeros(n_countries, n_types), 2,
                      variables$init_country_sd, FUN = "+")
  variables$init_country_level <- normal(country_mean, country_sd)

  c(variables, rho, floor, init_covariates, reversion, smooths)
}

# Reversion to susceptibility (#24). A fitness cost c of resistance, paid
# whatever the selection, makes the relative fitness of the resistant phenotype
# w = (1 - c)(1 + x' exp(beta)), so on the logit scale
#   logit q_t = logit q_0 - sum_{s <= t} log(1 + x_s' exp(beta)) - t kappa,
#   kappa = log(1 - c) <= 0,
# i.e. without selection, logit resistance changes by kappa per year and logit
# susceptibility by -kappa. The selection terms stay strictly positive. kappa
# is per class (kdr gives resistance to both DDT and the pyrethroids), expanded
# here to types: n_types, -reversion_rate when estimated, and NULL for none.
# For greta arrays or one draw in plain R.
reversion_kappa <- function(v, classes_index, options) {
  if (isFALSE(options$reversion)) {
    return(NULL)
  }
  -v$reversion_rate[classes_index]
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
#   logit_init_relative n_countries x n_types, the logit relative initial state
#                       (above init_frac_min) of each country, to which the
#                       initial-state covariates are added at each cell
#   init_coef           n_init_covs x n_types, their coefficients, or NULL
#   kappa_type          n_types, the per-year change in logit resistance from
#                       reversion (<= 0), or NULL for none (reversion_kappa())
#   smooth_weights_selection, smooth_weights_floor, selection_offset,
#   floor_log_hazard    the latent smooths and the floor (smooth_terms(),
#                       R/latent_smooth.R)
# and beta_class, which the figure scripts read. logit_init_country is the
# initial state without covariates, i.e. at a cell whose covariates are all 0
# (the mean).
#
# `v` is a named list of either greta arrays or plain R arrays for a single
# posterior draw (dimensions as in dynamical_variables(), vectors as vectors or
# one-column matrices).
dynamical_terms <- function(v, classes_index, types,
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

  # selection effects: doubly hierarchical, the centred rows (n_centred_rows)
  # variables in their own right, the others from their standard normal
  # deviations
  noncentred <- seq(n_centred_rows + 1, length(v$beta_overall))
  beta_class_sigma <- sweep(v$beta_class_raw, 1, v$sigma_overall[noncentred],
                            FUN = "*")
  class_noncentred <- sweep(beta_class_sigma, 1, v$beta_overall[noncentred],
                            FUN = "+")
  beta_type_sigma <- sweep(v$beta_type_raw, 1, v$sigma_class[noncentred],
                           FUN = "*")
  beta_class <- rbind(v$beta_class_centred, class_noncentred)
  beta_type <- rbind(v$beta_type_centred,
                     class_noncentred[, classes_index, drop = FALSE] +
                       beta_type_sigma)

  # initial state: the logit relative position above init_frac_min of each
  # country (#25), its level less the initial-state covariates' effect at the
  # country's mean covariates
  shift <- init_covariate_shift(v$init_coef, options)
  logit_init_relative <- if (is.null(shift)) v$init_country_level else
    v$init_country_level - shift

  init_min <- init_frac_constants(types)$min
  logit_init_country <- floored_logit(
    logit_init_relative,
    matrix(init_min, nrow(logit_init_relative), length(types), byrow = TRUE))

  # observation overdispersion per type
  rho_types <- inv_logit(v$logit_rho_type)

  c(list(beta_type = beta_type,
         logit_init_country = logit_init_country,
         rho_types = rho_types,
         logit_init_relative = logit_init_relative,
         init_coef = v$init_coef,
         kappa_type = reversion_kappa(v, classes_index, options),
         beta_class = beta_class),
    smooth_terms(v, options$smooth, classes_index))
}

# The logit of a + (1 - a) ilogit(l), for a floor `a` (conformable with l, or
# NULL for none): the initial state from its logit relative position l above
# the minimum init_frac_min, and bioassay mortality from logit q above the
# mortality floor. Computed without forming the probability, which rounds to 1
# in double precision when l is large: with 1 - p = (1 - a) ilogit(-l),
#   logit p = log(a + (1 - a) ilogit(l)) - log(1 - a) + softplus(l)
# For greta arrays or plain R.
floored_logit <- function(l, a) {
  if (is.null(a)) {
    return(l)
  }
  if (inherits(l, "greta_array")) {
    inv_logit <- greta::ilogit
    softplus <- greta::log1pe
  } else {
    inv_logit <- stats::plogis
    softplus <- function(x) -stats::plogis(-x, log.p = TRUE)
  }
  log(a + (1 - a) * inv_logit(l)) - log(1 - a) + softplus(l)
}


# Initial values for the model's `variables` from a cached set
# (dynamical_inits_file, posterior means from an earlier fit): the cached values
# of the variables the model has (where the dimensions match), and starts for
# the initial-state coefficients, the reversion rate, the floor and the latent
# smooths. Other variables start where greta puts them. The cached values are
# matched to the model by name: `levels` are the model's types, classes, regions
# and countries (build_dynamical_model()'s lookups$levels), matched to the
# cached fit's (attribute "levels"), and a level new to the model starts at the
# mean of the others; `columns` are the columns of the selection design matrix,
# matched to the cached fit's (attribute "columns").
dynamical_inits <- function(cached, variables, levels, columns = NULL) {
  cached_levels <- attr(cached, "levels")
  cached_columns <- attr(cached, "columns")
  if (is.null(cached_levels)) {
    stop("the cached initial values have no levels; remake them as in ",
         "fit_model.R")
  }
  cached <- unclass(cached)
  for (name in intersect(names(inits_levels), names(cached))) {
    x <- as.matrix(cached[[name]])
    for (d in which(!is.na(inits_levels[[name]]))) {
      level <- inits_levels[[name]][d]
      i <- match(levels[[level]], cached_levels[[level]])
      if (d == 1) {
        x <- x[i, , drop = FALSE]
        x[is.na(i), ] <- rep(colMeans(x, na.rm = TRUE), each = sum(is.na(i)))
      } else {
        x <- x[, i, drop = FALSE]
        x[, is.na(i)] <- rowMeans(x, na.rm = TRUE)
      }
    }
    cached[[name]] <- x
  }
  # the selection coefficients, by covariate (row): the cached values where
  # the column is in the cached fit's design, and weak selection (a log effect
  # of -4, no deviations) for new columns. The centred effects' rows are the
  # first n_centred_rows columns, and the deviations' the others. The cached
  # pop coefficients put on log population, which is near 0.8 where raw
  # population is near 0, drove p to 0 and stalled a short run at its initial
  # values
  if (!is.null(columns)) {
    stopifnot(!is.null(cached_columns))
    selection_starts <- c(beta_overall = -4, sigma_overall = 0.5,
                          sigma_class = 0.5, beta_class_centred = -4,
                          beta_type_centred = -4, beta_class_raw = 0,
                          beta_type_raw = 0)
    centred <- seq_len(n_centred_rows)
    rows_of <- function(design, name) {
      if (grepl("_centred$", name)) design[centred] else
        if (grepl("_raw$", name)) design[-centred] else design
    }
    for (name in intersect(names(selection_starts), names(cached))) {
      old <- as.matrix(cached[[name]])
      old_columns <- rows_of(cached_columns, name)
      if (nrow(old) != length(old_columns)) {
        stop("the cached ", name, " is not in this model's variables; ",
             "remake the cache from a fit of the model, as in fit_model.R")
      }
      new_columns <- rows_of(columns, name)
      new <- matrix(selection_starts[[name]], length(new_columns), ncol(old))
      shared <- match(new_columns, old_columns)
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
  # inside the constraint on the initial-state coefficients
  if (!is.null(out$init_coef)) {
    out$init_coef <- pmin(as.matrix(out$init_coef), -0.05)
  }
  # the other new terms start near the model without them. Left to greta, the
  # reversion rate starts around 1 per year, which drives p to 1 at most
  # assays and stalled a short run at its initial values. The floors start
  # low, and the latent smooths flat (raw weights 0), their sds inside their
  # priors
  starts <- c(init_coef = -0.05, reversion_rate = 0.01, floor_flat = 0.002,
              floor_sd_type = 0.3, floor_type_raw = 0,
              smooth_raw_selection = 0, smooth_sd_selection = 0.02,
              smooth_raw_floor = 0, smooth_sd_floor = 0.02)
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
# 30 to 60 leapfrog steps, redrawn every 10 iterations, target acceptance
# 0.65, 4 chains, 2,000 warmup and 1,500 samples.
#   Lmin, Lmax     range of the number of leapfrog steps, drawn afresh for each
#                  burst of iterations
#   accept_target  target acceptance of the step-size adaptation
#   pb_update      iterations per burst while sampling, so how often the
#                  number of leapfrog steps is redrawn. With it fixed for a
#                  burst, a parameter whose trajectory returns near its start
#                  hardly moves for the whole burst
dynamical_mcmc_settings <- function(n_chains = 4,
                                    warmup = 2000,
                                    n_samples = 1500,
                                    Lmin = 30,
                                    Lmax = 60,
                                    accept_target = 0.65,
                                    pb_update = 10) {
  list(n_chains = n_chains,
       warmup = warmup,
       n_samples = n_samples,
       Lmin = Lmin,
       Lmax = Lmax,
       accept_target = accept_target,
       pb_update = pb_update)
}

# Sample the model `m` with `settings` (dynamical_mcmc_settings()), with chain
# i starting from element i of `inits` (dynamical_chain_inits()). mcmc()
# matches the initial values to the variables by name, so `variables` (the
# model's variables, as from build_dynamical_model()) are put in the calling
# frame.
run_dynamical_mcmc <- function(m, variables, inits,
                               settings = dynamical_mcmc_settings()) {
  list2env(variables, environment())
  stopifnot(length(inits) == settings$n_chains)
  sampler <- windowed_hmc(Lmin = settings$Lmin, Lmax = settings$Lmax,
                          accept_target = settings$accept_target)
  mcmc(m,
       chains = settings$n_chains,
       initial_values = inits,
       warmup = settings$warmup,
       sampler = sampler,
       n_samples = settings$n_samples,
       pb_update = settings$pb_update)
}


# The levels each dimension of a cached variable is indexed by (NA: the
# selection or initial-state covariates), for dynamical_inits()
inits_levels <- list(
  init_region_sd = "types", init_country_sd = "types",
  logit_init_mean = "types", logit_rho_type = "types",
  rho_class_raw = "classes", reversion_rate = "classes",
  floor_flat = "classes", floor_type_raw = "types",
  init_region_raw = c("regions", "types"),
  init_country_level = c("countries", "types"),
  beta_class_raw = c(NA, "classes"), beta_type_raw = c(NA, "types"),
  beta_class_centred = c(NA, "classes"), beta_type_centred = c(NA, "types"),
  init_coef = c(NA, "types"))

# The cached initial values for the fits: posterior means of every variable
# of a 4-chain full fit with the centred hierarchy (October 2026), with the
# columns of its selection design as attribute "columns" and the names of its
# types, classes, regions and countries as "levels". Not in git: remake it
# from a fit's draws as in fit_model.R.
dynamical_inits_file <- "temporary/inits_refit.RDS"

# The cached initial values the fits start from: the files in IR_CUBE_INITS,
# comma-separated, or dynamical_inits_file if it is unset or empty. With
# several, the chains are split between them (dynamical_chain_inits()), e.g.
# to start chains from different points
dynamical_inits_files <- function() {
  files <- Sys.getenv("IR_CUBE_INITS")
  if (files == "") dynamical_inits_file else strsplit(files, ",")[[1]]
}

# Initial values for each of `n_chains` chains, as a list for
# run_dynamical_mcmc(): dynamical_inits() of each cached file in `files`, the
# chains split into equal consecutive groups, one per file, in order. With one
# file every chain starts from the same values
dynamical_chain_inits <- function(files, variables, levels, columns,
                                  n_chains) {
  stopifnot(length(files) >= 1, n_chains %% length(files) == 0)
  per_file <- lapply(files, function(file) {
    dynamical_inits(readRDS(file), variables, levels = levels,
                    columns = columns)
  })
  per_file[rep(seq_along(files), each = n_chains / length(files))]
}


# Bioassay mortality from the fraction susceptible q, with the floor f
# (conformable with q, or NULL for none): f + (1 - f) q. For greta arrays or
# plain R.
floored_mortality <- function(q, floor) {
  if (is.null(floor)) {
    return(q)
  }
  floor + (1 - floor) * q
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
# susceptible) with the mortality floor applied.
build_dynamical_model <- function(train_df,
                                  df,
                                  x_cell_years,
                                  cell_years_index,
                                  classes_index,
                                  types,
                                  options = dynamical_model_options(),
                                  x_cells_init = NULL) {

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
  if (!is.null(x_init)) {
    x_cells <- x_init[seq_len(n_unique_cells), , drop = FALSE]
    centre <- t(vapply(seq_len(n_countries), function(country) {
      rows <- lookups$cell_country_lookup == country
      if (any(rows)) colMeans(x_cells[rows, , drop = FALSE]) else
        colMeans(x_cells)
    }, numeric(ncol(x_cells))))
    colnames(centre) <- colnames(x_cells)
    options$init_covariate_centre <- centre
  }
  # the basis of the latent smooths, with the cells with bioassays of
  # smooth_classes (of the full data, whatever the fold), over which the
  # selection smooth is centred and the floor's mean hazard is taken
  # (smooth_basis()), recorded in the options for the plain-R predictions
  smooth_class_ids <- which(lookups$levels$classes %in% smooth_classes)
  options$smooth <- smooth_basis(
    class_cells = unique(df$cell[df$class_id %in% smooth_class_ids]),
    classes = lookups$levels$classes)

  variables <- dynamical_variables(n_covs = n_covs,
                                   n_classes = n_classes,
                                   n_types = n_types,
                                   n_regions = n_regions,
                                   n_countries = n_countries,
                                   types = types,
                                   classes_index = classes_index,
                                   options = options,
                                   country_region_index =
                                     lookups$country_region_index)
  terms <- dynamical_terms(variables,
                           classes_index = classes_index,
                           types = types,
                           options = options)
  # each latent smooth at every cell_id, n_unique_cells x 1, and whether each
  # type has them (1 or 0)
  basis_cells <- prediction_basis(options$smooth, lookups$cells)
  selection_cells <- basis_product(basis_cells,
                                   terms$smooth_weights_selection) -
    terms$selection_offset
  floor_cells <- basis_product(basis_cells, terms$smooth_weights_floor)
  type_smooth <- as.numeric(options$smooth$term_classes[classes_index])
  # the log multiplier of selection, u_s(x), at cell-type pairs with cell_ids
  # `cell` and type_ids `type` (0 for the types without it)
  pair_selection <- function(cell, type) {
    selection_cells[cell] * type_smooth[type]
  }
  # the floor at rows with cell_ids `cell` and type_ids `type`
  row_floor <- function(cell, type) {
    floor_value(terms$floor_log_hazard[type],
                floor_cells[cell] * type_smooth[type])
  }

  # predicted mortality (the fraction susceptible) at the (cell_id, type_id,
  # year_id) of `rows`, computing the states only for the cell-type pairs there
  mortality <- function(rows) {
    pairs <- distinct(tibble(cell_id = as.integer(rows$cell_id),
                             type_id = as.integer(rows$type_id)))
    states <- closed_form_states(terms, x_cell_years, pairs$cell_id,
                                 pairs$type_id, lookups$cell_country_lookup,
                                 n_times, types, x_init,
                                 log_multiplier = pair_selection(
                                   pairs$cell_id, pairs$type_id))
    pair_index <- match(paste(rows$cell_id, rows$type_id),
                        paste(pairs$cell_id, pairs$type_id))
    floored_mortality(states[cbind(pair_index, rows$year_id)],
                      row_floor(rows$cell_id, rows$type_id))
  }

  # the predicted mortality at every cell, type and year, as n_unique_cells x
  # n_types x n_times, for after sampling (created before model(), it would be
  # in the model's graph)
  all_states <- function() {
    pairs <- expand.grid(cell_id = seq_len(n_unique_cells),
                         type_id = seq_len(n_types))
    states <- closed_form_states(terms, x_cell_years, pairs$cell_id,
                                 pairs$type_id, lookups$cell_country_lookup,
                                 n_times, types, x_init,
                                 log_multiplier = pair_selection(
                                   pairs$cell_id, pairs$type_id))
    floor <- row_floor(pairs$cell_id, pairs$type_id) %*% t(rep(1, n_times))
    p <- floored_mortality(states, floor)
    dim(p) <- c(n_unique_cells, n_types, n_times)
    p
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
# greta.dynamics recorded it). This is computed by one greta op (#25). The
# cumulative log fitness of a cell is multiplied by exp(u_s(x)), the latent
# smooth of selection (#47), outside the log so that one cumulative sum serves
# every cell: where x' exp(beta) is small, exp(u) log(1 + x' exp(beta)) is
# close to log(1 + exp(u) x' exp(beta)), a multiplier on the selection
# effects.
#
# log w is computed from the linear predictor as
#   m + log(exp(-m) + x' exp(beta - m)),  m = max(0, max_k beta_k)
# per type, with the gradient stopped through m. That is exact for any m, and
# cannot overflow in float64. The covariates are all non-negative.

# The TensorFlow side. Arguments are tensors with a leading batch dimension B
# (greta's), then the constants:
#   beta_type            (B, n_covs, n_types)
#   logit_init           (B, J, 1) logit q_0 of each pair
#   log_multiplier       (B, J, 1) log multiplier of each pair's cumulative
#                        log fitness, u_s(x)
#   kappa_type           (B, n_types, 1) reversion kappa per type (optional)
#   x_pairs              (J, n_times, n_covs) covariates of each pair's cell
#   pair_type            (J) 0-based type of each pair
# Returns (B, J, n_times): the fraction susceptible for each pair and year.
tf_closed_form_states <- function(beta_type, logit_init, log_multiplier,
                                  kappa_type = NULL, x_pairs, pair_type) {
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

  logit_q <- logit_init - tf$exp(log_multiplier) * tf$cumsum(log_w, axis = 2L)

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
# none). log_multiplier is u_s(x) at each pair (a J x 1 greta array).
closed_form_states <- function(terms, x_cell_years, pair_cell, pair_type,
                               cell_country_lookup, n_times, types,
                               x_init = NULL, log_multiplier) {
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
  logit_init <- floored_logit(
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
         logit_init,
         log_multiplier),
    kappa,
    list(operation_args = list(
           x_pairs = x_pairs,
           pair_type = as.integer(pair_type - 1)),
         tf_operation = "tf_closed_form_states",
         tf_function_env = op_env,
         dim = c(length(pair_cell), n_times))))
}
