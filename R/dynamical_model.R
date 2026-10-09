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
# Source from the repo root, after R/packages.R and R/functions.R.

source("R/greta_setup.R")
source("R/model_covariates.R")
source("R/latent_smooth.R")
source("R/windowed_hmc.R")


# model options ------------------------------------------------------------

# Switches for model terms. The defaults are V5f (#47): population d_half 270
# in selection_design(); a mortality floor per insecticide class, shifted by
# a latent smooth for the pyrethroids and DDT, centred over the cells with
# their bioassays, the floor at a flat smooth with a half-normal prior of
# scale 0.05; a latent smooth of the strength of selection for every class,
# centred over every modelled cell; both smooths at a fixed range of 1,500
# km, each sd with a half-normal prior of scale 0.5 (smooth_options()); and
# the beta-binomial likelihood. V5h, the default before, centred both
# smooths over every modelled cell and had a logit-normal prior on the floor
# intercepts; the fits before V5h had no floor and no smooths (#37), and
# those before #37 d_half 50 and an estimated floor. A fit's own options are
# saved with it (model_options), with the smooths' basis and centring, and
# the scripts that use a fit take them from there (complete_model_options(),
# R/species_fit_helpers.R).
#   mortality_floor   TRUE (the default) for an estimated floor on bioassay
#                     mortality, the mortality of a fully resistant
#                     population (#14), or FALSE for none
#   floor_prior       the Beta shape parameters of the prior of the floor:
#                     Beta(1, 4) by default (Beta(1, 49) before V5h;
#                     dynamical_variables()); for the smooths' floor with
#                     smooth_options(floor_intercept_prior =
#                     "beta_moments"), the normal prior of each logit floor
#                     intercept has its logit's mean and sd. Unused by the default floor, whose prior
#                     is smooth_options(floor_intercept_prior = )
#   init_covariates   names of static covariates of the initial state, from
#                     init_covariate_names(selection_columns), or NULL for
#                     none (#19)
#   selection_columns how the selection design matrix is built
#                     (selection_design(), R/model_covariates.R; #23);
#                     build_dynamical_model() checks the matrix's columns
#                     against it
#   reversion         reversion to susceptibility (#24): "estimated" for one
#                     rate per class, or FALSE for none
#   smooth            smooth_options() (the default; R/latent_smooth.R) for
#                     latent spatial smooths of the strength of selection and
#                     of the mortality floor (V5, #47), or FALSE for none;
#                     with mortality_floor = TRUE, the floor is per class by
#                     default (smooth_options(floor_intercepts = )), and
#                     without the smooths it is one floor
#   centred           which hierarchy levels are sampled centred
#                     (centred_options()); by default the data-informed levels
#                     (centred_options_data_informed(); #48). The same model
#                     either way: it changes only the coordinates HMC moves in
dynamical_model_options <- function(mortality_floor = TRUE,
                                    floor_prior = c(1, 4),
                                    init_covariates =
                                      init_covariate_names(selection_columns),
                                    selection_columns = selection_design(),
                                    reversion = "estimated",
                                    smooth = smooth_options(),
                                    centred = centred_options_data_informed()) {
  list(mortality_floor = mortality_floor,
       floor_prior = floor_prior,
       init_covariates = init_covariates,
       selection_columns = selection_columns,
       reversion = reversion,
       smooth = smooth,
       centred = centred)
}

# Which levels of the selection and overdispersion hierarchies are sampled
# centred.
#
# A hierarchical effect b ~ N(m, s) can be sampled non-centred, as a standard
# normal deviation z with b = m + s z (as greta variable z), or centred, as b
# itself (variable b with prior N(m, s)). The model and its posterior are the
# same either way; only the coordinates HMC moves in differ. Where the data
# say little about b (posterior sd near the prior sd s), the centred form is a
# funnel: as s shrinks, b is squeezed towards m, and no one step size suits
# both ends. Where the data pin b down (posterior sd much less than s), the
# non-centred form is the problem: with b fixed by the data, z has to move
# as (b - m) / s whenever s or m moves, a curved ridge in (s, z) that HMC can
# only follow with small steps; centred, b sits still and s moves freely.
#   selection  names of selection design columns (selection_column_names())
#              whose class and type levels are both centred:
#                beta_class[c, ] ~ N(beta_overall[c], sigma_overall[c])
#                beta_type[c, ] ~ N(beta_class[c, class of type],
#                                   sigma_class[c])
#              The other columns stay non-centred at both levels. At least
#              one column must stay non-centred.
#   rho_type   TRUE to centre the type level of the overdispersion,
#              logit rho_type ~ N(logit rho_class, rho_sigma_type); the class
#              level stays non-centred
centred_options <- function(selection = NULL, rho_type = FALSE) {
  list(selection = selection,
       rho_type = rho_type)
}

# The centring of the data-informed levels (#48): both levels of the net use,
# IRS and population selection effects, and the overdispersion of the types. In
# five fits of October 2026 (posterior sd over the prior sd at that level):
#   - the type-level effects of these columns were 0.04-0.65 for the
#     pyrethroids and mostly 0.1-1.5 for the other types (up to 3.9 for net
#     use on the organophosphates in one fit), and their class-level effects
#     0.2-1.2;
#   - the types' overdispersion was 0.06-0.33 (the classes' 0.5-1.2);
#   - the crop columns' type-level effects were mostly 0.8-2.8, i.e.
#     prior-dominated, so they stay non-centred: centred, they would be
#     funnels.
# On those fits, the local curvature of the posterior predicted a step size
# about 4 times larger with this centring, with no new funnel in the tails.
centred_options_data_informed <- function(
    selection = c("nets", "irs", "pop_enc:g_dom")) {
  centred_options(selection = selection, rho_type = TRUE)
}

check_dynamical_model_options <- function(options) {
  reversion <- options$reversion
  init_covariates <- options$init_covariates
  stopifnot(
    # build_dynamical_model() sets init_covariate_centre
    setequal(setdiff(names(options), dynamical_built_options),
             names(dynamical_model_options())),
    is.null(options$init_covariate_centre) ||
      identical(colnames(options$init_covariate_centre), init_covariates),
    isFALSE(options$mortality_floor) || isTRUE(options$mortality_floor),
    is.numeric(options$floor_prior), length(options$floor_prior) == 2,
    all(options$floor_prior > 0),
    isFALSE(reversion) || identical(reversion, "estimated"),
    is.null(init_covariates) ||
      (is.character(init_covariates) && !anyDuplicated(init_covariates) &&
         all(init_covariates %in%
               init_covariate_names(options$selection_columns))))
  # errors on a design selection_design() does not build
  complete_selection_design(options$selection_columns)
  centred <- options$centred
  columns <- selection_column_names(options$selection_columns)
  stopifnot(
    setequal(names(centred), names(centred_options())),
    is.null(centred$selection) ||
      (is.character(centred$selection) && !anyDuplicated(centred$selection)),
    isFALSE(centred$rho_type) || isTRUE(centred$rho_type))
  if (!all(centred$selection %in% columns)) {
    stop("centred selection columns not in the design: ",
         toString(setdiff(centred$selection, columns)))
  }
  if (all(columns %in% centred$selection)) {
    stop("at least one selection column must stay non-centred")
  }
  check_smooth_options(options$smooth)
  if (smooth_on(options) && !isFALSE(options$smooth$floor) &&
      !isTRUE(options$mortality_floor)) {
    stop("the smooth of the floor (smooth_options(floor = )) needs ",
         "mortality_floor = TRUE")
  }
  invisible(options)
}

# Whether each selection effect (row of beta_class and beta_type, column of
# the selection design) is centred, as a logical vector in design order. All
# FALSE for options saved before centring was an option.
centred_selection_rows <- function(options) {
  columns <- selection_column_names(options$selection_columns)
  columns %in% options$centred$selection
}

# Whether the type level of the overdispersion is centred
centred_rho_type <- function(options) {
  isTRUE(options$centred$rho_type)
}

# The options build_dynamical_model() adds, for the plain-R predictions
dynamical_built_options <- "init_covariate_centre"

# Options of experiments of #47 since removed, with the value that is the
# model without them: the likelihood (the weighted binomial), the species
# model (species_options()) and the kdr covariate (kdr_options()); and of the
# latent smooths, the shear of the selection smooth and the smooth of the
# initial state (smooth_options(shear = , init = ))
removed_options <- list(likelihood = "beta_binomial", species = FALSE,
                        kdr = FALSE)
removed_smooth_options <- list(shear = FALSE, init = FALSE)

# `options` (a list, or a call of dynamical_model_options()) without the
# removed options (`removed`, by name), each of which must have the value
# that is the model without it, and so too its smooth options (a list, or a
# call of smooth_options()) without removed_smooth_options; a fit with
# another value needs the code that had it (branch weighted-binomial)
drop_removed_options <- function(options, removed = removed_options) {
  smooth <- options[["smooth"]]
  if (is.list(smooth) || (is.call(smooth) &&
                          identical(smooth[[1]], quote(smooth_options)))) {
    options[["smooth"]] <- drop_removed_options(smooth,
                                                removed_smooth_options)
  }
  for (name in intersect(names(options), names(removed))) {
    # in a call, the argument as written if it needs code since removed
    value <- if (is.call(options)) {
      tryCatch(eval(options[[name]]), error = function(e) options[[name]])
    } else {
      options[[name]]
    }
    if (!identical(value, removed[[name]])) {
      stop("the option ", name, " = ", deparse(value), " is not in the ",
           "current model, which is ", name, " = ",
           deparse(removed[[name]]))
    }
    options[[name]] <- NULL
  }
  options
}

# The options of an option string (IR_CUBE_MODEL_OPTIONS, a call of
# dynamical_model_options()), with the removed options dropped
# (drop_removed_options()), so that the option strings of fits made before
# they were removed (e.g. a job's job.options) give the same model
model_options_from_string <- function(expression) {
  call <- str2lang(expression)
  stopifnot(identical(call[[1]], quote(dynamical_model_options)))
  eval(drop_removed_options(call))
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

  # the names of the types, classes, regions and countries, in id order
  in_id_order <- function(name, id) df[[name]][match(seq_len(max(df[[id]])),
                                                     df[[id]])]
  levels <- list(types = in_id_order("insecticide_type", "type_id"),
                 classes = in_id_order("insecticide_class", "class_id"),
                 regions = in_id_order("region", "region_id"),
                 countries = in_id_order("country_name", "country_id"))

  list(country_region_index = country_region_index,
       cell_country_lookup = cell_country_lookup,
       levels = levels)
}


# parameters -----------------------------------------------------------------

# The model's variables (the free parameters), with their priors, as a named
# list of greta arrays. Every element is passed to model(), so every one has
# named columns in the draws.
#
# countries is "centred" for the model, or "noncentred" for the same posterior
# with init_country_raw ~ N(0, 1) in place of init_country_level, which is
# then mean + sd raw and returned as attribute "init_country_level": for
# R/floor_profile.R, which maximises the log density, unbounded in the
# centred form where a type's country levels can all meet their means as its
# init_country_sd goes to 0
dynamical_variables <- function(n_covs, n_classes, n_types, n_regions,
                                n_countries, types,
                                options = dynamical_model_options(),
                                country_region_index = NULL,
                                classes_index = NULL,
                                countries = c("centred", "noncentred")) {
  countries <- match.arg(countries)

  init <- init_frac_constants(types)
  # whether each selection effect (row) is centred (centred_options())
  centred <- centred_selection_rows(options)
  stopifnot(length(centred) == n_covs)
  n_noncentred <- sum(!centred)

  variables <- list(
    # initial fractions susceptible: a prior logit-mean per type, and IID
    # deviations by region and by country within region
    init_region_sd = normal(0, 1, truncation = c(0, Inf), dim = n_types),
    init_country_sd = normal(0, 1, truncation = c(0, Inf), dim = n_types),
    init_region_raw = normal(0, 1, dim = c(n_regions, n_types)),
    # hierarchical regression coefficients: overall -> class -> type, with
    # the standard normal deviations of the non-centred rows at each level
    # (all rows unless some are centred, below)
    beta_overall = normal(0, 1, dim = n_covs),
    beta_class_raw = normal(0, 1, dim = c(n_noncentred, n_classes)),
    beta_type_raw = normal(0, 1, dim = c(n_noncentred, n_types)),
    sigma_overall = normal(0, 1, dim = n_covs, truncation = c(0, Inf)),
    sigma_class = normal(0, 1, dim = n_covs, truncation = c(0, Inf)),
    logit_init_mean = normal(qlogis(init$relative_prior), 1, dim = n_types)
  )

  # The centred rows of the selection effects (centred_options()): the class
  # and type effects themselves, beta_class_centred (n_centred x n_classes)
  # and beta_type_centred (n_centred x n_types), with the priors the
  # non-centred rows imply,
  #   beta_class[c, ] ~ N(beta_overall[c], sigma_overall[c])
  #   beta_type[c, ] ~ N(beta_class[c, class of type], sigma_class[c])
  # so the model is the same. The mean and sd of each row are repeated
  # across its columns.
  if (any(centred)) {
    rows <- which(centred)
    n_centred <- length(rows)
    class_mean <- sweep(zeros(n_centred, n_classes), 1,
                        variables$beta_overall[rows], FUN = "+")
    class_sd <- sweep(zeros(n_centred, n_classes), 1,
                      variables$sigma_overall[rows], FUN = "+")
    variables$beta_class_centred <- normal(class_mean, class_sd)
    type_mean <- variables$beta_class_centred[, classes_index]
    type_sd <- sweep(zeros(n_centred, n_types), 1,
                     variables$sigma_class[rows], FUN = "+")
    variables$beta_type_centred <- normal(type_mean, type_sd)
  }

  # Observation overdispersion per type, nested in class, on the logit scale
  # and non-centred (the type level can be centred, below), as the
  # replicate-assay estimate in R/fig_illustrate_bioassay_variability.R (#20):
  #   logit rho_type = rho_mu + rho_sigma_class z_class + rho_sigma_type z_type
  # with the same priors. The prior centre is that of the replicate-assay rho
  # (0.15); rho here also absorbs misfit of the model, and the class-level rho
  # of the fits before the refit reached posterior means of 0.24-0.38 despite
  # a half-normal prior with sd 0.025, so no stronger prior is put on it.
  rho <- list(
    rho_mu = normal(qlogis(0.15), 1),
    rho_sigma_class = normal(0, 0.5, truncation = c(0, Inf)),
    rho_sigma_type = normal(0, 0.5, truncation = c(0, Inf)),
    rho_class_raw = normal(0, 1, dim = n_classes)
  )
  # The type level, non-centred, or centred (centred_options()): logit rho
  # of each type itself, with the prior the non-centred form implies,
  #   logit rho_type ~ N(logit rho_class, rho_sigma_type)
  if (centred_rho_type(options)) {
    logit_rho_class <- logit_rho_classes(rho)
    rho$logit_rho_type <- normal(logit_rho_class[classes_index],
                                 rho$rho_sigma_type)
  } else {
    rho$rho_type_raw <- normal(0, 1, dim = n_types)
  }

  # Floor on bioassay mortality (#14): predicted mortality f + (1 - f) q_t,
  # f the mortality of a fully resistant population (mechanisms with finite
  # protection at the discriminating dose, handling deaths). Beta(1, 49): mode
  # at 0 (no floor), mean 0.02, P(f > 0.1) = 0.006. WHO tests with control
  # mortality above 20% are discarded and those at 5-20% Abbott-corrected.
  # Beta(1, 9) left a second mode once the initial state was constrained
  # (#19): f near 0.27, about 59 lower in log posterior than f near 0.002,
  # which trapped whole chains and folds. options$floor_prior sets another
  # prior, e.g. Beta(1, 4) (#47).
  floor <- if (isTRUE(options$mortality_floor)) {
    list(mortality_floor = beta(options$floor_prior[1],
                                options$floor_prior[2]))
  }
  # The floor of the latent smooths (V5, #47), in its place:
  # plogis(floor_intercept + u_f(x)), the intercept per class or one for all
  # (smooth_options(floor_intercepts = )); the intercept is the logit floor
  # where u_f is 0, its mean over the cells with bioassays of the classes it
  # applies to (smooth_centre()). Its prior (smooth_options(
  # floor_intercept_prior = )): by default (V5f), on the floor itself at a
  # flat smooth, f0 = plogis(floor_intercept), sampled as floor_flat,
  #   f0 ~ half-normal(0, 0.05), truncated to [0, 1]:
  # P(f0 > 0.1) = 0.046, small but not negligible; floors up to ~0.1 are
  # nearly free and higher floors are penalised (log density 2.8 at 0, 0.8 at
  # 0.1, -5.2 at 0.2, -15.2 at 0.3), so high floors must come from the smooth
  # locally.
  # Or (V5 to V5h; "beta_moments") on the intercept, the normal prior with
  # the mean and variance of the logit of a Beta(floor_prior) variable,
  # digamma(a) - digamma(b) and trigamma(a) + trigamma(b) (for Beta(1, 4),
  # N(-1.83, 1.39^2); logit_beta_moments())
  if (smooth_floor_on(options)) {
    by_class <- identical(options$smooth$floor_intercepts, "class")
    n_intercepts <- if (by_class) n_classes else 1
    floor <- if (smooth_floor_flat_on(options)) {
      list(floor_flat = normal(0, options$smooth$floor_intercept_prior$scale,
                               dim = n_intercepts, truncation = c(0, 1)))
    } else {
      prior <- logit_beta_moments(options$floor_prior)
      list(floor_intercept = normal(prior$mean, prior$sd,
                                    dim = n_intercepts))
    }
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
  prior <- country_level_prior(c(variables, init_covariates),
                               country_region_index, options)
  init_country_level <- NULL
  if (countries == "centred") {
    variables$init_country_level <- normal(prior$mean, prior$sd)
  } else {
    variables$init_country_raw <- normal(0, 1, dim = c(n_countries, n_types))
    init_country_level <- prior$mean + prior$sd * variables$init_country_raw
  }

  # The latent smooths (V5, #47; R/latent_smooth.R), each with standard
  # normal raw weights of its basis functions (non-centred), and its marginal
  # sd and inverse range with the penalised-complexity priors of
  # smooth_prior_rates(): exponential, on the inverse range because in two
  # dimensions the prior of the range rho is the density of 1 / rho for an
  # exponential 1 / rho. With V5's PC priors (smooth_options(range = NULL,
  # sd_prior = c(1, 0.05))), P(rho < 1,500 km) = 0.05 and P(sd > 1) = 0.05;
  # with smooth_options(sd_prior = list(family = "half_normal", scale = s))
  # (the default, s = 0.5), sd ~ N(0, s^2) truncated to sd > 0 in place of
  # the sd's PC prior (smooth_sd_prior_distribution()). With a fixed range
  # (smooth_options(range = ), the default), there is no inverse range
  # variable
  smooths <- list()
  if (smooth_on(options)) {
    rates <- smooth_prior_rates(options$smooth)
    for (kind in smooth_kinds(options)) {
      names <- smooth_variable_names(kind)
      smooths[[names[["raw"]]]] <- normal(0, 1,
                                          dim = nrow(options$smooth$indices))
      smooths[[names[["sd"]]]] <- smooth_sd_prior_distribution(
        options$smooth)
      if (!smooth_range_fixed(options$smooth)) {
        smooths[[names[["inv_range"]]]] <- exponential(rates$range)
      }
    }
  }

  out <- c(variables, rho, floor, init_covariates, reversion, smooths)
  attr(out, "init_country_level") <- init_country_level
  out
}

# The mean and sd of logit(f) for f ~ Beta(shapes[1], shapes[2])
logit_beta_moments <- function(shapes) {
  list(mean = digamma(shapes[1]) - digamma(shapes[2]),
       sd = sqrt(trigamma(shapes[1]) + trigamma(shapes[2])))
}

# The prior of the countries' logit relative initial states (#25), N(mean,
# sd), as list(mean, sd), both countries x types: mean, the region's level,
# logit_init_mean plus init_region_sd times its deviation init_region_raw,
# plus the initial-state covariates' effect at the country's mean
# (init_covariate_shift()); sd, init_country_sd of the type. From the
# variables `v`, greta arrays or one draw in plain R (matrices, vectors as
# vectors or one-column matrices).
country_level_prior <- function(v, country_region_index, options) {
  n_countries <- length(country_region_index)
  is_greta <- inherits(v$init_region_raw, "greta_array")
  if (is_greta) {
    zero <- zeros(n_countries, length(v$init_country_sd))
  } else {
    v[c("init_region_sd", "logit_init_mean", "init_country_sd")] <- lapply(
      v[c("init_region_sd", "logit_init_mean", "init_country_sd")], c)
    zero <- matrix(0, n_countries, length(v$init_country_sd))
  }
  region_level <- sweep(sweep(v$init_region_raw, 2, v$init_region_sd,
                              FUN = "*"),
                        2, v$logit_init_mean, FUN = "+")
  country_mean <- if (is_greta) region_level[country_region_index, ] else
    region_level[country_region_index, , drop = FALSE]
  shift <- init_covariate_shift(v$init_coef, options)
  if (!is.null(shift)) {
    country_mean <- country_mean + shift
  }
  list(mean = country_mean,
       sd = sweep(zero, 2, v$init_country_sd, FUN = "+"))
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
#   mortality_floor     the floor on bioassay mortality, or NULL for none
#   logit_init_relative n_countries x n_types, the logit relative initial state
#                       (above init_frac_min) of each country, to which the
#                       initial-state covariates are added at each cell
#   init_coef           n_init_covs x n_types, their coefficients, or NULL
#   kappa_type          n_types, the per-year change in logit resistance from
#                       reversion (<= 0), or NULL for none (reversion_kappa())
# and beta_class, which the figure scripts read. With the latent smooths (V5),
# smooth_weights, a list of the weights of each smooth's basis functions
# (smooth_weight_terms()), and with a floor, floor_intercept.
# logit_init_country is the initial state without covariates, i.e. at a cell
# whose covariates are all 0 (the mean).
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
    v <- plain_variables(v)
    inv_logit <- stats::plogis
  }

  # selection effects: doubly hierarchical. The non-centred rows from their
  # standard normal deviations, with the centred rows (variables in their
  # own right; centred_options()) put back in their places
  centred <- centred_selection_rows(options)
  noncentred <- noncentred_selection_effects(v, !centred, classes_index)
  beta_class <- noncentred$beta_class
  beta_type <- noncentred$beta_type
  if (any(centred)) {
    beta_class <- stack_rows(v$beta_class_centred, beta_class, centred)
    beta_type <- stack_rows(v$beta_type_centred, beta_type, centred)
  }

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

  # observation overdispersion per type, its type level non-centred or
  # centred (centred_options())
  if (centred_rho_type(options)) {
    logit_rho_type <- v$logit_rho_type
  } else {
    logit_rho_class <- logit_rho_classes(v)
    logit_rho_type <- logit_rho_class[classes_index] +
      v$rho_sigma_type * v$rho_type_raw
  }
  rho_types <- inv_logit(logit_rho_type)

  terms <- list(beta_type = beta_type,
                logit_init_country = logit_init_country,
                rho_types = rho_types,
                mortality_floor = v$mortality_floor,
                logit_init_relative = logit_init_relative,
                init_coef = v$init_coef,
                kappa_type = reversion_kappa(v, classes_index, options),
                beta_class = beta_class)
  if (smooth_on(options)) {
    terms$smooth_weights <- smooth_weight_terms(v, options)
    terms$floor_intercept <- smooth_floor_intercept(v, options)
  }
  terms
}

# One draw of the variables in plain R (`v`, a named list of arrays) with
# vectors as plain vectors, one-column matrices included
plain_variables <- function(v) {
  lapply(v, function(x) {
    if (length(dim(x)) <= 1 || (is.matrix(x) && ncol(x) == 1)) c(x) else x
  })
}

# The selection effects of the non-centred rows `rows` (logical, in design
# order) from their standard normal deviations: beta_class (rows x n_classes)
# and beta_type (rows x n_types),
#   beta_class = beta_overall + sigma_overall beta_class_raw
#   beta_type = beta_class[, class of type] + sigma_class beta_type_raw
# The hyperparameters are subset only when some rows are centred, so that
# without centring the model is exactly as before. For greta arrays or plain
# R.
noncentred_selection_effects <- function(v, rows, classes_index) {
  beta_overall <- v$beta_overall
  sigma_overall <- v$sigma_overall
  sigma_class <- v$sigma_class
  if (!all(rows)) {
    beta_overall <- beta_overall[which(rows)]
    sigma_overall <- sigma_overall[which(rows)]
    sigma_class <- sigma_class[which(rows)]
  }
  beta_class_sigma <- sweep(v$beta_class_raw, 1, sigma_overall, FUN = "*")
  beta_class <- sweep(beta_class_sigma, 1, beta_overall, FUN = "+")
  beta_type_sigma <- sweep(v$beta_type_raw, 1, sigma_class, FUN = "*")
  beta_type <- beta_class[, classes_index, drop = FALSE] + beta_type_sigma
  list(beta_class = beta_class, beta_type = beta_type)
}

# The logit overdispersion of each class, non-centred:
#   logit rho_class = rho_mu + rho_sigma_class rho_class_raw
# For greta arrays or plain R.
logit_rho_classes <- function(v) {
  v$rho_mu + v$rho_sigma_class * v$rho_class_raw
}

# The rows of the centred and non-centred selection effects (`centred` and
# `noncentred`, matrices with the same columns) as one matrix, in design order
# (`is_centred`, logical over its rows). For greta arrays or plain R. In
# greta, rbind() is one concat; on the full data, with 4 chains, this and
# the centred priors added about 0.2 ms to the 150 ms of a gradient, against
# 1.2 ms filling a matrix of zeros by `[<-` (#48).
stack_rows <- function(centred, noncentred, is_centred) {
  stacked <- rbind(centred, noncentred)
  # the stacked rows' places in design order; the default centred columns
  # (nets, irs, population) are the first in the design, so need no reordering
  design_order <- order(c(which(is_centred), which(!is_centred)))
  if (identical(design_order, seq_along(design_order))) {
    return(stacked)
  }
  stacked[design_order, , drop = FALSE]
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
# the floor, the initial-state coefficients and the reversion rate. Other
# variables start where greta puts them. The cached values are matched to the
# model by name: `levels` are the model's types, classes, regions and countries
# (build_dynamical_model()'s lookups$levels), matched to the cached fit's
# (attribute "levels"), and a level new to the model starts at the mean of the
# others; `columns` are the columns of the selection design matrix, matched to
# the cached fit's (attribute "columns"). The cached values are those of the
# non-centred model; for a model with centred levels (`options`, the model's,
# as build_dynamical_model() returns them, and `classes_index`), they are
# moved to the centred variables at the same point (centre_variables()).
# Cached variables the model does not have are left out.
dynamical_inits <- function(cached, variables, levels, columns = NULL,
                            options = NULL, classes_index = NULL) {
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
  # of -4, no deviations) for new columns. The cached pop coefficients put on
  # log population, which is near 0.8 where raw population is near 0, drove p
  # to 0 and stalled a short run at its initial values
  if (!is.null(columns)) {
    stopifnot(!is.null(cached_columns))
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
  if (any(c("beta_class_centred", "logit_rho_type") %in% names(variables))) {
    stopifnot(!is.null(options), !is.null(classes_index))
    cached <- centre_variables(cached, classes_index, options)
  }
  # the floor where the floor smooth is 0 in the form the model samples it:
  # the floor itself, floor_flat, from a cached logit floor_intercept (V5 to
  # V5h), or the reverse
  # (smooth_floor_intercept())
  if ("floor_flat" %in% names(variables) && is.null(cached$floor_flat) &&
      !is.null(cached$floor_intercept)) {
    cached$floor_flat <- plogis(as.matrix(cached$floor_intercept))
  }
  if ("floor_intercept" %in% names(variables) &&
      is.null(cached$floor_intercept) && !is.null(cached$floor_flat)) {
    cached$floor_intercept <- qlogis(as.matrix(cached$floor_flat))
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
  # assays and stalled a short run at its initial values. The floors of the
  # latent smooths (#47) start at the cached mortality_floor if there is one
  # (so that cached values in either floor mode, R/floor_mode_inits.R, start
  # there), or else the cached floor of the first class at a flat smooth
  # (floor_flat, or plogis(floor_intercept)), or else at 0.02. The latent
  # smooths start flat (raw weights 0), with sd 0.3 and range 3,000 km (with
  # V5's PC priors, the means of sd and 1 / range are 0.33 and 1 / 4,500 km;
  # a fixed range is no variable, so has no start). These starts are for
  # variables the cache does not have: a cache made from one draw of a fit
  # with the smooths (R/draw_inits.R) gives them all
  cached_floor <- if (!is.null(cached$mortality_floor)) {
    c(cached$mortality_floor)[1]
  } else if (!is.null(cached$floor_flat)) {
    c(cached$floor_flat)[1]
  } else {
    0.02
  }
  starts <- c(mortality_floor = 0.02, init_coef = -0.05, reversion_rate = 0.01,
              floor_intercept = qlogis(cached_floor),
              floor_flat = cached_floor,
              smooth_raw_selection = 0, smooth_sd_selection = 0.3,
              smooth_inv_range_selection = 1 / 3,
              smooth_raw_floor = 0, smooth_sd_floor = 0.3,
              smooth_inv_range_floor = 1 / 3)
  for (name in intersect(names(starts), setdiff(names(variables),
                                                names(out)))) {
    out[[name]] <- array(starts[[name]], dim(variables[[name]]))
  }
  do.call(greta::initials, out)
}


# The same point in the variables of the model with centred levels
# (options$centred), from the non-centred model's variables `v` (one draw, or
# cached posterior means): the centred rows of the selection effects as the
# effects themselves, in place of their standard normal deviations,
#   beta_class_centred = beta_overall + sigma_overall beta_class_raw
#   beta_type_centred = beta_class_centred[, class of type] +
#                       sigma_class beta_type_raw
# and, with the type level of the overdispersion centred,
#   logit_rho_type = logit_rho_class[class of type] +
#                    rho_sigma_type rho_type_raw
# in place of rho_type_raw. The new variables are matrices (vectors as one
# column); a level whose variables are not all in `v` is left as it is.
# Plain R; the inverse of noncentre_variables().
centre_variables <- function(v, classes_index, options) {
  centred <- centred_selection_rows(options)
  p <- plain_variables(v)
  selection <- c("beta_overall", "sigma_overall", "sigma_class",
                 "beta_class_raw", "beta_type_raw")
  if (any(centred) && all(selection %in% names(v)) &&
      nrow(p$beta_class_raw) == length(centred)) {
    all_rows <- rep(TRUE, length(centred))
    effects <- noncentred_selection_effects(p, all_rows, classes_index)
    v$beta_class_centred <- effects$beta_class[centred, , drop = FALSE]
    v$beta_type_centred <- effects$beta_type[centred, , drop = FALSE]
    v$beta_class_raw <- p$beta_class_raw[!centred, , drop = FALSE]
    v$beta_type_raw <- p$beta_type_raw[!centred, , drop = FALSE]
  }
  rho <- c("rho_mu", "rho_sigma_class", "rho_class_raw", "rho_sigma_type",
           "rho_type_raw")
  if (centred_rho_type(options) && all(rho %in% names(v))) {
    logit_rho_class <- logit_rho_classes(p)
    logit_rho_type <- logit_rho_class[classes_index] +
      p$rho_sigma_type * p$rho_type_raw
    v$logit_rho_type <- as.matrix(logit_rho_type)
    v$rho_type_raw <- NULL
  }
  v
}

# The same point in the non-centred model's variables, from one draw `v` of
# the variables of the model with centred levels (options$centred): the
# standard normal deviations of the centred rows of the selection effects
# and of the overdispersion of the types,
#   beta_class_raw = (beta_class - beta_overall) / sigma_overall
#   beta_type_raw = (beta_type - beta_class[, class of type]) / sigma_class
#   rho_type_raw = (logit_rho_type - logit_rho_class[class of type]) /
#                  rho_sigma_type
# in place of the effects, as matrices (vectors as one column). Plain R; the
# inverse of centre_variables().
noncentre_variables <- function(v, classes_index, options) {
  centred <- centred_selection_rows(options)
  p <- plain_variables(v)
  if (any(centred)) {
    rows <- which(centred)
    class_deviation <- sweep(p$beta_class_centred, 1, p$beta_overall[rows],
                             FUN = "-")
    class_raw <- sweep(class_deviation, 1, p$sigma_overall[rows], FUN = "/")
    type_deviation <- p$beta_type_centred -
      p$beta_class_centred[, classes_index, drop = FALSE]
    type_raw <- sweep(type_deviation, 1, p$sigma_class[rows], FUN = "/")
    v$beta_class_raw <- stack_rows(class_raw, p$beta_class_raw, centred)
    v$beta_type_raw <- stack_rows(type_raw, p$beta_type_raw, centred)
    v$beta_class_centred <- NULL
    v$beta_type_centred <- NULL
  }
  if (centred_rho_type(options)) {
    logit_rho_class <- logit_rho_classes(p)
    rho_deviation <- p$logit_rho_type - logit_rho_class[classes_index]
    v$rho_type_raw <- as.matrix(rho_deviation / p$rho_sigma_type)
    v$logit_rho_type <- NULL
  }
  v
}

# One draw `i` of the variables `v` (a named list of draws x dim arrays), as
# arrays of dim(variable)
variables_at_draw <- function(v, i) {
  lapply(v, function(a) {
    d <- dim(a)[-1]
    array(a[i + (seq_len(prod(d)) - 1) * nrow(a)], d)
  })
}

# Draws of the variables of a model with centred levels (`draws`, a named
# list of draws x dim arrays, as from calculate() or extract_parameter()) as
# draws of the non-centred model's (noncentre_variables()), the form
# fit_model.R caches for dynamical_inits(). Unchanged without centring.
noncentred_draws <- function(draws, classes_index, options) {
  if (!any(centred_selection_rows(options)) && !centred_rho_type(options)) {
    return(draws)
  }
  n_draws <- nrow(draws[[1]])
  converted <- lapply(seq_len(n_draws), function(i) {
    noncentre_variables(variables_at_draw(draws, i), classes_index, options)
  })
  lapply(setNames(nm = names(converted[[1]])), function(name) {
    values <- lapply(converted, function(draw) draw[[name]])
    stacked <- matrix(unlist(values), n_draws, byrow = TRUE)
    array(stacked, c(n_draws, dim(as.array(values[[1]]))))
  })
}


# sampling -------------------------------------------------------------------

# The sampler settings for the dynamical model, used by fit_fold()
# (R/fit_validation_fold.R) and fit_model.R. The arguments override single
# settings, e.g. for a smoke test. The defaults, and the evidence for them, are
# in doc/cv_run_plan.md (section 3, sampling settings): windowed_hmc() with
# 30 to 60 leapfrog steps (60 to 120 before the centred selection hierarchy,
# #48), redrawn every 10 iterations, target acceptance 0.65, 4 chains, 2,000
# warmup and 1,500 samples (3,000 before #48).
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
  logit_init_mean = "types", rho_type_raw = "types",
  rho_class_raw = "classes", reversion_rate = "classes",
  init_region_raw = c("regions", "types"),
  init_country_level = c("countries", "types"),
  beta_class_raw = c(NA, "classes"), beta_type_raw = c(NA, "types"),
  init_coef = c(NA, "types"))

# The cached initial values for the fits: posterior means of every variable
# of the default model of the time, from a 4-chain fit to the interpolation
# fold (September 2026, before reversion was added; reversion_rate starts at
# 0.01), with the columns of its selection design as attribute "columns" and
# the names of its types, classes, regions and countries as "levels". Not in
# git: remake it from a fit's draws as in fit_model.R. The V5h fits (#47)
# started from temporary/inits_floor_low.RDS instead (IR_CUBE_INITS).
dynamical_inits_file <- "temporary/inits_refit.RDS"

# The cached initial values the fits start from: the files in IR_CUBE_INITS,
# comma-separated, or dynamical_inits_file if it is unset or empty. With
# several, the chains are split between them (dynamical_chain_inits()), e.g.
# to start chains in both mortality-floor modes (#37; R/floor_mode_inits.R),
# or each chain from its own draw of a full fit (one file per chain;
# R/draw_inits.R)
dynamical_inits_files <- function() {
  files <- Sys.getenv("IR_CUBE_INITS")
  if (files == "") dynamical_inits_file else strsplit(files, ",")[[1]]
}

# Initial values for each of `n_chains` chains, as a list for
# run_dynamical_mcmc(): dynamical_inits() of each cached file in `files`, the
# chains split into equal consecutive groups, one per file, in order. With one
# file every chain starts from the same values. `options` (the model's) and
# `classes_index` are needed for a model with centred levels
dynamical_chain_inits <- function(files, variables, levels, columns,
                                  n_chains, options = NULL,
                                  classes_index = NULL) {
  stopifnot(length(files) >= 1, n_chains %% length(files) == 0)
  per_file <- lapply(files, function(file) {
    dynamical_inits(readRDS(file), variables, levels = levels,
                    columns = columns, options = options,
                    classes_index = classes_index)
  })
  per_file[rep(seq_along(files), each = n_chains / length(files))]
}


# The floors in the columns of `x` (named), on the floor's scale: as they are
# for mortality_floor and the floor of the latent smooths where u_f is 0
# (floor_flat, V5f), and plogis() of the intercepts of the floor of the
# latent smooths with the other prior (floor_intercept, where u_f is 0; #47)
floor_values <- function(x) {
  intercept <- grepl("^floor_intercept", colnames(x))
  x[, intercept] <- plogis(x[, intercept])
  x
}

# Bioassay mortality from the fraction susceptible q, with the floor f (a
# scalar, or one per draw conformable with q; NULL for none): f + (1 - f) q.
# For greta arrays or plain R.
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
#   countries         "centred" (the model), or "noncentred" for the same
#                     posterior in other coordinates (dynamical_variables()),
#                     for R/floor_profile.R
#
# Returns a list with the model, its variables (as passed to model()), the
# derived terms, and two functions for predictions after sampling:
# mortality(rows) returns the greta array of predicted bioassay mortality at the
# (cell_id, type_id, year_id) of `rows`, and all_states() the predicted
# mortality at every cell, type and year. Mortality is the state (the fraction
# susceptible) with the mortality floor applied, if there is one. With the
# latent smooths (options$smooth, V5, #47), mortality is computed by
# outer_mortality(), with each cell's smooths.
build_dynamical_model <- function(train_df,
                                  df,
                                  x_cell_years,
                                  cell_years_index,
                                  classes_index,
                                  types,
                                  options = dynamical_model_options(),
                                  x_cells_init = NULL,
                                  countries = "centred") {

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
  # the outer form (outer_mortality()) for the latent smooths, and the
  # closed-form op otherwise
  outer <- smooth_on(options)
  # the mask cell of each cell_id
  cells <- df$cell[match(seq_len(n_unique_cells), df$cell_id)]
  # the basis of the latent smooths (V5), set up for these cells (all of them,
  # whatever the fold), each smooth centred over the cells of the bioassays it
  # applies to (of the full data, whatever the fold), and recorded in the
  # options for the plain-R predictions; and the basis at each cell_id
  # (smooth_basis_at(); NULL without them)
  if (smooth_on(options)) {
    term_class_ids <- which(lookups$levels$classes %in% smooth_classes)
    options$smooth <- smooth_box(
      options$smooth, cells, lookups$levels$classes,
      class_cells = unique(df$cell[df$class_id %in% term_class_ids]))
  }
  basis_cells <- prediction_basis(options, cells)
  # the latent smooth `kind` at rows with cell_id `cell` and class_id `class`
  # (0 where the class has none), or NULL without it; smooth_cells, each
  # smooth at every cell_id, is made with the terms below
  row_smooth <- function(kind, cell, class) {
    if (!kind %in% smooth_kinds(options)) {
      return(NULL)
    }
    u <- smooth_cells[[kind]][cell]
    if (isTRUE(options$smooth[[kind]])) u else
      u * smooth_class_weight(options, kind, class)
  }
  # the floor of the latent smooths at rows with cell_id `cell` and type_id
  # `type`, or NULL without it
  row_floor <- function(cell, type) {
    if (!smooth_floor_on(options)) {
      return(NULL)
    }
    class <- classes_index[type]
    smooth_floor_value(
      terms$floor_intercept[smooth_intercept_index(options, class)],
      row_smooth("floor", cell, class))
  }
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

  variables <- dynamical_variables(n_covs = n_covs,
                                   n_classes = n_classes,
                                   n_types = n_types,
                                   n_regions = n_regions,
                                   n_countries = n_countries,
                                   types = types,
                                   options = options,
                                   country_region_index =
                                     lookups$country_region_index,
                                   classes_index = classes_index,
                                   countries = countries)
  # the country levels, a variable when centred, and derived when not
  levels_noncentred <- attr(variables, "init_country_level")
  terms <- dynamical_terms(
    if (is.null(levels_noncentred)) variables else
      c(variables, list(init_country_level = levels_noncentred)),
    classes_index = classes_index,
    types = types,
    options = options)
  # each latent smooth at every cell_id, n_unique_cells x 1
  smooth_cells <- lapply(terms$smooth_weights, function(weights) {
    basis_product(basis_cells, weights)
  })

  # predicted mortality (the fraction susceptible) at the (cell_id, type_id,
  # year_id) of `rows`, computing the states only for the cell-type pairs
  # there
  mortality <- function(rows) {
    pairs <- distinct(tibble(cell_id = as.integer(rows$cell_id),
                             type_id = as.integer(rows$type_id)))
    pair_index <- match(paste(rows$cell_id, rows$type_id),
                        paste(pairs$cell_id, pairs$type_id))
    if (outer) {
      return(outer_mortality(terms, x_cell_years, pairs$cell_id,
                             pairs$type_id, lookups$cell_country_lookup,
                             n_times, types, x_init,
                             row_pair = pair_index, row_year = rows$year_id,
                             row_floor = row_floor(rows$cell_id,
                                                   rows$type_id),
                             row_selection = row_smooth(
                               "selection", rows$cell_id,
                               classes_index[rows$type_id])))
    }
    states <- closed_form_states(terms, x_cell_years, pairs$cell_id,
                                 pairs$type_id, lookups$cell_country_lookup,
                                 n_times, types, x_init)
    rows_states <- states[cbind(pair_index, rows$year_id)]
    floored_mortality(rows_states, terms$mortality_floor)
  }

  # the predicted mortality at every cell, type and year, as n_unique_cells x
  # n_types x n_times, for after sampling (created before model(), it would be
  # in the model's graph)
  all_states <- function() {
    pairs <- expand.grid(cell_id = seq_len(n_unique_cells),
                         type_id = seq_len(n_types))
    if (!outer) {
      states <- closed_form_states(terms, x_cell_years, pairs$cell_id,
                                   pairs$type_id, lookups$cell_country_lookup,
                                   n_times, types, x_init)
      dim(states) <- c(n_unique_cells, n_types, n_times)
      return(floored_mortality(states, terms$mortality_floor))
    }
    # every pair (fastest) and year
    n_pairs <- nrow(pairs)
    row_pair <- rep(seq_len(n_pairs), n_times)
    row_cell <- pairs$cell_id[row_pair]
    row_type <- pairs$type_id[row_pair]
    p <- outer_mortality(terms, x_cell_years, pairs$cell_id, pairs$type_id,
                         lookups$cell_country_lookup, n_times, types,
                         x_init, row_pair = row_pair,
                         row_year = rep(seq_len(n_times), each = n_pairs),
                         row_floor = row_floor(row_cell, row_type),
                         row_selection = row_smooth(
                           "selection", row_cell, classes_index[row_type]))
    dim(p) <- c(n_unique_cells, n_types, n_times)
    p
  }

  # likelihood: the beta-binomial, with the estimated rho of each type
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

# The outer form (outer_mortality(); #47): the cumulative log fitness of each
# pair and year, sum_{s <= t} log w_s, as (B, J, n_times), from the same
# selection term as tf_closed_form_states(). Arguments as there.
tf_cumulative_log_fitness <- function(beta_type, x_pairs, pair_type) {
  tf <- tensorflow::tf
  pair_type <- tf$reshape(tf$constant(pair_type, dtype = tf$int32), list(-1L))
  x_pairs <- tf$constant(x_pairs, dtype = beta_type$dtype)
  shifted <- tf_shifted_selection(beta_type, x_pairs, pair_type)
  tf$cumsum(tf_log_fitness(shifted), axis = 2L)
}

# The selection term x' exp(beta) of each pair and year, shifted for the
# log fitness (tf_log_fitness()): as list(m, selection), m the shift,
# (B, J, 1), and selection = x' exp(beta - m), (B, J, n_times). m is
# max(0, max_k beta_k) per type, with its gradient stopped: exp(beta - m) is
# then at most 1 and cannot overflow, the covariates are all non-negative, and
# the shift cancels in tf_log_fitness(), so it is exact for any m. Arguments
# as tf_closed_form_states(), x_pairs and pair_type as tensors.
tf_shifted_selection <- function(beta_type, x_pairs, pair_type) {
  tf <- tensorflow::tf
  dtype <- beta_type$dtype
  m <- tf$stop_gradient(tf$maximum(tf$reduce_max(beta_type, axis = 1L,
                                                  keepdims = TRUE),
                                   tf$constant(0, dtype = dtype)))
  effect <- tf$gather(tf$exp(beta_type - m), pair_type, axis = 2L)
  m <- tf$transpose(tf$gather(m, pair_type, axis = 2L), c(0L, 2L, 1L))
  list(m = m,
       selection = tf$einsum("jtp,bpj->bjt", x_pairs, effect))
}

# The log fitness log(1 + x' exp(beta)) from a shifted selection term
# (tf_shifted_selection()), as m + log(exp(-m) + selection), (B, J, n_times)
tf_log_fitness <- function(shifted) {
  tf <- tensorflow::tf
  shifted$m + tf$math$log(tf$exp(-shifted$m) + shifted$selection)
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
  inputs <- pair_inputs(terms, x_cell_years, pair_cell, pair_type,
                        cell_country_lookup, n_times, types, x_init)
  x_pairs <- inputs$x_pairs
  logit_init <- inputs$logit_init

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

# The product basis %*% weights of a data matrix `basis` (n x m) and a greta
# array `weights` (m x 1), as an n x 1 greta array, by one greta op. greta's
# %*% would copy the data matrix to every chain's slice of the batch at every
# evaluation of the density (for the latent smooths, 3,290 x 539 per chain);
# the einsum of tf_basis_product() broadcasts it instead
basis_product <- function(basis, weights) {
  op_env <- new.env(parent = globalenv())
  op_env$tf_basis_product <- tf_basis_product
  greta:::op("basis_product", weights,
             operation_args = list(basis = basis),
             tf_operation = "tf_basis_product",
             tf_function_env = op_env,
             dim = c(nrow(basis), 1L))
}

# The TensorFlow side of basis_product(): weights (B, m, 1), basis an R
# matrix (n, m); returns (B, n, 1)
tf_basis_product <- function(weights, basis) {
  tf <- tensorflow::tf
  basis <- tf$constant(basis, dtype = weights$dtype)
  tf$einsum("nm,bmk->bnk", basis, weights)
}

# The covariates and initial state of the cell-type pairs (pair_cell,
# pair_type), for closed_form_states() and outer_mortality(): x_pairs,
# J x n_times x n_covs, and logit_init, logit q_0 of each pair (a J x 1 greta
# array). Arguments as closed_form_states().
pair_inputs <- function(terms, x_cell_years, pair_cell, pair_type,
                        cell_country_lookup, n_times, types, x_init = NULL) {
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

  list(x_pairs = x_pairs, logit_init = logit_init)
}

# Predicted bioassay mortality by the outer form (#47), for the latent
# smooths, at pair row_pair (of the cell-type pairs pair_cell, pair_type)
# and year index row_year of each row, a greta array with one element per
# row. row_selection is the latent smooth of selection at each row (NULL
# without it). Other arguments as closed_form_states().
#
# The cumulative log fitness C_t of the closed form is multiplied by the
# multiplier of selection at the row's cell:
#   logit q_t = logit q_0 - exp(u_s(x)) C_t - t kappa.
# The multiplier is outside the log so that one cumulative sum C serves
# every cell, and costs only these few operations at each row; where
# x' exp(beta) is small, exp(u) log(1 + x' exp(beta)) is close to
# log(1 + exp(u) x' exp(beta)), a multiplier on the selection effects. With
# u_s 0 this is the closed form. The floor is mortality_floor, or row_floor,
# that of the latent smooths (smooth_floor_value()) at each row, if given.
outer_mortality <- function(terms, x_cell_years, pair_cell, pair_type,
                            cell_country_lookup, n_times, types,
                            x_init = NULL, row_pair, row_year,
                            row_floor = NULL, row_selection = NULL) {
  inputs <- pair_inputs(terms, x_cell_years, pair_cell, pair_type,
                        cell_country_lookup, n_times, types, x_init)

  # the TensorFlow function and the helpers it calls, found in this small
  # environment, which is saved with the node (see closed_form_states())
  op_env <- new.env(parent = globalenv())
  for (name in c("tf_cumulative_log_fitness", "tf_shifted_selection",
                 "tf_log_fitness")) {
    op_env[[name]] <- get(name)
    environment(op_env[[name]]) <- op_env
  }
  cumulative <- greta:::op("cumulative_log_fitness",
                           terms$beta_type,
                           operation_args = list(
                             x_pairs = inputs$x_pairs,
                             pair_type = as.integer(pair_type - 1)),
                           tf_operation = "tf_cumulative_log_fitness",
                           tf_function_env = op_env,
                           dim = c(length(pair_cell), n_times))

  rows <- list(logit_init = inputs$logit_init[row_pair],
               cumulative = cumulative[cbind(row_pair, row_year)],
               reversion = if (!is.null(terms$kappa_type)) {
                 terms$kappa_type[pair_type[row_pair]] * row_year
               })
  logit_q <- outer_logit(rows, log_selection = row_selection)
  floored_mortality(ilogit(logit_q),
                    if (is.null(row_floor)) terms$mortality_floor else
                      row_floor)
}

# The outer form's logit q: logit q_0 - exp(log_selection) C - reversion,
# from `rows`, a list of logit_init, cumulative and reversion (NULL for
# none), and the log multiplier of selection (NULL for none, a multiplier of
# 1). For greta arrays (one element per row) or plain R (draws x cells, the
# log multiplier draws x cells, and reversion one per draw).
outer_logit <- function(rows, log_selection = NULL) {
  selection <- rows$cumulative
  if (!is.null(log_selection)) {
    selection <- exp(log_selection) * selection
  }
  logit_q <- rows$logit_init - selection
  if (!is.null(rows$reversion)) {
    logit_q <- logit_q - rows$reversion
  }
  logit_q
}
