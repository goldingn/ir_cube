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
source("R/species.R")
source("R/kdr_covariate.R")
source("R/latent_smooth.R")
source("R/weighted_binomial.R")
source("R/windowed_hmc.R")


# model options ------------------------------------------------------------

# Switches for model terms. The defaults (population d_half 270 in
# selection_design(), no mortality floor; #37) are not those of the fits
# before #37 (d_half 50, an estimated floor). A fit's own options are saved
# with it (model_options), and the scripts that use a fit take them from there.
#   mortality_floor   TRUE for an estimated floor on bioassay mortality, the
#                     mortality of a fully resistant population (#14), or
#                     FALSE for none
#   floor_prior       the Beta shape parameters of the prior of the floor:
#                     Beta(1, 49) by default (dynamical_variables())
#   init_covariates   names of static covariates of the initial state, from
#                     init_covariate_names(selection_columns), or NULL for
#                     none (#19)
#   selection_columns how the selection design matrix is built
#                     (selection_design(), R/model_covariates.R; #23);
#                     build_dynamical_model() checks the matrix's columns
#                     against it
#   reversion         reversion to susceptibility (#24): "estimated" for one
#                     rate per class, or FALSE for none
#   species           FALSE (the default) for one trajectory for the whole
#                     complex; species_options() (R/species.R) for two, one
#                     for An. arabiensis and one for the other members, mixed
#                     at each bioassay by its arabiensis share (#47). With it,
#                     species_options() sets the floors of both species, and
#                     mortality_floor must be FALSE
#   kdr               FALSE (the default) for no kdr covariate;
#                     kdr_options() (R/kdr_covariate.R) for the map of total
#                     kdr as a covariate of the strength of selection and the
#                     fitness cost, with or without the species model (#47)
#   smooth            FALSE (the default) for no latent smooths;
#                     smooth_options() (R/latent_smooth.R) for latent spatial
#                     smooths of the strength of selection and of the
#                     mortality floor, in place of the kdr covariate (V5,
#                     #47). Not with the species model or the kdr covariate;
#                     with mortality_floor = TRUE, the floor is per class by
#                     default (smooth_options(floor_intercepts = )). With
#                     smooth_options(init = TRUE), a third smooth, of the
#                     initial state, with a loading per type, replaces the
#                     hierarchy of regions and countries (no init_region_*,
#                     init_country_* or init_country_level variables); the
#                     limits init_frac_min and the initial-state covariates
#                     stay
#   likelihood        the likelihood of the bioassays: "beta_binomial" (the
#                     default), with the overdispersion rho per type
#                     estimated, nested in class; or "weighted_binomial"
#                     (#47; R/weighted_binomial.R), the binomial log
#                     likelihood of each bioassay weighted by its design
#                     effect 1 / (1 + (n - 1) rho), with rho per type fixed at
#                     the replicate-bioassay estimates (replicate_rho()), and
#                     no overdispersion parameters
#   centred           which hierarchy levels are sampled centred
#                     (centred_options()); by default the data-informed levels
#                     (centred_options_data_informed(); #48). The same model
#                     either way: it changes only the coordinates HMC moves in
dynamical_model_options <- function(mortality_floor = FALSE,
                                    floor_prior = c(1, 49),
                                    init_covariates =
                                      init_covariate_names(selection_columns),
                                    selection_columns = selection_design(),
                                    reversion = "estimated",
                                    species = FALSE,
                                    kdr = FALSE,
                                    smooth = FALSE,
                                    likelihood = "beta_binomial",
                                    centred = centred_options_data_informed()) {
  list(mortality_floor = mortality_floor,
       floor_prior = floor_prior,
       init_covariates = init_covariates,
       selection_columns = selection_columns,
       reversion = reversion,
       species = species,
       kdr = kdr,
       smooth = smooth,
       likelihood = likelihood,
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
#              level stays non-centred. Ignored with the weighted binomial
#              likelihood, which has no overdispersion parameters
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
    # build_dynamical_model() sets init_covariate_centre, and replicate_rho
    # with the weighted binomial likelihood
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
  if (!(is.character(options$likelihood) && length(options$likelihood) == 1 &&
        options$likelihood %in% dynamical_likelihoods)) {
    stop("likelihood must be one of ", toString(dynamical_likelihoods))
  }
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
  check_species_options(options$species)
  check_kdr_options(options$kdr)
  if (species_on(options) && isTRUE(options$mortality_floor)) {
    stop("with the species model, the floors are set by ",
         "species_options(floors = ); leave mortality_floor FALSE")
  }
  if (!isFALSE(kdr_floor(options)) &&
      (!isTRUE(options$mortality_floor) || species_on(options))) {
    stop("the kdr-dependent floor (kdr_options(floor = )) replaces the ",
         "constant floor: it needs mortality_floor = TRUE and no species model")
  }
  check_smooth_options(options$smooth)
  if (smooth_on(options) && (species_on(options) || kdr_on(options))) {
    stop("the latent smooths (smooth_options()) replace the kdr covariate: ",
         "leave kdr and species FALSE")
  }
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

# Whether the type level of the overdispersion is centred: never with the
# weighted binomial likelihood, which has no overdispersion parameters
centred_rho_type <- function(options) {
  isTRUE(options$centred$rho_type) && rho_estimated(options)
}

# The options build_dynamical_model() adds, for the plain-R predictions
dynamical_built_options <- c("init_covariate_centre", "replicate_rho")

# The likelihoods of the bioassays (dynamical_model_options(likelihood = ))
dynamical_likelihoods <- c("beta_binomial", "weighted_binomial")

# A model's likelihood: options saved before it was an option have the
# beta-binomial
model_likelihood <- function(options) {
  if (is.null(options$likelihood)) "beta_binomial" else options$likelihood
}

# Whether the model estimates the overdispersion rho (the beta-binomial), or
# fixes it at the replicate-bioassay estimates (the weighted binomial)
rho_estimated <- function(options) {
  model_likelihood(options) == "beta_binomial"
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
  # the hierarchy of regions and countries of the initial state, or in its
  # place the smooth of the initial state (smooth_options(init = TRUE))
  hierarchy <- !smooth_init_on(options)

  variables <- c(
    # initial fractions susceptible: a prior logit-mean per type, and IID
    # deviations by region and by country within region
    if (hierarchy) {
      list(
        init_region_sd = normal(0, 1, truncation = c(0, Inf), dim = n_types),
        init_country_sd = normal(0, 1, truncation = c(0, Inf), dim = n_types),
        init_region_raw = normal(0, 1, dim = c(n_regions, n_types)))
    },
    list(
      # hierarchical regression coefficients: overall -> class -> type, with
      # the standard normal deviations of the non-centred rows at each level
      # (all rows unless some are centred, below)
      beta_overall = normal(0, 1, dim = n_covs),
      beta_class_raw = normal(0, 1, dim = c(n_noncentred, n_classes)),
      beta_type_raw = normal(0, 1, dim = c(n_noncentred, n_types)),
      sigma_overall = normal(0, 1, dim = n_covs, truncation = c(0, Inf)),
      sigma_class = normal(0, 1, dim = n_covs, truncation = c(0, Inf)),
      logit_init_mean = normal(qlogis(init$relative_prior), 1, dim = n_types)
    ))

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
  # With the weighted binomial likelihood rho is fixed (replicate_rho()), and
  # there are none of these variables.
  rho <- NULL
  if (rho_estimated(options)) {
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
  }

  # Floor on bioassay mortality (#14): predicted mortality f + (1 - f) q_t,
  # f the mortality of a fully resistant population (mechanisms with finite
  # protection at the discriminating dose, handling deaths). Beta(1, 49): mode
  # at 0 (no floor), mean 0.02, P(f > 0.1) = 0.006. WHO tests with control
  # mortality above 20% are discarded and those at 5-20% Abbott-corrected.
  # Beta(1, 9) left a second mode once the initial state was constrained
  # (#19): f near 0.27, about 59 lower in log posterior than f near 0.002,
  # which trapped whole chains and folds. options$floor_prior sets another
  # prior, e.g. Beta(1, 4) for a reference fit to the species model (#47).
  floor <- if (isTRUE(options$mortality_floor)) {
    list(mortality_floor = beta(options$floor_prior[1],
                                options$floor_prior[2]))
  }
  # The kdr-dependent floor (#47; kdr_options(floor = )), in its place:
  # plogis(floor_intercept + floor_kdr k(x)), with one intercept, or one per
  # class (floor = "class"). The intercept, the logit floor at the mean kdr,
  # has the normal prior with the mean and variance of the logit of a
  # Beta(floor_prior) variable, digamma(a) - digamma(b) and trigamma(a) +
  # trigamma(b) (for Beta(1, 4), N(-1.83, 1.39^2)); the slope N(0, 1)
  floor_intercept <- function(by_class) {
    prior <- logit_beta_moments(options$floor_prior)
    if (by_class) normal(prior$mean, prior$sd, dim = n_classes) else
      normal(prior$mean, prior$sd)
  }
  if (!isFALSE(kdr_floor(options))) {
    floor <- list(
      floor_intercept = floor_intercept(identical(kdr_floor(options), "class")),
      floor_kdr = normal(0, 1))
  }
  # The floor of the latent smooths (V5, #47), in its place too:
  # plogis(floor_intercept + u_f(x)), the intercept per class or one for all
  # (smooth_options(floor_intercepts = )), with the same prior; the
  # intercept is the floor where u_f is 0, its mean over the modelled cells
  if (smooth_floor_on(options)) {
    floor <- list(floor_intercept = floor_intercept(
      identical(options$smooth$floor_intercepts, "class")))
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
  # with init_region_sd, which 5 regions barely identify. With the smooth of
  # the initial state there are no country levels, and `countries` is moot.
  init_country_level <- NULL
  if (hierarchy) {
    stopifnot(length(country_region_index) == n_countries)
    prior <- country_level_prior(c(variables, init_covariates),
                                 country_region_index, options)
    if (countries == "centred") {
      variables$init_country_level <- normal(prior$mean, prior$sd)
    } else {
      variables$init_country_raw <- normal(0, 1, dim = c(n_countries, n_types))
      init_country_level <- prior$mean + prior$sd * variables$init_country_raw
    }
  }

  # The species model (#47): log multipliers on arabiensis's log fitness from
  # selection and on its fitness cost (outer_mortality()), N(0, 1) like the
  # log selection effects (beta_overall), centred on no difference between
  # the species and putting 95% of the prior mass of each multiplier between
  # 0.14 and 7.1; and with species_options(floors = TRUE), a mortality floor
  # for each species, estimated separately, each with the prior
  # species_options()$floor_prior (Beta(1, 4) by default: mean 0.2, P(f >
  # 0.5) = 0.06, which allows the plateaus of #37)
  species <- if (species_on(options)) {
    floor_prior <- options$species$floor_prior
    c(list(gamma_selection = normal(0, 1)),
      if (identical(options$reversion, "estimated")) {
        list(gamma_cost = normal(0, 1))
      },
      if (isTRUE(options$species$floors)) {
        list(other_floor = beta(floor_prior[1], floor_prior[2]),
             arabiensis_floor = beta(floor_prior[1], floor_prior[2]))
      })
  }

  # The kdr covariate (#47): slopes of the log multipliers of the cumulative
  # log fitness and of the fitness cost on a cell's standardised kdr
  # (outer_mortality()), N(0, 1) like gamma_selection: one pair, or with the
  # species model one pair per species, each on its own band
  kdr <- if (kdr_on(options)) {
    suffixes <- if (species_on(options)) c("_other", "_arabiensis") else ""
    slopes <- list()
    for (suffix in suffixes) {
      slopes[[paste0("delta_selection", suffix)]] <- normal(0, 1)
      if (identical(options$reversion, "estimated")) {
        slopes[[paste0("delta_cost", suffix)]] <- normal(0, 1)
      }
    }
    slopes
  }

  # The latent smooths (V5, #47; R/latent_smooth.R), each with standard
  # normal raw weights of its basis functions (non-centred), and its marginal
  # sd and inverse range with the penalised-complexity priors of
  # smooth_prior_rates(): exponential, on the inverse range because in two
  # dimensions the prior of the range rho is the density of 1 / rho for an
  # exponential 1 / rho. With the defaults, P(rho < 1,500 km) = 0.05 and
  # P(sd > 1) = 0.05; with smooth_options(sd_prior = list(family =
  # "half_normal", scale = s)), sd ~ N(0, s^2) truncated to sd > 0 in place
  # of the sd's PC prior (smooth_sd_prior_distribution()); and with the
  # shear, its loading b ~ N(1, 0.5) (smooth_shear_prior). With a fixed range
  # (smooth_options(range = )), there is no inverse range variable. The
  # smooth of the initial state (smooth_options(init = TRUE)) has its sd
  # fixed at 1, and in its place a loading lambda >= 0 per type, each with
  # the sd's prior, so each type's sd of the field
  smooths <- list()
  if (smooth_on(options)) {
    rates <- smooth_prior_rates(options$smooth)
    for (kind in smooth_kinds(options)) {
      names <- smooth_variable_names(kind)
      smooths[[names[["raw"]]]] <- normal(0, 1,
                                          dim = nrow(options$smooth$indices))
      if (kind == "init") {
        smooths[[names[["loading"]]]] <- smooth_sd_prior_distribution(
          options$smooth, dim = n_types)
      } else {
        smooths[[names[["sd"]]]] <- smooth_sd_prior_distribution(
          options$smooth)
      }
      if (!smooth_range_fixed(options$smooth)) {
        smooths[[names[["inv_range"]]]] <- exponential(rates$range)
      }
    }
    if (isTRUE(options$smooth$shear)) {
      smooths$smooth_shear <- normal(smooth_shear_prior$mean,
                                     smooth_shear_prior$sd)
    }
  }

  out <- c(variables, rho, floor, init_covariates, reversion, species, kdr,
           smooths)
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
#                       (with the weighted binomial likelihood the fixed
#                       replicate rho, a plain vector; fixed_rho_types())
#   mortality_floor     the floor on bioassay mortality, or NULL for none
#   logit_init_relative n_countries x n_types, the logit relative initial state
#                       (above init_frac_min) of each country, to which the
#                       initial-state covariates are added at each cell
#   init_coef           n_init_covs x n_types, their coefficients, or NULL
#   kappa_type          n_types, the per-year change in logit resistance from
#                       reversion (<= 0), or NULL for none (reversion_kappa())
# and beta_class, which the figure scripts read. With the species model
# (#47), also gamma_selection and gamma_cost (NULL without reversion), the log
# multipliers of arabiensis's log fitness from selection and of its fitness
# cost (outer_mortality()), and other_floor and arabiensis_floor, the floors
# of the other members of the complex and of arabiensis (NULL for none); the
# other terms are shared by both species. With the kdr covariate, also its
# slopes (kdr_slope_names, outer_mortality()), and with the kdr-dependent
# floor, floor_intercept and floor_kdr. With the latent smooths (V5),
# smooth_weights, a list of the weights of each smooth's basis functions
# (smooth_weight_terms()), with a floor, floor_intercept, and with the smooth
# of the initial state, init_loading, its loadings per type
# (smooth_init_loadings()); logit_init_relative is then the same in every
# country (init_field_level()), and a cell adds lambda u_init(x).
# logit_init_country is the initial state without covariates, i.e. at a cell
# whose covariates are all 0 (the mean), and without u_init.
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
  # country's mean covariates; or with the smooth of the initial state, the
  # same in every country, logit_init_mean less the covariates' effect at
  # their mean over the modelled cells (init_field_level())
  shift <- init_covariate_shift(v$init_coef, options)
  level <- if (smooth_init_on(options)) {
    init_field_level(v$logit_init_mean, options)
  } else {
    v$init_country_level
  }
  logit_init_relative <- if (is.null(shift)) level else level - shift

  init_min <- init_frac_constants(types)$min
  logit_init_country <- floored_logit(
    logit_init_relative,
    matrix(init_min, nrow(logit_init_relative), length(types), byrow = TRUE))

  # observation overdispersion per type, its type level non-centred or
  # centred (centred_options()); or with the weighted binomial likelihood,
  # fixed: the replicate rho, as plain R even in greta
  if (!rho_estimated(options)) {
    rho_types <- fixed_rho_types(options, types)
  } else {
    if (centred_rho_type(options)) {
      logit_rho_type <- v$logit_rho_type
    } else {
      logit_rho_class <- logit_rho_classes(v)
      logit_rho_type <- logit_rho_class[classes_index] +
        v$rho_sigma_type * v$rho_type_raw
    }
    rho_types <- inv_logit(logit_rho_type)
  }

  terms <- list(beta_type = beta_type,
                logit_init_country = logit_init_country,
                rho_types = rho_types,
                mortality_floor = v$mortality_floor,
                logit_init_relative = logit_init_relative,
                init_coef = v$init_coef,
                kappa_type = reversion_kappa(v, classes_index, options),
                beta_class = beta_class)
  if (species_on(options)) {
    terms <- c(terms, list(gamma_selection = v$gamma_selection,
                           gamma_cost = v$gamma_cost,
                           other_floor = v$other_floor,
                           arabiensis_floor = v$arabiensis_floor))
  }
  if (kdr_on(options)) {
    terms <- c(terms, v[intersect(c(kdr_slope_names, "floor_intercept",
                                    "floor_kdr"), names(v))])
  }
  if (smooth_on(options)) {
    terms$smooth_weights <- smooth_weight_terms(v, options)
    terms$floor_intercept <- v$floor_intercept
    terms$init_loading <- smooth_init_loadings(v, options)
  }
  terms
}

# The logit relative initial state of each country at covariates 0 before
# the smooth of the initial state (smooth_options(init = TRUE)), as
# n_countries x n_types: logit_init_mean of each type, the same in every
# country, the number of countries that of the rows of
# options$init_covariate_centre (build_dynamical_model()). One row per
# country, as with the hierarchy, so that the predictions index it by country
# either way; a cell adds lambda u_init(x) and its covariates' effects. For
# greta arrays or plain R
init_field_level <- function(logit_init_mean, options) {
  n_countries <- nrow(options$init_covariate_centre)
  stopifnot(!is.null(n_countries))
  if (inherits(logit_init_mean, "greta_array")) {
    zero <- zeros(n_countries, length(logit_init_mean))
  } else {
    logit_init_mean <- c(logit_init_mean)
    zero <- matrix(0, n_countries, length(logit_init_mean))
  }
  sweep(zero, 2, logit_init_mean, FUN = "+")
}

# The fixed rho of each of `types` with the weighted binomial likelihood, as
# a plain vector: that build_dynamical_model() recorded in the options, or if
# it has not, the replicate rho (replicate_rho(), R/weighted_binomial.R)
fixed_rho_types <- function(options, types) {
  rho <- options$replicate_rho
  if (is.null(rho)) {
    rho <- replicate_rho(types)
  }
  stopifnot(all(types %in% names(rho)))
  unname(rho[types])
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
# Cached variables the model does not have are left out: e.g. the
# overdispersion parameters of every cached file (inits_refit.RDS,
# inits_floor_low.RDS, inits_floor_high.RDS) for a model with the weighted
# binomial likelihood, which has none.
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
  # assays and stalled a short run at its initial values. The species
  # multipliers start at no difference between the species, and both species'
  # floors, and the kdr-dependent floor at the mean kdr, at the cached
  # mortality_floor if there is one (so that cached values in either floor
  # mode, R/floor_mode_inits.R, start there), or else at 0.02 (#47). The
  # latent smooths start flat (raw weights 0), with sd 0.3 and range 3,000
  # km (with the default priors, the means of sd and 1 / range are 0.33 and
  # 1 / 4,500 km; a fixed range is no variable, so has no start), the shear
  # loading at its prior mean, 1, and the floor where they are 0 at the
  # cached floor too; the smooth of the initial state flat too, with every
  # type's loading 0.3 (its sd of the field), and logit_init_mean at the
  # cached value (the mean of the regions' levels in the model with the
  # hierarchy, whose region and country variables are left out)
  species_floor <- if (!is.null(cached$mortality_floor)) {
    c(cached$mortality_floor)[1]
  } else {
    0.02
  }
  starts <- c(mortality_floor = 0.02, init_coef = -0.05, reversion_rate = 0.01,
              gamma_selection = 0, gamma_cost = 0,
              other_floor = species_floor, arabiensis_floor = species_floor,
              setNames(rep(0, length(kdr_slope_names)), kdr_slope_names),
              floor_intercept = qlogis(species_floor), floor_kdr = 0,
              smooth_raw_selection = 0, smooth_sd_selection = 0.3,
              smooth_inv_range_selection = 1 / 3,
              smooth_raw_floor = 0, smooth_sd_floor = 0.3,
              smooth_inv_range_floor = 1 / 3, smooth_shear = 1,
              smooth_raw_init = 0, smooth_loading_init = 0.3,
              smooth_inv_range_init = 1 / 3)
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
# of the default model, from a 4-chain fit to the interpolation fold
# (September 2026, before reversion was added; reversion_rate starts at
# 0.01), with the columns of its selection design as attribute "columns" and
# the names of its types, classes, regions and countries as "levels". Not in
# git: remake it from a fit's draws as in fit_model.R.
dynamical_inits_file <- "temporary/inits_refit.RDS"

# The cached initial values the fits start from: the files in IR_CUBE_INITS,
# comma-separated, or dynamical_inits_file if it is unset or empty. With
# several, the chains are split between them (dynamical_chain_inits()), e.g.
# to start chains in both mortality-floor modes (#37; R/floor_mode_inits.R)
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
# for mortality_floor, other_floor and arabiensis_floor, and plogis() of the
# intercepts of the kdr-dependent floor or the floor of the latent smooths
# (floor_intercept, the floor at the mean kdr, or where u_f is 0; #47)
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
# susceptible) with the mortality floor applied, if there is one.
#
# With the species model (options$species, #47), mortality is computed for
# arabiensis and for the other members of the complex (outer_mortality()),
# and mortality(rows) is their mixture a pA + (1 - a) pG at each row's
# arabiensis share a (arabiensis_share(), R/species.R; the rows need species
# and cell), which is the mean of the beta-binomial. all_states(species_mix)
# is the mixture at the arabiensis fraction r(x) of each cell ("complex", the
# default), or the mortality of "arabiensis" or "other" alone. With the kdr
# covariate (options$kdr), mortality is computed by outer_mortality() too, at
# each cell's standardised kdr, and with the latent smooths (options$smooth,
# V5), with each cell's smooths; with the smooth of the initial state, each
# cell-type pair's initial state is shifted by lambda[type] u_init(x).
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
  species <- species_on(options)
  # the outer form (outer_mortality()) for the species model, the kdr
  # covariate or the latent smooths, and the closed-form op otherwise
  outer <- species || kdr_on(options) || smooth_on(options)
  # the mask cell of each cell_id
  cells <- df$cell[match(seq_len(n_unique_cells), df$cell_id)]
  # the standardised kdr at each cell_id (NULL without it), standardised over
  # these cells, all of them whatever the fold, recorded in the options for
  # the plain-R predictions
  if (kdr_on(options)) {
    options$kdr <- standardise_kdr(options$kdr, cells)
  }
  kdr_cells <- prediction_kdr(options, cells)
  # with the kdr-dependent floor by class, whether each class has the kdr
  # term, recorded in the options for the plain-R predictions
  if (identical(kdr_floor(options), "class")) {
    options$kdr$floor_classes <- lookups$levels$classes %in% kdr_floor_classes
  }
  # the basis of the latent smooths (V5), set up for these cells (all of them,
  # whatever the fold) and recorded in the options for the plain-R
  # predictions, and the centred basis at each cell_id (NULL without them)
  if (smooth_on(options)) {
    options$smooth <- smooth_box(options$smooth, cells,
                                 lookups$levels$classes)
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
  # the floor at rows with cell_id `cell` and type_id `type`: the
  # kdr-dependent floor, or the floor of the latent smooths, or NULL without
  # either
  row_floor <- function(cell, type) {
    if (smooth_floor_on(options)) {
      class <- classes_index[type]
      return(smooth_floor_value(
        terms$floor_intercept[smooth_intercept_index(options, class)],
        row_smooth("floor", cell, class)))
    }
    if (isFALSE(kdr_floor(options))) {
      return(NULL)
    }
    k <- kdr_cells[cell, "complex"]
    class <- rep(1L, length(type))
    if (identical(kdr_floor(options), "class")) {
      class <- classes_index[type]
      k <- k * options$kdr$floor_classes[class]
    }
    kdr_floor_value(terms$floor_intercept[class], terms$floor_kdr, k)
  }
  # the smooth of the initial state at cell-type pairs with cell_ids `cell`
  # and type_ids `type`, lambda[type] u_init(cell), or NULL without it
  pair_init <- function(cell, type) {
    if (!smooth_init_on(options)) {
      return(NULL)
    }
    terms$init_loading[type] * smooth_cells$init[cell]
  }
  x_init <- select_init_covariates(x_cells_init, options, n_unique_cells)
  # the centred country levels are at each country's mean initial-state
  # covariates over its modelled cells (all of them, whatever the fold; the
  # overall mean for a country with none), recorded in the options for the
  # plain-R predictions. With the smooth of the initial state there are no
  # country levels, and the covariates are centred at their mean over all
  # the modelled cells, the same in every row (one per country, which also
  # gives init_field_level() the number of countries; no columns without
  # covariates): one centre, as a centre per country would put steps at the
  # borders into the initial state
  options$init_covariate_centre <- NULL
  if (!is.null(x_init)) {
    x_cells <- x_init[seq_len(n_unique_cells), , drop = FALSE]
    centre <- if (smooth_init_on(options)) {
      matrix(colMeans(x_cells), n_countries, ncol(x_cells), byrow = TRUE)
    } else {
      t(vapply(seq_len(n_countries), function(country) {
        rows <- lookups$cell_country_lookup == country
        if (any(rows)) colMeans(x_cells[rows, , drop = FALSE]) else
          colMeans(x_cells)
      }, numeric(ncol(x_cells))))
    }
    colnames(centre) <- colnames(x_cells)
    options$init_covariate_centre <- centre
  } else if (smooth_init_on(options)) {
    options$init_covariate_centre <- matrix(0, n_countries, 0)
  }
  # with the weighted binomial likelihood, the fixed rho of each type,
  # recorded in the options for the plain-R predictions
  options$replicate_rho <- NULL
  if (!rho_estimated(options)) {
    options$replicate_rho <- replicate_rho(types)
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
  # there; with scale "log", as list(log_p, log_not_p), log mortality and
  # log(1 - mortality) computed from the logit (floored_log_probs(),
  # mixture_log_probs()), for the weighted binomial likelihood
  mortality <- function(rows, scale = c("probability", "log")) {
    scale <- match.arg(scale)
    pairs <- distinct(tibble(cell_id = as.integer(rows$cell_id),
                             type_id = as.integer(rows$type_id)))
    pair_index <- match(paste(rows$cell_id, rows$type_id),
                        paste(pairs$cell_id, pairs$type_id))
    if (outer) {
      p <- outer_mortality(terms, x_cell_years, pairs$cell_id,
                           pairs$type_id, lookups$cell_country_lookup,
                           n_times, types, x_init,
                           row_pair = pair_index, row_year = rows$year_id,
                           row_kdr = kdr_cells[rows$cell_id, , drop = FALSE],
                           row_floor = row_floor(rows$cell_id, rows$type_id),
                           row_selection = row_smooth(
                             "selection", rows$cell_id,
                             classes_index[rows$type_id]),
                           pair_init = pair_init(pairs$cell_id,
                                                 pairs$type_id),
                           scale = scale)
      if (!species) {
        return(p)
      }
      share <- arabiensis_share(rows, options)
      if (scale == "log") {
        return(mixture_log_probs(share, p$arabiensis, p$other))
      }
      return(share * p$arabiensis + (1 - share) * p$other)
    }
    states <- closed_form_states(terms, x_cell_years, pairs$cell_id,
                                 pairs$type_id, lookups$cell_country_lookup,
                                 n_times, types, x_init,
                                 logit = scale == "log")
    rows_states <- states[cbind(pair_index, rows$year_id)]
    if (scale == "log") {
      return(floored_log_probs(rows_states, terms$mortality_floor))
    }
    floored_mortality(rows_states, terms$mortality_floor)
  }

  # the predicted mortality at every cell, type and year, as n_unique_cells x
  # n_types x n_times, for after sampling (created before model(), it would be
  # in the model's graph)
  all_states <- function(species_mix = c("complex", "arabiensis", "other")) {
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
    species_mix <- match.arg(species_mix)
    n_pairs <- nrow(pairs)
    row_pair <- rep(seq_len(n_pairs), n_times)
    row_cell <- pairs$cell_id[row_pair]
    row_type <- pairs$type_id[row_pair]
    p <- outer_mortality(terms, x_cell_years, pairs$cell_id, pairs$type_id,
                         lookups$cell_country_lookup, n_times, types,
                         x_init, row_pair = row_pair,
                         row_year = rep(seq_len(n_times), each = n_pairs),
                         row_kdr = kdr_cells[row_cell, , drop = FALSE],
                         row_floor = row_floor(row_cell, row_type),
                         row_selection = row_smooth(
                           "selection", row_cell, classes_index[row_type]),
                         pair_init = pair_init(pairs$cell_id, pairs$type_id))
    if (species) {
      p <- switch(species_mix,
                  arabiensis = p$arabiensis,
                  other = p$other,
                  complex = {
                    r <- prediction_share(options, cells)[row_cell]
                    r * p$arabiensis + (1 - r) * p$other
                  })
    }
    dim(p) <- c(n_unique_cells, n_types, n_times)
    p
  }

  # likelihood: the beta-binomial with the estimated rho, or the weighted
  # binomial (R/weighted_binomial.R) with each bioassay's design-effect
  # weight at the fixed replicate rho. population_mortality_vec, the
  # predicted mortality at each bioassay, is then exp(log p): one
  # computation of the states serves both
  if (rho_estimated(options)) {
    population_mortality_vec <- mortality(train_df)
    rho <- terms$rho_types[train_df$type_id]
    distribution(train_df$died) <- betabinomial_p_rho(
      N = train_df$mosquito_number,
      p = population_mortality_vec,
      rho = rho)
  } else {
    log_probs <- mortality(train_df, scale = "log")
    population_mortality_vec <- exp(log_probs$log_p)
    weight <- design_effect_weight(train_df$mosquito_number,
                                   terms$rho_types[train_df$type_id])
    distribution(train_df$died) <- weighted_binomial(
      size = train_df$mosquito_number,
      log_p = log_probs$log_p,
      log_not_p = log_probs$log_not_p,
      weight = weight)
  }

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
# effects of the type (without the smooth of the initial state, which
# pair_inputs() adds). For greta arrays (a column vector) or one draw in plain
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
#   logit                TRUE to return the logit of the fraction susceptible
# Returns (B, J, n_times): the fraction susceptible for each pair and year, or
# its logit.
tf_closed_form_states <- function(beta_type, logit_init, kappa_type = NULL,
                                  x_pairs, pair_type, logit = FALSE) {
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

  if (logit) logit_q else tf$sigmoid(logit_q)
}

# The species model (#47): the cumulative log fitness of each pair and year,
# sum_{s <= t} log w_s, as (B, J, n_times), from the same selection term as
# tf_closed_form_states(). Arguments as there.
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
# none). With logit = TRUE, the logit of the fraction susceptible, for the
# weighted binomial likelihood (floored_log_probs()).
closed_form_states <- function(terms, x_cell_years, pair_cell, pair_type,
                               cell_country_lookup, n_times, types,
                               x_init = NULL, logit = FALSE) {
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
           pair_type = as.integer(pair_type - 1),
           logit = logit),
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
# array). pair_init is the smooth of the initial state at each pair, lambda
# u_init (a J x 1 greta array; NULL without it), added to the logit relative
# initial state. Other arguments as closed_form_states().
pair_inputs <- function(terms, x_cell_years, pair_cell, pair_type,
                        cell_country_lookup, n_times, types, x_init = NULL,
                        pair_init = NULL) {
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
  if (!is.null(pair_init)) {
    l <- l + pair_init
  }
  logit_init <- floored_logit(
    l, init_frac_constants(types)$min[pair_type])

  list(x_pairs = x_pairs, logit_init = logit_init)
}

# Predicted bioassay mortality by the outer form (#47), for the species model,
# the kdr covariate or the latent smooths, at pair row_pair (of the cell-type
# pairs pair_cell, pair_type) and year index row_year of each row: with the
# species model, a list of two greta arrays, `other` (the other members of the
# complex) and `arabiensis`, and otherwise one, each with one element per row.
# row_kdr is the standardised kdr at each row's cell, rows x kdr_bands() (NULL
# without the kdr covariate), and row_selection the latent smooth of
# selection at each row (NULL without it); pair_init the smooth of the
# initial state at each pair (pair_inputs(); NULL without it). Other
# arguments as closed_form_states().
#
# Each trajectory multiplies the cumulative log fitness C_t and the reversion
# of the closed form by its own factors:
#   logit q_t = logit q_0 - exp(s) C_t - exp(c) t kappa,
# with log multipliers s and c of 0 for the other members of the complex, or
# the whole complex, and gamma_selection and gamma_cost for arabiensis, i.e.
# arabiensis's fitness is (1 + x' exp(beta))^exp(gamma_selection); plus, with
# the kdr covariate, delta_selection k(x) and delta_cost k(x) for the
# trajectory's kdr k(x) at the cell and its own slopes (delta_selection_other
# and so on with the species model). The multipliers are outside the log so
# that one cumulative sum C serves every trajectory and cell, and each costs
# only these few operations at each row; where x' exp(beta) is small,
# exp(s) log(1 + x' exp(beta)) is close to log(1 + exp(s) x' exp(beta)), a
# multiplier on the selection effects. With every log multiplier 0 this is
# the closed form. With the latent smooths (V5), the log multiplier of the
# one trajectory's cumulative log fitness is the smooth of selection u_s(x)
# at the row, row_selection. Each species has its own mortality floor,
# other_floor and arabiensis_floor (none without species_options(floors =
# TRUE)); one trajectory has mortality_floor, or row_floor, the kdr-dependent
# floor (kdr_floor_value()) or that of the latent smooths
# (smooth_floor_value()) at each row, if given. With scale "log", each
# trajectory is the list(log_p, log_not_p) of floored_log_probs(), log
# mortality and log(1 - mortality), for the weighted binomial likelihood.
outer_mortality <- function(terms, x_cell_years, pair_cell, pair_type,
                            cell_country_lookup, n_times, types,
                            x_init = NULL, row_pair, row_year,
                            row_kdr = NULL, row_floor = NULL,
                            row_selection = NULL, pair_init = NULL,
                            scale = c("probability", "log")) {
  scale <- match.arg(scale)
  inputs <- pair_inputs(terms, x_cell_years, pair_cell, pair_type,
                        cell_country_lookup, n_times, types, x_init,
                        pair_init)

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
  # a trajectory with species offsets gamma (NULL for none), the kdr of
  # `band` and the slopes named with `suffix`, and `floor`
  trajectory <- function(gamma_selection, gamma_cost, band, suffix, floor) {
    k <- if (!is.null(row_kdr)) row_kdr[, band]
    logit_q <- outer_logit(
      rows,
      log_selection = outer_log_multiplier(
        gamma_selection, terms[[paste0("delta_selection", suffix)]], k,
        row_selection),
      log_cost = outer_log_multiplier(
        gamma_cost, terms[[paste0("delta_cost", suffix)]], k))
    if (scale == "log") {
      return(floored_log_probs(logit_q, floor))
    }
    floored_mortality(ilogit(logit_q), floor)
  }
  if (is.null(terms$gamma_selection)) {
    return(trajectory(NULL, NULL, "complex", "",
                      if (is.null(row_floor)) terms$mortality_floor else
                        row_floor))
  }
  list(other = trajectory(NULL, NULL, "other", "_other", terms$other_floor),
       arabiensis = trajectory(terms$gamma_selection, terms$gamma_cost,
                               "arabiensis", "_arabiensis",
                               terms$arabiensis_floor))
}

# The outer form's logit q: logit q_0 - exp(log_selection) C -
# exp(log_cost) reversion, from `rows`, a list of logit_init, cumulative and
# reversion (NULL for none), and the log multipliers (NULL for none, a
# multiplier of 1). For greta arrays (one element per row) or plain R (draws
# x cells, log multipliers one per draw or draws x cells, and reversion one
# per draw).
outer_logit <- function(rows, log_selection = NULL, log_cost = NULL) {
  selection <- rows$cumulative
  if (!is.null(log_selection)) {
    selection <- exp(log_selection) * selection
  }
  logit_q <- rows$logit_init - selection
  if (!is.null(rows$reversion)) {
    reversion <- rows$reversion
    if (!is.null(log_cost)) {
      reversion <- exp(log_cost) * reversion
    }
    logit_q <- logit_q - reversion
  }
  logit_q
}

# A trajectory's log multiplier: its species offset gamma (NULL for none) plus
# the kdr slope delta times its kdr k (both NULL without the kdr covariate),
# plus the latent smooth (NULL without it), or NULL if there is none. For
# greta arrays (k and smooth one per row) or plain R (gamma and delta one per
# draw, k one per cell, smooth draws x cells, giving draws x cells)
outer_log_multiplier <- function(gamma, delta, k, smooth = NULL) {
  kdr <- if (!is.null(delta)) {
    if (inherits(delta, "greta_array")) delta * k else outer(delta, k)
  }
  out <- if (is.null(gamma)) kdr else if (is.null(kdr)) gamma else gamma + kdr
  if (is.null(smooth)) out else if (is.null(out)) smooth else out + smooth
}
