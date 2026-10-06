# functions for the joint model of the Vgsc-995F and 995S allele frequencies
# (#46), used by R/kdr_joint.R

# The alleles at codon 995 are L (wild type, the reference), F and S, with
# frequencies p = softmax(0, eta_F, eta_S), so that p_F + p_S <= 1. Each
# variant's linear predictor is
#
#   eta_v = alpha_v + s_v(longitude, latitude) + beta_v * year_c
#           [+ delta_v for An. arabiensis]
#           [+ an An. arabiensis deviation smooth]
#
# where s_v is the thin plate (k = 80) chosen by cross-validation for that
# variant (R/kdr_cv.R, R/kdr_cv_extrap.R), with its smoothing parameters
# fixed: 995F range space 10^1, null space 10^10; 995S range space 10^2, null
# space unpenalised. eta_v is the log ratio of the frequencies of v and L.
#
# Records typing one variant are binomial (alleles of that variant out of
# 2n); samples typing both are trinomial (F, S and L alleles). Samples of the
# complex not identified to species can be a mixture of An. arabiensis and the
# rest of the complex, in the proportions r(x) of the Vector Atlas arabiensis
# fraction layer. The model is fitted by penalised maximum likelihood (Fisher
# scoring), with the penalties on the deviance scale as in mgcv, so that with
# one variant it is the same fit as gam().

joint_variants <- c(F = "995F", S = "995S")

# the smoothing parameters of the spatial smooths: range space, null space
# (NA for unpenalised)
joint_sp <- list(F = c(10 ^ 1, 10 ^ 10), S = c(10 ^ 2, NA))

# observations -----------------------------------------------------------------

# One row per binomial record or trinomial sample, with allele counts y_F and
# y_S out of m = 2n (y_S NA for a 995F-only record and vice versa). Samples
# with one record of each variant are trinomial, unless their counts sum to
# more than m by more than one (rounding), when the records are kept as
# separate binomials and flagged; samples with more than one record of a
# variant are also kept as separate binomials.
joint_observations <- function(records) {
  alleles <- records %>%
    filter(variant %in% joint_variants, !is.na(n_tested), n_tested > 0) %>%
    mutate(
      code = names(joint_variants)[match(variant, joint_variants)],
      m = 2 * n_tested,
      y = round(frequency * m)
    ) %>%
    group_by(sample_id) %>%
    mutate(
      records_F = sum(code == "F"),
      records_S = sum(code == "S")
    ) %>%
    ungroup()

  paired <- alleles %>%
    filter(records_F == 1, records_S == 1, code == "F") %>%
    select(-y, -record_id, -code, -variant, -frequency) %>%
    left_join(alleles %>% filter(code == "F") %>%
                select(sample_id, y_F = y, record_id_F = record_id),
              by = "sample_id") %>%
    left_join(alleles %>% filter(code == "S") %>%
                select(sample_id, y_S = y, record_id_S = record_id),
              by = "sample_id") %>%
    mutate(
      excess = y_F + y_S - m,
      # a count over by one is rounding
      y_S = if_else(excess == 1 & y_S >= y_F, y_S - 1, y_S),
      y_F = if_else(excess == 1 & y_F > y_S + 1, y_F - 1, y_F),
      consistent = y_F + y_S <= m
    )

  trinomial <- paired %>%
    filter(consistent) %>%
    mutate(type = "trinomial",
           record_id = paste(record_id_F, record_id_S))
  binomial <- alleles %>%
    filter(!(sample_id %in% trinomial$sample_id)) %>%
    mutate(
      type = paste0("binomial_", code),
      y_F = if_else(code == "F", y, NA),
      y_S = if_else(code == "S", y, NA),
      flag_inconsistent = sample_id %in% paired$sample_id
    )

  bind_rows(trinomial, binomial) %>%
    mutate(
      inconsistent_pair = coalesce(flag_inconsistent, FALSE),
      types_F = type %in% c("trinomial", "binomial_F"),
      types_S = type %in% c("trinomial", "binomial_S"),
      species_type = case_when(
        species_group == "arabiensis" ~ "arabiensis",
        species_group == "gambiae_complex" ~ "complex",
        .default = "other"
      )
    ) %>%
    select(record_id, sample_id, type, y_F, y_S, m, longitude, latitude,
           year_start, year_c, species_group, species_type, study, country,
           types_F, types_S, inconsistent_pair)
}

# the An. arabiensis fraction of the complex at some locations, from the
# Vector Atlas layer, filling cells outside it from the layer averaged over
# half a degree, and then the median over the locations
arabiensis_fraction_at <- function(longitude, latitude) {
  layer <- rast("data/clean/arabiensis_fraction.tif")
  points <- cbind(longitude, latitude)
  r <- terra::extract(layer, points)[[1]]
  coarse <- terra::aggregate(layer, 12, mean, na.rm = TRUE)
  r <- coalesce(r, terra::extract(coarse, points)[[1]])
  coalesce(r, median(r, na.rm = TRUE))
}

# bases and design matrices ----------------------------------------------------

# the thin-plate smooth of a variant, set up on the locations of the
# observations typing it, as gam() does for the separate fits
joint_smooth <- function(obs, code) {
  locations <- obs[obs[[paste0("types_", code)]],
                   c("longitude", "latitude")]
  smoothCon(s(longitude, latitude, bs = "tp", k = 80),
            data = as.data.frame(locations),
            absorb.cons = TRUE,
            null.space.penalty = !is.na(joint_sp[[code]][2]))[[1]]
}

# A candidate model: whether and how the species enter. `species` is "none",
# "own" (offsets for An. arabiensis and for samples identified only to the
# complex) or "mixture" (an An. arabiensis offset, with samples identified
# only to the complex a mixture by r); `deviation` is NULL, or the factor by
# which an An. arabiensis deviation smooth's smoothing parameters exceed the
# variant's own.
joint_candidate <- function(species = "none", deviation = NULL) {
  list(species = species, deviation = deviation)
}

# The design matrices of a variant at some observations: for the rest of the
# complex (g) and for An. arabiensis (a), and the names of the coefficients.
joint_design <- function(obs, smooth, candidate) {
  basis <- PredictMat(smooth, as.data.frame(obs[, c("longitude",
                                                    "latitude")]))
  colnames(basis) <- sprintf("s%02i", seq_len(ncol(basis)))
  n <- nrow(obs)
  g <- cbind(intercept = 1, year = obs$year_c, basis)
  a <- g
  if (candidate$species != "none") {
    g <- cbind(g, arabiensis = 0)
    a <- cbind(a, arabiensis = 1)
  }
  if (candidate$species == "own") {
    complex <- as.numeric(obs$species_type == "complex")
    g <- cbind(g, complex = complex)
    a <- cbind(a, complex = 0)
  }
  if (!is.null(candidate$deviation)) {
    deviation <- basis
    colnames(deviation) <- sprintf("arabiensis_%s", colnames(basis))
    g <- cbind(g, deviation * 0)
    a <- cbind(a, deviation)
  }
  list(g = g, a = a)
}

# the penalty matrix of a variant's coefficients, given the columns of its
# design and its smooth
joint_penalty <- function(columns, smooth, code, candidate) {
  p <- length(columns)
  penalty <- matrix(0, p, p, dimnames = list(columns, columns))
  sp <- joint_sp[[code]]
  smooth_penalty <- sp[1] * smooth$S[[1]]
  if (!is.na(sp[2])) {
    smooth_penalty <- smooth_penalty + sp[2] * smooth$S[[2]]
  }
  basis_columns <- grep("^s[0-9]+$", columns)
  penalty[basis_columns, basis_columns] <- smooth_penalty
  if (!is.null(candidate$deviation)) {
    # the deviation has its null space penalised too, so it shrinks to zero
    deviation_penalty <- candidate$deviation * sp[1] * smooth$S[[1]]
    deviation_penalty <- deviation_penalty +
      if (length(smooth$S) > 1) {
        candidate$deviation * sp[1] * smooth$S[[2]]
      } else {
        candidate$deviation * sp[1] * null_space_penalty(smooth)
      }
    deviation_columns <- grep("^arabiensis_s", columns)
    penalty[deviation_columns, deviation_columns] <- deviation_penalty
  }
  penalty
}

# a penalty on the null space of a smooth's (single) penalty, scaled like it
null_space_penalty <- function(smooth) {
  eigen_s <- eigen(smooth$S[[1]], symmetric = TRUE)
  null <- eigen_s$values < max(eigen_s$values) * 1e-8
  u <- eigen_s$vectors[, null, drop = FALSE]
  penalty <- tcrossprod(u)
  penalty * norm(smooth$S[[1]], "F") / norm(penalty, "F")
}

# the mixture weights of the An. arabiensis component of each observation
arabiensis_weight <- function(obs, candidate) {
  case_when(
    candidate$species == "none" ~ 0,
    obs$species_type == "arabiensis" ~ 1,
    obs$species_type == "complex" & candidate$species == "mixture" ~ obs$r,
    .default = 0
  )
}

# fitting ----------------------------------------------------------------------

# allele frequencies (L, F, S) from linear predictors
softmax3 <- function(eta_F, eta_S) {
  top <- pmax(0, eta_F, eta_S)
  e_L <- exp(-top)
  e_F <- exp(eta_F - top)
  e_S <- exp(eta_S - top)
  total <- e_L + e_F + e_S
  list(L = e_L / total, F = e_F / total, S = e_S / total)
}

# the frequencies, and their derivatives with respect to the coefficients of
# each variant, at some observations
joint_frequencies <- function(theta, designs, weight_a,
                              derivatives = TRUE) {
  codes <- names(designs)
  eta <- function(code, part) {
    if (code %in% codes) {
      as.numeric(designs[[code]][[part]] %*% theta[[code]])
    } else {
      rep(-Inf, length(weight_a))
    }
  }
  components <- list(
    g = softmax3(eta("F", "g"), eta("S", "g")),
    a = softmax3(eta("F", "a"), eta("S", "a"))
  )
  weights <- list(g = 1 - weight_a, a = weight_a)
  p <- map(c(L = "L", F = "F", S = "S"), \(k) {
    weights$g * components$g[[k]] + weights$a * components$a[[k]]
  })
  if (!derivatives) {
    return(list(p = p))
  }
  # dp_k / d theta, as a matrix with a row per observation and the
  # coefficients of each variant in turn
  gradient <- map(c(L = "L", F = "F", S = "S"), \(k) {
    map(setNames(codes, codes), \(v) {
      Reduce(`+`, map(c("g", "a"), \(part) {
        q <- components[[part]]
        d <- if (k == v) q[[k]] * (1 - q[[v]]) else -q[[k]] * q[[v]]
        designs[[v]][[part]] * (weights[[part]] * d)
      }))
    }) %>%
      do.call(cbind, .)
  })
  list(p = p, gradient = gradient)
}

# the log-likelihood terms (u, the score weights; w, the information
# weights) of each observation, given the frequencies
joint_weights <- function(obs, p) {
  tri <- obs$type == "trinomial"
  bin_F <- obs$type == "binomial_F"
  bin_S <- obs$type == "binomial_S"
  m <- obs$m
  y_F <- coalesce(obs$y_F, 0)
  y_S <- coalesce(obs$y_S, 0)
  y_L <- m - y_F - y_S
  zero <- numeric(nrow(obs))
  u <- list(
    L = if_else(tri, y_L / p$L, zero),
    F = case_when(tri ~ y_F / p$F,
                  bin_F ~ y_F / p$F - (m - y_F) / (1 - p$F),
                  .default = 0),
    S = case_when(tri ~ y_S / p$S,
                  bin_S ~ y_S / p$S - (m - y_S) / (1 - p$S),
                  .default = 0)
  )
  w <- list(
    L = if_else(tri, m / p$L, zero),
    F = case_when(tri ~ m / p$F,
                  bin_F ~ m / (p$F * (1 - p$F)),
                  .default = 0),
    S = case_when(tri ~ m / p$S,
                  bin_S ~ m / (p$S * (1 - p$S)),
                  .default = 0)
  )
  list(u = u, w = w)
}

# the log-likelihood of each observation, without its constant
joint_loglik <- function(obs, p) {
  xlogp <- function(y, prob) if_else(y > 0, y * log(prob), 0)
  m <- obs$m
  y_F <- coalesce(obs$y_F, 0)
  y_S <- coalesce(obs$y_S, 0)
  case_when(
    obs$type == "trinomial" ~ xlogp(y_F, p$F) + xlogp(y_S, p$S) +
      xlogp(m - y_F - y_S, p$L),
    obs$type == "binomial_F" ~ xlogp(y_F, p$F) + xlogp(m - y_F, 1 - p$F),
    obs$type == "binomial_S" ~ xlogp(y_S, p$S) + xlogp(m - y_S, 1 - p$S)
  )
}

# the deviance of each observation: twice its log-likelihood under the
# observed frequencies less that under the model
joint_deviance <- function(obs, p) {
  m <- obs$m
  y_F <- coalesce(obs$y_F, 0)
  y_S <- coalesce(obs$y_S, 0)
  saturated <- list(L = (m - y_F - y_S) / m, F = y_F / m, S = y_S / m)
  p <- map(p, \(x) pmin(pmax(x, 1e-12), 1 - 1e-12))
  2 * (joint_loglik(obs, saturated) - joint_loglik(obs, p))
}

# Fit a candidate to some observations by penalised maximum likelihood
# (Fisher scoring with step halving). Returns the coefficients, their
# covariance (unscaled, and scaled by the Pearson estimate of
# overdispersion), the smooths and the edf.
fit_joint <- function(obs, candidate, codes = c("F", "S"),
                      max_iterations = 200, tolerance = 1e-4) {
  if (!"r" %in% names(obs)) {
    obs$r <- NA_real_
  }
  codes <- setNames(codes, codes)
  smooths <- map(codes, \(code) joint_smooth(obs, code))
  designs <- map(codes, \(code) {
    joint_design(obs, smooths[[code]], candidate)
  })
  columns <- map(designs, \(d) colnames(d$g))
  penalty <- Matrix::bdiag(
    imap(columns, \(cols, code) {
      joint_penalty(cols, smooths[[code]], code, candidate)
    })
  ) %>%
    as.matrix()
  weight_a <- arabiensis_weight(obs, candidate)
  sizes <- map_int(columns, length)
  split_theta <- function(theta) {
    split(theta, rep(names(sizes), sizes))[names(sizes)]
  }

  # start from the pooled frequencies
  pooled <- c(
    F = sum(obs$y_F, na.rm = TRUE) / sum(obs$m[obs$types_F]),
    S = sum(obs$y_S, na.rm = TRUE) / sum(obs$m[obs$types_S])
  )
  pooled_L <- 1 - sum(pooled[codes])
  starts <- cumsum(c(1, sizes))[seq_along(sizes)]
  theta <- numeric(sum(sizes))
  theta[starts] <- log(pooled[codes] / pooled_L)

  objective <- function(theta) {
    p <- joint_frequencies(split_theta(theta), designs, weight_a,
                           derivatives = FALSE)$p
    sum(joint_loglik(obs, p)) -
      0.5 * as.numeric(crossprod(theta, penalty %*% theta))
  }
  information_score <- function(theta) {
    frequencies <- joint_frequencies(split_theta(theta), designs, weight_a)
    weights <- joint_weights(obs, frequencies$p)
    list(
      score = Reduce(`+`, map(c("L", "F", "S"), \(k) {
        crossprod(frequencies$gradient[[k]], weights$u[[k]])
      })),
      information = Reduce(`+`, map(c("L", "F", "S"), \(k) {
        crossprod(frequencies$gradient[[k]] * sqrt(weights$w[[k]]))
      })),
      p = frequencies$p
    )
  }

  # Fisher scoring converges linearly here (the expected information is not
  # the Hessian for observations typing one variant), so we stop when the
  # largest step in a coefficient is below `tolerance`, or when no step
  # improves the penalised log-likelihood
  current <- objective(theta)
  converged <- FALSE
  for (iteration in seq_len(max_iterations)) {
    terms <- information_score(theta)
    step <- solve(terms$information + penalty,
                  as.numeric(terms$score) - as.numeric(penalty %*% theta))
    size <- 1
    repeat {
      value <- objective(theta + size * step)
      if (is.finite(value) && value >= current) break
      size <- size / 2
      if (size < 1e-8) break
    }
    if (size < 1e-8) {
      converged <- max(abs(step)) < 10 * tolerance
      break
    }
    theta <- theta + size * step
    current <- value
    if (max(abs(size * step)) < tolerance) {
      converged <- TRUE
      break
    }
  }

  # covariance, edf and Pearson overdispersion at the estimate
  terms <- information_score(theta)
  covariance <- solve(terms$information + penalty)
  edf <- sum(diag(covariance %*% terms$information))
  residual_df <- sum(if_else(obs$type == "trinomial", 2, 1)) - edf
  scale <- sum(joint_pearson(obs, terms$p)) / residual_df
  names(theta) <- unlist(imap(columns, \(cols, code) {
    paste(code, cols, sep = ":")
  }))
  dimnames(covariance) <- list(names(theta), names(theta))
  list(theta = theta, split_theta = split_theta, covariance = covariance,
       scaled_covariance = covariance * scale, scale = scale, edf = edf,
       iterations = iteration, converged = converged,
       smooths = smooths, candidate = candidate, columns = columns,
       sizes = sizes, deviance = sum(joint_deviance(obs, terms$p)))
}

# the frequencies of a fitted model at some observations (with their species
# types and, for mixtures, r)
predict_observations <- function(fit, obs) {
  if (!"r" %in% names(obs)) {
    obs$r <- NA_real_
  }
  designs <- map(fit$smooths, \(smooth) {
    joint_design(obs, smooth, fit$candidate)
  })
  joint_frequencies(fit$split_theta(fit$theta), designs,
                    arabiensis_weight(obs, fit$candidate),
                    derivatives = FALSE)$p
}

# Pearson chi-squared of each observation
joint_pearson <- function(obs, p) {
  m <- obs$m
  y_F <- coalesce(obs$y_F, 0)
  y_S <- coalesce(obs$y_S, 0)
  chi <- function(y, prob) (y - m * prob) ^ 2 / (m * prob)
  case_when(
    obs$type == "trinomial" ~ chi(y_F, p$F) + chi(y_S, p$S) +
      chi(m - y_F - y_S, p$L),
    obs$type == "binomial_F" ~ chi(y_F, p$F) / (1 - p$F),
    obs$type == "binomial_S" ~ chi(y_S, p$S) / (1 - p$S)
  )
}

# prediction -------------------------------------------------------------------

# The frequencies of a fitted model at new locations (columns longitude,
# latitude, year_c), for the rest of the complex (g), for An. arabiensis (a),
# and for the complex as a mixture by `r` (if given), with the standard error
# of logit(total kdr) for each, using the scaled covariance.
predict_joint <- function(fit, new, r = NULL) {
  new <- as.data.frame(new)
  new$species_type <- "other"
  designs <- map(fit$smooths, \(smooth) {
    joint_design(new, smooth, fit$candidate)
  })
  theta <- fit$split_theta(fit$theta)
  parts <- list(g = rep(0, nrow(new)), a = rep(1, nrow(new)))
  if (!is.null(r)) {
    parts$complex <- r
  }
  map(parts, \(weight_a) {
    frequencies <- joint_frequencies(theta, designs, weight_a)
    p <- frequencies$p
    total <- p$F + p$S
    # d logit(total) / d theta = (dp_F + dp_S) / (total * (1 - total))
    gradient <- (frequencies$gradient$F + frequencies$gradient$S) /
      (total * (1 - total))
    se <- sqrt(rowSums((gradient %*% fit$scaled_covariance) * gradient))
    tibble(p_F = p$F, p_S = p$S, total = total, total_logit_se = se)
  })
}
