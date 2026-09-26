# Do the static pixel effect p and the survey effect s take over the omega
# hotspots? (#21; doc/two_stage_plan.md, "Error structure")
#
#   Rscript R/two_stage_error_structure_hotspots.R [types]
#
# Fits omega_xi_u to all of a type's data (as R/two_stage_maps.R does, on the
# adopted omega5000_xi2500 meshes, t0 = 1995, T = the last data year), stage A
# then stage B, without and with the optional terms (+p, +p+s), and applies the
# hotspot definition of R/two_stage_hotspot_diagnostics.R to each fit:
#
#   omega hotspots: omega-mesh nodes that are strict local extrema of omega-hat
#   among their mesh neighbours, near the data (an assay within 150 km), with
#   |omega-hat| in the top 2% of the near-data nodes;
#   map hotspots: the same for the 2020 correction omega + xi at the omega
#   nodes, thresholded on its prominence over the 2-ring.
#
# Because the 2% threshold is relative, the count at it changes little by
# construction; the count above the threshold of the fit without p and s (a
# fixed |omega-hat|), and the size of the peaks, say whether the hotspots fall.
#
# m (the dynamical model's posterior mean logit at each assay) is read from the
# cache of R/two_stage_covariate_diagnostics.R, as in the hotspot diagnostics.
# Writes outputs/two_stage/error_structure_hotspots.csv.

suppressMessages({
  sink("/dev/null")
  source("R/validation_folds.R")
  sink()
})
source("R/two_stage_correction.R")
source("R/two_stage_pql.R")

arguments <- commandArgs(trailingOnly = TRUE)
focus_types <- if (length(arguments)) strsplit(arguments[1], ",")[[1]] else
  c("Deltamethrin", "Permethrin")
output_file <- "outputs/two_stage/error_structure_hotspots.csv"
hotspot_quantile <- 0.98
data_region_km <- 150
map_year <- 2020

report <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", sprintf(...), "\n")
  flush(stdout())
}

dynamical_file <- "outputs/two_stage/covariate_diagnostics/dynamical_components.rds"
m_ref <- readRDS(dynamical_file)$m_ref
stopifnot(length(m_ref) == nrow(df))
rho_table <- read.csv("outputs/bioassay_rho_hierarchical.csv")

# mesh neighbours and the 2-ring, as in R/two_stage_hotspot_diagnostics.R
mesh_neighbours <- function(mesh) {
  tv <- mesh$graph$tv
  edges <- rbind(tv[, 1:2], tv[, 2:3], tv[, c(1, 3)])
  n <- mesh$n
  adjacency <- Matrix::sparseMatrix(i = c(edges[, 1], edges[, 2]),
                                    j = c(edges[, 2], edges[, 1]),
                                    x = 1, dims = c(n, n))
  adjacency <- (adjacency > 0) * 1
  ring2 <- as((adjacency + adjacency %*% adjacency) > 0, "dMatrix")
  Matrix::diag(ring2) <- 0
  adjacency <- as(adjacency, "CsparseMatrix")
  list(neighbours = split(adjacency@i + 1, rep(seq_len(n), diff(adjacency@p))),
       ring2 = Matrix::drop0(ring2))
}
local_extremum <- function(f, geometry) {
  vapply(seq_along(f), function(k) {
    nb <- geometry$neighbours[[as.character(k)]]
    if (length(nb) == 0) return(FALSE)
    if (f[k] > 0) all(f[k] > f[nb]) else all(f[k] < f[nb])
  }, logical(1))
}
prominence <- function(f, geometry) {
  f - as.numeric(geometry$ring2 %*% f) / Matrix::rowSums(geometry$ring2)
}

rows_out <- list()
for (type in focus_types) {
  k <- match(type, types)
  rows <- which(df$type_id == k)
  rho <- rho_table$rho[rho_table$insecticide_type == type]
  train <- tibble(lon = df$longitude[rows], lat = df$latitude[rows],
                  year = df$year_start[rows], cell = df$cell[rows],
                  died = df$died[rows],
                  mosquito_number = df$mosquito_number[rows],
                  m = m_ref[rows], rho = rho)
  stage_a <- empirical_logit(train$died, train$mosquito_number, rho)
  train$z <- stage_a$z
  train$v <- stage_a$v
  coords <- coords_km(train)
  train$survey <- survey_id(df$citation[rows], df$country_name[rows],
                            train$year, coords)
  T_k <- max(train$year)
  meshes <- suppressMessages(build_correction_meshes(coords,
                                                     "omega5000_xi2500"))
  geometry <- mesh_neighbours(meshes$omega)
  loc <- meshes$omega$loc[, 1:2]
  near_data <- vapply(seq_len(meshes$omega$n), function(j) {
    any((coords[, 1] - loc[j, 1]) ^ 2 + (coords[, 2] - loc[j, 2]) ^ 2 <
          data_region_km ^ 2)
  }, logical(1))

  base_threshold <- NULL
  for (terms in list(character(0), "p", c("p", "s"))) {
    label <- paste(c("omega_xi_u", terms), collapse = "_")
    set.seed(2026 + k)
    time <- system.time({
      fit_a <- fit_correction(train, "omega_xi_u", t0 = baseline_year,
                              T = T_k, mesh = meshes$omega,
                              mesh_xi = meshes$xi,
                              pixel_effect = "p" %in% terms,
                              survey_effect = "s" %in% terms)
      fit_b <- fit_correction_pql(fit_a, train)
    })[["elapsed"]]
    for (stage in c("A", "B")) {
      fit <- if (stage == "A") fit_a else fit_b
      w <- fit$mode[fit$blocks$w_omega]
      threshold <- quantile(abs(w[near_data]), hotspot_quantile)
      if (stage == "B" && is.null(base_threshold)) base_threshold <- threshold
      peak <- local_extremum(w, geometry) & near_data
      hotspots <- which(peak & abs(w) >= threshold)
      # the 2020 correction omega + xi at the omega nodes
      year_col <- map_year - fit$t0
      x <- fit$mode[fit$blocks$x]
      xi_nodes <- if (year_col <= fit$n_years) {
        as.numeric(mesh_basis(fit$mesh_xi, loc) %*%
                     x[(year_col - 1) * fit$mesh_xi$n +
                         seq_len(fit$mesh_xi$n)])
      } else rep(0, meshes$omega$n)
      correction <- w + xi_nodes
      prom <- prominence(correction, geometry)
      prom_threshold <- quantile(abs(prom[near_data]), hotspot_quantile)
      map_peaks <- which(local_extremum(correction, geometry) & near_data &
                           abs(prom) >= prom_threshold)
      omega_at_data <- as.numeric(mesh_basis(fit$mesh, coords) %*% w)
      rows_out[[length(rows_out) + 1]] <- tibble(
        insecticide_type = type, model = label, stage = stage,
        convergence = fit$opt$convergence,
        range_omega = fit$hyper$range_omega,
        sigma_omega = fit$hyper$sigma_omega, tau = fit$hyper$tau,
        sigma_p = if (is.null(fit$hyper$sigma_p)) NA else fit$hyper$sigma_p,
        sigma_s = if (is.null(fit$hyper$sigma_s)) NA else fit$hyper$sigma_s,
        range_eta = fit$hyper$range_eta, sigma_eta = fit$hyper$sigma_eta,
        sd_omega_hat_at_data = sd(omega_at_data),
        threshold = threshold,
        n_hotspots = length(hotspots),
        # at the fixed |omega-hat| threshold of the stage-B fit without p, s
        n_above_base_threshold = sum(peak & abs(w) >=
                                       if (is.null(base_threshold)) threshold
                                     else base_threshold),
        median_abs_hotspot = median(abs(w[hotspots])),
        max_abs_omega = max(abs(w[near_data])),
        n_map_hotspots = length(map_peaks),
        map_prominence_threshold = prom_threshold,
        median_map_prominence = median(abs(prom[map_peaks])),
        seconds = time
      )
      report("%-13s %-16s %s conv=%i range_omega=%.0f sigma_omega=%.2f tau=%.3f sigma_p=%s sigma_s=%s | hotspots %i (threshold %.2f), above base %i, median |w| %.2f; map hotspots %i (prominence threshold %.2f) | %.0f s",
             type, label, stage, fit$opt$convergence, fit$hyper$range_omega,
             fit$hyper$sigma_omega, fit$hyper$tau,
             format(fit$hyper$sigma_p, digits = 3),
             format(fit$hyper$sigma_s, digits = 3), length(hotspots),
             threshold, tail(rows_out, 1)[[1]]$n_above_base_threshold,
             median(abs(w[hotspots])), length(map_peaks), prom_threshold,
             time)
    }
    rm(fit_a, fit_b, fit)
    invisible(gc())
  }
  write.csv(bind_rows(rows_out), output_file, row.names = FALSE)
}
report("wrote %s", output_file)
