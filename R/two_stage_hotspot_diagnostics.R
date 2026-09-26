# What makes the hotspots in the second-stage correction maps (#21)?
#
#   Rscript R/two_stage_hotspot_diagnostics.R
#
# The correction maps (figures/two_stage/*_correction_map.png) show compact
# hotspots. The stage-A model is
#
#   z = m + omega(s) + xi(s, t) + u(pixel-year) + e,   e ~ N(0, v)
#
# with no static per-pixel term. Three candidate causes:
#
#   (A) mesh artefact: omega's fitted range (~40-50 km for the pyrethroids and
#       DDT on the adopted omega5000_xi2500 mesh) is about one mesh edge, so
#       the field is barely resolved; a local offset in the data is fitted by
#       pushing one node to an extreme value, and linear interpolation then
#       spreads it over that node's whole star (a triangular footprint whose
#       size is set by the local edge length, not the range);
#   (B) spatial clustering: many nearby pixels share a strong residual of the
#       same sign, a genuine short-range spatial signal;
#   (C) one heavily sampled pixel persistently off the dynamical prediction
#       across years; with no static per-pixel term, omega must absorb it.
#
# The script reads the per-type fits of R/two_stage_maps.R (the latent mode is
# saved in fit.rds) for the adopted meshes and, if present, for the base meshes
# (omega mesh <= 2500 nodes; a backup of the maps made before the mesh change,
# see base_map_dir), and needs no refitting:
#
#   1. hotspots: nodes of the omega mesh that are local extrema of omega-hat
#      (among their mesh neighbours) with |omega-hat| in the top 2% of nodes
#      near the data. Separately, local extrema of the 2020 correction
#      omega + xi (at the omega nodes) with a large local prominence, with the
#      prominence split into its omega and xi parts, to say which field makes
#      the hotspots seen on the maps;
#   2. mesh fingerprint of each hotspot: half-maximum footprint (area of the
#      connected region where omega-hat has the peak's sign and at least half
#      its magnitude, evaluated on a 2 km grid) against the local edge length,
#      the footprint a single node's hat function would have, and the footprint
#      a point source would have if the Matern field were resolved; the node
#      value against the precision-weighted mean of the partial residual
#      z - m - xi - u of the assays in the node's star (the assays that load on
#      it); and the same hotspot on the base mesh;
#   3. data fingerprint: assays, pixels, years and sites loading on the node;
#      the share of the node's data "pull" (A' D r, the term that sets the
#      mode) from its single largest pixel; the sign agreement across pixels;
#      and whether the dominant pixel's residual is consistent across years;
#   4. a classification of each hotspot (rule in classify_hotspot() below) and
#      zoomed example maps;
#   5. what a static per-pixel nugget would absorb: variance components
#      (static pixel, pixel-year) of the residual before and after omega-hat,
#      the covariance of residuals between distinct nearby pixels against that
#      between years at the same pixel, and how much of omega-hat at the data
#      is explained by the pixel's own mean residual vs its neighbours'.
#
# Outputs: outputs/two_stage/hotspot_diagnostics.csv (one row per hotspot and
# mesh), outputs/two_stage/hotspot_summary.csv (per type and mesh),
# outputs/two_stage/hotspot_nugget.csv (step 5), figures
# figures/two_stage/hotspot_*.png.
#
# m (the dynamical model's posterior mean logit at each assay) is read from the
# cache of R/two_stage_covariate_diagnostics.R, which recomputes it exactly as
# R/two_stage_maps.R did; the script checks it against the saved fits (the
# mode of u must satisfy its own stationarity condition given m).

suppressMessages({
  sink("/dev/null")
  source("R/validation_folds.R")
  sink()
  library(ggplot2)
  library(patchwork)
})
source("R/two_stage_correction.R")

map_dir <- "outputs/two_stage/maps"
# the base-mesh maps, backed up before the mesh change. Optional: without them
# the mesh comparison is skipped
base_map_dir <- Sys.getenv(
  "BASE_MAP_DIR",
  paste0("/tmp/claude-1000/-home-nick-Dropbox-github-ir-cube/",
         "25b39fd7-3494-419a-929b-b732b45e3827/scratchpad/base/maps"))
dynamical_file <- "outputs/two_stage/covariate_diagnostics/dynamical_components.rds"
output_dir <- "outputs/two_stage"
figure_dir <- "figures/two_stage"

focus_types <- c("Deltamethrin", "Permethrin", "Lambda-cyhalothrin", "DDT",
                 "Bendiocarb", "Malathion", "Pirimiphos-methyl")
example_types <- c("Deltamethrin", "Permethrin", "Lambda-cyhalothrin", "DDT",
                   "Bendiocarb")

# hotspot threshold: top 2% of |omega-hat| among nodes near the data
hotspot_quantile <- 0.98
# "near the data": an assay within this distance of the node (km)
data_region_km <- 150
map_year <- 2020
# fixed neighbourhood for the mesh-independent data fingerprint (km)
neighbourhood_km <- 50
grid_step_km <- 2

report <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", sprintf(...), "\n")
  flush(stdout())
}

stopifnot(file.exists(dynamical_file))
m_ref <- readRDS(dynamical_file)$m_ref
stopifnot(length(m_ref) == nrow(df))
rho_table <- read.csv("outputs/bioassay_rho_hierarchical.csv")

meshes_available <- c(fine = map_dir,
                      base = if (dir.exists(base_map_dir)) base_map_dir)


# helpers ----------------------------------------------------------------------

# the half-maximum radius of a Matern (nu = 1) correlation, as a multiple of the
# practical range sqrt(8) / kappa: the footprint a single well-resolved point
# source would leave in the posterior mean
matern_half_radius <- local({
  f <- function(d) d * besselK(d, 1) - 0.5
  uniroot(f, c(1e-3, 5))$root / sqrt(8)
})

# mesh geometry: edges, adjacency, node areas and local edge length
mesh_geometry <- function(mesh) {
  tv <- mesh$graph$tv
  edges <- rbind(tv[, 1:2], tv[, 2:3], tv[, c(1, 3)])
  edges <- unique(cbind(pmin(edges[, 1], edges[, 2]),
                        pmax(edges[, 1], edges[, 2])))
  loc <- mesh$loc[, 1:2]
  length_km <- sqrt(rowSums((loc[edges[, 1], ] - loc[edges[, 2], ]) ^ 2))
  n <- mesh$n
  adjacency <- Matrix::sparseMatrix(i = c(edges[, 1], edges[, 2]),
                                    j = c(edges[, 2], edges[, 1]),
                                    x = 1, dims = c(n, n))
  # mean length of the edges at each node: the local resolution
  edge_sum <- Matrix::sparseMatrix(i = c(edges[, 1], edges[, 2]),
                                   j = c(edges[, 2], edges[, 1]),
                                   x = rep(length_km, 2), dims = c(n, n))
  local_edge <- Matrix::rowSums(edge_sum) / Matrix::rowSums(adjacency)
  # the 2-ring (neighbours and their neighbours, not the node) for a local
  # background
  ring2 <- (adjacency + adjacency %*% adjacency) > 0
  ring2 <- as(ring2, "dMatrix")
  Matrix::diag(ring2) <- 0
  ring2 <- Matrix::drop0(ring2)
  fem <- correction_fem(mesh)
  list(loc = loc, edges = edges, adjacency = adjacency,
       neighbours = split(adjacency@i + 1,
                          rep(seq_len(n), diff(adjacency@p))),
       local_edge = local_edge, ring2 = ring2,
       node_area = Matrix::diag(fem$M0))
}

# is each node a strict local extremum of f among its mesh neighbours, in the
# direction of its own sign?
local_extremum <- function(f, geometry) {
  vapply(seq_along(f), function(k) {
    nb <- geometry$neighbours[[as.character(k)]]
    if (length(nb) == 0) return(FALSE)
    if (f[k] > 0) all(f[k] > f[nb]) else all(f[k] < f[nb])
  }, logical(1))
}

# f at a node minus its mean over the node's 2-ring
prominence <- function(f, geometry) {
  f - as.numeric(geometry$ring2 %*% f) / Matrix::rowSums(geometry$ring2)
}

# the connected half-maximum region of omega-hat around a peak, on a grid:
# its area (km^2) and the number of mesh nodes inside it
half_max_footprint <- function(peak_xy, peak_value, mesh, w, geometry,
                               radius_km) {
  xs <- seq(peak_xy[1] - radius_km, peak_xy[1] + radius_km, by = grid_step_km)
  ys <- seq(peak_xy[2] - radius_km, peak_xy[2] + radius_km, by = grid_step_km)
  grid <- as.matrix(expand.grid(x = xs, y = ys))
  values <- as.numeric(mesh_basis(mesh, grid) %*% w)
  inside <- sign(values) == sign(peak_value) &
    abs(values) >= abs(peak_value) / 2
  r <- terra::rast(nrows = length(ys), ncols = length(xs),
                   xmin = min(xs) - grid_step_km / 2,
                   xmax = max(xs) + grid_step_km / 2,
                   ymin = min(ys) - grid_step_km / 2,
                   ymax = max(ys) + grid_step_km / 2, crs = "")
  # expand.grid runs x fastest from the bottom row; terra fills from the top
  m <- matrix(as.numeric(inside), nrow = length(ys), byrow = TRUE)
  terra::values(r) <- as.numeric(t(m[nrow(m):1, ]))
  r[r == 0] <- NA
  patches <- terra::patches(r, directions = 8)
  centre <- terra::extract(patches, matrix(peak_xy, 1))[1, 1]
  if (is.na(centre)) return(c(area = NA, n_nodes = NA, touches_edge = NA))
  in_patch <- terra::values(patches)[, 1] %in% centre
  area <- sum(in_patch) * grid_step_km ^ 2
  cell_nodes <- terra::cellFromXY(patches, geometry$loc)
  nodes_in <- sum(!is.na(cell_nodes) &
                    terra::values(patches)[cell_nodes, 1] %in% centre)
  # a footprint that reaches the box edge is truncated
  rows_cols <- terra::rowColFromCell(patches, which(in_patch))
  touches <- any(rows_cols[, 1] %in% c(1, nrow(patches)) |
                   rows_cols[, 2] %in% c(1, ncol(patches)))
  c(area = area, n_nodes = nodes_in, touches_edge = touches)
}

# consistency of a pixel's residual across years: per-year precision-weighted
# means, and a heterogeneity test against their sampling SEs (plus tau, the
# pixel-year SD the model allows)
year_consistency <- function(r, v, year, sign_ref, tau) {
  by_year <- tibble(r = r, w = 1 / v, year = year) %>%
    group_by(year) %>%
    summarise(mean = sum(w * r) / sum(w), var = 1 / sum(w) + tau ^ 2,
              .groups = "drop")
  n_years <- nrow(by_year)
  pooled <- sum(by_year$mean / by_year$var) / sum(1 / by_year$var)
  q <- sum((by_year$mean - pooled) ^ 2 / by_year$var)
  tibble(
    dominant_n_years = n_years,
    dominant_year_sign_agreement = mean(sign(by_year$mean) == sign_ref),
    dominant_pooled_mean = pooled,
    dominant_pooled_z = pooled * sqrt(sum(1 / by_year$var)),
    dominant_heterogeneity_p = if (n_years > 1) {
      pchisq(q, n_years - 1, lower.tail = FALSE)
    } else NA_real_
  )
}

# The rule. Data cause first, from where the data pull on the peak comes from,
# over the assays within 50 km of it (a fixed radius, so that both meshes are
# judged on the same data; the star-based shares are also saved):
#   C  one pixel supplies >= 50% of the pull in the hotspot's direction;
#      "C persistent" if that pixel has >= 2 years and >= 75% of its year means
#      have the hotspot's sign, "C one year" otherwise;
#   B  otherwise, >= 3 pixels loading on the node and >= 70% of them agree in
#      sign with the hotspot;
#   mixed  anything else (few pixels, or pixels of mixed sign).
# Mesh flag A (not exclusive of B/C: it is about the hotspot's shape, not its
# cause): the half-maximum footprint is >= 2x the area a resolved point source
# would have (pi (0.57 range)^2 ... see matern_half_radius) and contains at most
# 2 nodes, i.e. its extent is set by the node's star rather than the range
classify_hotspot <- function(h) {
  cause <- case_when(
    h$largest_pixel_share_50km >= 0.5 & h$dominant_n_years >= 2 &
      h$dominant_year_sign_agreement >= 0.75 ~ "C persistent",
    h$largest_pixel_share_50km >= 0.5 ~ "C one year",
    h$n_pixels_50km >= 3 & h$pixel_sign_agreement_50km >= 0.7 ~ "B",
    TRUE ~ "mixed"
  )
  mesh_flag <- !is.na(h$footprint_area_km2) &
    h$footprint_area_km2 >= 2 * h$point_source_area_km2 &
    h$footprint_nodes <= 2
  tibble(cause = cause, mesh_artefact = mesh_flag)
}

# per-pixel variance components by maximum likelihood: pixel-year means
# r_py = mu + a_p + b_py + noise, a ~ N(0, s_a^2) static per pixel,
# b ~ N(0, s_b^2) per pixel-year, noise with known variance se2. Within a pixel
# the covariance is diag(s_b^2 + se2) + s_a^2 11', whose determinant and
# inverse have closed forms
pixel_variance_components <- function(r, se2, pixel) {
  groups <- split(seq_along(r), pixel)
  nll <- function(par) {
    mu <- par[1]
    s_a2 <- exp(par[2])
    s_b2 <- exp(par[3])
    total <- 0
    d <- s_b2 + se2
    e <- r - mu
    sum_inv_d <- vapply(groups, function(g) sum(1 / d[g]), numeric(1))
    sum_e_d <- vapply(groups, function(g) sum(e[g] / d[g]), numeric(1))
    log_det <- sum(log(d)) + sum(log1p(s_a2 * sum_inv_d))
    quad <- sum(e ^ 2 / d) - sum(s_a2 * sum_e_d ^ 2 / (1 + s_a2 * sum_inv_d))
    0.5 * (log_det + quad)
  }
  opt <- optim(c(weighted.mean(r, 1 / se2), log(0.1), log(0.1)), nll,
               method = "BFGS")
  c(sigma2_pixel = exp(opt$par[2]), sigma2_pixel_year = exp(opt$par[3]))
}


# per type and mesh --------------------------------------------------------------

hotspot_rows <- list()
map_hotspot_rows <- list()
nugget_rows <- list()
covariance_rows <- list()
example_store <- list()

for (type in focus_types) {

  k <- match(type, types)
  rows <- which(df$type_id == k)
  rho <- rho_table$rho[rho_table$insecticide_type == type]
  stage_a <- empirical_logit(df$died[rows], df$mosquito_number[rows], rho)
  train <- tibble(lon = df$longitude[rows], lat = df$latitude[rows],
                  year = df$year_start[rows], cell = df$cell[rows],
                  z = stage_a$z, v = stage_a$v, m = m_ref[rows])
  coords <- coords_km(train)
  # a site is a distinct coordinate pair; a pixel a raster cell (~4.6 km)
  train$site <- match(paste(coords[, 1], coords[, 2]),
                      unique(paste(coords[, 1], coords[, 2])))

  for (mesh_name in names(meshes_available)) {

    fit_file <- file.path(meshes_available[[mesh_name]], type, "fit.rds")
    if (!file.exists(fit_file)) next
    fit <- readRDS(fit_file)
    mesh <- fit$mesh
    geometry <- mesh_geometry(mesh)
    range_km <- fit$hyper$range_omega
    tau <- fit$hyper$tau

    design <- correction_design(mesh, fit$mesh_xi, coords, train$year,
                                fit$t0, fit$T)
    w <- fit$mode[fit$blocks$w_omega]
    x <- fit$mode[fit$blocks$x]
    omega_i <- as.numeric(design$A_omega %*% w)
    xi_i <- as.numeric(design$A_xi %*% x)
    u_key <- match(paste(train$cell, train$year),
                   paste(fit$pixel_years$cell, fit$pixel_years$year))
    u_i <- fit$mode[fit$blocks$u][u_key]

    # m must be the m of the fit: at the mode, each u solves
    # u = sum((z - m - omega - xi) / v) / (sum(1 / v) + 1 / tau^2)
    u_check <- tapply((train$z - train$m - omega_i - xi_i) / train$v, u_key,
                      sum) /
      (tapply(1 / train$v, u_key, sum) + 1 / tau ^ 2)
    u_error <- max(abs(u_check - fit$mode[fit$blocks$u][
      as.integer(names(u_check))]))
    stopifnot(u_error < 1e-8)

    # the partial residual omega is fitting, and the residual after omega
    r_omega <- train$z - train$m - xi_i - u_i
    r0 <- train$z - train$m - xi_i
    r1 <- r0 - omega_i

    # nodes near the data define the threshold, so that the zeros of the
    # outer extension do not
    near_data <- vapply(seq_len(mesh$n), function(j) {
      any(abs(coords[, 1] - geometry$loc[j, 1]) < data_region_km &
            (coords[, 1] - geometry$loc[j, 1]) ^ 2 +
            (coords[, 2] - geometry$loc[j, 2]) ^ 2 < data_region_km ^ 2)
    }, logical(1))
    threshold <- quantile(abs(w[near_data]), hotspot_quantile)
    is_peak <- local_extremum(w, geometry)
    hotspots <- which(is_peak & abs(w) >= threshold & near_data)

    # 1b. which field makes the local structure of the 2020 map? xi at the
    # omega nodes, and the prominence of omega + xi split into its parts
    year_col <- map_year - fit$t0
    xi_nodes_year <- if (year_col <= fit$n_years) {
      x_year <- x[(year_col - 1) * fit$mesh_xi$n + seq_len(fit$mesh_xi$n)]
      as.numeric(mesh_basis(fit$mesh_xi, geometry$loc) %*% x_year)
    } else rep(0, mesh$n)
    correction_nodes <- w + xi_nodes_year
    prom_total <- prominence(correction_nodes, geometry)
    prom_omega <- prominence(w, geometry)
    prom_xi <- prominence(xi_nodes_year, geometry)
    prom_threshold <- quantile(abs(prom_total[near_data]), hotspot_quantile)
    map_peaks <- which(local_extremum(correction_nodes, geometry) &
                         abs(prom_total) >= prom_threshold & near_data)
    map_hotspot_rows[[length(map_hotspot_rows) + 1]] <- tibble(
      insecticide_type = type, mesh = mesh_name, node = map_peaks,
      value = correction_nodes[map_peaks], omega = w[map_peaks],
      xi = xi_nodes_year[map_peaks], prominence = prom_total[map_peaks],
      prominence_omega = prom_omega[map_peaks],
      prominence_xi = prom_xi[map_peaks],
      is_omega_hotspot = map_peaks %in% hotspots)

    report("%-18s %-4s nodes=%i range=%.0f km threshold |omega|=%.2f: %i omega hotspots, %i map (2020) hotspots",
           type, mesh_name, mesh$n, range_km, threshold, length(hotspots),
           length(map_peaks))

    # 2-3. each omega hotspot --------------------------------------------------
    A <- design$A_omega
    point_radius <- matern_half_radius * range_km
    for (h in hotspots) {
      peak_xy <- geometry$loc[h, ]
      edge <- geometry$local_edge[h]
      footprint <- half_max_footprint(peak_xy, w[h], mesh, w, geometry,
                                      radius_km = max(3 * edge, 4 * range_km))
      # assays that load on the node (its star), and those within one range
      star <- which(A[, h] > 0)
      distance <- sqrt((coords[, 1] - peak_xy[1]) ^ 2 +
                         (coords[, 2] - peak_xy[2]) ^ 2)
      within_range <- which(distance <= range_km)
      s <- sign(w[h])
      prec <- 1 / train$v
      star_mean_r <- if (length(star)) {
        sum(prec[star] * r_omega[star]) / sum(prec[star])
      } else NA_real_
      star_mean_omega <- if (length(star)) {
        sum(prec[star] * omega_i[star]) / sum(prec[star])
      } else NA_real_
      range_mean_r <- if (length(within_range)) {
        sum(prec[within_range] * r_omega[within_range]) /
          sum(prec[within_range])
      } else NA_real_

      # the node's data pull A[, h]' D r, by pixel
      pull <- tibble(cell = train$cell[star], year = train$year[star],
                     p = A[star, h] * prec[star] * r_omega[star]) %>%
        group_by(cell) %>%
        summarise(pull = s * sum(p), .groups = "drop") %>%
        arrange(desc(pull))
      positive_pull <- sum(pmax(pull$pull, 0))
      # the same pull without the mesh: all assays within a fixed 50 km of the
      # peak, unweighted by the basis, so the fine and base meshes (and nodes
      # placed at every site cluster) are judged on the same data. Its
      # largest pixel is the dominant pixel
      near <- which(distance <= neighbourhood_km)
      pull_near <- tibble(cell = train$cell[near],
                          p = prec[near] * r_omega[near]) %>%
        group_by(cell) %>%
        summarise(pull = s * sum(p), .groups = "drop") %>%
        arrange(desc(pull))
      dominant <- if (nrow(pull_near)) pull_near$cell[1] else
        if (nrow(pull)) pull$cell[1] else NA
      dominant_rows <- which(train$cell == dominant)
      consistency <- if (length(dominant_rows)) {
        year_consistency(r_omega[dominant_rows], train$v[dominant_rows],
                         train$year[dominant_rows], s, tau)
      } else {
        tibble(dominant_n_years = 0L, dominant_year_sign_agreement = NA,
               dominant_pooled_mean = NA, dominant_pooled_z = NA,
               dominant_heterogeneity_p = NA)
      }

      positive_near <- sum(pmax(pull_near$pull, 0))
      dominant_xy <- if (length(dominant_rows)) {
        colMeans(coords[dominant_rows, , drop = FALSE])
      } else c(NA, NA)

      hotspot_rows[[length(hotspot_rows) + 1]] <- tibble(
        insecticide_type = type, mesh = mesh_name, node = h,
        x_km = peak_xy[1], y_km = peak_xy[2],
        omega_node = w[h], threshold = threshold,
        range_km = range_km, sigma_omega = fit$hyper$sigma_omega,
        local_edge_km = edge,
        star_area_km2 = 3 * geometry$node_area[h],
        # a lone hat function is >= 1/2 on its star shrunk by half: 1/4 of it
        hat_area_km2 = 3 * geometry$node_area[h] / 4,
        point_source_area_km2 = pi * point_radius ^ 2,
        footprint_area_km2 = footprint[["area"]],
        footprint_radius_km = sqrt(footprint[["area"]] / pi),
        footprint_nodes = footprint[["n_nodes"]],
        footprint_truncated = as.logical(footprint[["touches_edge"]]),
        xi_2020_node = xi_nodes_year[h],
        prominence_2020 = prom_total[h], prominence_omega = prom_omega[h],
        prominence_xi = prom_xi[h],
        n_assays_star = length(star),
        n_pixels = n_distinct(train$cell[star]),
        n_sites = n_distinct(train$site[star]),
        n_years = n_distinct(train$year[star]),
        n_assays_range = length(within_range),
        n_pixels_range = n_distinct(train$cell[within_range]),
        n_sites_range = n_distinct(train$site[within_range]),
        n_years_range = n_distinct(train$year[within_range]),
        star_mean_residual = star_mean_r,
        range_mean_residual = range_mean_r,
        star_mean_omega = star_mean_omega,
        # node / data-level mean residual: > 1 is the node overshooting the
        # data it has to fit, compensating for interpolation dilution
        overshoot_ratio = w[h] / star_mean_r,
        # interpolated omega at the assays / the data-level mean residual: how
        # much of the local signal the field actually reproduces at the data
        fit_ratio = star_mean_omega / star_mean_r,
        largest_pixel = dominant,
        largest_pixel_share = if (positive_pull > 0) {
          max(pull$pull[1], 0) / positive_pull
        } else NA_real_,
        pixel_sign_agreement = mean(pull$pull > 0),
        n_pixels_50km = nrow(pull_near),
        largest_pixel_share_50km = if (positive_near > 0) {
          max(pull_near$pull[1], 0) / positive_near
        } else NA_real_,
        pixel_sign_agreement_50km = mean(pull_near$pull > 0),
        largest_pixel_50km = if (nrow(pull_near)) pull_near$cell[1] else NA,
        dominant_n_assays = length(dominant_rows),
        dominant_distance_km = sqrt(sum((dominant_xy - peak_xy) ^ 2)),
        dominant_omega = if (length(dominant_rows)) {
          mean(omega_i[dominant_rows])
        } else NA_real_
      ) %>%
        bind_cols(consistency)
    }

    # keep what the example figures need
    if (type %in% example_types) {
      example_store[[paste(type, mesh_name)]] <- list(
        mesh = mesh, geometry = geometry, w = w, coords = coords,
        pixel_residual = tibble(cell = train$cell, x = coords[, 1],
                                y = coords[, 2], r = r_omega,
                                prec = 1 / train$v) %>%
          group_by(cell) %>%
          summarise(x = mean(x), y = mean(y),
                    r = sum(prec * r) / sum(prec), prec = sum(prec),
                    .groups = "drop"))
    }

    # 5. what a static per-pixel nugget would absorb ------------------------------
    # pixel-year means of the residual before (r0) and after (r1) omega-hat
    pixel_year <- tibble(cell = train$cell, year = train$year, r0 = r0,
                         r1 = r1, omega = omega_i, prec = 1 / train$v,
                         x = coords[, 1], y = coords[, 2]) %>%
      group_by(cell, year) %>%
      summarise(r0 = sum(prec * r0) / sum(prec),
                r1 = sum(prec * r1) / sum(prec),
                omega = mean(omega), se2 = 1 / sum(prec),
                x = mean(x), y = mean(y), .groups = "drop")
    vc0 <- pixel_variance_components(pixel_year$r0, pixel_year$se2,
                                     pixel_year$cell)
    vc1 <- pixel_variance_components(pixel_year$r1, pixel_year$se2,
                                     pixel_year$cell)

    # omega-hat at each pixel against the pixel's own mean residual and the
    # mean residual of the other pixels within one range (leave own pixel out)
    pixel <- pixel_year %>%
      group_by(cell) %>%
      summarise(r0 = sum(r0 / se2) / sum(1 / se2), prec = sum(1 / se2),
                omega = mean(omega), n_years = n(), x = mean(x),
                y = mean(y), .groups = "drop")
    neighbour_radius <- max(range_km, 50)
    pixel$neighbour_r0 <- vapply(seq_len(nrow(pixel)), function(j) {
      d2 <- (pixel$x - pixel$x[j]) ^ 2 + (pixel$y - pixel$y[j]) ^ 2
      nb <- which(d2 <= neighbour_radius ^ 2 & seq_len(nrow(pixel)) != j)
      if (length(nb) == 0) return(0)
      # shrunk towards 0 with a prior precision of 1 / sigma_omega^2 ~ 4, as
      # a field would be
      sum(pixel$prec[nb] * pixel$r0[nb]) / (sum(pixel$prec[nb]) + 4)
    }, numeric(1))
    pixel$own_r0 <- pixel$prec * pixel$r0 / (pixel$prec + 4)
    r2 <- function(formula) {
      summary(lm(formula, data = pixel))$r.squared
    }

    # covariance of pixel-year residuals r0: same pixel in different years,
    # and distinct pixels by distance (in the same year and in different
    # years). xi-hat is already removed, so what remains is static structure
    # plus noise. Products are unbiased for the covariance whatever the noise;
    # weights 1 / (se2 + 0.5) keep the noisiest pixel-years from dominating
    bins <- c(-1, 0, 10, 25, 50, 100, 200)
    labels <- c("same pixel", "0-10 km", "10-25 km", "25-50 km", "50-100 km",
                "100-200 km")
    mu0 <- weighted.mean(pixel_year$r0, 1 / pixel_year$se2)
    e <- pixel_year$r0 - mu0
    wt <- 1 / (pixel_year$se2 + 0.5)
    # distinct pixels are split by whether the two pixel-years are in the same
    # year: nearby pixels sampled in the same year are often one survey, so
    # a shared same-year covariance would point to a survey-level effect
    # rather than a static site effect
    labels_all <- c(labels[1], as.vector(t(outer(
      labels[-1], c("same year", "other year"), paste))))
    acc <- matrix(0, length(labels_all), 3,
                  dimnames = list(labels_all, c("sum_wp", "sum_w", "n")))
    n_py <- nrow(pixel_year)
    for (j in seq_len(n_py - 1)) {
      o <- (j + 1):n_py
      d <- sqrt((pixel_year$x[o] - pixel_year$x[j]) ^ 2 +
                  (pixel_year$y[o] - pixel_year$y[j]) ^ 2)
      same <- pixel_year$cell[o] == pixel_year$cell[j]
      d[same] <- -0.5
      keep <- d <= 200
      if (!any(keep)) next
      b <- as.character(cut(d[keep], bins, labels = labels))
      distinct_pixel <- b != "same pixel"
      b[distinct_pixel] <- paste(
        b[distinct_pixel],
        ifelse(pixel_year$year[o[keep]][distinct_pixel] == pixel_year$year[j],
               "same year", "other year"))
      b <- factor(b, levels = rownames(acc))
      ww <- wt[j] * wt[o[keep]]
      acc[, "sum_wp"] <- acc[, "sum_wp"] +
        tapply(ww * e[j] * e[o[keep]], b, sum, default = 0)
      acc[, "sum_w"] <- acc[, "sum_w"] + tapply(ww, b, sum, default = 0)
      acc[, "n"] <- acc[, "n"] + tapply(ww, b, length, default = 0)
    }
    covariance_rows[[length(covariance_rows) + 1]] <- tibble(
      insecticide_type = type, mesh = mesh_name, bin = labels_all,
      covariance = acc[, "sum_wp"] / acc[, "sum_w"], n_pairs = acc[, "n"])

    nugget_rows[[length(nugget_rows) + 1]] <- tibble(
      insecticide_type = type, mesh = mesh_name,
      n_pixel_years = nrow(pixel_year), n_pixels = nrow(pixel),
      share_pixels_multi_year = mean(pixel$n_years > 1),
      sigma2_omega = fit$hyper$sigma_omega ^ 2, tau2 = tau ^ 2,
      var_omega_at_pixels = var(pixel$omega),
      sigma2_pixel_before_omega = vc0[["sigma2_pixel"]],
      sigma2_pixel_year_before_omega = vc0[["sigma2_pixel_year"]],
      sigma2_pixel_after_omega = vc1[["sigma2_pixel"]],
      sigma2_pixel_year_after_omega = vc1[["sigma2_pixel_year"]],
      r2_omega_own_pixel = r2(omega ~ own_r0),
      r2_omega_neighbours = r2(omega ~ neighbour_r0),
      r2_omega_both = r2(omega ~ own_r0 + neighbour_r0))
  }
}

hotspots <- bind_rows(hotspot_rows)
hotspots <- bind_cols(hotspots, classify_hotspot(hotspots))
map_hotspots <- bind_rows(map_hotspot_rows)
nugget <- bind_rows(nugget_rows)
covariances <- bind_rows(covariance_rows)


# 2c. the same hotspot on the base mesh -------------------------------------------

# match each fine-mesh hotspot to the nearest base-mesh hotspot of the same sign
# within max(150 km, 1.5 base edges)
if ("base" %in% hotspots$mesh) {
  fine <- filter(hotspots, mesh == "fine")
  base <- filter(hotspots, mesh == "base")
  matched <- lapply(seq_len(nrow(fine)), function(j) {
    b <- filter(base, insecticide_type == fine$insecticide_type[j],
                sign(omega_node) == sign(fine$omega_node[j]))
    if (nrow(b) == 0) return(NULL)
    d <- sqrt((b$x_km - fine$x_km[j]) ^ 2 + (b$y_km - fine$y_km[j]) ^ 2)
    i <- which.min(d)
    if (d[i] > max(150, 1.5 * b$local_edge_km[i])) return(NULL)
    tibble(insecticide_type = fine$insecticide_type[j], node = fine$node[j],
           mesh = "fine", base_node = b$node[i], base_shift_km = d[i],
           base_omega_node = b$omega_node[i],
           base_local_edge_km = b$local_edge_km[i],
           base_footprint_area_km2 = b$footprint_area_km2[i],
           base_dominant_omega = b$dominant_omega[i],
           base_largest_pixel = b$largest_pixel[i])
  })
  hotspots <- left_join(hotspots, bind_rows(matched),
                        by = c("insecticide_type", "node", "mesh"))
}

write.csv(hotspots, file.path(output_dir, "hotspot_diagnostics.csv"),
          row.names = FALSE)
write.csv(bind_rows(nugget), file.path(output_dir, "hotspot_nugget.csv"),
          row.names = FALSE)
write.csv(covariances, file.path(output_dir, "hotspot_covariance.csv"),
          row.names = FALSE)


# summaries ----------------------------------------------------------------------

summary_table <- hotspots %>%
  group_by(insecticide_type, mesh) %>%
  summarise(
    n_hotspots = n(),
    range_km = first(range_km),
    threshold = first(threshold),
    median_abs_omega = median(abs(omega_node)),
    n_C_persistent = sum(cause == "C persistent"),
    n_C_one_year = sum(cause == "C one year"),
    n_B = sum(cause == "B"),
    n_mixed = sum(cause == "mixed"),
    n_mesh_artefact = sum(mesh_artefact),
    median_edge_km = median(local_edge_km),
    median_footprint_radius_km = median(footprint_radius_km, na.rm = TRUE),
    median_point_source_radius_km = median(sqrt(point_source_area_km2 / pi)),
    median_footprint_over_hat = median(footprint_area_km2 / hat_area_km2,
                                       na.rm = TRUE),
    median_footprint_nodes = median(footprint_nodes, na.rm = TRUE),
    median_overshoot = median(overshoot_ratio, na.rm = TRUE),
    median_fit_ratio = median(fit_ratio, na.rm = TRUE),
    share_overshoot_above_1 = mean(overshoot_ratio > 1, na.rm = TRUE),
    median_node_over_dominant_omega = median(omega_node / dominant_omega,
                                             na.rm = TRUE),
    median_dominant_distance_km = median(dominant_distance_km, na.rm = TRUE),
    median_largest_pixel_share = median(largest_pixel_share, na.rm = TRUE),
    median_largest_pixel_share_50km = median(largest_pixel_share_50km,
                                             na.rm = TRUE),
    median_pixel_sign_agreement_50km = median(pixel_sign_agreement_50km,
                                              na.rm = TRUE),
    median_n_pixels_50km = median(n_pixels_50km),
    share_dominant_multi_year = mean(dominant_n_years >= 2),
    median_pixel_sign_agreement = median(pixel_sign_agreement, na.rm = TRUE),
    median_n_pixels = median(n_pixels),
    median_share_prominence_omega = median(prominence_omega / prominence_2020),
    matched_in_base = if ("base_shift_km" %in% names(hotspots)) {
      mean(!is.na(base_shift_km))
    } else NA_real_,
    median_base_shift_km = if ("base_shift_km" %in% names(hotspots)) {
      median(base_shift_km, na.rm = TRUE)
    } else NA_real_,
    median_peak_fine_over_base = if ("base_omega_node" %in% names(hotspots)) {
      median(omega_node / base_omega_node, na.rm = TRUE)
    } else NA_real_,
    median_area_fine_over_base = if ("base_footprint_area_km2" %in%
                                     names(hotspots)) {
      median(footprint_area_km2 / base_footprint_area_km2, na.rm = TRUE)
    } else NA_real_,
    median_edge2_fine_over_base = if ("base_local_edge_km" %in%
                                      names(hotspots)) {
      median((local_edge_km / base_local_edge_km) ^ 2, na.rm = TRUE)
    } else NA_real_,
    median_dominant_omega_fine_over_base = if ("base_dominant_omega" %in%
                                               names(hotspots)) {
      median(dominant_omega / base_dominant_omega, na.rm = TRUE)
    } else NA_real_,
    .groups = "drop")

map_summary <- map_hotspots %>%
  group_by(insecticide_type, mesh) %>%
  summarise(n_map_hotspots = n(),
            median_share_omega = median(prominence_omega / prominence),
            share_omega_dominated = mean(prominence_omega / prominence > 0.5),
            share_coincide_omega_hotspot = mean(is_omega_hotspot),
            .groups = "drop")
summary_table <- left_join(summary_table, map_summary,
                           by = c("insecticide_type", "mesh"))
write.csv(summary_table, file.path(output_dir, "hotspot_summary.csv"),
          row.names = FALSE)

options(width = 200)
print(as.data.frame(summary_table), digits = 3)
print(as.data.frame(nugget), digits = 3)
print(as.data.frame(tidyr::pivot_wider(covariances %>% select(-n_pairs),
                                       names_from = insecticide_type,
                                       values_from = covariance)),
      digits = 3)


# figures ------------------------------------------------------------------------

# zoomed maps of example hotspots: omega-hat on a fine grid, the mesh, and
# the pixels coloured by their precision-weighted mean partial residual
# z - m - xi - u, for the fine mesh and the base mesh side by side
zoom_panel <- function(store, centre, half_width, limit, title) {
  xs <- seq(centre[1] - half_width, centre[1] + half_width, length.out = 150)
  ys <- seq(centre[2] - half_width, centre[2] + half_width, length.out = 150)
  grid <- as.matrix(expand.grid(x = xs, y = ys))
  field <- tibble(x = grid[, 1], y = grid[, 2],
                  omega = as.numeric(mesh_basis(store$mesh, grid) %*% store$w))
  in_box <- function(x, y) {
    abs(x - centre[1]) <= half_width & abs(y - centre[2]) <= half_width
  }
  loc <- store$geometry$loc
  edges <- store$geometry$edges
  segments <- tibble(x = loc[edges[, 1], 1], y = loc[edges[, 1], 2],
                     xend = loc[edges[, 2], 1], yend = loc[edges[, 2], 2]) %>%
    filter(in_box(x, y) | in_box(xend, yend))
  nodes <- tibble(x = loc[, 1], y = loc[, 2], omega = store$w) %>%
    filter(in_box(x, y))
  pixels <- filter(store$pixel_residual, in_box(x, y))
  clamp <- function(v) pmax(pmin(v, limit), -limit)
  ggplot() +
    geom_raster(aes(x, y, fill = clamp(omega)), data = field) +
    geom_segment(aes(x, y, xend = xend, yend = yend), data = segments,
                 colour = "grey40", linewidth = 0.2) +
    geom_point(aes(x, y), data = nodes, size = 0.6, colour = "grey20") +
    geom_point(aes(x, y, fill = clamp(r), size = prec), data = pixels,
               shape = 21, colour = "black", stroke = 0.3) +
    scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac",
                         limits = c(-limit, limit),
                         name = "omega-hat /\npixel residual") +
    scale_size_area(max_size = 4, name = "pixel\nprecision") +
    coord_equal(xlim = centre[1] + c(-1, 1) * half_width,
                ylim = centre[2] + c(-1, 1) * half_width, expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(axis.text = element_blank(), panel.grid = element_blank(),
          plot.title = element_text(size = 7.5))
}

for (type in intersect(example_types, unique(hotspots$insecticide_type))) {
  fine_store <- example_store[[paste(type, "fine")]]
  base_store <- example_store[[paste(type, "base")]]
  h_type <- filter(hotspots, insecticide_type == type, mesh == "fine") %>%
    arrange(desc(abs(omega_node)))
  # the largest hotspot of each cause, then the largest remaining, up to 3
  pick <- bind_rows(h_type %>% group_by(cause) %>% slice(1) %>% ungroup(),
                    h_type) %>%
    distinct(node, .keep_all = TRUE) %>%
    arrange(desc(abs(omega_node))) %>%
    slice(1:min(3, n()))
  panels <- list()
  for (j in seq_len(nrow(pick))) {
    h <- pick[j, ]
    half_width <- max(2.5 * h$local_edge_km, 3 * h$range_km, 60)
    limit <- max(abs(h$omega_node), abs(h$star_mean_residual), 0.3,
                 na.rm = TRUE)
    title <- sprintf(paste("fine mesh: %s, omega-hat %.2f, edge %.0f km,",
                           "range %.0f km,\nlargest-pixel share (50 km) %.2f,",
                           "overshoot %.1f"),
                     h$cause, h$omega_node, h$local_edge_km, h$range_km,
                     h$largest_pixel_share_50km, h$overshoot_ratio)
    panels[[length(panels) + 1]] <- zoom_panel(
      fine_store, c(h$x_km, h$y_km), half_width, limit, title)
    if (!is.null(base_store)) {
      base_title <- if (!is.null(h$base_omega_node) &&
                        !is.na(h$base_omega_node)) {
        sprintf("base mesh: omega-hat %.2f, edge %.0f km, shift %.0f km",
                h$base_omega_node, h$base_local_edge_km, h$base_shift_km)
      } else "base mesh: no matching hotspot"
      panels[[length(panels) + 1]] <- zoom_panel(
        base_store, c(h$x_km, h$y_km), half_width, limit, base_title)
    }
  }
  figure <- wrap_plots(panels, ncol = if (is.null(base_store)) 1 else 2) +
    plot_annotation(
      title = sprintf("%s: largest omega hotspots", type),
      subtitle = paste0("omega-hat (background), mesh edges and nodes, and ",
                        "pixels coloured by their mean partial residual\n",
                        "z - m - xi - u (size: precision); left the adopted ",
                        "mesh, right the base mesh, same window"))
  ggsave(file.path(figure_dir, sprintf("hotspot_%s.png", type)), figure,
         width = 8, height = 3.3 * nrow(pick), dpi = 150, bg = "white")
}

# footprint radius against local edge length, with the resolved point-source
# radius for reference
footprint_plot <- hotspots %>%
  filter(!is.na(footprint_radius_km)) %>%
  ggplot(aes(local_edge_km, footprint_radius_km, colour = mesh)) +
  geom_abline(slope = sqrt(6 * sqrt(3) / (16 * pi)), linetype = 2,
              colour = "grey50") +
  geom_hline(aes(yintercept = matern_half_radius * range_km, colour = mesh),
             data = distinct(hotspots, insecticide_type, mesh, range_km),
             linetype = 3) +
  geom_point(size = 1, alpha = 0.7) +
  facet_wrap(~insecticide_type, scales = "free") +
  labs(x = "mean edge length at the peak node (km)",
       y = "half-maximum footprint radius (km)",
       title = "omega hotspots: footprint vs mesh resolution",
       subtitle = paste("dashed: a lone hat function on a regular mesh;",
                        "dotted: a resolved point source at the fitted range")) +
  theme_minimal(base_size = 9)
ggsave(file.path(figure_dir, "hotspot_footprint_vs_edge.png"), footprint_plot,
       width = 9, height = 6, dpi = 150, bg = "white")

report("done")
