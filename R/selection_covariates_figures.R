# Distributions of the selection covariates at the bioassays and across Africa
# (#23): population (min-max scaled, raw and log), net use (raw and hinged at
# the knots tested in selection_shape_analysis.R) and IRS, by year, with the
# number of bioassays per year.
#
# The net-use cube starts in 2000 and IRS in 1997; earlier years repeat the
# first layer, as in the model.
#
# Outputs:
#   figures/selection_covariates_*.png
#   outputs/selection_covariates_summary.csv

library(tidyverse)
library(terra)
library(patchwork)
source("R/functions.R")

# build df and the covariate cubes exactly as fit_model.R does
fit_model_exprs <- parse("R/fit_model.R")
fit_model_text <- vapply(fit_model_exprs,
                         function(e) paste(deparse(e), collapse = " "),
                         "")
last_expr <- which(startsWith(fit_model_text, "x_cell_years <-"))
stopifnot(length(last_expr) == 1)
for (i in seq_len(last_expr)) {
  if (grepl("^source\\(\"R/(packages|functions).R\"\\)", fit_model_text[i])) {
    next
  }
  eval(fit_model_exprs[[i]], envir = globalenv())
}

# knots used in selection_shape_analysis.R (quartiles of net use over the
# pair-years there)
nets_knots <- c(0.180, 0.433, 0.605)

# log population, min-max scaled over 2000-2030 as in selection_shape_analysis.R
pop_raw_cube <- c(rast("data/clean/pop_cube.tif"),
                  rast("data/clean/pop_cube_future.tif"))
log_pop_range <- range(log(global(pop_raw_cube, "range", na.rm = TRUE)))
pop_raw_cube <- pre_pad_cube(pop_raw_cube[[paste0("pop_", 2000:final_data_year)]],
                             baseline_year)

model_years <- baseline_year:final_data_year

# covariate values at a set of cells for every year, long format
extract_cells <- function(cells) {
  layers <- list(nets = nets_cube, irs = irs_cube, pop = pop_cube,
                 pop_raw = pop_raw_cube)
  imap(layers, function(cube, name) {
    cube <- cube[[paste0(name %>% str_remove("_raw"), "_", model_years)]]
    terra::extract(cube, cells) %>%
      as_tibble() %>%
      mutate(cell = cells) %>%
      pivot_longer(-cell, names_to = "year", values_to = name) %>%
      mutate(year = as.numeric(str_sub(year, -4)))
  }) %>%
    reduce(left_join, by = c("cell", "year")) %>%
    mutate(log_pop = (log(pop_raw) - log_pop_range[1]) / diff(log_pop_range))
}

# bioassays: one row per bioassay, with the covariates of its pixel-year. The
# classes have near-identical covariate distributions (see the summary table),
# so the figures pool them.
bioassay_covs <- df %>%
  select(cell, year = year_start, insecticide_class) %>%
  left_join(extract_cells(unique(df$cell)), by = c("cell", "year")) %>%
  mutate(class_group = if_else(insecticide_class == "Pyrethroids",
                               "pyrethroid", "other classes"),
         group = "bioassay pixel-years")
stopifnot(!anyNA(bioassay_covs$log_pop))

# all pixels: a random sample of mask pixels, every year
set.seed(2026)
mask_cells <- which(!is.na(values(mask)[, 1]))
sample_cells <- sample(mask_cells, 20000)
africa_covs <- extract_cells(sample_cells) %>%
  filter(!is.na(nets), !is.na(irs), !is.na(pop), !is.na(log_pop)) %>%
  mutate(group = "all pixels (20,000 sampled)")

group_colours <- c("bioassay pixel-years" = "#c0392b",
                   "all pixels (20,000 sampled)" = "#34495e")

# add the hinged net-use versions
add_hinges <- function(data) {
  for (k in nets_knots) {
    data[[sprintf("min(nets, %.2f)", k)]] <- pmin(data$nets, k)
  }
  data
}
bioassay_covs <- add_hinges(bioassay_covs)
africa_covs <- add_hinges(africa_covs)
all_covs <- bind_rows(bioassay_covs, africa_covs) %>%
  mutate(group = factor(group, names(group_colours)))

versions <- c("pop", "log_pop", "nets", sprintf("min(nets, %.2f)", nets_knots),
              "irs")
version_labels <- c(pop = "population, min-max scaled (model)",
                    log_pop = "log population, min-max scaled",
                    nets = "net use",
                    irs = "IRS coverage, scaled")
version_labels <- c(version_labels,
                    setNames(sprintf("min(net use, %.2f)", nets_knots),
                             sprintf("min(nets, %.2f)", nets_knots)))
pop_versions <- version_labels[c("pop", "log_pop")]
net_versions <- version_labels[c("nets", sprintf("min(nets, %.2f)", nets_knots))]

long_covs <- all_covs %>%
  select(group, class_group, year, all_of(versions)) %>%
  pivot_longer(all_of(versions), names_to = "version") %>%
  mutate(version = factor(version, versions, version_labels[versions]))

# median and 50% / 90% bands by year
year_bands <- long_covs %>%
  group_by(group, version, year) %>%
  summarise(q05 = quantile(value, 0.05), q25 = quantile(value, 0.25),
            q50 = median(value), q75 = quantile(value, 0.75),
            q95 = quantile(value, 0.95), n = n(),
            .groups = "drop")

# where the bioassays sit in the all-pixel range: the share of all pixels below
# the bioassays' 5% quantile and above their 95% quantile, by year
outside_share <- year_bands %>%
  filter(group == "bioassay pixel-years") %>%
  select(version, year, bioassay_q05 = q05, bioassay_q95 = q95) %>%
  right_join(long_covs %>% filter(group == "all pixels (20,000 sampled)"),
             by = c("version", "year")) %>%
  group_by(version, year) %>%
  summarise(`below bioassay 5%` = mean(value < bioassay_q05),
            `above bioassay 95%` = mean(value > bioassay_q95),
            .groups = "drop") %>%
  pivot_longer(c(`below bioassay 5%`, `above bioassay 95%`),
               names_to = "side", values_to = "share")

# population on a pseudo-log scale, so that zero, near-zero and populated
# pixels are all distinguishable
pop_trans <- scales::pseudo_log_trans(sigma = 1e-7, base = 10)
pop_breaks <- c(0, 1e-6, 1e-4, 1e-2, 1)

theme_panels <- function() {
  theme_minimal(base_size = 9) +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 8))
}

# bands by year, bioassays and all pixels in side-by-side columns with a shared
# y axis per row
plot_bands <- function(data) {
  ggplot(data, aes(year, q50, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.2, colour = NA) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.35, colour = NA) +
    geom_line(linewidth = 0.8) +
    facet_grid(version ~ group, scales = "free_y",
               labeller = label_wrap_gen(22)) +
    scale_colour_manual(values = group_colours, guide = "none") +
    scale_fill_manual(values = group_colours, guide = "none") +
    labs(x = NULL, y = "median, 50% and 90% bands") +
    theme_panels()
}

plot_outside <- function(data) {
  ggplot(data, aes(year, share, linetype = side)) +
    geom_line(linewidth = 0.7, colour = group_colours[2]) +
    facet_wrap(~ version, nrow = 1, labeller = label_wrap_gen(22)) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "share of all pixels", linetype = NULL,
         subtitle = "All pixels outside the bioassays' 5-95% range that year") +
    theme_panels()
}

plot_density <- function(data) {
  ggplot(data, aes(value, colour = group)) +
    geom_density(linewidth = 0.8, adjust = 0.7) +
    scale_colour_manual(values = group_colours, name = NULL) +
    labs(x = NULL, y = "density, all years pooled") +
    theme_panels()
}

dir.create("figures", showWarnings = FALSE)
save_figure <- function(plot, name, height = 5) {
  ggsave(file.path("figures", paste0("selection_covariates_", name, ".png")),
         plot, width = 8, height = height, dpi = 200, bg = "white")
}


# 1. population, raw and log min-max scaled

pop_density_raw <- plot_density(long_covs %>%
                                  filter(version == pop_versions[1])) +
  scale_x_continuous(trans = pop_trans, breaks = pop_breaks,
                     labels = c("0", "1e-6", "1e-4", "0.01", "1")) +
  labs(subtitle = paste(pop_versions[1], "(pseudo-log axis)"))
pop_density_log <- plot_density(long_covs %>%
                                  filter(version == pop_versions[2])) +
  labs(subtitle = pop_versions[2], y = NULL)
population_distribution <- (pop_density_raw | pop_density_log) /
  plot_outside(outside_share %>% filter(version %in% pop_versions)) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Population covariate: bioassays and all pixels") &
  theme(legend.position = "bottom")
save_figure(population_distribution, "population_distribution", height = 6)

population_by_year <-
  (plot_bands(year_bands %>% filter(version == pop_versions[1])) +
     scale_y_continuous(trans = pop_trans, breaks = pop_breaks,
                        labels = c("0", "1e-6", "1e-4", "0.01", "1")) +
     labs(y = "pseudo-log axis")) /
  (plot_bands(year_bands %>% filter(version == pop_versions[2])) +
     labs(y = NULL)) +
  plot_annotation(title = "Population covariate by year: median, 50% and 90% bands")
save_figure(population_by_year, "population_by_year", height = 6)


# 2. net use, raw and hinged

nets_by_year <- plot_bands(year_bands %>% filter(version %in% net_versions)) +
  labs(title = "Net use, raw and hinged, by year: median, 50% and 90% bands")
save_figure(nets_by_year, "nets_by_year", height = 8)

nets_low <- all_covs %>%
  group_by(group, year) %>%
  summarise(`net use < 0.05` = mean(nets < 0.05),
            `net use < 0.18` = mean(nets < 0.18),
            .groups = "drop") %>%
  pivot_longer(starts_with("net use"), names_to = "threshold",
               values_to = "share") %>%
  ggplot(aes(year, share, colour = group, linetype = threshold)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = group_colours, name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "share", linetype = NULL,
       subtitle = "Share of bioassays / pixels with low net use") +
  theme_panels() +
  theme(legend.box = "vertical")
nets_outside <- plot_outside(outside_share %>%
                               filter(version %in% net_versions))
save_figure(nets_low / nets_outside +
              plot_annotation(title = "Net use: low values, and coverage of the all-pixel range"),
            "nets_coverage", height = 6.5)

period_labels <- c("1995-2004", "2005-2012", "2013-2024")
nets_histogram <- all_covs %>%
  mutate(period = cut(year, c(-Inf, 2004, 2012, Inf), labels = period_labels)) %>%
  ggplot(aes(nets, after_stat(density), fill = group)) +
  geom_histogram(binwidth = 0.025, boundary = 0) +
  geom_vline(xintercept = nets_knots, linetype = 2, colour = "grey30") +
  facet_grid(group ~ period, scales = "free_y",
             labeller = label_wrap_gen(20)) +
  scale_fill_manual(values = group_colours, guide = "none") +
  labs(x = "net use (dashed: hinge knots)", y = "density",
       title = "Net use by period: bioassay pixel-years and all pixels") +
  theme_panels()
save_figure(nets_histogram, "nets_histogram")


# 3. IRS, compactly

irs_label <- version_labels["irs"]
irs_any <- all_covs %>%
  group_by(group, year) %>%
  summarise(share = mean(irs > 0), .groups = "drop") %>%
  ggplot(aes(year, share, colour = group)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = group_colours, name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "share with IRS > 0", subtitle = "Any IRS") +
  theme_panels()
irs_plot <- plot_bands(year_bands %>% filter(version == irs_label)) /
  (irs_any | plot_outside(outside_share %>% filter(version == irs_label))) +
  plot_annotation(title = "IRS coverage (scaled): bioassays and all pixels")
save_figure(irs_plot, "irs", height = 6)


# 4. bioassays per year

count_plot <- bioassay_covs %>%
  count(class_group, year) %>%
  ggplot(aes(year, n, fill = class_group)) +
  geom_col() +
  scale_fill_manual(values = c(pyrethroid = "#c0392b",
                               `other classes` = "#e6b0aa"),
                    name = NULL) +
  labs(x = NULL, y = "bioassays", title = "Bioassays per year") +
  theme_panels()
save_figure(count_plot, "bioassays_per_year", height = 3)


# numeric summary: quantiles of each version at the bioassays (pooled and by
# class) and across all pixels (all years pooled), the share of all pixels
# outside the bioassays' 5-95% range, and the correlation with year at the
# bioassays

summary_groups <- bind_rows(
  long_covs %>% mutate(group = as.character(group)),
  long_covs %>%
    filter(group == "bioassay pixel-years") %>%
    mutate(group = paste("bioassays:", class_group))
)
bioassay_range <- long_covs %>%
  filter(group == "bioassay pixel-years") %>%
  group_by(version) %>%
  summarise(bq05 = quantile(value, 0.05), bq95 = quantile(value, 0.95))
covariate_summary <- summary_groups %>%
  left_join(bioassay_range, by = "version") %>%
  group_by(version, group) %>%
  summarise(n = n(),
            q05 = quantile(value, 0.05),
            q25 = quantile(value, 0.25),
            q50 = median(value),
            q75 = quantile(value, 0.75),
            q95 = quantile(value, 0.95),
            sd = sd(value),
            share_below_bioassay_q05 = mean(value < bq05),
            share_above_bioassay_q95 = mean(value > bq95),
            cor_with_year = cor(value, year),
            .groups = "drop") %>%
  arrange(version, group)
print(covariate_summary, n = Inf, width = Inf)
write_csv(covariate_summary, "outputs/selection_covariates_summary.csv")
