# Diagnostics of the initial-state population covariate (pop_2000) against
# where and when the bioassays were done: whether few early bioassays come
# from populated pixels, so that the initial state there is weakly pinned.
#
#   Rscript R/diag_init_pop_data.R
#
# Needs no fitted model: the bioassays are the fit's subset
# (subset_modelled_bioassays()), and the covariates are built by the shared
# functions of R/model_covariates.R. Writes
#   figures/diag_pop_transform_maps.png   encounter transform at d_half 5, 50
#                                         and 200, 2000 and 2020, and log10
#                                         density
#   figures/diag_pop_init_covariate.png   standardised pop_2000 as used, at
#                                         all cells and at the bioassays
#   figures/diag_bioassay_timing_pop.png  bioassays per year by density class
#   figures/diag_bioassay_pop_period.png  pop_2000 at the bioassays by period
#   outputs/diag_bioassay_pop_counts.csv  bioassays by period, density class
#                                         and type

source("R/packages.R")
source("R/functions.R")
source("R/bioassay_subset.R")
source("R/model_covariates.R")

baseline_year <- 1995
final_data_year <- 2024

mask <- rast("data/clean/raster_mask.tif")
borders <- readRDS("data/clean/country_borders.RDS")

# the fit's design: encounter transform, d_half = 50
design <- selection_design()

# bioassays, as in fit_model.R
ir_africa <- readRDS("data/clean/all_gambiae_complex_data.RDS")
df <- subset_modelled_bioassays(ir_africa, mask,
                                baseline_year = baseline_year,
                                final_data_year = final_data_year)

# population density (people per km2, empty pixels 0) in 2000 and 2020
area <- terra::cellSize(mask, unit = "km")
pop <- pop_zero_filled(rast("data/clean/pop_cube.tif")[[c("pop_2000",
                                                          "pop_2020")]])
density <- terra::mask(pop / area, mask)
names(density) <- c("2000", "2020")

# standardised initial-state covariates as the fit uses them
init_layers <- init_covariate_layers(design)

cells <- unique(df$cell)
df <- df %>%
  mutate(
    density_2000 = terra::extract(density[["2000"]], cell)[, 1],
    pop_2000_std = terra::extract(init_layers[["pop_2000"]], cell)[, 1],
    density_class = cut(density_2000,
                        c(-Inf, 10, 100, 1000, Inf),
                        labels = c("<10", "10-100", "100-1,000", ">1,000"),
                        right = FALSE),
    period = cut(year_start, c(-Inf, 2005, 2013, Inf),
                 labels = c("1995-2004", "2005-2012", "2013-2024"),
                 right = FALSE),
    insecticide_class = factor(insecticide_class),
    insecticide_type = factor(insecticide_type,
                              levels = insecticides_plot_order)
  )
stopifnot(!anyNA(df$pop_2000_std))

density_colours <- setNames(
  c("#c6dbef", "#6baed6", "#2171b5", "#08306b"),
  levels(df$density_class))

# coarser rasters for plotting (mean of 3 x 3 pixels)
plot_agg <- function(r) terra::aggregate(r, 3, mean, na.rm = TRUE)

map_theme <- theme_minimal(base_size = 9) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 9))

country_lines <- geom_sf(data = borders, fill = NA, colour = grey(0.5),
                         linewidth = 0.1)


# figure 1a: the encounter transform at three half-saturation densities ------

d_halves <- c(5, 50, 200)
encounter <- lapply(d_halves, function(d_half) {
  out <- 1 - exp(-density * log(2) / d_half)
  names(out) <- sprintf("d_half = %i, %s", d_half, names(density))
  out
})
encounter <- plot_agg(do.call(c, encounter))

p_encounter <- ggplot() +
  geom_spatraster(data = encounter) +
  country_lines +
  facet_wrap(~lyr, ncol = 2) +
  scale_fill_viridis_c(name = "Encounter\ntransform\n(0-1)",
                       limits = c(0, 1), na.value = "transparent") +
  map_theme

log_density <- plot_agg(log10(density + 0.1))
names(log_density) <- paste("log10 density,", names(density))
p_log <- ggplot() +
  geom_spatraster(data = log_density) +
  country_lines +
  facet_wrap(~lyr, ncol = 2) +
  scale_fill_viridis_c(name = "log10 people\nper km²\n(+0.1)",
                       option = "magma", na.value = "transparent") +
  map_theme

p1a <- p_encounter / p_log +
  plot_layout(heights = c(3, 1)) +
  plot_annotation(
    caption = paste(
      "Encounter transform 1 - exp(-d log 2 / d_half) of population density d",
      "(people per km², WorldPop, empty pixels 0),\non a common 0-1 scale.",
      "The fit uses d_half = 50 people per km². Bottom: log10(d + 0.1).",
      "Pixels averaged over 3 x 3 for plotting."))
ggsave("figures/diag_pop_transform_maps.png", p1a, width = 7.5, height = 12,
       dpi = 150, bg = "white")
rm(encounter, log_density, p_encounter, p_log, p1a)
invisible(gc())


# figure 1b: the standardised initial-state covariate, and its values at the
# bioassays ------------------------------------------------------------------

pop_std <- plot_agg(init_layers[["pop_2000"]])
names(pop_std) <- "Standardised pop_2000 (d_half = 50), all mask cells"
sites <- df %>%
  distinct(longitude, latitude, cell, pop_2000_std)

p_std_map <- ggplot() +
  geom_spatraster(data = pop_std) +
  country_lines +
  geom_point(aes(longitude, latitude), data = sites, shape = 21,
             fill = NA, colour = "red", size = 0.5, stroke = 0.2) +
  facet_wrap(~lyr) +
  scale_fill_viridis_c(name = "Standardised\npop_2000\n(sd units)",
                       na.value = "transparent") +
  map_theme

mask_values <- terra::values(init_layers[["pop_2000"]], mat = FALSE)
mask_values <- mask_values[!is.na(mask_values)]
std_compare <- bind_rows(
  tibble(set = "All mask cells", value = mask_values),
  tibble(set = "Bioassay sites (unique cells)",
         value = unique(select(df, cell, pop_2000_std))$pop_2000_std),
  tibble(set = "Bioassays (records)", value = df$pop_2000_std)
)
rm(mask_values)
std_summary <- std_compare %>%
  group_by(set) %>%
  summarise(n = n(), mean = mean(value), median = median(value),
            share_above_1sd = mean(value > 1), .groups = "drop")
print(std_summary)

p_std_hist <- ggplot(std_compare, aes(value, after_stat(density),
                                      fill = set)) +
  geom_histogram(binwidth = 0.1, position = "identity", alpha = 0.45,
                 boundary = 0) +
  scale_fill_manual(values = c("grey40", "red", "orange"), name = NULL) +
  labs(x = "Standardised pop_2000 (sd units over all mask cells)",
       y = "Density") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom")

p1b <- p_std_map / p_std_hist +
  plot_layout(heights = c(2.2, 1)) +
  plot_annotation(caption = paste(
    "Top: the initial-state population covariate as used in the fit,",
    "standardised over all mask cells; red circles are bioassay sites.\n",
    "Bottom: its distribution over all mask cells, unique bioassay cells and",
    "bioassay records."))
ggsave("figures/diag_pop_init_covariate.png", p1b, width = 7.5, height = 10,
       dpi = 150, bg = "white")
rm(pop_std, p_std_map, p1b)
invisible(gc())


# figure 2a: bioassays per year by density class -----------------------------

counts_year <- df %>%
  count(insecticide_class, year_start, density_class)

p2a_counts <- ggplot(counts_year, aes(year_start, n, fill = density_class)) +
  geom_col(width = 0.9) +
  facet_wrap(~insecticide_class, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = density_colours,
                    name = "2000 density\n(people per km²)") +
  labs(x = "Year", y = "Number of bioassays") +
  theme_minimal(base_size = 9)

p2a_share <- ggplot(counts_year, aes(year_start, n, fill = density_class)) +
  geom_col(width = 0.9, position = "fill") +
  facet_wrap(~insecticide_class, ncol = 2) +
  scale_fill_manual(values = density_colours,
                    name = "2000 density\n(people per km²)") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Year", y = "Share of bioassays in the year") +
  theme_minimal(base_size = 9)

p2a <- p2a_counts / p2a_share +
  plot_layout(guides = "collect") +
  plot_annotation(caption = paste(
    "Bioassays (records in the fitted data) per year, by insecticide class,",
    "stratified by the 2000 population density of the bioassay pixel.",
    "\nTop: counts; bottom: share of each year's bioassays."))
ggsave("figures/diag_bioassay_timing_pop.png", p2a, width = 8, height = 9,
       dpi = 150, bg = "white")


# figure 2b: pop_2000 at the bioassays by period ------------------------------

p2b <- ggplot(df, aes(period, pop_2000_std, fill = insecticide_type)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3, varwidth = FALSE) +
  facet_wrap(~insecticide_type, ncol = 3) +
  scale_fill_manual(values = insecticide_colours(), guide = "none") +
  geom_text(aes(x = period, label = n, y = -1.6),
            data = count(df, insecticide_type, period),
            size = 2.5, inherit.aes = FALSE) +
  labs(x = "Period of the bioassay",
       y = "Standardised pop_2000 at the bioassay (sd units)") +
  theme_minimal(base_size = 9) +
  labs(caption = paste(
    "Boxes: median and interquartile range of the standardised pop_2000",
    "covariate at the bioassay pixels, by period;\nwhiskers 1.5 IQR;",
    "numbers are bioassays. Over all mask cells the covariate has mean 0 and",
    "sd 1."))
ggsave("figures/diag_bioassay_pop_period.png", p2b, width = 8, height = 7,
       dpi = 150, bg = "white")


# counts table ---------------------------------------------------------------

counts <- df %>%
  count(insecticide_type, period, density_class) %>%
  complete(insecticide_type, period, density_class, fill = list(n = 0)) %>%
  arrange(insecticide_type, period, density_class)
write_csv(counts, "outputs/diag_bioassay_pop_counts.csv")

df %>%
  count(period, density_class) %>%
  pivot_wider(names_from = density_class, values_from = n) %>%
  print()
df %>%
  group_by(period) %>%
  summarise(n = n(), mean_std = mean(pop_2000_std),
            median_std = median(pop_2000_std),
            share_100plus = mean(density_2000 >= 100)) %>%
  print()
df %>%
  filter(insecticide_type %in% c("Deltamethrin", "Permethrin", "DDT",
                                 "Bendiocarb")) %>%
  count(insecticide_type, period, density_class) %>%
  pivot_wider(names_from = density_class, values_from = n,
              values_fill = 0) %>%
  print(n = 30)
