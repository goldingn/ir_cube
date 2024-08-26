# make plots of the available data by time, space, and insecticide type

source("R/packages.R")
source("R/functions.R")

# load the mask
mask <- rast("data/clean/raster_mask.tif")

ir_africa <- readRDS("data/clean/all_gambiae_complex_data.RDS")

types <- unique(ir_africa$insecticide_type)
classes <- unique(ir_africa$insecticide_class)

ir_africa$class_id <- match(ir_africa$insecticide_class, classes)
ir_africa$type_id <- match(ir_africa$insecticide_type, types)

# index to the classes for each type
classes_index <- ir_africa %>%
  distinct(type_id, class_id) %>%
  arrange(type_id) %>%
  pull(class_id)

insecticides_plot <- tibble(
  insecticide = types,
  class = classes[classes_index]
) %>%
  arrange(desc(class), insecticide) %>%
  pull(insecticide)

colour_types <- scales::hue_pal(direction = -1)(9)
types_plot_id <- match(insecticides_plot, types)
colours_plot <- colour_types[types_plot_id]

df_plot <- ir_africa %>%
  rename(
    year = year_start,
    country = country_name,
    insecticide = insecticide_type
  ) %>%
  group_by(
    country, year, insecticide
  ) %>%
  filter(
    year >= 2000
  ) %>%
  summarise(
    mosquito_number = sum(mosquito_number),
    .groups = "drop"
  ) %>%
  mutate(
    type_id = match(insecticide, types),
    type_id = factor(type_id),
    country = factor(country)
  )

for (this_insecticide in types) {
  df_plot %>%
    filter(
      insecticide == this_insecticide
    ) %>%
    ggplot(
      aes(
        x = year,
        y = country,
        size = mosquito_number,
      )
    ) +
    geom_point(
      fill = colours_plot[match(this_insecticide, types)],
      shape = 21,
      alpha = 1
    ) +
    ylab("") +
    xlab("") +
    scale_size_continuous(
      labels = scales::number_format(accuracy = 1000, big.mark = ","),
      limits = range(df_plot$mosquito_number)
    ) +
    scale_x_continuous(
      limits = range(df_plot$year)
    ) +
    scale_y_discrete(
      labels = rev(sort(unique(df_plot$country))),
      drop = FALSE
    ) +
    guides(
      size = guide_legend(title = "No. tested")
    ) +
    theme_minimal() +
    ggtitle(
      sprintf("Data availability for %s",
              this_insecticide)
    )
  
  ggsave(sprintf("figures/data_time_%s.png",
                 this_insecticide),
         bg = "white",
         scale = 0.8,
         width = 8,
         height = 8)
  
}


