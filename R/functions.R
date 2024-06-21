# misc functions

# list of the countries in MTM data that are in malaria-endemic Africa
mtm_africa_countries <- function() {
  africa_countries <- c("Burkina Faso",
                        "Botswana",
                        "Mali",
                        "Kenya",
                        "Malawi",
                        "Cameroon",
                        "Niger",
                        "Benin",
                        "Côte d’Ivoire",
                        "Togo",
                        "Nigeria",
                        "United Republic of Tanzania",
                        "Ethiopia",
                        "Mozambique",
                        "Uganda",
                        "Ghana",
                        "Senegal",
                        "Democratic Republic of the Congo",
                        "Somalia",
                        "South Africa",
                        "Zambia",
                        "Eritrea",
                        "Sudan",
                        "Sierra Leone",
                        "Angola",
                        "Namibia",
                        "Djibouti",
                        "Rwanda",
                        "Chad",
                        "Liberia",
                        "Zimbabwe",
                        "Eswatini",
                        "Guinea",
                        "Equatorial Guinea",
                        "Congo",
                        "Central African Republic",
                        "Madagascar",
                        "Comoros",
                        "Sao Tome and Principe",
                        "Mauritania",
                        "Cabo Verde",
                        "Gabon",
                        "Burundi",
                        "Guinea-Bissau",
                        "Gambia")
}


# given the percent mortality, infer the sample size as the integer closest to
# 100 that implies an integer number that died. If one can't be found, within
# the limits, set to the default instead
infer_sample_size <- function(mortality_percent,
                              limits = c(50, 200),
                              default = 100,
                              tol = 1e-12) {
  
  # try integer values between the limits
  possible_integers <- seq(limits[1], limits[2])
  distance <- abs(possible_integers - default)
  
  # for all of these, see if the mortality rate gives an integer
  died <- possible_integers %*% t(mortality_percent / 100)
  is_integer <- abs(died - round(died)) < tol
  
  # find the value closest to 100 that gives an integer
  valid_distance <- sweep(ifelse(is_integer, 0, Inf),
                          1,
                          distance,
                          FUN = "+")
  
  # for any that don't have a match, put in the index for 100
  missing_index <- which(possible_integers == default)
  
  # function to find the smallest, without dropping values
  closest <- function(distances) {
    closest_index <- which.min(distances)
    if (length(closest_index) == 0) {
      closest_index <- missing_index
    } 
    closest_index
  }

  closest_index <- apply(valid_distance, 2, closest)
  
  # missing <- !apply(is_integer, 2, any)
  # closest_index[missing] <- missing_index
  
  # look up the corresponding integers
  possible_integers[closest_index]
  
}

# compute the mode of a sample
sample_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# given a vector of sample sizes N, predicted proportions p, and overdispersion
# parameter rho (on unit interval, where 0 means no overdispersioan and 1 means
# maximum), returna betatbinomial-distributed variable greta array
betabinomial_p_rho <- function(N, p, rho) {
  
  # model the observation (betabinomial) sd as a multiplier on the binomial sd,
  # accounting for additional error due to nonindependent sampling of individuals
  # from the population. This is based on the INLA parameterisation
  
  # solve for a and b:
  #   p = a / (a + b)
  #   rho = 1 / (a + b + 1)
  a <- p * (1 / rho - 1)
  b <- a * (1 - p) / p

  # define betabinomial according to the greta interface
  beta_binomial(size = N, alpha = a, beta = b)
  
}

nearish <- function(x, y) {
  dplyr::near(x, y, tol = 1e-2)
}

# given a raster file, force it to be written to disk and not in memory (save to
# disk and reload). filename is the name fo the file to write it to, dots is any
# other arguments for writeraster
force_to_disk <- function(raster, filename = tempfile(fileext = ".tif"), ...) {
  terra::writeRaster(x = raster,
                     filename = filename,
                     overwrite = TRUE,
                     ...)
  terra::rast(filename)
}

# random RNG seed (copied from greta, used to fix downstream seeds, conditional
# on top-level randomness)
get_seed <- function ()  {
  sample.int(n = 2^30, size = 1)
}

# get rid of all the craps in plotting the IR maps
theme_ir_maps <- function() {
  theme_minimal() +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          # legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
}

