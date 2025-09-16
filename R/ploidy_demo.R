# Compare the implications of ploidy model detail on temporal trend in IR.

# Compute the proportion with the resistant allele (and phenotype) at the next
# timestep under haploidy, given the proportion from the previous timestep and
# the selection coefficient 's' (the fitness advantage 'w' minus 1)
haploid_next <- function(p, s = 0.1) {
  q <- 1 - p
  w <- 1 + s
  p * w / (p * w + q)
}

# Compute the proportion with the resistance allele at the next timestep under
# diploidy, given the proportion from the previous timestep, the selection
# coefficient 's' (the fitness advantage 'w' minus 1), and the degree of
# heterozygote dominance 'h' (whether heterozygotes are more like the homozygote
# resistant phenotype or the homozygote susceptible phenotype regarding
# fitness), where h = 1 means heterozygotes have identical phenotype to
# homozygote resistant, h = 0 means heterozygotes have identical phenotype to
# homozygote susceptible, and h = 0.5 menans they are intermediate. Ref here for
# a nice example notation and discussion re. mossies:
# https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0001387
diploid_next <- function(p, s = 0.2, h = 0.5) {
  q <- 1 - p
  numerator <- p ^ 2 * (1 + s) + p * q * (1 + h * s)
  denominator <- 1 + s * (p ^ 2 + 2 * h * ( p * q))
  numerator / denominator
}

# fit haploid parameters that most closely approximate the output of the
# different diploid models

# given a timeseries of proportions outputs, compute the selection coefficient
# and initial state of the haploidy model that best fits, using the logit
# approximation
estimate_haploid_s <- function(proportions) {
  logits <- qlogis(proportions)
  keep <- proportions > 0.25 & proportions < 0.75
  y <- logits[keep]
  x <- seq_along(y)
  slope <- coef(lm(y ~ x))[2]
  w <- exp(slope)
  s <- w - 1
  
  # work back to the initial state using this logistic approximation, using the value closest to 0.5
  t_midpoint <- which.min(abs(logits))
  logit_init <- logits[t_midpoint] - t_midpoint * slope
  init <- plogis(logit_init)
  list(s = s,
       init = init)
}

# iterate these recursion functions to solve for the proportion with the
# resistant allele, over time. 'next_state_function' is a function with the
# first argument giving the proportion with the resistant allele at the previous
# timestep and returning the proportion at this timestep, optionally with other
# arguments. 'p_init' is the initial fraction with the allele. 'n_times' is the
# number of timesteps to iterate. Anything provided via the dots argument is
# passed as an additional argument to next_state_function
solve_selection <- function(next_state_function,
                            p_init = 0.001,
                            n_times = 100,
                            ...) {
  
  result <- rep(NA, n_times)
  result[1] <- p <- p_init
  for (time in 2:n_times) {
    result[time] <- p <- next_state_function(p, ...)
  }
  
  result
    
}

# convert proportions with the resistant allele to the proportions expressing
# the resistant phenotype, given the specified level of heterozygote dominance.
# If h = NA, assume a haploid model where the proportion with the allele *is*
# the proportion with the phenotype
apply_phenotype <- function(proportion_resistant_allele, h = NA) {
  
  if (is.na(h)) {
    proportion_resistant_phenotype <- proportion_resistant_allele
  } else {
    # otherwise, compute the proportions homozygotic resistant (express
    # phenotype), the proportion homozygotic sysceptible (don't express the
    # phenotype), and the proportion heterozygotic (fraction h express the
    # phenotype)
    proportion_susceptible_allele <- 1 - proportion_resistant_allele
    proportion_hom_resistant <- proportion_resistant_allele ^ 2
    proportion_hom_susceptible <- proportion_susceptible_allele ^ 2
    proportion_het <- 2* proportion_resistant_allele * proportion_susceptible_allele
    proportion_resistant_phenotype <- proportion_hom_resistant * 1 +
      proportion_hom_susceptible * 0 +
      proportion_het * h
  }
  
  proportion_resistant_phenotype
  
}


# compute resistant allele proportions under diploidy for a range of values for
# heterozygote dominance
diploid_p_h01 <- solve_selection(diploid_next, h = 0.1)
diploid_p_h025 <- solve_selection(diploid_next, h = 0.25)
diploid_p_h05 <- solve_selection(diploid_next, h = 0.5)
diploid_p_h075 <- solve_selection(diploid_next, h = 0.75)
diploid_p_h1 <- solve_selection(diploid_next, h = 1)

diploid_pheno_h01 <- apply_phenotype(diploid_p_h01, h = 0.1)
diploid_pheno_h025 <- apply_phenotype(diploid_p_h025, h = 0.25)
diploid_pheno_h05 <- apply_phenotype(diploid_p_h05, h = 0.5)
diploid_pheno_h075 <- apply_phenotype(diploid_p_h075, h = 0.75)
diploid_pheno_h1 <- apply_phenotype(diploid_p_h1, h = 1)

# compute proportions under haploidy
haploid_p <- solve_selection(haploid_next)
haploid_pheno <- apply_phenotype(haploid_p)

hs <- c(0.25, 0.5, 0.75, 1)
diploid_s <- 0.35

png("figures/ploidy_demo.png",
    width = 960,
    height = 960,
    pointsize = 20)
par(mfrow = n2mfrow(length(hs)))
for (h in hs) {
  # solve diploid model for phenotype
  diploid_p <- solve_selection(diploid_next, s = diploid_s, h = h)
  diploid_pheno <- apply_phenotype(diploid_p, h = h)
  # approximate phenotype with haploidy
  approx <- estimate_haploid_s(diploid_pheno)
  haploid_p <- solve_selection(haploid_next,
                               s = approx$s,
                               p_init = approx$init)
  # plot both
  plot(diploid_p,
       type = "l",
       lwd = 3,
       col = "grey",
       main = sprintf("Diploid het. dominance: h=%s\nblack = resistant phenotype,\ngrey = resistant genotype,\ndashed = haploid approximation", h),
       ylim = c(0, 1),
       xlab = "generations",
       ylab = "frequency")
  lines(diploid_pheno,
        lwd = 3)
  lines(haploid_p,
        lwd = 3,
        lty = 2)
}
dev.off()
