# Compare the implications of ploidy model detail on temporal trend in IR.

# Compute the proportion with the resistance phenotype at the next timestep under haploidy,
# given the proportion from the previous timestep and the selection coefficient 's'
# (the fitness advantage 'w' minus 1)
haploid_next <- function(p, s = 0.1) {
  q <- 1 - p
  w <- 1 + s
  p * w / (p * w + q)
}

# Compute the proportion with the resistance phenotype at the next timestep
# under diploidy, given the proportion from the previous timestep, the selection
# coefficient 's' (the fitness advantage 'w' minus 1), and the degree of
# heterozygote dominance 'h' (whether heterozygotes are more like the homozygote
# resistant phenotype or the homozygote susceptible phenotype), where h = 1
# means heterozygotes have identical phenotype to homozygote resistant, h = 0
# means heterozygotes have identical phenotype to homozygote susceptible, and h
# = 0.5 menans they are intermediate. Ref here for a nice example notation and
# discussion re. mossies:
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
# of the haploidy model that best fits, using the logit approximation
estimate_haploid_s <- function(proportions) {
  logits <- qlogis(proportions)
  intercept <- logits[1]
  times <- seq_along(proportions) - 1
  not_first <- times != times[1]
  slope <- mean((logits[not_first] - intercept) / times[not_first])
  w <- exp(slope)
  s <- w - 1
  s
}

# iterate these recursion functions to solve for the proportion expressing the
# phenotype, over time. 'next_state_function' is a function with the first
# argument giving the proportion expressing the phenotype at the previous
# timestep and returning the proportion at this timestep, optionally with other
# arguments. 'p_init' is the initial fraction expressing the phenotype.
# 'n_times' is the number of timesteps to iterate. Anything provided via the
# dots argument is passed as an additional argument to next_state_function
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

# compute phenotype proportions under diploidy for a range of values for
# heterozygote dominance
diploid_p_h01 <- solve_selection(diploid_next, h = 0.1)
diploid_p_h025 <- solve_selection(diploid_next, h = 0.25)
diploid_p_h05 <- solve_selection(diploid_next, h = 0.5)
diploid_p_h075 <- solve_selection(diploid_next, h = 0.75)
diploid_p_h1 <- solve_selection(diploid_next, h = 1)

# compute phenotype proportions under haploidy
haploid_p <- solve_selection(haploid_next)

hs <- c(0.25, 0.5, 0.75, 1)
diploid_s <- 0.35

png("figures/ploidy_demo.png",
    width = 960,
    height = 960,
    pointsize = 20)
par(mfrow = n2mfrow(length(hs)))
for (h in hs) {
  # solve diploid model
  diploid_p <- solve_selection(diploid_next, s = diploid_s, h = h)
  # approximate with haploidy
  similar_haploid_s <- estimate_haploid_s(diploid_p)
  haploid_p <- solve_selection(haploid_next, s = similar_haploid_s)
  # plot both
  plot(diploid_p,
       type = "l",
       lwd = 3,
       main = sprintf("Diploid het. dominance: h=%s\n(dashed = haploid approx.)", h),
       ylim = c(0, 1),
       xlab = "generations",
       ylab = "proportion resistant")
  lines(haploid_p,
        lwd = 3,
        lty = 2)
}
dev.off()
