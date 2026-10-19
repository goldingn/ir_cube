# Helpers shared by the two-stage scripts (#21). Functions only; sources
# R/validation_scoring.R (functions and settings) for the rules the scoring
# uses: thin_draws(), rho_lookup() and string_seed(); and read_assay_draws().

source("R/validation_scoring.R")

# timestamped progress line
report <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", sprintf(...), "\n")
  flush(stdout())
}

# peak resident memory of this process so far, from the kernel
peak_memory_gb <- function() {
  status <- readLines("/proc/self/status")
  as.numeric(gsub("\\D", "", grep("^VmHWM", status, value = TRUE))) / 1024 ^ 2
}
