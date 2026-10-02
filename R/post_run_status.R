# Summarise an unattended post-processing run (docker/post_pod.sh in the
# RunPod image repository) as STATUS.md: each step's status and duration, and
# the key tables, from whatever outputs exist. Any table whose inputs are
# missing is noted and skipped.
#
#   Rscript R/post_run_status.R [run dir]
#
# Reads, under the run dir (default "."): status.tsv (written by the driver),
# outputs/mode_check/*.csv (R/chain_mode_check.R), outputs/posterior_summary.csv
# (R/summarise_posterior.R), outputs/sensitivity/sensitivity_*.csv
# (R/sensitivity_compare.R), outputs/cv_summary.csv and
# outputs/cv_variance_explained.csv (R/validation_metrics.R,
# R/variance_explained.R), outputs/two_stage/cv_headline_two_stage.csv
# (R/two_stage_metrics.R). Plain R only.

arguments <- commandArgs(trailingOnly = TRUE)
run_dir <- if (length(arguments) >= 1) arguments[1] else "."
path <- function(...) file.path(run_dir, ...)

md_table <- function(x, digits = 3) {
  if (is.null(x) || nrow(x) == 0) return("(none)")
  x <- as.data.frame(x)
  for (j in seq_along(x)) {
    if (is.numeric(x[[j]]) && all(x[[j]] == round(x[[j]]), na.rm = TRUE)) {
      x[[j]] <- ifelse(is.na(x[[j]]), "", format(x[[j]], big.mark = ""))
    } else if (is.numeric(x[[j]])) {
      x[[j]] <- ifelse(is.na(x[[j]]), "",
                       formatC(signif(x[[j]], digits), format = "fg",
                               digits = digits, flag = "#"))
      x[[j]] <- sub("\\.$", "", x[[j]])
    } else {
      x[[j]] <- ifelse(is.na(x[[j]]), "", as.character(x[[j]]))
    }
  }
  c(paste("|", paste(names(x), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
    apply(x, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}

read <- function(file) {
  if (!file.exists(file)) return(NULL)
  tryCatch(read.csv(file, stringsAsFactors = FALSE, check.names = FALSE,
                    encoding = "UTF-8"),
           error = function(e) NULL)
}

section <- function(title, body) {
  lines <- tryCatch(body, error = function(e) {
    paste("(failed to summarise:", conditionMessage(e), ")")
  })
  c(paste("##", title), "", lines, "")
}

out <- c("# Post-processing status", "",
         sprintf("Written %s, from %s.",
                 format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                 normalizePath(run_dir)), "")

out <- c(out, section("Steps", {
  status_file <- path("status.tsv")
  if (!file.exists(status_file)) {
    "(no status.tsv)"
  } else {
    s <- if (file.size(status_file) > 0) {
      read.delim(status_file, header = FALSE, stringsAsFactors = FALSE,
                 col.names = c("step", "status", "exit", "start", "end",
                               "minutes", "note"),
                 colClasses = "character", fill = TRUE, quote = "")
    }
    running <- list.files(path("running"))
    if (length(running)) {
      started <- format(file.mtime(file.path(path("running"), running)),
                        "%H:%M", tz = "UTC")
      s <- rbind(s, data.frame(step = running, status = "RUNNING", exit = "",
                               start = started, end = "", minutes = "",
                               note = ""))
    }
    c(md_table(s), "",
      sprintf("%i OK, %i failed, %i skipped, %i running.",
              sum(s$status == "OK"), sum(s$status == "FAILED"),
              sum(s$status == "SKIPPED"), sum(s$status == "RUNNING")),
      "Logs: logs/<step>.log. Memory over time: logs/memory.tsv.")
  }
}))

mode_files <- Sys.glob(path("outputs/mode_check/*.csv"))
mode_files <- mode_files[!grepl("_spread\\.csv$", mode_files)]
mode <- if (length(mode_files)) {
  do.call(rbind, lapply(mode_files, read))
}

out <- c(out, section("Mode check", {
  if (is.null(mode)) {
    "(no mode check results)"
  } else {
    c("Per-chain mean mortality floor; chains whose mean floor is more than",
      "0.05 above the lowest, and stuck chains, are dropped",
      "(R/drop_stuck_chains.R). Rank Rhat and ESS over every parameter.", "",
      md_table(mode[, c("fit", "chains", "floor_per_chain", "drop",
                        "rhat_all", "rhat_kept", "n_rhat_105_kept",
                        "ess_bulk_min_kept", "ess_tail_min_kept",
                        "top_spread")]))
  }
}))

out <- c(out, section("Main-fit diagnostics", {
  main <- mode[mode$fit == "r2_main", , drop = FALSE]
  if (is.null(mode) || nrow(main) == 0) {
    "(no mode check of the main fit)"
  } else {
    md_table(data.frame(
      quantity = c("chains kept", "floor (kept chains)", "worst rank Rhat",
                   "worst parameter", "params > 1.01", "params > 1.05",
                   "bulk ESS min", "bulk ESS median", "tail ESS min"),
      value = c(sprintf("%i of %i", main$chains -
                          lengths(strsplit(as.character(main$drop), ",")),
                        main$chains),
                sprintf("%.4f", main$floor_kept),
                sprintf("%.3f", main$rhat_kept), main$worst_kept,
                main$n_rhat_101_kept, main$n_rhat_105_kept,
                sprintf("%.0f", main$ess_bulk_min_kept),
                sprintf("%.0f", main$ess_bulk_median_kept),
                sprintf("%.0f", main$ess_tail_min_kept))))
  }
}))

out <- c(out, section("Posterior summary (main fit)", {
  x <- read(path("outputs/posterior_summary.csv"))
  if (is.null(x)) {
    "(no outputs/posterior_summary.csv)"
  } else {
    ci <- function(m, l, u, d = 3) sprintf("%.*f [%.*f, %.*f]", d, m, d, l, d, u)
    scalar <- x[x$parameter %in% c("mortality_floor", "reversion_rate"), ]
    t1 <- data.frame(parameter = scalar$element,
                     `mean [95% CI]` = ci(scalar$mean, scalar$lower,
                                          scalar$upper, 4),
                     `prior sd` = scalar$prior_sd,
                     contraction = scalar$contraction,
                     `rank Rhat` = scalar$rhat, `bulk ESS` = scalar$ess_bulk,
                     check.names = FALSE)
    rho <- x[x$parameter == "rho_types", ]
    t2 <- data.frame(parameter = rho$element,
                     model = ci(rho$mean, rho$lower, rho$upper, 2),
                     external = if ("external_rho" %in% names(rho))
                       ci(rho$external_rho, rho$external_lower,
                          rho$external_upper, 2) else NA)
    init <- x[x$parameter == "init_coef", ]
    t3 <- data.frame(parameter = init$element,
                     `mean [95% CI]` = ci(init$mean, init$lower, init$upper, 2),
                     contraction = init$contraction, check.names = FALSE)
    beta <- x[x$parameter == "beta_type", ]
    beta$covariate <- sub("^beta_type\\[(.*), [^,]*\\]$", "\\1", beta$element)
    t4 <- do.call(rbind, lapply(split(beta, factor(beta$covariate,
                                                   unique(beta$covariate))),
                                function(b) data.frame(
      covariate = b$covariate[1],
      `type means, range` = sprintf("%.2f to %.2f", min(b$mean), max(b$mean)),
      `95% CI below 0` = sprintf("%i of %i", sum(b$upper < 0), nrow(b)),
      `median contraction` = median(b$contraction), check.names = FALSE)))
    c(md_table(t1, 3), "", "rho per type:", "", md_table(t2), "",
      "Initial-state coefficients:", "", md_table(t3, 2), "",
      "Selection coefficients (log effect, beta_type), per covariate:", "",
      md_table(t4, 2), "",
      sprintf(paste("All %i parameters: worst rank Rhat %.3f, minimum bulk",
                    "ESS %.0f, minimum tail ESS %.0f."),
              nrow(x), max(x$rhat), min(x$ess_bulk), min(x$ess_tail)))
  }
}))

out <- c(out, section("Predicted against observed mortality (main fit)", {
  x <- read(path("outputs/sensitivity/sensitivity_in_sample.csv"))
  if (is.null(x)) {
    "(no outputs/sensitivity/sensitivity_in_sample.csv)"
  } else {
    c(paste("Posterior mean mortality (500 draws, predict.R's path) at the",
            "modelled assays' pixel-years against observed, per type."), "",
      md_table(x[x$fit == "main", setdiff(names(x), "fit")]), "",
      "Every fit, correlation per type:", "",
      md_table(reshape(x[, c("fit", "type", "r")], idvar = "type",
                       timevar = "fit", direction = "wide")))
  }
}))

out <- c(out, section("Cross-validation", {
  v <- read(path("outputs/cv_variance_explained.csv"))
  s <- read(path("outputs/cv_summary.csv"))
  h <- read(path("outputs/two_stage/cv_headline_two_stage.csv"))
  lines <- character(0)
  if (!is.null(v)) {
    lines <- c(lines, "Variance explained (%), outputs/cv_variance_explained.csv:",
               "", md_table(v[, intersect(c("experiment", "quantity", "kind",
                                            "estimate", "lower", "upper"),
                                          names(v))]), "")
  }
  if (!is.null(s)) {
    lines <- c(lines, "Scores, outputs/cv_summary.csv:", "",
               md_table(s[, intersect(c("model", "experiment", "n", "elpd",
                                        "crps", "mse", "skill", "coverage_50",
                                        "coverage_95", "mean_pit"),
                                      names(s))]), "")
  }
  if (!is.null(h)) {
    lines <- c(lines,
               "Two-stage, outputs/two_stage/cv_headline_two_stage.csv:", "",
               md_table(h[, intersect(c("experiment", "fold", "model", "n",
                                        "elpd", "diff_elpd", "crps",
                                        "explained", "explained_lower",
                                        "explained_upper", "diff_explained",
                                        "cover50", "cover95"), names(h))]),
               "")
  }
  if (length(lines) == 0) "(no cross-validation tables)" else lines
}))

out <- c(out, section("Sensitivity", {
  s <- read(path("outputs/sensitivity/sensitivity_summary.csv"))
  p <- read(path("outputs/sensitivity/sensitivity_parameters.csv"))
  lines <- character(0)
  if (!is.null(s)) {
    lines <- c(lines,
      paste("Against the main fit: mean |difference| in 2025 mortality and in",
            "the 2014-2024 change at the bioassay pixels (max over types, and",
            "that type), and over all map cells in 2025; map_noise_mean is",
            "the Monte Carlo noise expected in the map mean |difference|.",
            "cv_trigger: a bioassay-pixel difference above 0.03."), "",
      md_table(s), "")
  }
  if (!is.null(p)) {
    wide <- reshape(p[, c("fit", "parameter", "mean")], idvar = "parameter",
                    timevar = "fit", direction = "wide")
    names(wide) <- sub("^mean\\.", "", names(wide))
    lines <- c(lines, "Posterior means:", "", md_table(wide), "")
  }
  if (length(lines) == 0) "(no sensitivity tables)" else lines
}))

writeLines(out, path("STATUS.md.tmp"))
file.rename(path("STATUS.md.tmp"), path("STATUS.md"))
