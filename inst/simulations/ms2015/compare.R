# Compare the stored gPIPE results with pipe.design 0.5.1 and with the
# printed Tables III and IV of Mander & Sweeting (2015).
# Usage: Rscript compare.R [results folder = results]
#
# z is the difference in Monte Carlo standard errors (SE).
# - Against pipe.design: SE of the difference of two independent runs.
# - Against the paper: the paper gives no SE, so its SE is taken from the
#   gPIPE run, scaled to the paper's number of trials, plus the rounding of
#   the printed value (to 1 point, or 0.1 doses).
ms2015_compare <- function(dir = "results") {
  key <- c("table", "scenario", "admis", "selection", "metric", "band")
  read <- function(f) utils::read.csv(file.path(dir, f), colClasses = c(scenario = "character", band = "character"))
  g <- read("ms2015_gpipe.csv")
  p <- read("ms2015_pipe_design.csv")
  paper <- read("ms2015_paper.csv")
  z <- function(diff, se) ifelse(diff == 0, 0, diff / se)

  a <- merge(g, p, by = key, suffixes = c("_gpipe", "_pipe_design"))
  a$diff <- a$estimate_gpipe - a$estimate_pipe_design
  a$z <- z(a$diff, sqrt(a$se_gpipe^2 + a$se_pipe_design^2))

  b <- merge(g, paper, by = key, suffixes = c("_gpipe", "_paper"))
  round_sd <- ifelse(b$metric == "mean_n_rec", 0.1, 1) / sqrt(12)
  b$diff <- b$estimate_gpipe - b$estimate_paper
  b$z <- z(b$diff, sqrt(b$se^2 * (1 + b$nsim_gpipe / b$nsim_paper) + round_sd^2))
  list(pipe_design = a, paper = b)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  cmp <- ms2015_compare(if (length(args)) args[1] else "results")
  for (ref in names(cmp)) {
    x <- cmp[[ref]]
    cat(sprintf("\ngPIPE against %s: %d cells\n", ref, nrow(x)))
    cat(sprintf("  largest |difference|: %.2f points (percentages), %.3f doses (mean number recommended)\n",
                max(abs(x$diff[x$metric != "mean_n_rec"])), max(abs(x$diff[x$metric == "mean_n_rec"]))))
    cat(sprintf("  |z| > 2: %d cells (%.1f expected by chance), |z| > 3: %d, largest |z|: %.2f\n",
                sum(abs(x$z) > 2), 0.0455 * sum(x$diff != 0), sum(abs(x$z) > 3), max(abs(x$z))))
    cat(sprintf("  chi-square of z over %d non-zero cells: %.1f\n", sum(x$diff != 0), sum(x$z^2)))
    worst <- x[order(-abs(x$z)), ][1:8, ]
    print(worst[c("table", "scenario", "admis", "selection", "metric", "band",
                  grep("^estimate", names(x), value = TRUE), "diff", "z")], row.names = FALSE, digits = 3)
  }
}
