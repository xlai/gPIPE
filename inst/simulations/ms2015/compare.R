# Compare the stored gPIPE results with pipe.design 0.5.1 and with the
# printed Tables III and IV of Mander & Sweeting (2015).
# Usage: Rscript compare.R [results folder = results]
#
# The pass criterion is one chi-square test per comparison. For each
# configuration, d is the vector of differences over its cells and V the
# Monte Carlo covariance of d. The statistic d' V^-1 d is chi-square with
# one degree of freedom per cell when the two sources agree. The bands of
# one metric sum to 100, so one band per metric is left out, and so are
# cells that are 0 in every trial. Configurations are independent runs, so
# their statistics and degrees of freedom add up to one test.
# - Against pipe.design: V is the sum of the two runs' covariances.
# - Against the paper: the paper gives no covariance, so V is the gPIPE
#   covariance scaled to the paper's number of trials, plus the gPIPE run's
#   own, plus the rounding of the printed value (to 1 point, or 0.1 doses).
#   This test is approximate.
#
# z, the difference of one cell in standard errors, is kept for reading the
# tables. It is not the criterion: with 239 cells some are large by chance.
ms2015_compare <- function(dir = "results") {
  key <- c("table", "scenario", "admis", "selection")
  read <- function(f) {
    utils::read.csv(file.path(dir, f), colClasses = c(scenario = "character"), na.strings = character())
  }
  id <- function(x) paste(x$table, x$scenario, x$admis, x$selection)
  cov_matrix <- function(v, k) {
    m <- matrix(0, k, k)
    m[cbind(v$cell_i, v$cell_j)] <- v$cov
    m[cbind(v$cell_j, v$cell_i)] <- v$cov
    m
  }
  # Cells that enter the chi-square of one configuration. A cell with no
  # Monte Carlo variance (0 in every trial of both runs) carries no
  # information and is left out. With `drop_last`, so is the last remaining
  # band of each percentage metric, because the bands of a metric sum to 100.
  free_cells <- function(metric, variance, drop_last) {
    keep <- variance > 0
    if (drop_last) {
      for (m in c("exp", "rec")) {
        i <- which(keep & metric == m)
        if (length(i)) keep[i[length(i)]] <- FALSE
      }
    }
    keep
  }
  wald <- function(d, v, keep) {
    stopifnot(all(d[!keep & diag(v) == 0] == 0))
    d <- d[keep]
    c(chisq = sum(d * solve(v[keep, keep, drop = FALSE], d)), df = sum(keep))
  }
  g <- read("ms2015_gpipe.csv")
  p <- read("ms2015_pipe_design.csv")
  paper <- read("ms2015_paper.csv")
  gv <- read("ms2015_gpipe_cov.csv")
  pv <- read("ms2015_pipe_design_cov.csv")
  stopifnot(identical(g[c(key, "metric", "band")], p[c(key, "metric", "band")]))

  cells_pd <- cbind(g[c(key, "metric", "band")], estimate_gpipe = g$estimate,
                    estimate_pipe_design = p$estimate, diff = g$estimate - p$estimate)
  se <- sqrt(g$se^2 + p$se^2)
  cells_pd$z <- ifelse(cells_pd$diff == 0, 0, cells_pd$diff / se)

  cell_id <- function(x) paste(id(x), x$metric, x$band)
  in_paper <- match(cell_id(g), cell_id(paper))
  cells_paper <- cbind(g[c(key, "metric", "band")], estimate_gpipe = g$estimate,
                       estimate_paper = paper$estimate[in_paper],
                       nsim_paper = paper$nsim[in_paper])[!is.na(in_paper), ]
  cells_paper$diff <- cells_paper$estimate_gpipe - cells_paper$estimate_paper
  round_var <- ifelse(cells_paper$metric == "mean_n_rec", 0.1, 1)^2 / 12
  g_paper <- g[!is.na(in_paper), ]
  se <- sqrt(g_paper$se^2 * (1 + g_paper$nsim / cells_paper$nsim_paper) + round_var)
  cells_paper$z <- cells_paper$diff / se

  by_config <- lapply(unique(id(g)), function(cfg) {
    rows <- id(g) == cfg
    vg <- cov_matrix(gv[id(gv) == cfg, ], sum(rows))
    vp <- cov_matrix(pv[id(pv) == cfg, ], sum(rows))
    metric <- g$metric[rows]
    v <- vg + vp
    a <- wald(g$estimate[rows] - p$estimate[rows], v, free_cells(metric, diag(v), drop_last = TRUE))
    # The printed values are rounded one by one, so they need not sum to 100
    # and every band stays in.
    has <- !is.na(in_paper[rows])
    scale <- 1 + g$nsim[rows][has] / paper$nsim[in_paper[rows][has]]
    rv <- ifelse(metric[has] == "mean_n_rec", 0.1, 1)^2 / 12
    vb <- vg[has, has] * scale + diag(rv)
    b <- wald(g$estimate[rows][has] - paper$estimate[in_paper[rows][has]], vb,
              free_cells(metric[has], diag(vg)[has], drop_last = FALSE))
    data.frame(g[which(rows)[1], key], chisq_pipe_design = a[["chisq"]], df_pipe_design = a[["df"]],
               chisq_paper = b[["chisq"]], df_paper = b[["df"]],
               row.names = NULL)
  })
  by_config <- do.call(rbind, by_config)
  total <- function(chisq, df) {
    c(chisq = sum(chisq), df = sum(df), p = stats::pchisq(sum(chisq), sum(df), lower.tail = FALSE))
  }
  list(pipe_design = cells_pd, paper = cells_paper, by_config = by_config,
       test_pipe_design = total(by_config$chisq_pipe_design, by_config$df_pipe_design),
       test_paper = total(by_config$chisq_paper, by_config$df_paper))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  cmp <- ms2015_compare(if (length(args)) args[1] else "results")
  for (ref in c("pipe_design", "paper")) {
    x <- cmp[[ref]]
    test <- cmp[[paste0("test_", ref)]]
    cat(sprintf("\ngPIPE against %s: %d cells\n", ref, nrow(x)))
    cat(sprintf("  chi-square %.1f on %d degrees of freedom, p = %.3f\n",
                test[["chisq"]], test[["df"]], test[["p"]]))
    cat(sprintf("  largest |difference|: %.2f points (percentages), %.3f doses (mean number recommended)\n",
                max(abs(x$diff[x$metric != "mean_n_rec"])), max(abs(x$diff[x$metric == "mean_n_rec"]))))
    cat(sprintf("  cells with |z| > 2: %d (%.1f expected by chance), largest |z|: %.2f\n",
                sum(abs(x$z) > 2), 0.0455 * sum(x$diff != 0), max(abs(x$z))))
  }
  cat("\nBy configuration:\n")
  print(cmp$by_config, row.names = FALSE, digits = 3)
}
