#' Plot simulated trials and calibrations
#'
#' `plot()` for a `PipeSimulation` shows, by dose combination, the percentage
#' of trials that selected it (`type = "selected"`, the final dose set) or the
#' mean number of patients treated (`type = "patients"`). Two-drug grids are
#' drawn as a tile map with the values printed; grids with other numbers of
#' drugs as a bar per combination. For a `PipeCalibration` it plots FDR and TPR
#' against `eps`, with the target `alpha` and the chosen `eps` marked.
#'
#' These methods need the ggplot2 package (in Suggests) and return a ggplot
#' object, which can be changed or saved in the usual way.
#'
#' @param x A `PipeSimulation` from [simulate()] or a `PipeCalibration` from
#'   [calibrate_eps()].
#' @param type `"selected"` or `"patients"` (simulations only).
#' @param ... Not used.
#' @return A ggplot object.
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   g <- DoseGrid(levels = c(A = 3, B = 3))
#'   des <- ms_design(g, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 9, N = 15)
#'   sim <- simulate(des, nsim = 10, seed = 1, truth = seq(0.1, 0.5, length.out = 9))
#'   plot(sim)
#'   plot(sim, type = "patients")
#' }
#' }
#' @name plot.PipeSimulation
#' @aliases plot.PipeCalibration
NULL

utils::globalVariables(".data")

need_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Plotting needs the ggplot2 package; install it first.", call. = FALSE)
  }
}

method(plot, PipeSimulation) <- function(x, type = c("selected", "patients"), ...) {
  need_ggplot2()
  type <- match.arg(type)
  des <- x@design
  grid <- des@grid
  val <- if (type == "selected") 100 * colMeans(x@final_set) else colMeans(x@n)
  lab <- if (type == "selected") "% of trials selecting" else "Mean patients"
  df <- data.frame(combo = combo_labels(grid), value = val)
  if (length(grid@levels) == 2L) {
    drugs <- names(grid@levels)
    df$a <- factor(grid@combos[, 1], levels = sort(unique(grid@combos[, 1])))
    df$b <- factor(grid@combos[, 2], levels = sort(unique(grid@combos[, 2])))
    ggplot2::ggplot(df, ggplot2::aes(x = .data$b, y = .data$a, fill = .data$value)) +
      ggplot2::geom_tile(colour = "white") +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", .data$value))) +
      ggplot2::scale_fill_gradient(low = "white", high = "steelblue", name = lab) +
      ggplot2::labs(x = drugs[2], y = drugs[1]) +
      ggplot2::theme_minimal()
  } else {
    df$combo <- factor(df$combo, levels = df$combo)
    ggplot2::ggplot(df, ggplot2::aes(x = .data$combo, y = .data$value)) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::labs(x = "Combination", y = lab) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
}

method(plot, PipeCalibration) <- function(x, ...) {
  need_ggplot2()
  tr <- x@trace
  long <- rbind(data.frame(eps = tr$eps, metric = "FDR", value = tr$fdr),
                data.frame(eps = tr$eps, metric = "TPR", value = tr$tpr))
  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$eps, y = .data$value, colour = .data$metric)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = x@alpha, linetype = "dashed") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "eps", y = NULL, colour = NULL)
  if (!is.na(x@eps)) p <- p + ggplot2::geom_vline(xintercept = x@eps, linetype = "dotted")
  p
}
