#' PIPE design
#'
#' A `PipeDesign` holds every setting of a PIPE dose-finding design for drug
#' combinations. Its validator checks the settings together, so an invalid
#' design cannot be built. Most users create one with [ms_design()] (Mander &
#' Sweeting 2015) or [cheung_design()] (Cheung et al. 2022), and then call
#' [estimate()], [next_dose()] or `simulate()`.
#'
#' Coding: a combination is "above" when its parameter exceeds `theta`
#' (toxic for `target = "MTD"`, effective for `target = "MED"`).
#'
#' @param grid A [DoseGrid].
#' @param target `"MTD"` or `"MED"`.
#' @param theta Target threshold.
#' @param outcome `"binary"` (Beta priors `prior_a`, `prior_b`) or
#'   `"normal"` (pooled-variance t posterior, with `sigma0` while there is no
#'   variance estimate).
#' @param prior_a,prior_b Beta prior parameters, one per combination or one
#'   for all.
#' @param sigma0 Known outcome SD used before the variance can be estimated.
#' @param rule `"cheung"` (epsilon-tapering, sampling from S1 u S0) or `"ms"`
#'   (Mander & Sweeting: admissible set, safety rule, RPII).
#' @param eps,taper Asymmetry of the gain and whether it tapers from 0.5 to
#'   `eps` over the trial (`"cheung"` only; `"ms"` uses 0.5).
#' @param weights `"n+1"` or `"1"`; `NULL` for the rule's default.
#' @param selection How the next dose is drawn: `"equal"` or `"prob"` for
#'   `"cheung"`; `"min_S"`, `"inverse_S"`, `"equal"` or `"inverse_sqrt_S"`
#'   for `"ms"`. `NULL` for the rule's default.
#' @param admis,constraint,safety M&S options: `"closest"` or `"adjacent"`;
#'   `"neighbouring"` or `"non-neighbouring"`; the safety threshold, or `NULL`
#'   for no safety rule.
#' @param N,cohort Maximum number of patients and cohort size.
#' @param start Index of the first combination; `NULL` for the default (top
#'   dose for Cheung MED designs, combination 1 otherwise).
#' @param delta,within_set Riviere's single pick: the interval half-width and
#'   whether to pick within the final set or among all tried doses.
#' @return A `PipeDesign` object.
#' @examples
#' g <- DoseGrid(levels = c(A = 3, B = 3))
#' PipeDesign(grid = g, target = "MTD", theta = 0.3, prior_a = 0.5, prior_b = 0.5,
#'            rule = "cheung", eps = 0.25, N = 30)
#' @name PipeDesign
#' @export
PipeDesign := new_class(
  package = "gPIPE",
  properties = list(
    grid = DoseGrid,
    target = new_property(class_character, default = "MTD"),
    theta = class_double,
    outcome = new_property(class_character, default = "binary"),
    prior_a = new_property(class_double, default = 1),
    prior_b = new_property(class_double, default = 1),
    sigma0 = new_property(class_double, default = 1),
    rule = new_property(class_character, default = "cheung"),
    eps = new_property(class_double, default = 0.5),
    taper = new_property(class_logical, default = TRUE),
    weights = new_property(NULL | class_character, default = NULL),
    selection = new_property(NULL | class_character, default = NULL),
    admis = new_property(class_character, default = "closest"),
    constraint = new_property(class_character, default = "neighbouring"),
    safety = new_property(NULL | class_double, default = 0.8),
    N = class_numeric,
    cohort = new_property(class_numeric, default = 1),
    start = new_property(NULL | class_numeric, default = NULL),
    delta = new_property(class_double, default = 0.1),
    within_set = new_property(class_logical, default = TRUE)
  ),
  validator = function(self) {
    tryCatch({
      design_list(self)
      NULL
    }, error = conditionMessage)
  }
)

# The plain list the trial loop runs on, built once from the S7 object.
design_list <- function(design) {
  K <- n_combos(design@grid)
  outcome <- match.arg(design@outcome, c("binary", "normal"))
  if (length(design@theta) != 1L) {
    stop("`theta` must be a single number.", call. = FALSE)
  }
  model <- if (outcome == "binary") {
    for (prior in list(design@prior_a, design@prior_b)) {
      if (!length(prior) %in% c(1L, K) || any(prior <= 0)) {
        stop("Beta priors must be positive, one per combination or one for all.", call. = FALSE)
      }
    }
    outcome_beta(design@theta, design@prior_a, design@prior_b)
  } else {
    outcome_normal(design@theta, design@sigma0)
  }
  safety <- if (design@rule == "ms") design@safety else NULL
  pipe_trial_design(
    grid = design@grid, model = model, target = design@target, N = design@N,
    cohort = design@cohort, start = design@start, rule = design@rule, eps = design@eps,
    taper = design@taper, weights = design@weights, selection = design@selection,
    admis = design@admis, constraint = design@constraint, safety = safety,
    delta = design@delta, within_set = design@within_set
  )
}

#' Convenience constructors for the two PIPE designs
#'
#' `ms_design()` sets up the original PIPE design of Mander & Sweeting (2015)
#' for binary toxicity (MTD), with priors given by their medians and
#' strengths. `cheung_design()` sets up PIPE with epsilon-tapering of Cheung
#' et al. (2022), for MED (default) or MTD, with normal or binary outcomes.
#'
#' @inheritParams PipeDesign
#' @param prior_median,prior_strength Prior median of each combination's
#'   toxicity probability and prior strength a + b (one per combination or one
#'   for all). M&S call strength 1 "strong" and 1 / K "weak".
#' @param prior_method `"bisection"` reproduces pipe.design; `"exact"` solves
#'   the median exactly.
#' @param ... Other [PipeDesign] properties.
#' @return A [PipeDesign].
#' @examples
#' g <- DoseGrid(levels = c(A = 4, B = 4))
#' ms_design(g, theta = 0.2, prior_median = 0.2, prior_strength = 1 / 16, N = 50)
#' cheung_design(g, theta = 0.5, sigma0 = 1, eps = 0.05, N = 60)
#' @name pipe-designs
#' @export
ms_design <- function(grid, theta, prior_median, prior_strength, N, cohort = 1,
                      admis = "closest", selection = "min_S", constraint = "neighbouring",
                      safety = 0.8, prior_method = "bisection", ...) {
  pr <- beta_prior_from_median(rep_len(c(prior_median), n_combos(grid)), prior_strength,
                               method = prior_method)
  PipeDesign(grid = grid, target = "MTD", theta = theta, outcome = "binary",
             prior_a = pr$a, prior_b = pr$b, rule = "ms", N = N, cohort = cohort,
             admis = admis, selection = selection, constraint = constraint, safety = safety,
             ...)
}

#' @rdname pipe-designs
#' @export
cheung_design <- function(grid, theta, eps, N, target = "MED", outcome = "normal",
                          sigma0 = 1, prior_a = 1, prior_b = 1, taper = TRUE,
                          selection = "equal", cohort = 1, ...) {
  PipeDesign(grid = grid, target = target, theta = theta, outcome = outcome, sigma0 = sigma0,
             prior_a = prior_a, prior_b = prior_b, rule = "cheung", eps = eps, taper = taper,
             selection = selection, N = N, cohort = cohort, safety = NULL, ...)
}

# Trial data as sufficient statistics. `data` is a data frame with one row
# per patient, in the order treated, a column per drug (named as the grid's
# drugs) holding levels, and the outcome in `y`.
design_data <- function(design, data) {
  grid <- design@grid
  K <- n_combos(grid)
  if (is.null(data) || !nrow(data)) {
    return(list(n = integer(K), s1 = numeric(K), s2 = numeric(K), history = integer()))
  }
  drugs <- names(grid@levels)
  if (!all(c(drugs, "y") %in% names(data))) {
    stop(sprintf("`data` needs columns %s and `y`.", paste(drugs, collapse = ", ")), call. = FALSE)
  }
  levels <- as.matrix(data[drugs])
  if (!is.numeric(levels) || anyNA(levels) || any(levels != round(levels))) {
    stop("Dose levels in `data` must be whole numbers.", call. = FALSE)
  }
  if (!is.numeric(data$y) || anyNA(data$y)) {
    stop("`y` must be numeric with no missing values.", call. = FALSE)
  }
  if (design@outcome == "binary" && !all(data$y %in% c(0, 1))) {
    stop("`y` must be 0 or 1 for a binary outcome.", call. = FALSE)
  }
  k <- combo_index(grid, levels)
  if (anyNA(k)) {
    stop("`data` has dose combinations that are not in the grid.", call. = FALSE)
  }
  y <- data$y
  list(n = tabulate(k, K), s1 = vapply(seq_len(K), function(j) sum(y[k == j]), numeric(1)),
       s2 = vapply(seq_len(K), function(j) sum(y[k == j]^2), numeric(1)), history = k)
}

#' PIPE estimate and next dose
#'
#' `estimate()` gives the PIPE estimate from the data so far: which
#' combinations exceed `theta`, and the design's final dose set (MTD or MED
#' set for Cheung designs, recommended phase II doses for M&S designs).
#' `next_dose()` gives the doses the design may give next, with their
#' probabilities, and draws one.
#'
#' @param design A [PipeDesign].
#' @param data A data frame with one row per patient in the order treated: a
#'   column per drug (named as in the grid) with the dose level, and `y` with
#'   the outcome. `NULL` before any patient.
#' @param ... Not used.
#' @return `estimate()`: a `PipeEstimate` with properties `p` (posterior
#'   probability of exceeding `theta`), `gamma` (0/1 estimate), `q` (M&S
#'   safety probabilities, or `NULL`), `set` and `pick` (Riviere's single dose,
#'   binary MTD designs only), indices into the grid's combinations.
#'   `next_dose()`: a list with `dose` (index, `NA` when the trial stops),
#'   `label`, `stopped`, and `candidates`, a data frame of the candidate doses
#'   and their probabilities.
#' @examples
#' g <- DoseGrid(levels = c(A = 3, B = 3))
#' des <- ms_design(g, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 9, N = 30)
#' dat <- data.frame(A = c(1, 1, 2, 2), B = c(1, 1, 2, 2), y = c(0, 0, 0, 1))
#' estimate(des, dat)
#' next_dose(des, dat)
#' @usage estimate(design, data = NULL, ...)
#' next_dose(design, data = NULL, ...)
#' @name estimate
#' @export estimate
estimate := new_generic(
  dispatch_args = "design",
  fun = function(design, data = NULL, ...) S7_dispatch()
)

#' @rdname estimate
#' @name next_dose
#' @export next_dose
next_dose := new_generic(
  dispatch_args = "design",
  fun = function(design, data = NULL, ...) S7_dispatch()
)

PipeEstimate := new_class(
  package = "gPIPE",
  properties = list(
    design = PipeDesign,
    p = class_double,
    gamma = class_integer,
    q = NULL | class_double,
    set = class_integer,
    pick = class_integer
  )
)

method(estimate, PipeDesign) <- function(design, data = NULL) {
  d <- design_list(design)
  x <- design_data(design, data)
  fin <- trial_final(d, x$n, x$s1, x$s2, x$history)
  PipeEstimate(design = design, p = fin$p, gamma = as.integer(fin$gamma), q = fin$q,
               set = as.integer(fin$final_set), pick = as.integer(fin$pick))
}

method(next_dose, PipeDesign) <- function(design, data = NULL) {
  d <- design_list(design)
  x <- design_data(design, data)
  labels <- combo_labels(design@grid)
  if (sum(x$n) >= d$N) {
    stop("The trial has already enrolled `N` patients.", call. = FALSE)
  }
  prob <- numeric(n_combos(design@grid))
  if (!sum(x$n)) {
    prob[d$start] <- 1
    stopped <- FALSE
  } else {
    step <- trial_step(d, x$n, x$s1, x$s2, x$history)
    prob <- step$prob
    stopped <- step$stopped
  }
  cand <- which(prob > 0)
  dose <- if (stopped) NA_integer_ else draw_dose(cand, prob[cand])
  list(dose = dose, label = labels[dose], stopped = stopped,
       candidates = data.frame(dose = cand, label = labels[cand], prob = prob[cand]))
}

method(print, PipeEstimate) <- function(x, ...) {
  grid <- x@design@grid
  labels <- combo_labels(grid)
  cat(sprintf("PIPE estimate (%s, theta = %g)\n", x@design@target, x@design@theta))
  if (length(grid@levels) == 2L) {
    cat("gamma (1 = above theta), rows = ", names(grid@levels)[1], " levels:\n", sep = "")
    print(grid_as_array(grid, x@gamma))
  }
  what <- if (x@design@rule == "ms") "Recommended phase II doses" else paste(x@design@target, "set")
  cat(what, ": ", if (length(x@set)) paste(labels[x@set], collapse = ", ") else "none", "\n", sep = "")
  if (!is.na(x@pick)) cat("Single pick: ", labels[x@pick], "\n", sep = "")
  invisible(x)
}

PipeSimulation := new_class(
  package = "gPIPE",
  properties = list(
    design = PipeDesign,
    summary = class_data.frame,
    n = class_double,
    s1 = class_double,
    gamma = class_double,
    final_set = class_logical,
    trials = class_list
  )
)

#' Simulate PIPE trials
#'
#' Runs `nsim` trials of a [PipeDesign] under a true dose-outcome scenario.
#' This is a method for [stats::simulate()], which gPIPE re-exports.
#'
#' @param object A [PipeDesign].
#' @param nsim Number of trials.
#' @param seed Optional seed, passed to [set.seed()].
#' @param ... For a [PipeDesign]:
#'   * `truth` (required): the true scenario. For binary outcomes the
#'     probabilities, as a vector over the combinations or an array over the
#'     full grid; for normal outcomes a list with `mean` (vector or array) and
#'     `sd`.
#'   * `keep_trials`: keep each trial's full result in `@trials` (default
#'     `FALSE`).
#' @return A `PipeSimulation` with `@summary` (one row per trial: patients
#'   enrolled, early stop, single pick), and trials x combinations matrices
#'   `@n` (patients), `@s1` (events or outcome sums), `@gamma` (final
#'   estimate) and `@final_set` (membership of the final dose set).
#' @examples
#' g <- DoseGrid(levels = c(A = 2, B = 2))
#' des <- ms_design(g, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 4, N = 12)
#' sim <- simulate(des, nsim = 5, seed = 1, truth = c(0.1, 0.2, 0.3, 0.5))
#' sim
#' @usage simulate(object, nsim = 1, seed = NULL, ...)
#' @name simulate.PipeDesign
#' @aliases simulate
#' @importFrom stats simulate
#' @rawNamespace export(simulate)
method(simulate, PipeDesign) <- function(object, nsim = 1, seed = NULL, ..., truth,
                                         keep_trials = FALSE) {
  grid <- object@grid
  as_vec <- function(x) if (length(dim(x)) > 1L) grid_from_array(grid, x) else c(x)
  truth <- if (is.list(truth)) list(mean = as_vec(truth$mean), sd = truth$sd) else as_vec(truth)
  sim <- pipe_simulate(design_list(object), truth, n_trials = nsim, seed = seed,
                       keep_trials = keep_trials)
  storage.mode(sim$n) <- storage.mode(sim$gamma) <- "double"
  PipeSimulation(design = object, summary = sim$summary, n = sim$n, s1 = sim$s1,
                 gamma = sim$gamma, final_set = sim$final_set, trials = sim$trials)
}

method(print, PipeSimulation) <- function(x, ...) {
  grid <- x@design@grid
  cat(sprintf("%d simulated PIPE trials (%s, rule = %s)\n", nrow(x@summary), x@design@target,
              x@design@rule))
  cat(sprintf("Mean patients: %.1f; stopped early: %.1f%%\n", mean(x@summary$n_enrolled),
              100 * mean(x@summary$stopped)))
  sel <- 100 * colMeans(x@final_set)
  if (length(grid@levels) == 2L) {
    cat("% of trials with each combination in the final set (rows = ", names(grid@levels)[1],
        " levels):\n", sep = "")
    print(round(grid_as_array(grid, sel), 1))
  }
  invisible(x)
}
