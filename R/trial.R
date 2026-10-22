# Trial loop and simulation (#18, step 5).
#
# The loop runs on plain vectors: a design is a list read once at the start,
# and the data are integer and numeric vectors over the combinations. Two
# rules choose the next dose:
#
# - "cheung": Cheung et al. (2022) Algorithm 1 with optional eps-tapering.
#   After each cohort the PIPE estimate uses eps_n = taper_eps(n, N, eps),
#   with n the patients enrolled so far (or eps itself without tapering), and
#   the next dose is drawn from S1 u S0 with equal probability ("equal",
#   Algorithm 1-EQ) or proportional to the probability of the favourable
#   state ("prob", Algorithm 1-PP). The final estimate uses eps.
# - "ms": Mander & Sweeting (2015) as pipe.design 0.5.1 runs it, through
#   ms_decide(). The trial stops early only through the safety rule, and then
#   gives no RPII.

# Design as a plain list, with defaults filled in and checked.
#   grid, model:   a DoseGrid and an outcome model (outcome_beta(), ...);
#   target:        "MTD" or "MED";
#   N, cohort:     maximum number of patients and cohort size;
#   start:         combination index of the first dose;
#   rule:          "cheung" or "ms";
#   eps, taper, weights, selection: as described above; for "ms" the
#     estimate always uses eps = 0.5 and w = 1, and selection is one of
#     ms_selection_weights()'s ("min_S" by default);
#   admis, constraint, safety: M&S options, see ms_decide();
#   delta, within_set: Riviere single pick, see riviere_pick().
pipe_trial_design <- function(grid, model, target = c("MTD", "MED"), N, cohort = 1L,
                              start = 1L, rule = c("cheung", "ms"), eps = 0.5, taper = TRUE,
                              weights = NULL, selection = NULL, admis = "closest",
                              constraint = "neighbouring", safety = 0.8, delta = 0.1,
                              within_set = TRUE) {
  target <- match.arg(target)
  rule <- match.arg(rule)
  if (rule == "ms" && model$family != "binary") {
    stop("The M&S rule needs a binary outcome model.", call. = FALSE)
  }
  if (rule == "ms" && target != "MTD") {
    stop("The M&S rule is for MTD only.", call. = FALSE)
  }
  if (!(start %in% seq_len(n_combos(grid)))) {
    stop("`start` must be a combination index of the grid.", call. = FALSE)
  }
  if (!(N >= 1 && cohort >= 1)) {
    stop("`N` and `cohort` must be at least 1.", call. = FALSE)
  }
  check_eps(eps)
  weights <- weights %||% if (rule == "cheung") "n+1" else "1"
  selection <- selection %||% if (rule == "cheung") "equal" else "min_S"
  if (rule == "cheung") {
    selection <- match.arg(selection, c("equal", "prob"))
  } else {
    selection <- match.arg(selection, c("min_S", "inverse_S", "equal", "inverse_sqrt_S"))
  }
  list(grid = grid, model = model, target = target, N = as.integer(N),
       cohort = as.integer(cohort), start = as.integer(start), rule = rule, eps = eps,
       taper = isTRUE(taper), weights = weights, selection = selection, admis = admis,
       constraint = constraint, safety = safety, delta = delta, within_set = within_set)
}

# One trial under `truth` (see the outcome models). Returns
#   n, s1, s2:      data per combination at the end;
#   doses:          the combination given to each cohort;
#   n_enrolled:     patients enrolled;
#   stopped:        TRUE when the trial stopped before N (M&S safety rule);
#   gamma:          the final PIPE estimate (gamma = 1: exceeds theta);
#   final_set:      maximal safe doses (MTD) or minimal effective doses (MED)
#                   for "cheung", the RPII for "ms" (empty when stopped);
#   pick:           single dose by Riviere's rule (binary MTD only, else NA).
pipe_trial <- function(design, truth) {
  grid <- design$grid
  model <- design$model
  K <- n_combos(grid)
  N <- design$N
  n <- integer(K)
  s1 <- s2 <- numeric(K)
  doses <- integer()
  dose <- design$start
  stopped <- FALSE

  repeat {
    size <- min(design$cohort, N - sum(n))
    y <- model$draw(dose, size, truth)
    n[dose] <- n[dose] + size
    s1[dose] <- s1[dose] + sum(y)
    s2[dose] <- s2[dose] + sum(y^2)
    doses <- c(doses, dose)
    if (sum(n) >= N) break
    p <- model$exceedance(n, s1, s2)
    if (design$rule == "cheung") {
      eps_n <- if (design$taper) taper_eps(sum(n), N, design$eps) else design$eps
      gamma <- trial_estimate(design, p, n, eps_n)
      cand <- cheung_candidates(grid, gamma)
      fav <- if (design$target == "MED") p[cand] else 1 - p[cand]
      prob <- if (design$selection == "equal") rep(1, length(cand)) else fav
      dose <- draw_dose(cand, prob)
    } else {
      dec <- ms_decide(grid, p, n, S = model$a + model$b + n, history = doses,
                       admis = design$admis, constraint = design$constraint,
                       safety = design$safety, selection = design$selection,
                       start = design$start)
      if (dec$stopped) {
        stopped <- TRUE
        break
      }
      dose <- draw_dose(seq_len(K), dec$prob)
    }
  }

  p <- model$exceedance(n, s1, s2)
  if (design$rule == "cheung") {
    gamma <- trial_estimate(design, p, n, design$eps)
    final_set <- pipe_final_set(grid, gamma, design$target)
  } else {
    dec <- ms_decide(grid, p, n, S = model$a + model$b + n, history = doses,
                     admis = design$admis, constraint = design$constraint,
                     safety = design$safety, start = design$start, final = TRUE)
    gamma <- dec$mode
    final_set <- if (stopped) integer() else dec$rpII
  }
  pick <- NA_integer_
  if (model$family == "binary" && design$target == "MTD" && length(final_set)) {
    pick <- riviere_pick(final_set, model$theta, model$a, model$b, n, s1,
                         delta = design$delta, within_set = design$within_set)
  }
  list(n = n, s1 = s1, s2 = s2, doses = doses, n_enrolled = sum(n), stopped = stopped,
       gamma = gamma, final_set = final_set, pick = pick)
}

trial_estimate <- function(design, p, n, eps) {
  w <- pipe_weights(n, design$weights)
  pipe_estimate(design$grid, pipe_gain(p, eps, w, design$target), design$target)
}

# One draw from `x` with weights `prob` (safe when only one entry is left;
# equal weights when every weight is 0).
draw_dose <- function(x, prob) {
  if (!any(prob > 0)) prob <- rep(1, length(x))
  keep <- prob > 0
  x <- x[keep]
  if (length(x) == 1L) x else sample(x, 1L, prob = prob[keep])
}

# Many trials. Returns the per-trial results (`trials`) and, over trials x
# combinations, the patients (`n`), the events or outcome sums (`s1`), the
# final estimate (`gamma`) and the final-set membership (`final_set`); and a
# data frame with n_enrolled, stopped and pick per trial.
pipe_simulate <- function(design, truth, n_trials, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  K <- n_combos(design$grid)
  trials <- lapply(seq_len(n_trials), function(i) pipe_trial(design, truth))
  by_combo <- function(f) matrix(unlist(lapply(trials, f)), n_trials, K, byrow = TRUE)
  list(
    trials = trials,
    n = by_combo(function(t) t$n),
    s1 = by_combo(function(t) t$s1),
    gamma = by_combo(function(t) t$gamma),
    final_set = by_combo(function(t) seq_len(K) %in% t$final_set),
    summary = data.frame(
      trial = seq_len(n_trials),
      n_enrolled = vapply(trials, function(t) t$n_enrolled, numeric(1)),
      stopped = vapply(trials, function(t) t$stopped, logical(1)),
      pick = vapply(trials, function(t) t$pick, integer(1))
    )
  )
}
