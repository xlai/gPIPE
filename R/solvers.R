# PIPE solvers (#18, step 3; #31).
#
# The PIPE estimate is the up-set U (gamma = 1 set) that maximises
# sum(d[U]), with d from pipe_gain(). One exact solver finds it on every
# grid: max-weight closure by min-cut. It only reads the cover relation
# (grid@edges), so it does not depend on the number of drugs.
#
# Tie rule (decided 22 Oct): the smallest favourable set, that is the largest
# U for MTD (largest toxic set) and the smallest U for MED. The largest and
# smallest optimal up-sets are unique, being the union and intersection of all
# optimal ones.
#
# The Mander & Sweeting safety probabilities q are a different quantity: an
# average over all up-sets, which a min-cut cannot give. upset_exceedance()
# computes them exactly with a sweep over the combinations, without listing
# the up-sets.

# ---- Min-cut ----------------------------------------------------------------

# Max-weight closure. Network: s -> k with capacity d_k when d_k > 0,
# k -> t with capacity -d_k when d_k < 0, and lower -> upper with infinite
# capacity for every cover edge, so a cut with finite capacity keeps U closed
# upwards. After a max flow (Edmonds-Karp), the nodes reachable from s in the
# residual graph form the smallest optimal U, and the nodes that cannot reach
# t form the largest one. A residual capacity counts only when it exceeds
# tol = 1e-10 (1 + max |d|). The tolerance follows the largest single gain,
# not the total, so that weights n + 1 on a large grid do not turn small real
# gains into ties.
solve_mincut <- function(grid, d, largest = TRUE) {
  grid <- grid_plain(grid)
  K <- grid$K
  s <- K + 1L
  t <- K + 2L
  pos <- which(d > 0)
  neg <- which(d < 0)
  from <- c(rep(s, length(pos)), neg, grid$lower)
  to <- c(pos, rep(t, length(neg)), grid$upper)
  cap <- c(d[pos], -d[neg], rep(Inf, length(grid$lower)))
  tol <- 1e-10 * (1 + max(abs(d)))

  # Edge i and its reverse are stored at 2i - 1 and 2i.
  m <- length(from)
  head <- as.integer(rbind(to, from))
  tail <- as.integer(rbind(from, to))
  res <- as.numeric(rbind(cap, 0))
  out_edges <- split(seq_len(2 * m), factor(tail, levels = seq_len(t)))
  in_edges <- split(seq_len(2 * m), factor(head, levels = seq_len(t)))
  rev_edge <- function(e) e + ifelse(e %% 2 == 1, 1L, -1L)

  # Breadth-first search from s; returns the edge used to reach each node.
  search <- function() {
    via <- integer(t)
    seen <- logical(t)
    seen[s] <- TRUE
    queue <- s
    while (length(queue)) {
      u <- queue[1]
      queue <- queue[-1]
      for (e in out_edges[[u]]) {
        v <- head[e]
        if (!seen[v] && res[e] > tol) {
          seen[v] <- TRUE
          via[v] <- e
          queue <- c(queue, v)
        }
      }
    }
    list(seen = seen, via = via)
  }

  repeat {
    bfs <- search()
    if (!bfs$seen[t]) break
    path <- integer()
    v <- t
    while (v != s) {
      path <- c(path, bfs$via[v])
      v <- tail[bfs$via[v]]
    }
    flow <- min(res[path])
    res[path] <- res[path] - flow
    back <- rev_edge(path)
    res[back] <- res[back] + flow
  }

  if (!largest) {
    return(as.integer(bfs$seen[seq_len(K)]))
  }
  # Nodes that can reach t: search backwards along edges with residual > tol.
  reach <- logical(t)
  reach[t] <- TRUE
  queue <- t
  while (length(queue)) {
    v <- queue[1]
    queue <- queue[-1]
    for (e in in_edges[[v]]) {
      u <- tail[e]
      if (!reach[u] && res[e] > tol) {
        reach[u] <- TRUE
        queue <- c(queue, u)
      }
    }
  }
  as.integer(!reach[seq_len(K)])
}

# ---- Estimate ---------------------------------------------------------------

# PIPE estimate gamma (0/1 over the combinations) from score differences d.
pipe_estimate <- function(grid, d, target = c("MTD", "MED")) {
  target <- arg_one(target, c("MTD", "MED"))
  grid <- grid_plain(grid)
  if (length(d) != grid$K) {
    stop("`d` must have one entry per combination.", call. = FALSE)
  }
  if (anyNA(d)) {
    stop("`d` must not have missing values.", call. = FALSE)
  }
  solve_mincut(grid, d, largest = target == "MTD")
}

# ---- Safety probabilities ---------------------------------------------------

# Sweep plan of a grid. The sweep visits the combinations in list order, which
# puts every lower cover before its upper covers. After combination k it only
# needs the 0/1 values of the "open" combinations: those visited that still
# have an unvisited upper cover. A state is one 0/1 pattern on the open
# combinations that an up-set can show. Step k lists the transitions:
#   parent: state before k;
#   choice: gamma_k (1 always allowed; 0 only when no lower cover is 1);
#   child:  state after k;
#   n:      number of states after k.
# The plan depends on the grid only and is cached. The number of states is
# small for two drugs (16 at most on a 6 x 6 grid, against 924 up-sets) and
# grows with the width of the grid (2,392 on 5 x 5 x 5); above `max_states`
# the sweep stops with an error.
upset_sweep_plan <- function(grid, max_states = 1e5) {
  grid <- grid_plain(grid)
  key <- grid$key %||% grid_key(grid)
  cached <- .solver_cache[[key]]
  if (!is.null(cached) && cached$peak <= max_states) {
    return(cached)
  }
  e <- grid$edges
  K <- grid$K
  if (any(e[, "lower"] >= e[, "upper"])) {
    stop("The sweep needs every lower cover to come before its upper cover in the list.",
         call. = FALSE)
  }
  # last_use[k]: the last upper cover of k in list order (0 when it has none)
  last_use <- integer(K)
  if (nrow(e)) {
    last <- tapply(e[, "upper"], e[, "lower"], max)
    last_use[as.integer(names(last))] <- last
  }
  open <- integer()
  states <- matrix(0L, nrow = 1, ncol = 0)
  steps <- vector("list", K)
  peak <- 1L
  for (k in seq_len(K)) {
    below <- e[e[, "upper"] == k, "lower"]
    forced <- if (length(below)) {
      rowSums(states[, match(below, open), drop = FALSE]) > 0L
    } else {
      rep(FALSE, nrow(states))
    }
    parent <- c(seq_len(nrow(states)), which(!forced))
    choice <- rep(1:0, c(nrow(states), sum(!forced)))
    keep <- last_use[c(open, k)] > k
    new <- cbind(states[parent, , drop = FALSE], choice)[, keep, drop = FALSE]
    open <- c(open, k)[keep]
    id <- if (ncol(new)) do.call(paste0, as.data.frame(new)) else rep("", nrow(new))
    first <- !duplicated(id)
    states <- new[first, , drop = FALSE]
    peak <- max(peak, nrow(states))
    if (peak > max_states) {
      stop("The safety probabilities need more than `max_states` sweep states on this grid.",
           call. = FALSE)
    }
    steps[[k]] <- list(parent = parent, choice = choice, child = match(id, id[first]),
                       n = nrow(states))
  }
  plan <- list(steps = steps, peak = peak)
  .solver_cache[[key]] <- plan
  plan
}

# Cache key of a grid: its levels and the combinations it includes.
grid_key <- function(grid) {
  paste(c(grid$levels, "|", which(!is.na(grid$lookup))), collapse = ",")
}

.solver_cache <- new.env(parent = emptyenv())

# Mander & Sweeting (2015) safety: the posterior probability that each
# combination is toxic, averaged over all up-sets with weights proportional
# to the contour probabilities prod_U p_k prod_notU (1 - p_k). It takes p, not
# the design's d: q never uses the design's eps or weights. p is clipped as
# in pipe_gain().
#
# The weight of an up-set is exp(sum of qlogis(p) over U) up to a constant, so
# q is a ratio of sums over up-sets. A forward pass over the sweep plan gives
# log_sum_exp of the weights of all partial up-sets that end in each state
# (alpha), a backward pass the same for all completions (beta), and
# q_k = sum over transitions with gamma_k = 1 of alpha x weight x beta / total.
upset_exceedance <- function(grid, p, max_states = 1e5) {
  grid <- grid_plain(grid)
  K <- grid$K
  if (length(p) != K) {
    stop("`p` must have one entry per combination.", call. = FALSE)
  }
  if (anyNA(p) || any(p < 0 | p > 1)) {
    stop("`p` must be probabilities, with no missing values.", call. = FALSE)
  }
  steps <- upset_sweep_plan(grid, max_states)$steps
  logit <- stats::qlogis(pmin(pmax(p, 1e-12), 1 - 1e-12))
  alpha <- vector("list", K + 1L)
  alpha[[1L]] <- 0
  for (k in seq_len(K)) {
    s <- steps[[k]]
    alpha[[k + 1L]] <- log_sum_exp_by(alpha[[k]][s$parent] + s$choice * logit[k], s$child)
  }
  total <- alpha[[K + 1L]]
  beta <- 0
  q <- numeric(K)
  for (k in rev(seq_len(K))) {
    s <- steps[[k]]
    after <- beta[s$child] + s$choice * logit[k]
    on <- s$choice == 1L
    q[k] <- sum(exp(alpha[[k]][s$parent[on]] + after[on] - total))
    beta <- log_sum_exp_by(after, s$parent)
  }
  # rounding can leave a sum a few ulp outside [0, 1]
  pmin(pmax(q, 0), 1)
}

# log(sum(exp(x))) within each group; groups are 1, 2, ... and all occur.
# All groups are shifted by the overall maximum first. A group more than 700
# below it would underflow (exp(-745) is 0), and such a group can become the
# largest again later in the sweep when p is close to 0 or 1. Those groups are
# redone with their own maximum. A group of -Inf only stays -Inf.
log_sum_exp_by <- function(x, group) {
  top <- max(x)
  out <- log(drop(rowsum(exp(x - top), group))) + top
  redo <- (out < top - 700)[group] & x > -Inf
  if (any(redo)) {
    g <- group[redo]
    own <- tapply(x[redo], g, max)
    part <- log(drop(rowsum(exp(x[redo] - own[as.character(g)]), g))) + own
    out[as.integer(names(own))] <- part
  }
  out
}
