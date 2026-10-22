# PIPE solvers (#18, step 3).
#
# The PIPE estimate is the up-set U (gamma = 1 set) that maximises
# sum(d[U]), with d from pipe_gain(). Two exact solvers:
#
# - enumeration of every up-set as one 0/1 matrix C (configurations x
#   combinations), scored by C %*% d. Used for full grids of one or two drugs,
#   and the only route that gives the M&S safety probabilities q;
# - max-weight closure by min-cut, for any other grid.
#
# Tie rule (decided 22 Oct): the smallest favourable set, that is the largest
# U for MTD (largest toxic set) and the smallest U for MED. The largest and
# smallest optimal up-sets are unique, being the union and intersection of all
# optimal ones.

# ---- Enumeration ------------------------------------------------------------

# All up-sets of the grid's partial order as a 0/1 integer matrix, one row per
# configuration, one column per combination. Combinations are in array order,
# so every lower cover of k has a smaller index and each column only depends on
# columns already filled. Results are cached by grid.
grid_upsets <- function(grid, max_configs = 1e6) {
  key <- paste(c(grid@levels, "|", which(!is.na(grid@lookup))), collapse = ",")
  cached <- .solver_cache[[key]]
  if (!is.null(cached) && nrow(cached) <= max_configs) {
    return(cached)
  }
  e <- grid@edges
  configs <- matrix(0L, nrow = 1, ncol = 0)
  for (k in seq_len(n_combos(grid))) {
    below <- e[e[, "upper"] == k, "lower"]
    forced <- if (length(below)) {
      rowSums(configs[, below, drop = FALSE]) > 0L
    } else {
      rep(FALSE, nrow(configs))
    }
    free <- configs[!forced, , drop = FALSE]
    if (nrow(configs) + nrow(free) > max_configs) {
      stop("The grid has more than `max_configs` up-sets; use the min-cut solver.",
           call. = FALSE)
    }
    configs <- rbind(cbind(configs, 1L), cbind(free, rep(0L, nrow(free))))
  }
  .solver_cache[[key]] <- configs
  configs
}

.solver_cache <- new.env(parent = emptyenv())

# Best up-set among the rows of `configs` under the tie rule. Scores within
# tol = 1e-10 (1 + max |score|) of the best count as ties.
solve_enumeration <- function(configs, d, largest = TRUE) {
  score <- drop(configs %*% d)
  tol <- 1e-10 * (1 + max(abs(score)))
  best <- which(score >= max(score) - tol)
  size <- rowSums(configs[best, , drop = FALSE])
  pick <- best[if (largest) which.max(size) else which.min(size)]
  configs[pick, ]
}

# Mander & Sweeting (2015) safety: the posterior probability that each
# combination is toxic, averaged over all up-sets with weights
# softmax(C %*% d). With eps = 0.5 and w = 1, d = qlogis(p) and the weights
# are the contour probabilities prod_U p_k prod_notU (1 - p_k), normalised.
upset_exceedance <- function(configs, d) {
  score <- drop(configs %*% d)
  wts <- exp(score - max(score))
  drop(crossprod(configs, wts / sum(wts)))
}

# ---- Min-cut ----------------------------------------------------------------

# Max-weight closure. Network: s -> k with capacity d_k when d_k > 0,
# k -> t with capacity -d_k when d_k < 0, and lower -> upper with infinite
# capacity for every cover edge, so a cut with finite capacity keeps U closed
# upwards. After a max flow (Edmonds-Karp), the nodes reachable from s in the
# residual graph form the smallest optimal U, and the nodes that cannot reach
# t form the largest one. A residual capacity counts only when it exceeds
# tol = 1e-10 (1 + sum |d|).
solve_mincut <- function(grid, d, largest = TRUE) {
  K <- n_combos(grid)
  s <- K + 1L
  t <- K + 2L
  pos <- which(d > 0)
  neg <- which(d < 0)
  from <- c(rep(s, length(pos)), neg, grid@edges[, "lower"])
  to <- c(pos, rep(t, length(neg)), grid@edges[, "upper"])
  cap <- c(d[pos], -d[neg], rep(Inf, nrow(grid@edges)))
  tol <- 1e-10 * (1 + sum(abs(d)))

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
# method = "auto" enumerates on full grids of one or two drugs and uses
# min-cut otherwise.
pipe_estimate <- function(grid, d, target = c("MTD", "MED"),
                          method = c("auto", "enumerate", "mincut")) {
  target <- match.arg(target)
  method <- match.arg(method)
  if (length(d) != n_combos(grid)) {
    stop("`d` must have one entry per combination.", call. = FALSE)
  }
  if (method == "auto") {
    full <- !anyNA(grid@lookup)
    method <- if (full && length(grid@levels) <= 2L) "enumerate" else "mincut"
  }
  largest <- target == "MTD"
  if (method == "enumerate") {
    solve_enumeration(grid_upsets(grid), d, largest)
  } else {
    solve_mincut(grid, d, largest)
  }
}
