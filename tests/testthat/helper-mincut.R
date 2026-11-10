# The min-cut of solve_mincut() in plain R. This was the package's solver
# before the max flow moved to C (src/mincut.c); the tests compare the two.
ref_mincut <- function(grid, d, largest = TRUE) {
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
