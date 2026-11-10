/* Max-weight closure by min-cut, for solve_mincut() in R/solvers.R.
 *
 * Network: s -> k with capacity d_k when d_k > 0, k -> t with capacity -d_k
 * when d_k < 0, and lower -> upper with infinite capacity for every cover
 * edge. The maximum flow is found by shortest augmenting paths
 * (Edmonds-Karp). A residual capacity counts only when it exceeds tol.
 *
 * Returns gamma (0/1 over the K combinations): with largest = TRUE the
 * largest optimal up-set (the nodes that cannot reach t in the residual
 * graph), otherwise the smallest one (the nodes reachable from s).
 */
#include <R.h>
#include <Rinternals.h>

/* Breadth-first search over arcs with residual > tol. Forward (from `root`
 * along arcs) when backward = 0; against the arcs when backward = 1. Arc a
 * and its reverse are a and a ^ 1, so the arcs entering v are the reverses
 * of the arcs leaving v. `via[v]` is the arc used to reach v. */
static void bfs(int root, int backward, int n_nodes, const int *first, const int *arcs,
                const int *head, const double *res, double tol, int *seen, int *via,
                int *queue) {
  for (int v = 0; v < n_nodes; v++) seen[v] = 0;
  int qh = 0, qt = 0;
  seen[root] = 1;
  queue[qt++] = root;
  while (qh < qt) {
    int u = queue[qh++];
    for (int i = first[u]; i < first[u + 1]; i++) {
      int a = arcs[i];
      int v = head[a];
      double r = backward ? res[a ^ 1] : res[a];
      if (!seen[v] && r > tol) {
        seen[v] = 1;
        via[v] = a;
        queue[qt++] = v;
      }
    }
  }
}

SEXP gpipe_mincut(SEXP lower_, SEXP upper_, SEXP d_, SEXP tol_, SEXP largest_) {
  if (TYPEOF(lower_) != INTSXP || TYPEOF(upper_) != INTSXP || TYPEOF(d_) != REALSXP) {
    error("`lower` and `upper` must be integer and `d` numeric.");
  }
  const int K = LENGTH(d_);
  const int E = LENGTH(lower_);
  const int *lower = INTEGER(lower_);
  const int *upper = INTEGER(upper_);
  const double *d = REAL(d_);
  const double tol = asReal(tol_);
  const int largest = asLogical(largest_);
  const int s = K, t = K + 1, n_nodes = K + 2;

  /* Check the input here: a bad index or gain would write outside the arrays. */
  if (LENGTH(upper_) != E) error("`lower` and `upper` must have the same length.");
  for (int e = 0; e < E; e++) {
    if (lower[e] < 1 || lower[e] > K || upper[e] < 1 || upper[e] > K) {
      error("The cover edges must be combination indices between 1 and length(d).");
    }
  }
  int n_term = 0;
  for (int k = 0; k < K; k++) {
    if (!R_FINITE(d[k])) error("`d` must be finite.");
    if (d[k] != 0) n_term++;
  }
  const int n_arcs = 2 * (n_term + E);

  int *head = (int *) R_alloc(n_arcs, sizeof(int));
  int *tail = (int *) R_alloc(n_arcs, sizeof(int));
  double *res = (double *) R_alloc(n_arcs, sizeof(double));
  int a = 0;
  for (int k = 0; k < K; k++) {
    if (d[k] > 0) {
      tail[a] = s; head[a] = k; res[a] = d[k];
    } else if (d[k] < 0) {
      tail[a] = k; head[a] = t; res[a] = -d[k];
    } else {
      continue;
    }
    tail[a + 1] = head[a]; head[a + 1] = tail[a]; res[a + 1] = 0;
    a += 2;
  }
  for (int e = 0; e < E; e++) {
    tail[a] = lower[e] - 1; head[a] = upper[e] - 1; res[a] = R_PosInf;
    tail[a + 1] = head[a]; head[a + 1] = tail[a]; res[a + 1] = 0;
    a += 2;
  }

  /* arcs leaving each node: arcs[first[v]] ... arcs[first[v + 1] - 1] */
  int *first = (int *) R_alloc(n_nodes + 1, sizeof(int));
  int *fill = (int *) R_alloc(n_nodes, sizeof(int));
  int *arcs = (int *) R_alloc(n_arcs, sizeof(int));
  for (int v = 0; v <= n_nodes; v++) first[v] = 0;
  for (int i = 0; i < n_arcs; i++) first[tail[i] + 1]++;
  for (int v = 0; v < n_nodes; v++) {
    first[v + 1] += first[v];
    fill[v] = first[v];
  }
  for (int i = 0; i < n_arcs; i++) arcs[fill[tail[i]]++] = i;

  int *seen = (int *) R_alloc(n_nodes, sizeof(int));
  int *via = (int *) R_alloc(n_nodes, sizeof(int));
  int *queue = (int *) R_alloc(n_nodes, sizeof(int));

  for (;;) {
    bfs(s, 0, n_nodes, first, arcs, head, res, tol, seen, via, queue);
    if (!seen[t]) break;
    double flow = R_PosInf;
    for (int v = t; v != s; v = tail[via[v]]) {
      if (res[via[v]] < flow) flow = res[via[v]];
    }
    for (int v = t; v != s; v = tail[via[v]]) {
      res[via[v]] -= flow;
      res[via[v] ^ 1] += flow;
    }
    R_CheckUserInterrupt();
  }

  SEXP out = PROTECT(allocVector(INTSXP, K));
  int *gamma = INTEGER(out);
  if (largest) {
    bfs(t, 1, n_nodes, first, arcs, head, res, tol, seen, via, queue);
    for (int k = 0; k < K; k++) gamma[k] = !seen[k];
  } else {
    for (int k = 0; k < K; k++) gamma[k] = seen[k];
  }
  UNPROTECT(1);
  return out;
}
