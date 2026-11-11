// Max-weight closure by min-cut, for solve_mincut() in R/solvers.R.
//
// The same algorithm in plain R is ref_mincut() in
// tests/testthat/helper-mincut.R; the steps below follow it one for one.
//
// Network. Nodes 0 .. K-1 are the combinations, then the source s and the
// sink t.
//   s -> k      with capacity d_k    when d_k > 0   (gain from putting k in U)
//   k -> t      with capacity -d_k   when d_k < 0   (gain from leaving k out)
//   lower -> upper with infinite capacity for every cover edge, so a cut of
//               finite capacity never has `lower` in U without `upper`.
// After a maximum flow, the nodes that s can still reach form the smallest
// optimal up-set U, and the nodes that cannot reach t form the largest one.

#include <Rcpp.h>
#include <algorithm>
#include <vector>

namespace {

struct Arc {
  int from;
  int to;
  double residual;  // capacity still unused
  int reverse;      // index of the arc to -> from that carries flow back
};

struct Network {
  std::vector<Arc> arcs;
  std::vector<std::vector<int>> leaving;  // leaving[v]: arcs with from == v

  explicit Network(int n_nodes) : leaving(n_nodes) {}

  // An arc with the given capacity, and its reverse with capacity 0.
  void add_arc(int from, int to, double capacity) {
    int forward = arcs.size();
    int backward = forward + 1;
    arcs.push_back({from, to, capacity, backward});
    arcs.push_back({to, from, 0.0, forward});
    leaving[from].push_back(forward);
    leaving[to].push_back(backward);
  }
};

// Breadth-first search over arcs whose residual exceeds tol.
//   towards_root = false: nodes that `root` can reach;
//   towards_root = true:  nodes that can reach `root`.
// reached[v] says whether v was found, and via[v] is the arc that found it.
void search(const Network& net, int root, bool towards_root, double tol,
            std::vector<char>& reached, std::vector<int>& via) {
  std::fill(reached.begin(), reached.end(), 0);
  // The nodes found so far, in the order found; `next` is the first one whose
  // arcs have not been looked at yet.
  std::vector<int> to_visit;
  to_visit.reserve(reached.size());
  reached[root] = 1;
  to_visit.push_back(root);
  for (std::size_t next = 0; next < to_visit.size(); next++) {
    int v = to_visit[next];
    for (int a : net.leaving[v]) {
      const Arc& arc = net.arcs[a];
      // Going out from v we use the arc itself. Going towards the root we
      // need the arc that points at v, which is the reverse of this one.
      double residual = towards_root ? net.arcs[arc.reverse].residual : arc.residual;
      if (!reached[arc.to] && residual > tol) {
        reached[arc.to] = 1;
        via[arc.to] = a;
        to_visit.push_back(arc.to);
      }
    }
  }
}

}  // namespace

// gamma (0/1 over the combinations) of the optimal up-set: the largest one
// when `largest` is true, otherwise the smallest. `lower` and `upper` are the
// cover edges as 1-based combination indices; a residual capacity counts only
// when it exceeds `tol`.
// [[Rcpp::export]]
Rcpp::IntegerVector mincut_cpp(Rcpp::IntegerVector lower, Rcpp::IntegerVector upper,
                               Rcpp::NumericVector d, double tol, bool largest) {
  const int K = d.size();
  const int source = K;
  const int sink = K + 1;

  if (lower.size() != upper.size()) {
    Rcpp::stop("`lower` and `upper` must have the same length.");
  }

  // 1. Build the network.
  Network net(K + 2);
  for (int k = 0; k < K; k++) {
    if (!R_finite(d[k])) Rcpp::stop("`d` must be finite.");
    if (d[k] > 0) net.add_arc(source, k, d[k]);
    if (d[k] < 0) net.add_arc(k, sink, -d[k]);
  }
  for (int e = 0; e < lower.size(); e++) {
    if (lower[e] < 1 || lower[e] > K || upper[e] < 1 || upper[e] > K) {
      Rcpp::stop("The cover edges must be combination indices between 1 and length(d).");
    }
    net.add_arc(lower[e] - 1, upper[e] - 1, R_PosInf);
  }

  // 2. Maximum flow (Edmonds-Karp): while a path from source to sink has
  //    unused capacity, push as much flow along it as its tightest arc allows.
  std::vector<char> reached(K + 2);
  std::vector<int> via(K + 2);
  while (true) {
    search(net, source, false, tol, reached, via);
    if (!reached[sink]) break;

    double flow = R_PosInf;
    for (int v = sink; v != source; v = net.arcs[via[v]].from) {
      flow = std::min(flow, net.arcs[via[v]].residual);
    }
    for (int v = sink; v != source; v = net.arcs[via[v]].from) {
      Arc& arc = net.arcs[via[v]];
      arc.residual -= flow;
      net.arcs[arc.reverse].residual += flow;
    }
    Rcpp::checkUserInterrupt();
  }

  // 3. Read the up-set off the residual network.
  Rcpp::IntegerVector gamma(K);
  if (largest) {
    // largest U: every combination that cannot reach the sink
    search(net, sink, true, tol, reached, via);
    for (int k = 0; k < K; k++) gamma[k] = !reached[k];
  } else {
    // smallest U: every combination the source can still reach (`reached`
    // holds the last search from the source)
    for (int k = 0; k < K; k++) gamma[k] = reached[k];
  }
  return gamma;
}
