# Implementations under golden test (#17).
#
# Every golden test runs once per implementation: "reference" (the brute-force
# rules in helper-reference.R) and "gPIPE" (the package). An implementation is
# a list of plain functions over combination vectors (array order, drug 1
# fastest). A test skips an implementation that lacks the function it needs,
# so each refactor step (#18) makes its tests run by filling in the matching
# entry of gpipe_impl():
#
# - enumerate: takes levels; returns a 0/1 matrix of all up-sets, one row
#   each, or NULL when the grid is not supported.
# - estimate: takes levels, p, n, eps, weights ("1" or "n+1") and target
#   ("MTD" or "MED"); returns the 0/1 vector gamma. p = P(pi > theta | data).
# - cheung_candidates: takes levels and gamma; returns the combination
#   indices of S1 u S0.
# - ms_state: takes a spec; returns one M&S decision. See ref_ms_state() for
#   the spec and the returned fields.

reference_impl <- function() {
  list(
    enumerate = ref_upsets,
    estimate = ref_estimate,
    cheung_candidates = ref_cheung_candidates,
    ms_state = ref_ms_state
  )
}

gpipe_impl <- function() {
  list(
    enumerate = function(levels) grid_upsets(DoseGrid(levels = levels)),
    estimate = function(levels, p, n, eps, weights = c("1", "n+1"), target = c("MTD", "MED")) {
      target <- match.arg(target)
      w <- pipe_weights(n, match.arg(weights))
      pipe_estimate(DoseGrid(levels = levels), pipe_gain(p, eps, w, target), target)
    }
  )
}

golden_impls <- function() {
  list(reference = reference_impl(), gPIPE = gpipe_impl())
}

# Run `code(impl)` as one test per implementation, skipping implementations
# without `fn`.
test_each_impl <- function(desc, fn, code) {
  impls <- golden_impls()
  for (name in names(impls)) {
    impl <- impls[[name]]
    test_that(sprintf("[%s] %s", name, desc), {
      if (is.null(impl[[fn]])) skip(sprintf("no `%s` in the %s adapter yet (#18)", fn, name))
      code(impl)
    })
  }
}
