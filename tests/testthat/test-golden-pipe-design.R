# Golden tests (#17): Mander & Sweeting (2015) decisions against pipe.design
# 0.5.1 snapshots (fixtures/pipe_design_V5.json, fixtures/pipe_design_V6.json).
#
# Each state is one decision on fixed data with the stored prior a, b. No
# state has a near tie (best vs second log-score gap >= 7.6e-6), so the tie
# rule does not matter here and tolerances stay far below that.

expect_ms_state <- function(got, expected, id) {
  expect_equal(got$p, expected$p, tolerance = 1e-10, info = paste(id, "p"))
  expect_identical(as.integer(got$mode), expected$mode, info = paste(id, "mode contour"))
  expect_equal(got$q, expected$q, tolerance = 1e-8, info = paste(id, "q"))
  sets <- c("safe", "neighbourhood", "admissible")
  if (!expected$stopped) sets <- c(sets, "dominant", "candidates", "min_S")
  if (!is.null(expected$rp2d)) sets <- c(sets, "rp2d")
  for (field in sets) {
    expect_equal(sort(as.numeric(got[[field]])), expected[[field]], info = paste(id, field))
  }
  expect_identical(got$stopped, expected$stopped, info = paste(id, "stopped"))
  if (!expected$stopped) {
    expect_equal(got$prob, expected$prob, tolerance = 1e-10, info = paste(id, "1/S weights"))
    if (!is.null(expected$next_dose)) {
      expect_true(expected$next_dose %in% got$candidates, info = paste(id, "observed next dose"))
    }
  }
}

for (setting in c("V5", "V6")) {
  desc <- sprintf("decisions match pipe.design snapshots (%s)", setting)
  test_each_impl(desc, "ms_state", function(impl) {
    for (s in read_pipe_design_states(setting)) {
      expect_ms_state(impl$ms_state(s$spec), s$expected, s$id)
    }
  })
}

desc <- "trial stops when no neighbourhood dose is safe, even if d11 is"
test_each_impl(desc, "ms_state", function(impl) {
  for (setting in c("V5", "V6")) {
    s <- Filter(function(x) grepl("edge-safety-neighbourhood-empty", x$id),
                read_pipe_design_states(setting))[[1]]
    got <- impl$ms_state(s$spec)
    expect_true(1 %in% got$safe, info = s$id)
    expect_true(got$stopped, info = s$id)
  }
})

# Final states from 23,000 replayed pipe.design trials (#31) where a dose that
# is safe on the mode contour fails the safety rule. pipe.design finds the
# closest doses among the doses that pass the rule, so the dose below the
# excluded one is recommended.
desc <- "RP2D is closest within the doses that pass the safety rule"
test_each_impl(desc, "ms_state", function(impl) {
  skip_if_not_installed("jsonlite")
  v5 <- jsonlite::read_json(fixture_path("pipe_design_V5.json"))$meta
  states <- jsonlite::read_json(fixture_path("pipe_design_rpii_V5.json"))$states
  expect_length(states, 36)
  for (s in states) {
    spec <- list(levels = c(v5$V5$I, v5$V5$J), theta = v5$V5$theta,
                 a = c(json_matrix(v5$V5$prior_a)), b = c(json_matrix(v5$V5$prior_b)),
                 n = unlist(s$n), r = unlist(s$r), history = s$last_dose,
                 admis = s$admis, final = TRUE, safety = v5$epsilon)
    got <- impl$ms_state(spec)
    expect_false(got$stopped, info = s$id)
    expect_equal(as.numeric(got$rp2d), as.numeric(unlist(s$rpII)), info = s$id)
  }
})
