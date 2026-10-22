# Fixture readers for the golden tests (#17). See fixtures/README.md.

fixture_path <- function(...) test_path("fixtures", ...)

json_matrix <- function(x) do.call(rbind, lapply(x, unlist))

# 0-based [doseA, doseB] pairs to combination indices (array order).
json_doses <- function(x, I) {
  vapply(x, function(d) d[[1]] + 1 + d[[2]] * I, numeric(1))
}

# pipe.design snapshots as a list of states, each with `spec` (input for
# ms_state) and `expected` (outputs in combination-index form).
read_pipe_design_states <- function(setting) {
  f <- jsonlite::read_json(fixture_path(sprintf("pipe_design_%s.json", setting)))
  st <- f$meta[[setting]]
  I <- st$I
  lapply(f$states, function(x) {
    y <- x$outputs
    doses <- function(cells) sort(json_doses(cells, I))
    expected <- list(
      p = 1 - c(json_matrix(y$p)),
      mode = as.integer(c(json_matrix(y$mode_contour))),
      q = c(json_matrix(y$q)),
      safe = doses(y$safety_admissible),
      neighbourhood = doses(y$neighbourhood),
      admissible = doses(y$admissible),
      stopped = y$stopped
    )
    if (!y$stopped) {
      expected$dominant <- doses(y$dominant)
      expected$candidates <- doses(y$candidates)
      expected$min_S <- doses(y$min_S_ties)
      expected$prob <- c(json_matrix(y$prob_1_over_S))
      if (!is.null(y$next_dose_observed)) {
        expected$next_dose <- json_doses(list(y$next_dose_observed), I)
      }
    }
    if (x$final) expected$rpII <- doses(y$rpII)
    list(
      id = x$id,
      spec = list(
        levels = c(st$I, st$J), theta = st$theta,
        a = c(json_matrix(st$prior_a)), b = c(json_matrix(st$prior_b)),
        n = c(json_matrix(x$inputs$n)), r = c(json_matrix(x$inputs$r)),
        history = json_doses(x$inputs$history, I),
        admis = x$admis, final = x$final, safety = f$meta$epsilon
      ),
      expected = expected
    )
  })
}

read_cheung_table1 <- function() {
  utils::read.csv(fixture_path("cheung2022_table1.csv"), stringsAsFactors = FALSE)
}

# "A4.B5;A5.B4" to combination indices on a 2D grid with I levels of drug A.
parse_dose_labels <- function(x, I) {
  if (is.na(x) || x == "") return(integer())
  parts <- strsplit(strsplit(x, ";", fixed = TRUE)[[1]], ".", fixed = TRUE)
  sort(vapply(parts, function(l) {
    as.integer(sub("A", "", l[1])) + as.integer(I) * (as.integer(sub("B", "", l[2])) - 1L)
  }, integer(1)))
}
