# Golden tests (#17): the space of monotone configurations.

# Known counts: binomial(I + J, I) in 2D; 980 monotone 3 x 3 x 3 arrays.
enumeration_cases <- list(
  list(levels = c(6, 6), count = 924),   # pipe.design V5
  list(levels = c(4, 4), count = 70),    # pipe.design V6
  list(levels = c(5, 5), count = 252),   # Cheung et al. (2022)
  list(levels = c(3, 5), count = 56),
  list(levels = c(5, 3), count = 56),
  list(levels = c(1, 4), count = 5),
  list(levels = c(3, 3, 3), count = 980)
)

is_upset_row <- function(levels, gamma) {
  e <- ref_covers(levels)
  all(gamma[e[, "lower"]] <= gamma[e[, "upper"]])
}

test_each_impl("enumeration gives every up-set exactly once", "enumerate", function(impl) {
  for (case in enumeration_cases) {
    configs <- impl$enumerate(case$levels)
    if (is.null(configs)) next
    label <- paste(case$levels, collapse = "x")
    expect_equal(nrow(configs), case$count, info = label)
    expect_equal(ncol(configs), prod(case$levels), info = label)
    expect_true(all(apply(configs, 1, is_upset_row, levels = case$levels)), info = label)
    expect_false(anyDuplicated(apply(configs, 1, paste, collapse = "")) > 0, info = label)
  }
})

test_each_impl("pipe.design mode contours are in the enumeration", "enumerate", function(impl) {
  for (setting in c("V5", "V6")) {
    states <- read_pipe_design_states(setting)
    configs <- impl$enumerate(states[[1]]$spec$levels)
    keys <- apply(configs, 1, paste, collapse = "")
    for (s in states) {
      expect_true(paste(s$expected$mode, collapse = "") %in% keys, info = s$id)
    }
  }
})
