# Cheung et al. (2022), Table 2: stored simulation results (inst/simulations).
# The simulation itself is not rerun here; see inst/simulations/README.md.

sim_dir <- function() system.file("simulations", package = "gPIPE")

load_sim <- function() {
  skip_on_cran()
  dir <- sim_dir()
  skip_if(!nzchar(dir) || !file.exists(file.path(dir, "cheung-table2.R")), "simulation scripts not installed")
  old <- options(gpipe.sim.source_only = TRUE)
  on.exit(options(old))
  env <- new.env(parent = globalenv())
  env$here <- dir
  sys.source(file.path(dir, "cheung-table2.R"), envir = env)
  env$dir <- dir
  env
}

test_that("the ensemble is the 252 up-sets of the 5 x 5 grid", {
  e <- load_sim()
  sc <- e$cheung_upsets()
  expect_length(sc, 252L)
  expect_equal(nrow(unique(do.call(rbind, lapply(sc, c)))), 252L)
  g <- e$cheung_grid()
  expect_true(all(vapply(sc, function(m) is_upset(g, grid_from_array(g, m)), logical(1))))
})

test_that("the script's FDR and TPR agree with summary()", {
  e <- load_sim()
  sc <- e$cheung_upsets()
  des <- S7::set_props(e$cheung_designs()[["1-EQ"]], eps = 0.05)
  for (i in c(1L, 100L, 252L)) {
    sim <- simulate(des, nsim = 15, seed = i, truth = e$truth_for(sc[[i]]))
    s <- grid_from_array(des@grid, sc[[i]])
    m <- e$cheung_metrics(s, sim@n, sim@gamma, sim@final_set)
    met <- summary(sim)@metrics
    expect_equal(mean(m$fdr), met$estimate[met$metric == "fdr"])
    expect_equal(mean(m$tpr), met$estimate[met$metric == "tpr"])
  }
})

test_that("stored results have 400 trials in each of 252 scenarios per design", {
  e <- load_sim()
  for (f in c("cheung_table2_scenarios.csv", "cheung_table2_scenarios_cell.csv")) {
    x <- utils::read.csv(file.path(e$dir, "results", f))
    expect_setequal(x$design, c("non-adaptive", "1-EQ", "1-PP", "2-EQ", "2-PP"))
    expect_true(all(table(x$design) == 252L))
    expect_true(all(x$n_trials == 400L))
  }
})

test_that("the 15-cell tables are the ones the stored results give", {
  e <- load_sim()
  for (suffix in c("", "_cell")) {
    scen <- utils::read.csv(file.path(e$dir, "results", sprintf("cheung_table2_scenarios%s.csv", suffix)))
    cells <- utils::read.csv(file.path(e$dir, "results", sprintf("cheung_table2_cells%s.csv", suffix)))
    now <- e$cheung_table2_report(scen)
    expect_equal(now$gpipe, cells$gpipe, tolerance = 1e-8)
    expect_equal(now$z, cells$z, tolerance = 1e-6)
    expect_equal(now$pass, cells$pass)
  }
})

test_that("gPIPE agrees with Cheung's Table 2 in the ways #31 records", {
  e <- load_sim()
  read <- function(suffix) utils::read.csv(file.path(e$dir, "results", sprintf("cheung_table2_cells%s.csv", suffix)))
  pooled <- read("")
  cell <- read("_cell")
  cell_of <- function(x, metric) setNames(x$gpipe[x$metric == metric], x$design[x$metric == metric])
  tpr_order <- c("non-adaptive", "1-PP", "1-EQ", "2-EQ", "2-PP")
  # Not all 15 cells pass the 3-SE rule (see README); these are the recorded outcomes.
  expect_equal(sum(pooled$pass), 6L)
  expect_true(all(pooled$pass[pooled$metric == "fdr"]))
  expect_equal(sum(cell$pass), 6L)
  # Qualitative agreement: the TPR ranking of the designs is Cheung's, under both posteriors.
  for (x in list(pooled, cell)) {
    expect_equal(names(sort(cell_of(x, "tpr")))[1:5], tpr_order)
  }
  # Fixed design, non-adaptive: n_eff and FDR within 3 SE under both posteriors.
  for (x in list(pooled, cell)) {
    expect_true(all(x$pass[x$design == "non-adaptive" & x$metric %in% c("fdr", "n_eff")]))
  }
  # The per-cell posterior closes most of the TPR gap of the fixed design.
  gap <- function(x) abs(x$diff[x$design == "non-adaptive" & x$metric == "tpr"])
  expect_lt(gap(cell), gap(pooled) / 4)
})
