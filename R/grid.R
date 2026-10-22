#' Dose grid for combination dose-finding
#'
#' A `DoseGrid` describes the dose combinations of a trial with any number of
#' drugs. It owns every conversion between the array view of the grid
#' (one dimension per drug) and the vector view used by the estimators (one
#' entry per combination), and it stores the partial order of the
#' combinations.
#'
#' Combinations are numbered in R's array order: the first drug's level
#' changes fastest. For two drugs this is `c(m)` of an `I x J` matrix whose
#' row `i` is drug 1 at level `i` and column `j` is drug 2 at level `j`.
#'
#' Combination `x` lies below combination `y` when every drug's level in `x`
#' is at most its level in `y`. `@edges` holds the cover relation: the pairs
#' `(lower, upper)` with `upper` directly above `lower`. On a full grid these
#' are the unit steps in one drug.
#'
#' @param levels Integer vector, the number of dose levels of each drug. Names,
#'   if given, are used as drug names.
#' @param drugs Character vector of drug names. Defaults to the names of
#'   `levels`, or `"A"`, `"B"`, ... .
#' @param keep Optional logical vector or array over the full grid (array
#'   order) marking the combinations included in the trial. Defaults to all.
#'
#' @return A `DoseGrid` object with properties
#'   * `levels`: number of levels per drug (named by drug).
#'   * `combos`: integer matrix, one row per included combination, one column
#'     per drug, holding 1-based levels.
#'   * `lookup`: integer array over the full grid giving each combination's
#'     index, or `NA` when it is not included.
#'   * `edges`: integer matrix with columns `lower` and `upper`, the cover
#'     relation of the partial order.
#'
#' @examples
#' g <- DoseGrid(levels = c(A = 3, B = 2))
#' g@combos
#' grid_as_array(g, seq_len(n_combos(g)))
#'
#' # a grid with one combination left out
#' keep <- matrix(TRUE, 3, 2)
#' keep[3, 2] <- FALSE
#' DoseGrid(levels = c(3, 2), keep = keep)@edges
#' @name DoseGrid
#' @export
DoseGrid := new_class(
  package = "gPIPE",
  properties = list(
    levels = class_integer,
    combos = class_integer,
    lookup = class_integer,
    edges = class_integer
  ),
  constructor = function(levels = integer(), drugs = NULL, keep = NULL) {
    if (!is.numeric(levels) || any(is.na(levels)) || any(levels != round(levels))) {
      stop("`levels` must be whole numbers.", call. = FALSE)
    }
    levels <- as.integer(levels)
    if (is.null(drugs)) {
      drugs <- names(levels)
    }
    if (is.null(drugs) && length(levels) > 0) {
      drugs <- LETTERS[seq_along(levels)]
    }
    names(levels) <- drugs

    full <- as.matrix(expand.grid(lapply(levels, seq_len), KEEP.OUT.ATTRS = FALSE))
    storage.mode(full) <- "integer"
    if (is.null(keep)) {
      keep <- rep(TRUE, nrow(full))
    }
    keep <- as.logical(keep)
    if (length(keep) != nrow(full)) {
      stop("`keep` must have one entry per combination of the full grid.", call. = FALSE)
    }

    combos <- full[keep, , drop = FALSE]
    dimnames(combos) <- list(NULL, drugs)
    lookup <- array(NA_integer_, dim = if (length(levels)) unname(levels) else 0L)
    lookup[keep] <- seq_len(sum(keep))

    new_object(
      S7_object(),
      levels = levels,
      combos = combos,
      lookup = lookup,
      edges = grid_cover_edges(combos, lookup, all(keep))
    )
  },
  validator = function(self) {
    if (length(self@levels) == 0) {
      return(NULL)
    }
    if (any(is.na(self@levels)) || any(self@levels < 1L)) {
      "@levels must be positive integers"
    } else if (nrow(self@combos) == 0) {
      "the grid must include at least one combination"
    } else if (ncol(self@combos) != length(self@levels)) {
      "@combos must have one column per drug"
    } else if (!identical(as.integer(dim(self@lookup)), unname(self@levels))) {
      "@lookup must have the dimensions of the full grid"
    } else if (!identical(sort(self@lookup[!is.na(self@lookup)]), seq_len(nrow(self@combos)))) {
      "@lookup must index each row of @combos exactly once"
    } else if (ncol(self@edges) != 2L || any(self@edges < 1L | self@edges > nrow(self@combos))) {
      "@edges must be pairs of combination indices"
    }
  }
)

# Cover relation of the product order. On a full grid the covers are the unit
# steps; on a subset they are found from the comparability matrix.
grid_cover_edges <- function(combos, lookup, full_grid) {
  n <- nrow(combos)
  d <- ncol(combos)
  if (n == 0) {
    return(matrix(integer(), ncol = 2, dimnames = list(NULL, c("lower", "upper"))))
  }
  if (full_grid) {
    edges <- lapply(seq_len(d), function(dim) {
      lower <- which(combos[, dim] < dim(lookup)[dim])
      up <- combos[lower, , drop = FALSE]
      up[, dim] <- up[, dim] + 1L
      cbind(lower, lookup[up])
    })
    edges <- do.call(rbind, edges)
  } else {
    below <- matrix(TRUE, n, n)
    for (dim in seq_len(d)) {
      below <- below & outer(combos[, dim], combos[, dim], `<=`)
    }
    diag(below) <- FALSE
    # x covered by y: x < y and no z with x < z < y
    between <- (below + 0) %*% (below + 0) > 0
    cover <- below & !between
    edges <- which(cover, arr.ind = TRUE)
  }
  storage.mode(edges) <- "integer"
  edges <- edges[order(edges[, 1], edges[, 2]), , drop = FALSE]
  dimnames(edges) <- list(NULL, c("lower", "upper"))
  edges
}

#' Grid helpers
#'
#' Small helpers to move between the vector view and the array view of a
#' [DoseGrid] and to look up combinations.
#'
#' @param grid A [DoseGrid].
#' @param x A vector with one entry per combination.
#' @param a An array over the full grid.
#' @param fill Value used for combinations that are not in the grid.
#' @param levels An integer vector of levels (one per drug), or a matrix with
#'   one row per combination.
#'
#' @return `n_combos()` the number of combinations; `grid_as_array()` an array
#'   with one dimension per drug; `grid_from_array()` a vector over the
#'   combinations; `combo_index()` integer indices (`NA` when absent);
#'   `combo_labels()` labels such as `"A1.B2"`.
#'
#' @examples
#' g <- DoseGrid(levels = c(A = 3, B = 2))
#' n_combos(g)
#' a <- grid_as_array(g, seq_len(n_combos(g)))
#' grid_from_array(g, a)
#' combo_index(g, c(2, 1))
#' combo_labels(g)
#' @name grid-helpers
NULL

#' @rdname grid-helpers
#' @export
n_combos <- function(grid) {
  nrow(grid@combos)
}

#' @rdname grid-helpers
#' @export
grid_as_array <- function(grid, x, fill = NA) {
  if (length(x) != n_combos(grid)) {
    stop("`x` must have one entry per combination.", call. = FALSE)
  }
  out <- array(fill, dim = unname(grid@levels))
  out[!is.na(grid@lookup)] <- x[grid@lookup[!is.na(grid@lookup)]]
  out
}

#' @rdname grid-helpers
#' @export
grid_from_array <- function(grid, a) {
  if (!identical(as.integer(dim(as.array(a))), unname(grid@levels))) {
    stop("`a` must have the dimensions of the full grid.", call. = FALSE)
  }
  out <- a[!is.na(grid@lookup)]
  out[grid@lookup[!is.na(grid@lookup)]] <- out
  out
}

#' @rdname grid-helpers
#' @export
combo_index <- function(grid, levels) {
  levels <- matrix(as.integer(levels), ncol = length(grid@levels))
  inside <- apply(levels, 1, function(l) all(l >= 1L & l <= grid@levels))
  out <- rep(NA_integer_, nrow(levels))
  out[inside] <- grid@lookup[levels[inside, , drop = FALSE]]
  out
}

#' @rdname grid-helpers
#' @export
combo_labels <- function(grid) {
  parts <- lapply(seq_along(grid@levels), function(d) {
    paste0(names(grid@levels)[d], grid@combos[, d])
  })
  do.call(paste, c(parts, sep = "."))
}

# Indices of the combinations covering k from above / below.
upper_covers <- function(grid, k) {
  grid@edges[grid@edges[, "lower"] == k, "upper"]
}

lower_covers <- function(grid, k) {
  grid@edges[grid@edges[, "upper"] == k, "lower"]
}

# Combinations whose levels differ from combination k by at most `radius` in
# every drug (k included).
grid_neighbourhood <- function(grid, k, radius = 1L) {
  delta <- abs(sweep(grid@combos, 2, grid@combos[k, ]))
  which(apply(delta <= radius, 1, all))
}

# TRUE when gamma (0/1 over the combinations) is an up-set: gamma never drops
# when moving up the partial order.
is_upset <- function(grid, gamma) {
  e <- grid@edges
  all(gamma[e[, "lower"]] <= gamma[e[, "upper"]])
}

# Minimal elements of the set `inside` (logical over combinations) and maximal
# elements, both within the partial order of the grid.
minimal_elements <- function(grid, inside) {
  e <- grid@edges
  blocked <- unique(e[inside[e[, "lower"]] & inside[e[, "upper"]], "upper"])
  setdiff(which(inside), blocked)
}

maximal_elements <- function(grid, inside) {
  e <- grid@edges
  blocked <- unique(e[inside[e[, "lower"]] & inside[e[, "upper"]], "lower"])
  setdiff(which(inside), blocked)
}
