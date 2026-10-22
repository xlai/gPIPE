#' gPIPE: PIPE designs for drug-combination dose finding
#'
#' gPIPE implements the PIPE design (product of independent beta
#' probabilities escalation) for dose-finding trials of drug combinations,
#' for any number of drugs.
#'
#' - [DoseGrid()] sets up the dose combinations and their partial order.
#' - [ms_design()] sets up the original PIPE design of Mander and Sweeting
#'   (2015) for the maximum tolerated dose (MTD) contour, and
#'   [cheung_design()] the PIPE design with epsilon-tapering of Cheung et al.
#'   (2022) for the minimum effective dose (MED). [PipeDesign()] gives full
#'   control over the settings.
#' - [estimate()] and [next_dose()] work on the data of a running trial, and
#'   [simulate()][simulate.PipeDesign] gives operating characteristics.
#'
#' See `vignette("gPIPE")` for a worked example.
#'
#' @references
#' Mander, A. P. and Sweeting, M. J. (2015). A product of independent beta
#' probabilities dose escalation design for dual-agent phase I trials.
#' *Statistics in Medicine*, 34, 1261--1276. \doi{10.1002/sim.6434}
#'
#' Cheung, Y. K., Chandereng, T. and Diaz, K. M. (2022). A novel framework to
#' estimate multidimensional minimum effective doses using asymmetric
#' posterior gain and epsilon-tapering. *The Annals of Applied Statistics*.
#' \doi{10.1214/21-AOAS1603}
#'
#' @keywords internal
"_PACKAGE"
