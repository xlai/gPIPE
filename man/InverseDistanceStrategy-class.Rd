% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strategies-selection-inverse_distance.R
\docType{class}
\name{InverseDistanceStrategy-class}
\alias{InverseDistanceStrategy-class}
\alias{InverseDistanceStrategy}
\title{Inverse Distance Selection Strategy}
\arguments{
\item{admissibleDoses}{Vector of admissible doses}

\item{summaryStats}{Summary statistics for each dose combination}

\item{drugCombi}{DrugCombinationModel object containing posterior probabilities and target threshold}

\item{distance_type}{Type of distance metric to use: "absolute" or "squared"}
}
\value{
The selected dose combination
}
\description{
Select the next dose based on inverse distance to target threshold
}
\details{
A reference class that implements a selection strategy for choosing the next dose
based on the inverse of the distance between posterior probabilities and the target
toxicity threshold. This strategy favors doses whose posterior probabilities are
closest to the target.

The InverseDistanceStrategy class implements the selectDose method by calculating
probabilities inversely proportional to the distance between each dose's posterior
probability and the target toxicity threshold (theta). The distance can be calculated
using absolute difference or squared difference.
}

\examples{
\dontrun{
# Create base admissible rule class first if needed
library(gPIPE)

# Create an InverseDistanceStrategy using absolute distance
inverse_dist_strategy <- InverseDistanceStrategy$new()
}

}
\seealso{
\code{\link{SelectionStrategy}} for the base class
\code{\link{SmallestSampleSizeStrategy}}, \code{\link{PosteriorProbabilityStrategy}},
\code{\link{InverseSqrtSampleSizeStrategy}}, \code{\link{EqualRandomisationStrategy}} for other selection strategies
}
