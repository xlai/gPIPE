% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strategies-selection-posterior_probability.R
\docType{class}
\name{PosteriorProbabilityStrategy-class}
\alias{PosteriorProbabilityStrategy-class}
\alias{PosteriorProbabilityStrategy}
\title{Posterior Probability Selection Strategy}
\arguments{
\item{admissibleDoses}{Vector of admissible doses}

\item{summaryStats}{Summary statistics for each dose combination}

\item{drugCombi}{DrugCombinationModel object containing posterior probabilities}
}
\value{
The selected dose combination
}
\description{
Select the next dose based on posterior probabilities
}
\details{
A reference class that implements a selection strategy for choosing the next dose
based on posterior probabilities. This strategy selects doses with probabilities
proportional to their posterior probabilities of being the target dose.

The PosteriorProbabilityStrategy class implements the selectDose method by using
the posterior probabilities from the drug combination model to weight the selection
of admissible doses. Doses with higher posterior probabilities have a greater chance
of being selected.
}

\examples{
\dontrun{
# Create base admissible rule class first if needed
library(gPIPE)

# Create a PosteriorProbabilityStrategy
posterior_prob_strategy <- PosteriorProbabilityStrategy$new()
}

}
\seealso{
\code{\link{SelectionStrategy}} for the base class
\code{\link{SmallestSampleSizeStrategy}}, \code{\link{InverseSqrtSampleSizeStrategy}},
\code{\link{InverseDistanceStrategy}}, \code{\link{EqualRandomisationStrategy}} for other selection strategies
}
