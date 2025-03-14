% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TrialConstructor.R
\docType{class}
\name{PatientDataModel-class}
\alias{PatientDataModel-class}
\alias{PatientDataModel}
\title{Patient Data Model Class}
\arguments{
\item{drugCombiObject}{An object of class DrugCombination}

\item{admissibleRuleList}{A list of admissibility rules}

\item{selectionStrategyList}{A list of selection strategies}

\item{startingDoseLevel}{The starting dose level (character)}

\item{cohortSize}{The number of patients per cohort}

\item{maxCohorts}{The maximum number of cohorts}

\item{outcome}{The patient outcome (0 = no toxicity, 1 = toxicity)}

\item{cohort}{The cohort number}

\item{doseCombination}{The dose combination to validate}

\item{includeAllCombi}{Logical indicating whether to include all combinations or only those with patient data}

\item{outcomeProb}{Probability of toxicity for outcome generation}

\item{doseLevel}{Optional dose level (defaults to current dose level)}

\item{numPatients}{Optional number of patients (defaults to cohort size)}

\item{currentLevel}{The current dose level}

\item{valid_dose_config}{Configuration of valid doses}

\item{doseConfig}{Dose configuration}

\item{drugCombiModel}{Drug combination model}
}
\value{
Logical indicating if the dose combination is valid

List of summary statistics for each dose combination

Character string representing the next dose level

Logical indicating if stopping criteria have been met

Character string with stopping reason or NA if trial should continue

Numeric vector of MTD indices or NA if none found

Numeric vector of RP2D indices or NA if none found
}
\description{
Initialize a new PatientDataModel object

Add patient data to the model

Reset the trial to its initial state

Check if a dose combination is valid

Set the current dose level

Get summary statistics for each dose combination

Generate random patient data for simulation

Determine the next dose level based on admissibility rules and selection strategy

Check if any stopping criteria have been met

Get the reason for stopping the trial

Get the Maximum Tolerated Dose (MTD) indices

Get the Recommended Phase 2 Dose (RP2D) indices
}
\details{
A reference class that manages patient data and trial progression in dose-finding studies.
This class works with the \code{\link{DrugCombination-class}} class to track patient outcomes
at different dose combinations, apply admissibility rules, and implement dose selection strategies.
}
\section{Fields}{

\describe{
\item{\code{drugCombi}}{An object of class \code{\link{DrugCombination-class}}}

\item{\code{patientData}}{A data frame storing patient-level data}

\item{\code{admissibleRule}}{A list of admissibility rule objects}

\item{\code{selectionStrategy}}{A list of dose selection strategy objects}

\item{\code{currentCohort}}{The current cohort number}

\item{\code{currentDoseLevel}}{The current dose level being tested}

\item{\code{startingDoseLevel}}{The starting dose level for the trial}

\item{\code{cohortSize}}{The number of patients per cohort}

\item{\code{maxCohorts}}{The maximum number of cohorts allowed in the trial}

\item{\code{maxSampleSize}}{The maximum sample size (derived from cohortSize * maxCohorts)}
}}


\seealso{
\code{\link{DrugCombination-class}} for the drug combination class that this model works with
}
