#' Patient Data Model Class
#'
#' A reference class that manages patient data and trial progression in dose-finding studies.
#' This class works with the \code{\link{DrugCombination-class}} class to track patient outcomes
#' at different dose combinations, apply admissibility rules, and implement dose selection strategies.
#'
#' @field drugCombi An object of class \code{\link{DrugCombination-class}}
#' @field patientData A data frame storing patient-level data
#' @field admissibleRule A list of admissibility rule objects
#' @field selectionStrategy A list of dose selection strategy objects
#' @field currentCohort The current cohort number
#' @field currentDoseLevel The current dose level being tested
#' @field startingDoseLevel The starting dose level for the trial
#' @field cohortSize The number of patients per cohort
#' @field maxCohorts The maximum number of cohorts allowed in the trial
#' @field maxSampleSize The maximum sample size (derived from cohortSize * maxCohorts)
#'
#' @importFrom methods setRefClass
#' @seealso \code{\link{DrugCombination-class}} for the drug combination class that this model works with
#' @export
PatientDataModel <- setRefClass("PatientDataModel",
    fields = list(
        drugCombi = "DrugCombination",   # The DrugDose object
        patientData = "data.frame",  # Data frame to store patient level data        
        admissibleRule = "list",
        selectionStrategy = "list",
        currentCohort = "numeric",
        currentDoseLevel = "character",
        startingDoseLevel = "character",  # Starting dose level
        cohortSize = "numeric",           # Number of patients per cohort
        maxCohorts = "numeric",           # Maximum number of cohorts
        maxSampleSize = "numeric"         # Derived from cohortSize * maxCohorts
    ),
    methods = list(
        #' @description
        #' Initialize a new PatientDataModel object
        #'
        #' @param drugCombiObject An object of class DrugCombination
        #' @param admissibleRuleList A list of admissibility rules
        #' @param selectionStrategyList A list of selection strategies
        #' @param startingDoseLevel The starting dose level (character)
        #' @param cohortSize The number of patients per cohort
        #' @param maxCohorts The maximum number of cohorts
        initialize = function(drugCombiObject, 
                            admissibleRuleList = list(), 
                            selectionStrategyList = list(),
                            startingDoseLevel = NULL,
                            cohortSize = 3,
                            maxCohorts = 20) {
            drugCombi <<- drugCombiObject
            patientData <<- data.frame(
                doseCombination = I(list()),
                outcome = numeric(),
                stringsAsFactors = FALSE
            )
            admissibleRule <<- admissibleRuleList
            selectionStrategy <<- selectionStrategyList
            currentCohort <<- 1
            # Initialize trial design parameters
            startingDoseLevel <<- if (is.null(startingDoseLevel)) {
                names(drugCombi$getDoseCombinationsLevel())[1]
            } else {
                startingDoseLevel
            }
            cohortSize <<- cohortSize
            maxCohorts <<- maxCohorts
            maxSampleSize <<- cohortSize * maxCohorts
            # Validate starting dose level
            if (!isDoseCombinationValid(startingDoseLevel)) {
                stop("Invalid starting dose level.")
            }
        },
        
        #' @description
        #' Add patient data to the model
        #'
        #' @param doseCombination The dose combination administered to the patient
        #' @param outcome The patient outcome (0 = no toxicity, 1 = toxicity)
        #' @param cohort The cohort number
        addPatientData = function(doseCombination, outcome, cohort) {
            # Validate doseCombination
            if (!isDoseCombinationValid(doseCombination)) {
                stop("Invalid dose combination.")
            }

            # Adds patient data to the model
            patientData <<- rbind(patientData, data.frame(doseCombination = I(list(doseCombination)), outcome = outcome, cohort))
        },
        
        #' @description
        #' Reset the trial to its initial state
        resetTrial = function() {
            patientData <<- data.frame(
                doseCombination = I(list()),
                outcome = numeric(),
                cohort = numeric(),
                stringsAsFactors = FALSE
            )
            currentCohort <<- 1
            currentDoseLevel <<- startingDoseLevel
        },
        
        #' @description
        #' Check if a dose combination is valid
        #'
        #' @param doseCombination The dose combination to validate
        #' @return Logical indicating if the dose combination is valid
        isDoseCombinationValid = function(doseCombination) {
            # Check if each dose is valid for its corresponding drug
            if (!(doseCombination %in% names(drugCombi$doseCombinations))) {
                return(FALSE)
            }
            return(TRUE)
        },
        
        #' @description
        #' Set the current dose level
        #'
        #' @param doseLevel The dose level to set as current
        setCurrentDoseLevel = function(doseLevel) {
            # First, validate the dose level
            if (!isDoseCombinationValid(doseLevel)) {
                stop("Invalid dose combination.")
            }
            
            # Update the current dose level
            currentDoseLevel <<- doseLevel
        },
        
        #' @description
        #' Get summary statistics for each dose combination
        #'
        #' @param includeAllCombi Logical indicating whether to include all combinations or only those with patient data
        #' @return List of summary statistics for each dose combination
        getSummaryStats = function(includeAllCombi = TRUE) {
            # Computes summary statistics for each dose level combination
            summaryStats <- tapply(patientData$outcome, sapply(patientData$doseCombination, toString), function(x) {
                list(numPatients = length(x), numAdverseEvents = sum(x))
            })
            if (includeAllCombi) {
                # Initialize all combinations with (numPatients = 0, numAdverseEvents = 0)
                allCombinations <- setNames(lapply(names(drugCombi$doseCombinations), function(x) list(numPatients = 0, numAdverseEvents = 0)), names(drugCombi$doseCombinations))
                # Update the initialized combinations with actual data where available
                summaryStats <- modifyList(allCombinations, summaryStats)
            } 
            return(summaryStats)
        },
        
        #' @description
        #' Generate random patient data for simulation
        #'
        #' @param outcomeProb Probability of toxicity for outcome generation
        #' @param doseLevel Optional dose level (defaults to current dose level)
        #' @param numPatients Optional number of patients (defaults to cohort size)
        generateRandomPatientData = function(outcomeProb, doseLevel = NULL, numPatients = NULL) {
            # If only outcomeProb provided, use object state
            if (is.null(doseLevel) && is.null(numPatients)) {
                doseLevel <- currentDoseLevel
                numPatients <- cohortSize
            }
            # Validate dose level
            if (!isDoseCombinationValid(doseLevel)) {
                stop("Invalid dose combination.")
            }
            # Store current state if using custom dose
            if (doseLevel != currentDoseLevel) {
                previousDose <- currentDoseLevel
                currentDoseLevel <<- doseLevel
            }
            # Generate data
            for (i in 1:numPatients) {
                outcome <- rbinom(1, 1, outcomeProb)
                addPatientData(
                    doseCombination = doseLevel,
                    outcome = outcome,
                    cohort = currentCohort
                )
            }
            # Update cohort
            currentCohort <<- currentCohort + 1
        },
        
        #' @description
        #' Determine the next dose level based on admissibility rules and selection strategy
        #'
        #' @param currentLevel The current dose level
        #' @param valid_dose_config Configuration of valid doses
        #' @param drugCombiModel Drug combination model
        #' @return Character string representing the next dose level
        getNextDoseLevel = function(currentLevel, valid_dose_config, drugCombiModel) {
            admissibleDoses <- lapply(admissibleRule, function(rule) rule$isAdmissible(valid_dose_config, currentLevel, drugCombiModel))

            if (length(admissibleDoses) == 1) {
                # If the length of admissibleDoses is 1, directly assign it to admissibleDoses_final
                admissibleDoses_final <- admissibleDoses[[1]]
            } else {
            #### Manually coded section, needs update later ###  
                admissibleDoses_neighbour_safe <- Reduce(intersect, admissibleDoses[-1])
                if (length(intersect(admissibleDoses[[1]], admissibleDoses_neighbour_safe)) == 0) {
                    admissibleDoses_final <- NULL
                } else {
                    admissibleDoses_final <- intersect(admissibleDoses[[1]], admissibleDoses_neighbour_safe)
                }
            }
            # Apply the selection strategy to choose the next dose
            if (length(admissibleDoses_final) > 0) {
                summaryStats <- .self$getSummaryStats(includeAllCombi = TRUE)
                nextDose <- selectionStrategy[[1]]$selectDose(admissibleDoses_final, summaryStats, drugCombiModel)
            } else {
                nextDose <- NA
                warning("No admissible doses found based on the current rules.")
            }
            return(nextDose)
        },
        
        #' @description
        #' Check if any stopping criteria have been met
        #'
        #' @return Logical indicating if stopping criteria have been met
        isStoppingCriteriaMet = function() {
            # For now, just check if we've exceeded maximum cohorts
            # Returns TRUE if trial should stop
            return(currentCohort > maxCohorts)
        },
        
        #' @description
        #' Get the reason for stopping the trial
        #'
        #' @return Character string with stopping reason or NA if trial should continue
        getStoppingReason = function() {
            if (isStoppingCriteriaMet()) {
                return("Maximum number of cohorts reached")
            }
            return(NA)  # No stopping reason if criteria not met
        },
        
        #' @description
        #' Get the Maximum Tolerated Dose (MTD) indices
        #'
        #' @param doseConfig Dose configuration
        #' @param drugCombiModel Drug combination model
        #' @return Numeric vector of MTD indices or NA if none found
        getMTD = function(doseConfig, drugCombiModel){
            n_dose_level <- drugCombiModel$getNumberOfDoseLevels()
            # Calculate the neighbour sum using the provided function
            neighbour_sum <- calculateNeighbourSum(doseConfig, n_dose_level)
            mtd_indices <- which(neighbour_sum == 2 & doseConfig == 1)
            if (length(mtd_indices) == 0){
                return(NA)
            }
            else{
                return(mtd_indices)
            }            
        },
        
        #' @description
        #' Get the Recommended Phase 2 Dose (RP2D) indices
        #'
        #' @param doseConfig Dose configuration
        #' @param drugCombiModel Drug combination model
        #' @return Numeric vector of RP2D indices or NA if none found
        getRP2D = function(doseConfig, drugCombiModel){
            n_dose_level <- drugCombiModel$getNumberOfDoseLevels()
            # Calculate the neighbour sum using the provided function
            neighbour_sum <- calculateNeighbourSum(doseConfig, n_dose_level)
            rp2d_indices <- which(neighbour_sum == 2 & doseConfig == 1)
            if (length(rp2d_indices) == 0){
                cat('No RP2D found.\n')
                return(NA)
            }
            else{
                return(rp2d_indices)
            }                    
        }
    )
)

#' Create a new PatientDataModel object
#'
#' This function creates a new PatientDataModel object, which manages patient data and trial progression
#' in dose-finding studies.
#'
#' @param drugDoseObject A DrugCombination object
#' @param admissibleRuleList A list of admissibility rules
#' @param selectionStrategyList A list of selection strategies
#' @param startingDoseLevel The starting dose level (character)
#' @param cohortSize The number of patients per cohort
#' @param maxCohorts The maximum number of cohorts
#'
#' @return A PatientDataModel object
#' @export
createPatientDataModel <- function(drugDoseObject, 
                                 admissibleRuleList = list(),
                                 selectionStrategyList = list(),
                                 startingDoseLevel = NULL,
                                 cohortSize = 3,
                                 maxCohorts = 20) {
    
    return(PatientDataModel$new(
        drugCombiObject = drugDoseObject,
        admissibleRuleList = admissibleRuleList,
        selectionStrategyList = selectionStrategyList,
        startingDoseLevel = startingDoseLevel,
        cohortSize = cohortSize,
        maxCohorts = maxCohorts
    ))
}

#' Find indices closest to a boundary
#'
#' This helper function finds indices of 0s that are adjacent to 1s in a matrix.
#'
#' @param vectorized_matrix A vectorized matrix representation
#' @param nrow Number of rows in the matrix
#' @param ncol Number of columns in the matrix
#'
#' @return Vector of indices closest to the boundary
#' @keywords internal
find_closest_to_boundary <- function(vectorized_matrix, nrow, ncol) {
  # Reshape the vectorized matrix into a 2D matrix
  matrix_2d <- matrix(vectorized_matrix, nrow = nrow, ncol = ncol, byrow = TRUE)
  
  # Function to check if an element is close to boundary
  is_close_to_boundary <- function(i, j) {
    neighbours <- c(
      if (i > 1) matrix_2d[i - 1, j] else NA,    # Up
      if (i < nrow) matrix_2d[i + 1, j] else NA, # Down
      if (j > 1) matrix_2d[i, j - 1] else NA,    # Left
      if (j < ncol) matrix_2d[i, j + 1] else NA  # Right
    )
    any(neighbours == 1, na.rm = TRUE)
  }
  
  # Find linear indices of 0s close to the boundary
  closest_indices <- which(matrix_2d == 0 & apply(expand.grid(1:nrow, 1:ncol), 1, function(idx) {
    is_close_to_boundary(idx[1], idx[2])
  }))
  
  return(closest_indices)
}