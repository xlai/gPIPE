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
            startingDoseLevel <<- startingDoseLevel %||% names(drugCombi$getDoseCombinationsLevel())[1]
            cohortSize <<- cohortSize
            maxCohorts <<- maxCohorts
            maxSampleSize <<- cohortSize * maxCohorts
            # Validate starting dose level
            if (!isDoseCombinationValid(startingDoseLevel)) {
                stop("Invalid starting dose level.")
            }
        },
        addPatientData = function(doseCombination, outcome, cohort) {
            # Validate doseCombination
            if (!isDoseCombinationValid(doseCombination)) {
                stop("Invalid dose combination.")
            }

            # Adds patient data to the model
            patientData <<- rbind(patientData, data.frame(doseCombination = I(list(doseCombination)), outcome = outcome, cohort))
        },
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
        isDoseCombinationValid = function(doseCombination) {
            # Check if each dose is valid for its corresponding drug
            if (!(doseCombination %in% names(drugCombi$doseCombinations))) {
                return(FALSE)
            }
            return(TRUE)
        },
        setCurrentDoseLevel = function(doseLevel) {
            # First, validate the dose level
            if (!isDoseCombinationValid(doseLevel)) {
                stop("Invalid dose combination.")
            }
            
            # Update the current dose level
            currentDoseLevel <<- doseLevel
        },        
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
        isStoppingCriteriaMet = function() {
            # For now, just check if we've exceeded maximum cohorts
            # Returns TRUE if trial should stop
            return(currentCohort > maxCohorts)
        },
        # Can also add a method to get stopping reason if needed
        getStoppingReason = function() {
            if (isStoppingCriteriaMet()) {
                return("Maximum number of cohorts reached")
            }
            return(NA)  # No stopping reason if criteria not met
        },
        getMTD = function(doseConfig, drugCombiModel){
            n_dose_level <- drugCombiModel$getNumberOfDoseLevels()
            # Calculate the neighbour sum using the provided function
            neighbour_sum <- calculateNeighbourSum(doseConfig, n_dose_level)
            mtd_indices <- which(neighbour_sum == 2 & doseConfig == 1)
            if (length(mtd_indices) == 0){
#                cat('No MTD found.\n')
                return(NA)
            }
            else{
                return(mtd_indices)
            }            
        },
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

# Function to create a new PatientDataModel object
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