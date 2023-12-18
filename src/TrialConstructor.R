PatientDataModel <- setRefClass("PatientDataModel",
    fields = list(
        drugCombi = "DrugCombination",   # The DrugDose object
        patientData = "data.frame"  # Data frame to store patient level data
    ),
    methods = list(
        initialize = function(drugCombiObject) {
            drugCombi <<- drugCombiObject
            patientData <<- data.frame(doseCombination = I(list()), outcome = numeric(), stringsAsFactors = FALSE)
        },
        addPatientData = function(doseCombination, outcome) {
            # Validate doseCombination
            if (!isDoseCombinationValid(doseCombination)) {
                stop("Invalid dose combination.")
            }

            # Adds patient data to the model
            patientData <<- rbind(patientData, data.frame(doseCombination = I(list(doseCombination)), outcome = outcome))
        },
        isDoseCombinationValid = function(doseCombination) {
            # Check if each dose is valid for its corresponding drug
            if (!(doseCombination %in% names(drugCombi$doseCombinations))) {
                return(FALSE)
            }
            return(TRUE)
        },
        getSummaryStats = function(includeAllCombi = TRUE) {
            # Computes summary statistics for each dose level combination
            summaryStats <- tapply(patientData$outcome, sapply(patientData$doseCombination, toString), function(x) {
                list(numPatients = length(x), numAdverseEvents = sum(x))
            })
            if (includeAllCombi){
                missingCombinations <- setdiff(names(drugCombi$doseCombinations), names(summaryStats))
                summaryStats[missingCombinations] <- lapply(missingCombinations, function(x) {
                    list(numPatients = 0, numAdverseEvents = 0)
                })
            }
            return(summaryStats)
        },
        generateRandomPatientData = function(maxDoseLevel, numPatients, outcomeProb = 0.5) {
            # maxDoseLevels: A maximum dose labels for the drug combination
            # numPatients: Number of random patients to generate
            # outcomeProb: Probability of adverse event (default 0.5)
            for (i in 1:numPatients) {
                if (.self$isDoseCombinationValid(maxDoseLevel)){
                    validLevels_numeric <- drugCombi$getDoseCombinationsLevel(maxDoseLevel)
                    validLevels <- names(drugCombi$doseCombinations[drugCombi$getDoseCombinationsLevel() <= validLevels_numeric])
                        if (length(validLevels) == 0) {
                            stop("Invalid max dose level for drug: ", drugName)
                        }
                        doseCombinationSelected <- sample(validLevels, 1)                
                }
                # Determine outcome based on specified probability
                outcome <- rbinom(1, 1, outcomeProb)

                # Add generated patient data
                addPatientData(doseCombination = doseCombinationSelected, outcome = outcome)
            }
        }       
    )
)

# Function to create a new PatientDataModel object
createPatientDataModel <- function(drugDoseObject) {
    return(PatientDataModel$new(drugDoseObject))
}