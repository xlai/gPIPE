PatientDataModel <- setRefClass("PatientDataModel",
    fields = list(
        drugDose = "DrugDose",   # The DrugDose object
        patientData = "data.frame"  # Data frame to store patient level data
    ),
    methods = list(
        initialize = function(drugDoseObject) {
            drugDose <<- drugDoseObject
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
            for (drug in names(doseCombination)) {
                if (!(drug %in% names(drugDose$drugs)) || 
                    !(doseCombination[[drug]] %in% names(drugDose$drugs[[drug]]$DoseLevels))) {
                    return(FALSE)
                }
            }
            return(TRUE)
        },
        getSummaryStats = function() {
            # Computes summary statistics for each dose level combination
            summaryStats <- tapply(patientData$outcome, sapply(patientData$doseCombination, toString), function(x) {
                list(numPatients = length(x), numAdverseEvents = sum(x))
            })
            return(summaryStats)
        },
        generateRandomPatientData = function(maxDoseLevels, numPatients, outcomeProb = 0.5) {
            # maxDoseLevels: A list of maximum dose labels for each drug
            # numPatients: Number of random patients to generate
            # outcomeProb: Probability of adverse event (default 0.5)

            for (i in 1:numPatients) {
                doseCombination <- lapply(names(maxDoseLevels), function(drugName) {
                    maxLevel <- maxDoseLevels[[drugName]]
                    validLevels <- names(drugDose$drugs[[drugName]]$DoseLevels)
                    validLevels_numeric <- drugDose$drugs[[drugName]]$DoseLevels
                    validLevels <- validLevels[validLevels_numeric <= drugDose$drugs[[drugName]]$DoseLevels[[maxLevel]]]
                    if (length(validLevels) == 0) {
                        stop("Invalid max dose level for drug: ", drugName)
                    }
                    sample(validLevels, 1)
                })

                # Determine outcome based on specified probability
                outcome <- rbinom(1, 1, outcomeProb)

                # Add generated patient data
                addPatientData(doseCombination = doseCombination, outcome = outcome)
            }
        }       
    )
)

# Function to create a new PatientDataModel object
createPatientDataModel <- function(drugDoseObject) {
    return(PatientDataModel$new(drugDoseObject))
}