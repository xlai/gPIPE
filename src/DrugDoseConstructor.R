DrugDose <- setRefClass("DrugDose",
    fields = list(
        drugs = "list"  # A list of drugs with their dose levels and numeric order
    ),
    methods = list(
        initialize = function(dose_info = NULL, num_doses = NULL) {
            if (!is.null(dose_info)) {
                # Initialize from YAML file
                drugs <<- dose_info
            } else if (!is.null(num_doses)) {
                # Initialize from provided number of doses
                drugs <<- lapply(num_doses, function(dose_num) {
                    list(DoseLevels = setNames(1:dose_num, 1:dose_num))
                })
            } else {
                # Handle NULL scenario - set drugs to an empty list or a default value
                drugs <<- list()                
#                stop("Either dose_info (YAML content) or num_doses must be provided.")
            }
        },
        getDoseLabels = function(drug_name) {
            # Returns the dose labels for a given drug
            return(names(drugs[[drug_name]]$DoseLevels))
        },
        getNumericOrder = function(drug_name, dose_label = NULL) {
            # Validate input
            if (!drug_name %in% names(drugs)) {
                stop("Invalid drug name")
            }
            if (is.null(dose_label)) {
                # Return all numeric levels for the drug
                return(drugs[[drug_name]]$DoseLevels)
            } else {
                # Return the numeric level for the specified dose label
                if (!dose_label %in% names(drugs[[drug_name]]$DoseLevels)) {
                    stop("Invalid dose label")
                }
                return(drugs[[drug_name]]$DoseLevels[[dose_label]])
            }
        },
        addDoseLabel = function(drug_name, new_label = NULL, new_order) {
            # Validate input
            if (!drug_name %in% names(drugs)) {
                stop("Invalid drug name")
            }
            if (is.null(new_label)) {
                # Generate default label if not provided
                existing_labels <- names(drugs[[drug_name]]$DoseLevels)
                new_label_count <- sum(grepl("newlabel", existing_labels))
                new_label <- paste("newlabel", new_label_count + 1, sep="")
            }

            # Adjust existing labels' numeric order
            for (label in names(drugs[[drug_name]]$DoseLevels)) {
                if (drugs[[drug_name]]$DoseLevels[[label]] >= new_order) {
                    drugs[[drug_name]]$DoseLevels[[label]] <<- drugs[[drug_name]]$DoseLevels[[label]] + 1
                }
            }

            # Add the new label with its order
            drugs[[drug_name]]$DoseLevels[[new_label]] <<- new_order
        }        
        # ... (other methods)
    )
)

# Function to create a new DrugDose object
createDrugDose <- function(dose_info = NULL, num_doses = NULL) {
    return(DrugDose$new(dose_info = dose_info, num_doses = num_doses))
}
