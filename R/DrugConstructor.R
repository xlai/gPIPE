Drug <- setRefClass("Drug",
    fields = list(
        name = "character",
        doses = "list"
    ),
    methods = list(
        initialize = function(drugData = NULL, name = NULL, doseCount = NULL) {
            if (!is.null(drugData)) {
                # Initialize from YAML file
                name <<- drugData$name
                doses <<- drugData$doses
            } else if (!is.null(name) && !is.null(doseCount)) {
                # Initialize with provided name and a default dose level structure
                name <<- name
                doseNames <- paste0("dose", seq_len(doseCount))
                # Creating a list with these names, initializing with NULL or any default value
                doses <<- setNames(as.list(seq_len(doseCount)), doseNames)
            } else {
                # NULL scenario
                name <<- "Unnamed Drug"
                doses <<- list()
            }
        },
        getDoseLabels = function() {
            # Returns the ordered dose labels for a given drug
            return(names(sort(unlist(doses))))
        },
        print = function() {
            cat("Drug Name:", name, "\nDose Levels:", .self$getDoseLabels(), "\n")
        },
        getNumberOfDoseLevels = function() {
            return(length(doses))
        },
        getDoseLevels = function(doseLabel = NULL) {
            # Validate input
            if (is.null(doseLabel)) {
                # Return all numeric levels for the drug
                return(as.numeric(sort(unlist(doses))))
            } else {
                # Return the numeric level for the specified dose labels
                if (!doseLabel %in% names(doses)) {
                    stop("Invalid dose label")
                }
                return(as.numeric(doses[doseLabel]))
            }
        },
        addDose = function(new_label = NULL, new_level) {
            # Validate input
            if (is.null(new_label)) {
                # Generate default label if not provided
                existing_labels <- names(doses)
                new_label_count <- sum(grepl("newlabel", existing_labels))
                new_label <- paste("newlabel", new_label_count + 1, sep="")
            }

            # Check if the new level is higher than all existing levels
            if (all(unlist(doses) < new_level)) {
                # If so, just add the new dose
                doses[[new_label]] <<- new_level
            } else {
                # Adjust existing labels' numeric order
                for (label in names(doses)) {
                    if (doses[[label]] >= new_level) {
                        doses[[label]] <<- doses[[label]] + 1
                    }
                }
            # Add the new label with its order
            doses[[new_label]] <<- new_level
            }
        }        
        # ... (other methods)
    )
)

# Function to create a new Drug object
createDrug <- function(drugData = NULL, name = NULL, doseCount = NULL) {
    if (is.null(drugData)){
        return( Drug$new(name = name, doseCount = doseCount))
    }
    else{
        return(
            lapply(
                drugData,
                function(d) Drug$new(drugData = d, name = name, doseCount = doseCount))
        )
    }
}