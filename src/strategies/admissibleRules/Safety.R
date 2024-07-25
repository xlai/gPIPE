SafetyAdmissible <- setRefClass(
  "SafetyAdmissible",
  contains = "AdmissibleCombinationRule", 
  methods = list(
    isAdmissible = function(doseConfig, currentDoseLevel, drugCombiModel, safety_threshold = 0.8) {
      # Custom implementation
      doseConfig_posterior <- doseConfig$updatePipeEstimator(drugCombiModel$p_posterior)
      sum_posterior <- Reduce("+", doseConfig_posterior)
      normalized_posterior <- lapply(doseConfig_posterior, function(x) x / sum_posterior)

      validConfig_list <- lapply(doseConfig$validConfigs, function(x) x$currentConfig)

      # Use Map to multiply element-wise
      result <- Map(function(x, y) x * y, normalized_posterior, validConfig_list)
      result_collapsed <-  Reduce("+", result)

      return(which(result_collapsed < safety_threshold))
    }
  )
)