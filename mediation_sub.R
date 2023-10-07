# install 'meditaion' package
install.packages('mediation')

library(limma)
library(MEAL)
library(mediation)

# Assuming the probe IDs are listed in the 'probe_id' column of the dataset
probe_ids_ars0.7 <- unique(result_df_4_ars0.7$probe_id)


# Create an empty data frame to store the results
mediation_results_ars07 <- data.frame(probe_id = character(0),
                                ACME_average_estimate = numeric(0),
                                ACME_average_pvalue = numeric(0),
                                ADE_average_estimate = numeric(0),
                                ADE_average_pvalue = numeric(0),
                                Total_effect_estimate = numeric(0),
                                proportions_mediated_control = numeric(0),
                                proportions_mediated_treated = numeric(0),
                                proportions_mediated_average_estimate = numeric(0),
                                proportions_mediated_control_pvalue = numeric(0),
                                proportions_mediated_treated_pvalue = numeric(0),
                                proportions_mediated_average_pvalue = numeric(0))


# Assuming the probe IDs are listed in the 'probe_id' column of the dataset
probe_ids_ars0.5 <- unique(result_df_4_ars0.5$probe_id)

# Create an empty data frame to store the results
mediation_results_ars05 <- data.frame(probe_id = character(0),
                                      ACME_average_estimate = numeric(0),
                                      ACME_average_pvalue = numeric(0),
                                      ADE_average_estimate = numeric(0),
                                      ADE_average_pvalue = numeric(0),
                                      Total_effect_estimate = numeric(0),
                                      proportions_mediated_control = numeric(0),
                                      proportions_mediated_treated = numeric(0),
                                      proportions_mediated_average_estimate = numeric(0),
                                      proportions_mediated_control_pvalue = numeric(0),
                                      proportions_mediated_treated_pvalue = numeric(0),
                                      proportions_mediated_average_pvalue = numeric(0))


# Define a function to fit the mediation model and extract relevant information
fit_mediation <- function(probe_id) {
  # Select data for the current probe
  probe_data <- subset(gse40279_matrix, fData(gse40279_matrix)$ID == probe_id)
  
  # Convert "F" and "M" into 0 and 1
  probe_data$gender_binary <- as.numeric(probe_data$`gender:ch1` == "M")
  
  # Extract the age, mediator, and beta value as the predictor, mediator, and response variables
  x1 <- probe_data$age
  gender <- probe_data$gender_binary
  y <- assayData(probe_data)$exprs[1, ]
  
  # Fit the mediation model
  fit_M <-lm(gender ~ x1)
  fit_Y <- lm(y ~ x1 + gender + x1:gender)
  
  # Run the mediation analysis
  med_result <- mediate(fit_M, fit_Y, treat = 'x1', mediator = 'gender')
  
  # Summary of mediation analysis for a probe
  summary_mediation <- summary(med_result)
  
  # Extract relevant information from the mediation result
  ACME_average_estimate <- summary_mediation$d.avg
  ACME_average_pvalue <- summary_mediation$d.avg.p
  ADE_average_estimate <- summary_mediation$z.avg
  ADE_average_pvalue <- summary_mediation$z.avg.p
  Total_effect_estimate <- summary_mediation$tau.coef
  
  proportions_mediated_control <- summary_mediation$n0
  proportions_mediated_control_pvalue <- summary_mediation$n0.p
  proportions_mediated_treated <- summary_mediation$n1
  proportions_mediated_treated_pvalue <- summary_mediation$n1.p
  proportions_mediated_average_estimate <- summary_mediation$n.avg
  proportions_mediated_average_pvalue <- summary_mediation$n.avg.p
  
  return(data.frame(probe_id = probe_id,
                    ACME_average_estimate = ACME_average_estimate,
                    ACME_average_pvalue = ACME_average_pvalue,
                    ADE_average_estimate = ADE_average_estimate,
                    ADE_average_pvalue = ADE_average_pvalue,
                    Total_effect_estimate = Total_effect_estimate,
                    proportions_mediated_control = proportions_mediated_control,
                    proportions_mediated_treated = proportions_mediated_treated,
                    proportions_mediated_average_estimate = proportions_mediated_average_estimate,
                    proportions_mediated_control_pvalue = proportions_mediated_control_pvalue,
                    proportions_mediated_treated_pvalue = proportions_mediated_treated_pvalue,
                    proportions_mediated_average_pvalue = proportions_mediated_average_pvalue))
  
}


# Apply the function to each probe ID and rbind the results to the data frame
mediation_results_ars07 <- do.call(rbind, lapply(probe_ids_ars0.7, fit_mediation))

# Apply the function to each probe ID and rbind the results to the data frame
mediation_results_ars05 <- do.call(rbind, lapply(probe_ids_ars0.5, fit_mediation))

# Adjust the p.values from above
mediation_results$adj.p_value_indirect <- p.adjust(mediation_results$p_value_indirect, method = "BH")

# Select probes based on the significance of the indirect effect
significant_probes <- mediation_results[mediation_results$adj.p_value_indirect < 0.05, ]

# Probe ids of interest
probe_id_of_interest <- significant_probes$probe_id

# Select data for training
probe_data_of_interest <- subset(gse40279_matrix, fData(gse40279_matrix)$ID %in% probe_id_of_interest)











## Debug

# Create an empty data frame to store the results
mediation_results <- data.frame(probe_id = character(0),
                                ACME_average_estimate = numeric(0),
                                ACME_average_pvalue = numeric(0),
                                ADE_average_estimate = numeric(0),
                                ADE_average_pvalue = numeric(0),
                                Total_effect_estimate = numeric(0),
                                proportions_mediated_control = numeric(0),
                                proportions_mediated_control_pvalue = numeric(0),
                                proportions_mediated_treated = numeric(0),
                                proportions_mediated_treated_pvalue = numeric(0),
                                proportions_mediated_average_estimate = numeric(0),
                                proportions_mediated_average_pvalue = numeric(0))

# Run the mediation analysis for the specified probe
probe_id <- 'cg24724428'

# Select data for the selected probe
probe_data <- subset(gse40279_matrix, fData(gse40279_matrix)$ID == probe_id)

# Assuming "F" is Female and "M" is Male
probe_data$gender_binary <- as.numeric(probe_data$`gender:ch1` == "M")

# Extract the age, mediator, and beta value as the predictor, mediator, and response variables
x1 <- probe_data$age
gender <- probe_data$gender_binary
y <- assayData(probe_data)$exprs[1, ]

# Fit the mediation model
fit_M <- lm(gender ~ x1)
fit_Y <- lm(y ~ x1 + gender + x1:gender)

# Run the mediation analysis
med_result <- mediate(fit_M, fit_Y, treat = 'x1', mediator = 'gender')

# Summary of mediation analysis for a probe
summary_mediation <- summary(med_result)

# Extract relevant information from the mediation result
ACME_average_estimate <- summary_mediation$d.avg
ACME_average_pvalue <- summary_mediation$d.avg.p
ADE_average_estimate <- summary_mediation$z.avg
ADE_average_pvalue <- summary_mediation$z.avg.p
Total_effect_estimate <- summary_mediation$tau.coef

proportions_mediated_control <- summary_mediation$n0
proportions_mediated_control_pvalue <- summary_mediation$n0.p
proportions_mediated_treated <- summary_mediation$n1
proportions_mediated_treated_pvalue <- summary_mediation$n1.p
proportions_mediated_average_estimate <- summary_mediation$n.avg
proportions_mediated_average_pvalue <- summary_mediation$n.avg.p

mediation_test_df <- data.frame(probe_id = probe_id,
                                ACME_average_estimate = ACME_average_estimate,
                                ACME_average_pvalue = ACME_average_pvalue,
                                ADE_average_estimate = ADE_average_estimate,
                                ADE_average_pvalue = ADE_average_pvalue,
                                Total_effect_estimate = Total_effect_estimate,
                                proportions_mediated_control = proportions_mediated_control,
                                proportions_mediated_treated = proportions_mediated_treated,
                                proportions_mediated_average_estimate = proportions_mediated_average_estimate,
                                proportions_mediated_control_pvalue = proportions_mediated_control_pvalue,
                                proportions_mediated_treated_pvalue = proportions_mediated_treated_pvalue,
                                proportions_mediated_average_pvalue = proportions_mediated_average_pvalue)




ACME_estimate <- summary_mediation$d.avg
ACME_pvalue <- summary_mediation$d.avg.p
ADE_estimate <- summary_mediation$z.avg
ACME_pvalue <- summary_mediation$z.avg.p

Total_effect_estimate <- summary_mediation$tau.coef
proportions_mediated_estimate <- summary_mediation$n.avg
proportions_mediated_estimate

# Print the summary
print(summary_mediation)

# Check ACME_estimate (average)
ACME_control_estimate <- summary_mediation$d0
ACME_treated_estimate <- summary_mediation$d1
(ACME_control_estimate + ACME_treated_estimate)/2 == ACME_estimate

# Check ADE_estimate (average)
ADE_control_estimate <- summary_mediation$z0
ADE_treated_estimate <- summary_mediation$z1
(ADE_control_estimate + ADE_treated_estimate)/2 == ADE_estimate

# Check proportions_mediated_estimate (average)
proportions_mediated_control <- summary_mediation$n0
proportions_mediated_treated <- summary_mediation$n1
(proportions_mediated_control + proportions_mediated_treated)/2 == proportions_mediated_estimate

# Check Total_effect_estimate
total_effect_control <- ACME_control_estimate + ADE_control_estimate
total_effect_treated <- ACME_treated_estimate + ADE_treated_estimate
(total_effect_control + total_effect_treated)/2 == Total_effect_estimate

# Check proportions_mediated_control
ACME_control_estimate / proportions_mediated_control

total_effect_control
proportions_mediated_control






