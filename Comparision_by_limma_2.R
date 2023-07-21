## author: yao

library(limma)

# Extract gender from phenodata of gse40279_matrix
gender_40279 <- gse40279_matrix@phenoData@data[["characteristics_ch1.3"]]

# Function to extract gender from gender_40279 format.
extract_gender <- function(gender) {
  gender <- sub("^gender: ", "", gender)  # Remove "age: " prefix
  
  "# Convert 'M' to 1, and convert 'F' to 2.
  if (gender = 'M') {
    age <- 1
  }else{
    age <- 2
  }"
  
  return(gender)
}

# Apply the function to each element of the age vector
gender_40279_2 <- sapply(gender_40279, extract_gender)

# Assign the prcessed gender to a new column in pData of gse40279_matrix
pData(gse40279_matrix)$gender <- gender_40279_2



# Extract age from phenodata of gse40279_matrix
age <- pData(gse40279_matrix)$characteristics_ch1

# Remove "age (y):" and convert to numeric
age <- sub("^\\s*age \\(y\\): ", "", age)
age <- as.numeric(age)

# Assign age values to a new column in pData of gse40279_matrix
pData(gse40279_matrix)$age <- age


# Extract age from phenodata of gse40279
age_30870 <- pData(gse30870_matrix)$characteristics_ch1

# Function to extract numeric age from pData(gse30870_matrix)$characteristics_ch1 format.
extract_numeric_age <- function(age) {
  age <- sub("^age: ", "", age)  # Remove "age: " prefix
  age <- sub(" years$", "", age) # Remove " years" suffix
  
  # Convert numeric ages to numeric, and convert "Newborn" to 1.
  if (age != "Newborn") {
    age <- as.numeric(age)
  }else{
    age <- 1
  }
  
  return(age)
}

# Apply the function to each element of the age vector
age_30870 <- sapply(age_30870, extract_numeric_age)

# Remove the names from the resulting vector (Unnecessary)
# age_30870 <- unname(age_30870) 

# Assign age values to a new column in pData of gse40279_matrix
pData(gse30870_matrix)$age <- age_30870

# Assuming you have loaded the expression data into variables "gse30870_data" and "gse40279_data"
# Get the common probe IDs from both ExpressionSet objects
common_probes <- intersect(row.names(exprs(gse30870_matrix)), row.names(exprs(gse40279_matrix)))

probe_ID_30870 <- fData(gse30870_matrix)$ID
probe_ID_40279 <- fData(gse40279_matrix)$ID
# Subset the ExpressionSet objects to include only the common probes
subset_gse30870 <- gse30870_matrix[probe_ID_30870 %in% common_probes, ]
subset_gse40279 <- gse40279_matrix[probe_ID_40279 %in% common_probes, ]

# Load the phenotype data using pData()
phenotype_data_gse30870 <- pData(gse30870_matrix)
phenotype_data_gse40279 <- pData(gse40279_matrix)$gender

# Combine the phenotype data from GSE30870 and GSE40279
phenotype_data_combined <- rbind(phenotype_data_gse30870, phenotype_data_gse40279)

# Define the factor variables
factor_age <- pData(gse40279_matrix)$age

# factor_sample_size <- phenotype_data_combined$sample_size

# Define the mediator variable
mediator_gender <- pData(gse40279_matrix)$gender

# Define the confounder variables (assuming you have other phenotype variables in your dataset)
confounders <- phenotype_data_combined[, c("other_phenotype_variable_1", "other_phenotype_variable_2", ...)]

# Perform differential expression analysis and mediation analysis in chunks
chunk_size <- 1000  # Set the desired chunk size

# Split the common probes into chunks
probe_chunks <- split(common_probes, ceiling(seq_along(common_probes) / chunk_size))

# Split gse40279 probes into chunks
probe_chunks <- split(probe_ID_40279, ceiling(seq_along(probe_ID_40279) / chunk_size))


# Iterate over the probe chunks
results_list <- list()
for (chunk in probe_chunks) {
  # Subset the expression data for the current chunk
  subset_gse40279_chunk <- subset_gse40279[, chunk]
  
  # Merge the datasets based on the common probes
  exprs_gse40279_chunk <- exprs(subset_gse40279_chunk)
  
  # Perform differential expression analysis
  design <- model.matrix(~ factor_age + mediator_gender)
  fit <- lmFit(exprs_gse40279_chunk, design)
  fit <- eBayes(fit)
  results <- decideTests(fit, method = "global", adjust.method = "BH")
  DEGs_chunk <- rownames(subset_gse40279_chunk)[results$adj.P.Val < 0.05]
  
  # Perform mediation analysis
  model <- lm(exprs_gse40279_chunk[, DEGs_chunk] ~ factor_age + mediator_gender, data = phenotype_data_40279)
  mediation_results <- mediate(model, mediator_gender, exprs_gse40279_chunk[, DEGs_chunk])
  med_results <- summary(mediation_results)
  
  # Store the results for the current chunk
  results_list[[chunk]] <- med_results
}

# Combine the results from all chunks
final_results <- do.call(rbind, results_list)


results_list <- list()
for (chunk in probe_chunks) {
  
  # Subset the expression data for the current chunk
  subset_gse40279_chunk <- subset_gse40279[,chunk]
  
  # Error handling to handle out-of-bounds error
  if (chunk > ncol(subset_gse40279)) {
    print(paste("Skipping chunk", chunk, "- Invalid chunk index"))
    next
  }
  
  # Merge the datasets based on the common probes
  exprs_gse40279_chunk <- exprs(subset_gse40279_chunk)
  
  # Perform differential expression analysis
  design <- model.matrix(~ factor_age + mediator_gender)
  fit <- lmFit(exprs_gse40279_chunk, design)
  fit <- eBayes(fit)
  results <- decideTests(fit, method = "global", adjust.method = "BH")
  DEGs_chunk <- rownames(subset_gse40279_chunk)[results$adj.P.Val < 0.05]
  
  # Perform mediation analysis
  model <- lm(exprs_gse40279_chunk[, DEGs_chunk] ~ factor_age + mediator_gender, data = phenotype_data_40279)
  mediation_results <- mediate(model, mediator_gender, exprs_gse40279_chunk[, DEGs_chunk])
  med_results <- summary(mediation_results)
  
  # Store the results for the current chunk
  results_list[[chunk]] <- med_results
}

# Combine the results from all chunks
final_results <- do.call(rbind, results_list)



# Perform differential expression analysis
design <- model.matrix(~ factor_age + mediator_gender)
fit <- lmFit(gse40279_data, design)
fit <- eBayes(fit)
results <- decideTests(fit, method = "global", adjust.method = "BH")
results_2 <- topTable(fit, adjust.method = "BH")
# DEGs <- probe_ID_40279[results$adj.P.Val < 0.05]


# Perform mediation analysis
model <- lm(gse40279_data ~ factor_age + mediator_gender, data = phenotype_data_gse40279)
mediation_results <- mediate(model, mediator_gender, gse40279_data)
med_results <- summary(mediation_results)