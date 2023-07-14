library(limma)

# Assuming you have loaded the expression data into variables "exprs_gse30870" and "exprs_gse40279"

# Get the common probe IDs from both ExpressionSet objects
common_probes <- intersect(row.names(exprs_gse30870), row.names(exprs_gse40279))

# Subset the ExpressionSet objects to include only the common probes
subset_gse30870 <- exprs_gse30870[row.names(exprs_gse30870) %in% common_probes, ]
subset_gse40279 <- exprs_gse40279[row.names(exprs_gse40279) %in% common_probes, ]

# Load the phenotype data using pData()
phenotype_data_gse30870 <- pData(gse30870_matrix)[c("source_name_ch1")]
colnames(phenotype_data_gse30870) <- "age"
phenotype_data_gse40279 <- pData(gse40279_matrix)[c("age")]


# Convert "Newborns" to 1 and "Nonagenarians" to 95 in phenotype_data_gse30870
phenotype_data_gse30870$age <- ifelse(phenotype_data_gse30870$age == "Newborns", 1,
                                      ifelse(phenotype_data_gse30870$age == "Nonagenarians", 95,
                                             phenotype_data_gse30870$age))

# Combine the phenotype data from GSE30870 and GSE40279
phenotype_data_combined <- rbind(phenotype_data_gse30870, phenotype_data_gse40279)

# Define the factor variables
factor_age <- phenotype_data_combined$age
factor_age <- as.numeric(factor_age)
# factor_sample_size <- phenotype_data_combined$sample_size

# Define the mediator variable
mediator_gender <- phenotype_data_combined$gender

# Define the confounder variables (assuming you have other phenotype variables in your dataset)
confounders <- phenotype_data_combined[, c("other_phenotype_variable_1", "other_phenotype_variable_2", ...)]

# Perform differential expression analysis and mediation analysis in chunks
chunk_size <- 1000  # Set the desired chunk size

# Split the common probes into chunks
probe_chunks <- split(common_probes, ceiling(seq_along(common_probes) / chunk_size))

# Iterate over the probe chunks
results_list <- list()
for (chunk in probe_chunks) {
  # Subset the expression data for the current chunk
  subset_gse30870_chunk <- subset_gse30870[, chunk]
  subset_gse40279_chunk <- subset_gse40279[, chunk]
  
  # Merge the datasets based on the common probes
  exprs_combined_chunk <- merge(subset_gse40279_chunk, subset_gse30870_chunk)
  
  # Perform differential expression analysis
  design <- model.matrix(~ factor_age + factor_sample_size + confounders)
  fit <- lmFit(exprs_combined_chunk, design)
  fit <- eBayes(fit)
  results <- decideTests(fit, method = "global", adjust.method = "BH")
  DEGs_chunk <- rownames(exprs_combined_chunk)[results$adj.P.Val < 0.05]
  
  # Perform mediation analysis
  model <- lm(exprs_combined_chunk[, DEGs_chunk] ~ factor_age + factor_sample_size + mediator_gender + confounders, data = phenotype_data_combined)
  mediation_results <- mediate(model, mediator_gender, exprs_combined_chunk[, DEGs_chunk], confounders)
  med_results <- summary(mediation_results)
  
  # Store the results for the current chunk
  results_list[[chunk]] <- med_results
}

# Combine the results from all chunks
final_results <- do.call(rbind, results_list)
