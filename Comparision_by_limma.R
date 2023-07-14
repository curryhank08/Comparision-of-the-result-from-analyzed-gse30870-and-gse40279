
library(limma)

# Load the expression data
# Assuming you have loaded the expression data into variables "exprs_gse30870" and "exprs_gse40279"

# Get the common probe IDs from both ExpressionSet objects
common_probes <- intersect(row.names(exprs(gse30870_matrix)), row.names(exprs(gse40279_matrix)))

# Subset the ExpressionSet objects to include only the common probes
subset_gse30870 <- gse30870_matrix[, fData(gse30870_matrix) %in% common_probes]
subset_gse40279 <- gse40279_matrix[, fData(gse40279_matrix) %in% common_probes]

# Create a common identifier in both datasets (e.g., probe or gene ID)
subset_gse30870_data <- exprs(subset_gse30870)
subset_gse40279_data <- exprs(subset_gse40279)

# Merge the datasets based on the common identifier
exprs_combined <- merge(subset_gse40279_data, subset_gse30870_data, all = TRUE)

# Load the phenotype data using pData()
phenotype_data_gse30870 <- pData(gse30870_matrix)
phenotype_data_gse40279 <- pData(gse40279_matrix)

# Combine the phenotype data from GSE30870 and GSE40279
phenotype_data_combined <- rbind(phenotype_data_gse30870, phenotype_data_gse40279)

# Define the factor variables
factor_age <- phenotype_data_combined$age
factor_sample_size <- phenotype_data_combined$sample_size

# Define the mediator variable
mediator_gender <- phenotype_data_combined$gender

# Define the confounder variables (assuming you have other phenotype variables in your dataset)
confounders <- phenotype_data_combined[, c("other_phenotype_variable_1", "other_phenotype_variable_2", ...)]

# Perform differential expression analysis
design <- model.matrix(~ factor_age + factor_sample_size + confounders)
fit <- lmFit(exprs_combined, design)
fit <- eBayes(fit)
results <- decideTests(fit, method = "global", adjust.method = "BH")
DEGs <- rownames(exprs_combined)[results$adj.P.Val < 0.05]

# Perform mediation analysis
model <- lm(exprs_combined[, DEGs] ~ factor_age + factor_sample_size + mediator_gender + confounders, data = phenotype_data_combined)
mediation_results <- mediate(model, mediator_gender, exprs_combined[, DEGs], confounders)
med_results <- summary(mediation_results)
