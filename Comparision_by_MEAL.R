#### author: 姚博瀚

### Part 1
## Download gse40279 and gse30870 and run analysis by using "MEAL" package.
library(MEAL)
library(minfi)
library(limma)
library(ggplot2)

# Install remotes from CRAN (Optional)
install.packages("remotes")
# Download modified GEOquery package from my github 
# by using the function install_github() from 'remotes' package. 
library(remotes)
install_github("curryhank08/GEOquery_with_modifiable_timeout_seconds", force = TRUE)
# Load modified GEOquery
library(GEOquery)
# Setting the max timeout_seconds
options(timeout=100000)
# Check the input timeout_seconds
getOption("timeout")

# Download GSE40279 by the fuction getGEO() from modified GEOquery package.
gse40279 <- getGEO("GSE40279", GSEMatrix = TRUE, AnnotGPL = TRUE)
gse40279_matrix <- gse40279[[1]]
gse40279_data <- exprs(gse40279_matrix)

# Download GSE30870 by the fuction getGEO() from modified GEOquery package.
gse30870 <- getGEO("GSE30870", GSEMatrix = TRUE, AnnotGPL = TRUE)
gse30870_matrix <- gse30870[[1]]
gse30870_data <- exprs(gse30870_matrix)

# Create a new age from 'characteristics_ch1' and assign into pheno data of gse40279_matrix after processed.
age_40279 <- pData(gse40279_matrix)$characteristics_ch1

# Remove "age (y):" and convert to numeric
age_40279 <- sub("^\\s*age \\(y\\): ", "", age_40279)
age_40279 <- as.numeric(age_40279)

# The ^ character denotes the start of the string,
# \\s* matches any number of leading whitespace characters,
# and "age \\(y\\): " matches the exact string "age (y): ". 

# Count the numbers of 90above and 30below (30 and 90 are included)
gse40279_30below <- table(cut(age_40279, breaks = 0:30))
gse40279_30below <- data.frame(gse40279_30below)
colnames(gse40279_30below) = c("age", "Freq")

gse40279_90above <- table(cut(age_40279, breaks = 89:110))
gse40279_90above <- data.frame(gse40279_90above)

gse40279_age30_90 <- data.frame(below30 = sum(gse40279_30below), up90 = sum(gse40279_90above), row.names = "numbers")

age_values <- 1:30
freq_values <- table(factor(age_40279, levels = age_values))
gse40279_30below <- data.frame(freq_values)
colnames(gse40279_30below) = c("age", "Freq")

# Create the histogram plot
ggplot(gse40279_30below, aes(x = age, y = Freq)) +
  geom_col(fill = "skyblue", color = "black") +
  labs(title = "Age Distribution (gse40279)", x = "Age", y = "Count")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

a <- data.frame(age_40279)
ggplot(a, aes(x = age_40279)) + 
  geom_histogram()

# Create the histogram plot
ggplot(a, aes(x = factor(age_40279)) +
  geom_bar(fill = "skyblue", color = "black", stat = "count") +
  scale_x_continuous(breaks = unique(a$age_40279)) +
  labs(title = "Age Distribution", x = "Age", y = "Count") +
  theme_minimal()

'
# Assign age values to a new column in pData of gse40279_matrix
pData(gse40279_matrix)$age <- age_40279

# Define age categories based on specific age ranges
age_categories_40279 <- cut(age_40279,
                      breaks = c(0, 30, 89, Inf),
                      labels = c("Young", "Middle", "Old"),
                      include.lowest = TRUE)

# Assign age categories to the pData of gse40279_matrix
pData(gse40279_matrix)$age_category <- age_categories_40279

# Assign seqnames to the pData of gse40279_matrix from CHR
# fData(gse40279_matrix)$seqnames <- as.numeric(fData(gse40279_matrix)$CHR)
'
# Run MEAL pipeline on the categorized data
res_40279 <- runPipeline(set = gse40279_matrix,
                   variable_names = "age_category",
                   betas = TRUE,
                   analyses = c("DiffMean", "DiffVar"))

# Extract the result of the DiffMean analysis
result_40279 <- getProbeResults(res_40279, rid = 1, 
                               fNames = c("UCSC_RefGene_Name", "RANGE_START", "CHR", "ID"))

# Extract samples belonging to "Young" and "Old" age categories
gse40279_matrix_YO <- gse40279_matrix[, age_categories_40279 %in% c("Young", "Old")]
gse40279_matrix_YO_data <- exprs(gse40279_matrix_YO)
gse40279_matrix_YO$age_category <- factor(gse40279_matrix_YO$age_category, levels = c("Young", "Old"))

# subeset samples pheno info
gse40279_matrix_YO_info <- pData(gse40279_matrix_YO)

# Run MEAL pipeline on the subset YO
res_40279_YO <- runPipeline(set = gse40279_matrix_YO,
                          variable_names = "age_category",
                          betas = TRUE,
                          analyses = c("DiffMean", "DiffVar"))

# Extract the result of the DiffMean analysis
result_40279_YO <- getProbeResults(res_40279_YO, rid = 1, 
                                      fNames = c("UCSC_RefGene_Name", "RANGE_START", "CHR", "ID"))

# Remove rows with missing values
result_40279_YO_clean <- na.omit(result_40279_YO)


### Part 2 (30870)

# Create age categories
age_categories_30870 <- pData(gse30870_matrix)$source_name_ch1
age_categories_30870 <- as.factor(age_categories_30870)

# Assign age categories to the pData of gse40279_matrix
pData(gse30870_matrix)$age_category <- age_categories_30870

# Run MEAL pipeline on the categorized data
res_30870 <- runPipeline(set = gse30870_matrix,
                   variable_names = "age_category",
                   betas = TRUE,
                   analyses = c("DiffMean", "DiffVar"))

# Extract the result of the DiffMean analysis
result_30870 <- getProbeResults(res_30870, rid = 1, 
                               fNames = c("UCSC_RefGene_Name", "RANGE_START", "CHR", "ID"))

# Create a subset of gse30870_matrix to match gse40279_matrix_YO s probes
gse30870_matrix_sub <- gse30870_matrix[fData(gse30870_matrix)$ID %in% fData(gse40279_matrix_YO)$ID]

# Run MEAL pipeline on the categorized subset data
res_30870_sub <- runPipeline(set = gse30870_matrix_sub,
                         variable_names = "age_category",
                         betas = TRUE,
                         analyses = c("DiffMean", "DiffVar"))

# Extract the result of the DiffMean analysis
result_30870_sub <- getProbeResults(res_30870_sub, rid = 1, 
                                fNames = c("UCSC_RefGene_Name", "RANGE_START", "CHR", "ID"))


### Part 3
## Merge the results of analysis from gse30870 and gse40279, and then create a scatter plot.
library(ggplot2)

# Merge results
res_30870_YO_p <- data.frame(row.names = row.names(result_30870_sub), p_30870 = result_30870_sub$P.Value, UCSC_RefGene_name = result_30870_sub$UCSC_RefGene_Name, ProbeID = result_30870_sub$ID)
res_40279_YO_p <- data.frame(row.names = row.names(result_40279_YO), p_40279 = result_40279_YO$P.Value, UCSC_RefGene_name = result_40279_YO$UCSC_RefGene_Name, ProbeID = result_40279_YO$ID)

merged_30870_40279_p <- merge(res_30870_YO_p,
                          res_40279_YO_p,
                         by.x = "ProbeID",
                         by.y = "ProbeID")

plot(log(merged_30870_40279_p$p_30870), log(merged_30870_40279_p$p_40279), 
       xlab = "P.value from analysis of gse30870", ylab = "P.value from analysis of gse40279", 
       main = "Comparision of gse30870 and gse40279", 
       pch = 20, col = "#8bc34a", cex = 1)

condition1 <- (merged_30870_40279_p$p_30870 < 1e-5)
condition2 <- (merged_30870_40279_p$p_40279 < 1e-5)

filerd_data <- merged_30870_40279_p[condition1&condition2, ]

# Calculate the condition
condition <- abs(log10(matching_rows$p_30870 / matching_rows$p_40279)) > 25

# Subset the matching rows based on the condition
different_probe <- matching_rows[condition, c("probe_name_30870", "p_30870", "p_40279", "UCSC_RefGene_name.x")]

# Create a scatter plot with x-axis: p-value from limma and y-axis: p-value from MEAL.
plot(log(YO_30870_40279_p$p_30870), log(YO_30870_40279_p$p_40279), 
     xlab = "P.value from analysis of gse30870", ylab = "P.value from analysis of gse40279", 
     main = "Comparision of gse30870 and gse40279", 
     pch = 20, col = "#8bc34a", cex = 1)

'
different_probe <- data.frame(row.names = res_30870_YO_p$probe_name_30870)
colnames(different_probe) = c(p_30870, p_40279, UCSC_RefGene_name)

condition <- abs((res_30870_YO_p$p_30870)/(res_40279_YO_p$p_40279)) > 25
for x in res_30870_YO_p$probe_name_30870{
  if condition {
    rbind()
  }
} 
'


'
# unefficient way 
different_probe <- data.frame(matrix(ncol = 3))
colnames(different_probe) <- c("p_30870", "p_40279", "UCSC_RefGene_name")

for (x in res_30870_YO_p$probe_name_30870) {
  condition <- abs((res_30870_YO_p[res_30870_YO_p$probe_name_30870 == x, "p_30870"]) /
                     (res_40279_YO_p[res_40279_YO_p$probe_name_40279 == x, "p_40279"])) > 25
  
  row_to_add <- ifelse(condition,
                       c(res_30870_YO_p[res_30870_YO_p$probe_name_30870 == x, "p_30870"]),
                       NA)
  
  if (length(row_to_add) > 0 && !is.na(row_to_add)) {
    row_to_add <- c(row_to_add,
                    res_40279_YO_p[res_40279_YO_p$probe_name_40279 == x, "p_40279"],
                    res_30870_YO_p[res_30870_YO_p$probe_name_30870 == x, "UCSC_RefGene_name"])
    
    different_probe <- rbind(different_probe, row_to_add)
  }
}
'
