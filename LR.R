## author: haku yao
install.packages("ggplot2")
install.packages("devtools")
library(devtools)
install_github('seandavi/GEOquery')
library(limma)
library(MEAL)

# Assuming the probe IDs are listed in the 'probe_id' column of the dataset
probe_ids <- unique(fData(gse40279_matrix)$ID)

# Create an empty list to store the regression results for each probe
probe_regression_results <- list()

for (probe_id in probe_ids) {
  # Select data for the current probe
  probe_data <- subset(gse40279_matrix, fData(gse40279_matrix)$ID == probe_id)
  
  # Extract the age and beta value as the predictor and response variables
  x <- probe_data$age
  y <- assayData(probe_data)$exprs[1, ]
  
  # Create and fit the linear regression model
  model <- lm(y ~ x)
  
  # Store the probe ID and regression model in the results list
  probe_regression_results[[probe_id]] <- model
}


# Assuming the probe IDs are listed in the 'probe_id' column of the dataset
probe_ids <- unique(fData(gse40279_matrix)$ID)

# Define a function to fit the linear regression model
fit_regression <- function(probe_id) {
  # Select data for the current probe
  probe_data <- subset(gse40279_matrix, fData(gse40279_matrix)$ID == probe_id)
  
  # Extract the age and beta value as the predictor and response variables
  x <- probe_data$age
  y <- assayData(probe_data)$exprs[1, ]
  
  # Create and fit the linear regression model
  model <- lm(y ~ x)
  
  return(model)
}

# Apply the function to each probe ID and store the results in a list
probe_regression_results <- lapply(probe_ids, fit_regression)



# Assuming the probe IDs are listed in the 'probe_id' column of the dataset
probe_ids <- unique(fData(gse40279_matrix)$ID)

# Define a function to fit the linear regression model and extract relevant information
fit_regression <- function(probe_id) {
  # Select data for the current probe
  probe_data <- subset(gse40279_matrix, fData(gse40279_matrix)$ID == probe_id)
  
  # Extract the age and beta value as the predictor and response variables
  x <- probe_data$age
  y <- assayData(probe_data)$exprs[1, ]
  
  # Create and fit the linear regression model
  model <- lm(y ~ x)
  
  # Extract and return relevant information from the model
  model_summary <- summary(model)
  coefficients <- coef(model)
  
  return(list(probe_id = probe_id, model_summary = model_summary, coefficients = coefficients))
}

# Apply the function to each probe ID and store the results in a list
probe_regression_results <- lapply(probe_ids, fit_regression)



# Assuming the probe IDs are listed in the 'probe_id' column of the dataset
probe_ids <- unique(fData(gse40279_matrix)$ID)

# Create an empty data frame to store the results
result_df <- data.frame(probe_id = character(0),
                        coefficients = numeric(0),
                        r_squared = numeric(0),
                        p_value = numeriv(0),
                        # Add more columns as needed
                        stringsAsFactors = FALSE)

# Define a function to fit the linear regression model and extract relevant information
fit_regression <- function(probe_id) {
  # Select data for the current probe
  probe_data <- subset(gse40279_matrix, fData(gse40279_matrix)$ID == probe_id)
  
  # Extract the age and beta value as the predictor and response variables
  x <- probe_data$age
  y <- assayData(probe_data)$exprs[1, ]
  
  # Create and fit the linear regression model
  model <- lm(y ~ x)
  
  # Extract relevant information from the model
  model_summary <- summary(model)
  coefficients <- coef(model)
  r_squared <- summary(model)$r.squared
  p_value <- summary(model)$coefficients["X", "Pr(>|t|)"]
  
  return(data.frame(probe_id = probe_id,
                    coefficients = coefficients["x"],  # Extract coefficient for predictor 'x'
                    r_squared = r_squared,
                    p_value = p_value))
}

# Apply the function to each probe ID and rbind the results to the data frame
result_df <- do.call(rbind, lapply(probe_ids, fit_regression))




# Assuming the probe IDs are listed in the 'probe_id' column of the dataset
probe_ids <- unique(fData(gse40279_matrix)$ID)

# Create an empty data frame to store the results
result_df_2 <- data.frame(probe_id = character(0),
                        coefficients = numeric(0),
                        p_value = numeric(0),
                        r_squared = numeric(0))

# Define a function to fit the linear regression model and extract relevant information
fit_regression_2 <- function(probe_id) {
  # Select data for the current probe
  probe_data <- subset(gse40279_matrix, fData(gse40279_matrix)$ID == probe_id)
  
  # Extract the age and beta value as the predictor and response variables
  x <- probe_data$age
  y <- assayData(probe_data)$exprs[1, ]
  
  # Create and fit the linear regression model
  model <- lm(y ~ x)
  
  # Extract relevant information from the model
  model_summary <- summary(model)
  coefficients <- coef(model)
  p_value <- summary(model)$coefficients[2, "Pr(>|t|)"]
  r_squared <- summary(model)$r.squared
  
  return(data.frame(probe_id = probe_id,
                    coefficients = coefficients["x"],  # Extract coefficient for predictor 'x'
                    p_value = p_value,
                    r_squared = r_squared))
}

# Apply the function to each probe ID and rbind the results to the data frame
result_df_2 <- do.call(rbind, lapply(probe_ids, fit_regression_2))



# Assuming the probe IDs are listed in the 'probe_id' column of the dataset
probe_ids <- unique(fData(gse40279_matrix)$ID)

# Create an empty data frame to store the results
result_df_3 <- data.frame(probe_id = character(0),
                          coefficient_intercept = numeric(0),
                          coefficient_x1 = numeric(0),
                          coefficient_x2 = numeric(0),
                          p_value_intercept = numeric(0),
                          p_value_x1 = numeric(0),
                          p_value_x2 = numeric(0),
                          r_squared = numeric(0),
                          adjusted_r_squared = numeric(0),
                          p_value_f_statistic = numeric(0))

# Define a function to fit the linear regression model and extract relevant information
fit_regression_3 <- function(probe_id) {
  # Select data for the current probe
  probe_data <- subset(gse40279_matrix, fData(gse40279_matrix)$ID == probe_id)
  
  # Extract the age and beta value as the predictor and response variables
  x1 <- probe_data$age
  x2 <- probe_data$`gender:ch1`
  y <- assayData(probe_data)$exprs[1, ]
  
  # Create and fit the linear regression model
  model <- lm(y ~ x1 + x2)
  
  # Extract relevant information from the model
  model_summary <- summary(model)
  coefficients <- coef(model)
  p_value_x1 <- summary(model)$coefficients[2, "Pr(>|t|)"]
  p_value_x2 <- summary(model)$coefficients[3, "Pr(>|t|)"]
  r_squared <- summary(model)$r.squared
  adjusted_r_squared <- summary(model)$adj.r.squared
  f_statistic <- model_summary$fstatistic
  p_value_f_statistic <- pf(f_statistic[1], f_statistic[2], f_statistic[3], lower.tail = FALSE)
  
  return(data.frame(probe_id = probe_id,
                    coefficient_intercept = coefficients["(Intercept)"],
                    coefficient_x1 = coefficients["x1"],
                    coefficient_x2 = coefficients["x2M"],
                    p_value_intercept = coefficients["(Intercept)"],
                    p_value_x1 = p_value_x1,
                    p_value_x2 = p_value_x2,
                    r_squared = r_squared,
                    adjusted_r_squared = adjusted_r_squared,
                    p_value_f_statistic = p_value_f_statistic))
}

# Apply the function to each probe ID and rbind the results to the data frame
result_df_3 <- do.call(rbind, lapply(probe_ids, fit_regression_3))




## By limma (Not yet solved)
conditions <- gse40279_matrix$age
#f <- factor(conditions)
design <- model.matrix(~0+conditions)
colnames(design) <- "age"
fit <- lmFit(gse40279_matrix, design)
fit <- eBayes(fit)
result_limma <- topTable(fit, number = Inf, adjust.method = "BH", sort.by = "P")
result_limma_2 <- result_limma[, c("UCSC_RefGene_Name", "P.Value", "logFC")]



result_df_2$"-log10(p_value)" <- -log10(result_df_2$p_value)
library(ggplot2)
# Calculate the total number of p-values within the limitation
total_amount <- sum(-log10(result_df_2$p_value) > 30 & -log10(result_df_2$p_value) < 200)
ggplot(result_df_2, aes(x = -log10(p_value))) +
  geom_histogram(binwidth = 0.5, color = "black") +
  xlim(c(30, 200))+
  ylim(c(0, 75))+
  ggtitle("Distribution of log-transformed p-values") +
  geom_text(
    x = 100, y = 65,  # Position the label below the title
    vjust = 0,         # Adjust the vertical position below the title
    hjust = 0,         # Align the label to the center
    label = paste("Total amount of the probes with p.value < 1e-30 :", total_amount),
    color = "black",
    size = 3
  )

ggplot(result_df_2, aes(x = -log10(p_value))) +
  geom_histogram(binwidth = 0.5, color = "black") +
  xlim(c(30, 200))+
  ylim(c(0, 70))+
  ggtitle("Distribution of log-transformed p-values") +
  geom_text(
    x = 110, y = 65,  # Position the label below the title
    vjust = 0,         # Adjust the vertical position below the title
    hjust = 0,         # Align the label to the center
    label = paste("Total amount of the probes with p.value < 1e-30 :", total_amount),
    color = "black",
    size = 3
  ) +
  scale_x_continuous(
    breaks = seq(30, 200, by = 10),
    labels = seq(30, 200, by = 10),
    limits = c(30, 200),  # Set limits to match xlim
    expand = c(0, 0)      # Adjust the expand parameter to remove extra space
  ) +
  scale_y_continuous(
    breaks = seq(0, 70, by = 10),     # Specify the desired breaks for y-axis
    labels = seq(0, 70, by = 10),     # Specify the corresponding labels
    limits = c(0, 70),
    expand = c(0, 0)                  # Adjust the expand parameter for y-axis
  )

ggplot(result_df_2, aes(x = -log10(p_value))) +
  geom_histogram(binwidth = 0.5)+
  xlim(c(30, 200))+
  ylim(c(0, 75))+
  ggtitle("Distribution of log-transformed p_value")

ggplot(result_df_2, aes(x = p_value)) +
  geom_histogram(binwidth = 0.01)+
  ggtitle("Distribution of log-transformed p_value")

  
# Select Probe
probe_id <- 'cg11643285'

# Select data for the current probe
probe_data <- subset(gse40279_matrix, fData(gse40279_matrix)$ID == probe_id)

# Extract the age and beta value as the predictor and response variables
x1 <- probe_data$age
x2 <- probe_data$`gender:ch1`
y <- assayData(probe_data)$exprs[1, ]

# Create and fit the linear regression model
model <- lm(y ~ x1 + x2)

# Extract relevant information from the model
model_summary <- summary(model)
coefficients <- coef(model)
p_value_x1 <- summary(model)$coefficients[2, "Pr(>|t|)"]
p_value_x2 <- summary(model)$coefficients[3, "Pr(>|t|)"]
r_squared <- summary(model)$r.squared
adjusted_r_squared <- summary(model)$adj.r.squared
f_statistic <- model_summary$fstatistic
p_value_f_statistic <- pf(f_statistic[1], f_statistic[2], f_statistic[3], lower.tail = FALSE)

test_reuslt <- data.frame(probe_id = probe_id,
           coefficient_intercept = coefficients["(Intercept)"],
           coefficient_x1 = coefficients["x1"],
           coefficient_x2 = coefficients["x2M"],
           p_value_intercept = coefficients["(Intercept)"],
           p_value_x1 = p_value_x1,
           p_value_x2 = p_value_x2,
           r_squared = r_squared,
           adjusted_r_squared = adjusted_r_squared,
           p_value_f_statistic = p_value_f_statistic)


# Select data for the current probe
probe_data <- subset(gse40279_matrix, fData(gse40279_matrix)$ID == 'cg01403239')

# Extract the age and beta value as the predictor and response variables
X1 <- probe_data$age
X2 <- probe_data$`gender:ch1`
y <- assayData(probe_data)$exprs[1, ]

# Create and fit the linear regression model
model <- lm(y ~ X1+X2)
coef(model)["(Intercept)"]
summary_model <- summary.lm(model)
f_statistic <- summary_model$fstatistic
p_value_f_statistic <- pf(f_statistic[1], f_statistic[2], f_statistic[3], lower.tail = FALSE)



# Get the regression model for the specific probe
model_probe_cg16867657 <- probe_regression_results[['cg16867657']]

model_probe_cg22454769 <- probe_regression_results[['cg22454769']]

model_probe_cg23201812 <- probe_regression_results[['cg23201812']]

model_probe_cg04875128 <- probe_regression_results[['cg04875128']]

cor_cg16867657 <- cor(model_probe_cg16867657$model$x, y = model_probe_cg16867657$model$y, use = "everything", 
                      method = "pearson")
cor_cg22454769 <- cor(model_probe_cg22454769$model$x, y = model_probe_cg22454769$model$y, use = "everything",
                      method = "pearson")
cor_cg23201812 <- cor(model_probe_cg23201812$model$x, y = model_probe_cg23201812$model$y, use = "everything",
                      method = "pearson")
cor_cg04875128 <- cor(model_probe_cg04875128$model$x, y = model_probe_cg04875128$model$y, use = "everything",
                      method = "pearson")

## Create a data frame with all probes corresponding pearson correlation coefficients  
# Define a function to compute the correlation coefficient for a given probe's model
get_correlation_coefficient <- function(model) {
  correlation_coefficient <- cor(model$x, y = model$y, use = "everything", method = "pearson")
  return(correlation_coefficient)
}

# Create a vector of probe IDs
probe_ids <- names(probe_regression_results)

# Use lapply to apply the function to each probe's model
correlation_coefficients <- lapply(probe_ids, function(probe_id) {
  model <- probe_regression_results[[probe_id]]$model
  get_correlation_coefficient(model)
})

# Combine the results into a data frame
correlation_df <- data.frame(probe_id = probe_ids, correlation_coefficient = unlist(correlation_coefficients))
##




## Plot a specific probe's linear regression model with correlation coefficient
# Replace 'cg16867657' with the specific probe ID you want to plot
probe_id_of_interest <- 'cg12803060'

# Select data for the current probe
probe_data <- subset(gse40279_matrix, fData(gse40279_matrix)$ID == probe_id_of_interest)

# Extract the age and beta value as the predictor and response variables
X1 <- probe_data$age
X2 <- probe_data$`gender:ch1`
y <- assayData(probe_data)$exprs[1, ]

# Create and fit the linear regression model
model <- lm(y ~ X1+X2)
model_of_interest <- model

# Assuming you have the 'X' and 'y' data for the specific probe
# You can also create a data frame if you have the original data
data_of_interest <- data.frame(
  X1 = model_of_interest$model$X1,
  X2 = model_of_interest$model$X2,
  y = model_of_interest$model$y
)

# Calculate the correlation coefficient between 'X' and 'y'
cor_coefficient_X1 <- cor(data_of_interest$X1, data_of_interest$y)

# Calculate the total number of points
total_points <- nrow(data_of_interest)

library(ggplot2)
# Plot the data points and the linear regression line
plot_lm <- ggplot(data_of_interest, aes(x = X1, y = y, color = X2)) +
  geom_point(size = 1) +  # Scatter plot of data points
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Linear regression line
  geom_text(x = min(data_of_interest$X1), y = max(data_of_interest$y),
            label = paste("Correlation_X1:", round(cor_coefficient, 3)),
            hjust = 0, vjust = -1, color = "black", size = 4) +  # Add correlation coefficient
  geom_text(x = min(data_of_interest$X1), y = max(data_of_interest$y) - 0.03,  # Adjust y value to place it below the correlation
            label = paste("Total numbers of points:", total_points),
            hjust = 0, vjust = 0, color = "black", size = 4) +  # Add total number of points
  labs(title = paste("Linear Model for Probe:", probe_id_of_interest),
       x = "Age",
       y = "Beta Value") +
  theme(plot.margin = margin(t = 50, r = 10, b = 20, l = 10, unit = "pt"))  # Adjust the bottom margin

# Plot the data points and the linear regression line
plot_lm <- ggplot(data_of_interest, aes(x = X1, y = y, color = X2)) +
  geom_point(size = 1) +  # Scatter plot of data points with color mapping
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Linear regression line
  geom_text(x = min(data_of_interest$X1), y = max(data_of_interest$y),
            label = paste("Correlation_X1:", round(cor_coefficient_X1, 3)),
            hjust = 0, vjust = -0.6, color = "black", size = 4) +  # Adjust vjust to place it below the title
  geom_text(x = min(data_of_interest$X1), y = max(data_of_interest$y) - 0.03,  # Adjust y value to place it below the correlation
            label = paste("Total numbers of points:", total_points),
            hjust = 0, vjust = -0.8, color = "black", size = 4) +  # Adjust vjust to place it below the title
  labs(title = paste("Linear Model for Probe:", probe_id_of_interest),
       x = "Age",
       y = "Beta Value") +
  theme(
    plot.title = element_text(hjust = 0, vjust = 15, margin = margin(b = 20, unit = "pt")),  # Adjust title positioning
    plot.margin = margin(t = 50, r = 10, b = 20, l = 10, unit = "pt")  # Adjust the bottom margin
  )

# Plot the data points and the linear regression line
plot_lm <- ggplot(data_of_interest, aes(x = X1, y = y, color = X2)) +
  geom_point(size = 1) +  # Scatter plot of data points with color mapping
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Linear regression line
  labs(title = paste("Linear Model for Probe:", probe_id_of_interest),
       x = "Age",
       y = "Beta Value") +
  annotate("text", x = min(data_of_interest$X1), y = max(data_of_interest$y) + 0.01,
           label = paste("Correlation_X1:", round(cor_coefficient_X1, 3)),
           hjust = 0, color = "black", size = 4) +  # Adjust y to place it above the plot area
  annotate("text", x = min(data_of_interest$X1), y = min(data_of_interest$y) - 0.01,
           label = paste("Total numbers of points:", total_points),
           hjust = 0, color = "black", size = 4) +  # Adjust y to place it below the plot area
  theme(
    plot.title = element_text(hjust = 0, vjust = 10, margin = margin(b = 20, unit = "pt")),  # Adjust title positioning
    plot.margin = margin(t = 40, r = 10, b = 10, l = 10, unit = "pt")  # Adjust the bottom margin
  )

# Plot the data points and the linear regression line
plot_lm <- ggplot(data_of_interest, aes(x = X1, y = y, color = X2)) +
  geom_point(size = 1) +  # Scatter plot of data points with color mapping
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Linear regression line
  labs(title = paste("Linear Model for Probe:", probe_id_of_interest),
       x = "Age",
       y = "Beta Value",
       subtitle = c(paste("Correlation_X1:", round(cor_coefficient_X1, 3)), paste("Total numbers of points:", total_points))) +
  theme(
    plot.title = element_text(hjust = 0, vjust = 10, margin = margin(b = 20, unit = "pt")),  # Adjust title positioning
    plot.subtitle = element_text(hjust = 0, vjust = 15, margin = margin(b = 20, unit = "pt")),  # Adjust subtitle positioning
    plot.margin = margin(t = 40, r = 10, b = 10, l = 10, unit = "pt")  # Adjust the bottom margin
  )


# Create two subtitle lines
subtitle_lines <- c(
  paste("Correlation_X1:", round(cor_coefficient_X1, 3)),
  paste("Total numbers of points:", total_points)
)

# Plot the data points and the linear regression line
plot_lm <- ggplot(data_of_interest, aes(x = X1, y = y, color = X2)) +
  geom_point(size = 1) +  # Scatter plot of data points with color mapping
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Linear regression line
  labs(title = paste("Linear Model for Probe:", probe_id_of_interest),
       x = "Age",
       y = "Beta Value",
       subtitle = subtitle_lines) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 20, unit = "pt")),  # Adjust title positioning
    plot.subtitle = element_text(hjust = 0.5, margin = margin(b = 10, unit = "pt")),  # Adjust subtitle positioning
    plot.margin = margin(t = 40, r = 10, b = 10, l = 10, unit = "pt")  # Adjust the bottom margin
  )

# Plot the data points and the linear regression line
plot_lm <- ggplot(data_of_interest, aes(x = X1, y = y, color = X2)) +
  geom_point(size = 1) +  # Scatter plot of data points with color mapping
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Linear regression line
  labs(title = paste("Linear Model for Probe:", probe_id_of_interest),
       x = "Age",
       y = "Beta Value",
       subtitle = paste("Correlation_X1:", round(cor_coefficient_X1, 3)),
       caption = paste("Total numbers of points:", total_points)) +
  theme(
    plot.title = element_text(hjust = 0, vjust = 5),  # Adjust title positioning
    plot.subtitle = element_text(hjust = 0, vjust = 0),  # Adjust subtitle positioning
    plot.caption = element_text(hjust = 0, vjust = 185),  # Adjust caption positioning
    plot.margin = margin(t = 15, r = 0, b = 0, l = 5, unit = "pt")  # Adjust the bottom margin
  )

# Show the plot
print(plot_lm)
##

# Coefficients
coefficients(model_probe_cg16867657)

# R-squared value
summary(model_probe_cg16867657)$r.squared


# Merge results
res_30870_YO_p <- data.frame(row.names = row.names(result_30870_sub), p_30870 = result_30870_sub$P.Value, UCSC_RefGene_name = result_30870_sub$UCSC_RefGene_Name, ProbeID = result_30870_sub$ID)
res_40279_LR_p <- data.frame(row.names = result_df_2$probe_id, p_40279 = result_df_2$p_value, ProbeID = result_df_2$probe_id)

merged_30870_40279_p_2 <- merge(res_30870_YO_p,
                              res_40279_LR_p,
                              by.x = "ProbeID",
                              by.y = "ProbeID")

plot(log10(merged_30870_40279_p_2$p_30870), log10(merged_30870_40279_p_2$p_40279), 
     xlab = "log10(P.value) from analysis of gse30870", ylab = "log10(P.value) from LR analysis of gse40279", 
     main = "Scatter plot for gse30870 and gse40279 (8/8)", 
     pch = 20, col = "#8bc34a", cex = 1)

correlation_for_merge_p2 <- cor(log10(merged_30870_40279_p_2$p_30870), log10(merged_30870_40279_p_2$p_40279))
text(-45, -180, sprintf("Correlation: %.4f", correlation_for_merge_p2), adj = 0)


# Merge results with ars > 0.6
res_30870_YO_p <- data.frame(row.names = row.names(result_30870_sub), p_30870 = result_30870_sub$P.Value, UCSC_RefGene_name = result_30870_sub$UCSC_RefGene_Name, ProbeID = result_30870_sub$ID)
ars_06up_40279 <- result_df_3[result_df_3$adjusted_r_squared > 0.6, ]
res_40279_LR3_p <- data.frame(row.names = ars_06up_40279$probe_id, p_40279 = ars_06up_40279$p_value_x1, ProbeID = ars_06up_40279$probe_id)

merged_30870_40279_p_3 <- merge(res_30870_YO_p,
                                res_40279_LR3_p,
                                by.x = "ProbeID",
                                by.y = "ProbeID")
# Assuming merged_30870_40279_p_3$p_30870 and merged_30870_40279_p_3$p_40279 are your data vectors
x <- log10(merged_30870_40279_p_3$p_30870)
y <- log10(merged_30870_40279_p_3$p_40279)

# Fit a linear model
model <- lm(y ~ x)

# Calculate the R-squared and correlation values
ar_squared <- summary.lm(model)$adj.r.squared
correlation <- cor(x, y)

plot(log10(merged_30870_40279_p_3$p_30870), log10(merged_30870_40279_p_3$p_40279), 
     xlab = "log10(P.value) from analysis of gse30870", ylab = "log10(P.value) from LR analysis of gse40279", 
     main = "Scatter plot for gse30870 and gse40279 (probes with ars >0.6)", 
     pch = 20, col = "#8bc34a", cex = 1, xlim = c(-20, 0), ylim = c(-20, 0))

# Add R-squared and correlation values to the plot
text(-20, 0, sprintf("R-squared: %.4f", ar_squared), adj = 0)
text(-20, -1, sprintf("Correlation: %.4f", correlation), adj = 0)



# Merge results with ars > 0.5
res_30870_YO_p <- data.frame(row.names = row.names(result_30870_sub), p_30870 = result_30870_sub$P.Value, UCSC_RefGene_name = result_30870_sub$UCSC_RefGene_Name, ProbeID = result_30870_sub$ID)
ars_05up_40279 <- result_df_3[result_df_3$adjusted_r_squared > 0.5, ]
res_40279_LR3_p <- data.frame(row.names = ars_05up_40279$probe_id, p_40279 = ars_05up_40279$p_value_x1, ProbeID = ars_05up_40279$probe_id)

merged_30870_40279_p_4 <- merge(res_30870_YO_p,
                                res_40279_LR3_p,
                                by.x = "ProbeID",
                                by.y = "ProbeID")
# Assuming merged_30870_40279_p_3$p_30870 and merged_30870_40279_p_3$p_40279 are your data vectors
x <- log10(merged_30870_40279_p_4$p_30870)
y <- log10(merged_30870_40279_p_4$p_40279)

# Fit a linear model
model <- lm(y ~ x)

# Calculate the R-squared and correlation values
ar_squared <- summary.lm(model)$adj.r.squared
correlation <- cor(x, y)

plot(log10(merged_30870_40279_p_4$p_30870), log10(merged_30870_40279_p_4$p_40279), 
     xlab = "log10(P.value) from analysis of gse30870", ylab = "log10(P.value) from LR analysis of gse40279", 
     main = "Scatter plot for gse30870 and gse40279 (probes with ars > 0.5)", 
     pch = 20, col = "#8bc34a", cex = 1, xlim = c(-20, 0), ylim = c(-20, 0))

# Add R-squared and correlation values to the plot
text(-20, 0, sprintf("Adjusted r-squared: %.4f", ar_squared), adj = 0)
text(-20, -1, sprintf("Correlation: %.4f", correlation), adj = 0)



condition1 <- (merged_30870_40279_p_2$p_30870 < 1e-10)
condition2 <- (merged_30870_40279_p_2$p_40279 < 1e-10)
filterd_data <- merged_30870_40279_p_2[condition1&condition2, ]

# Create a scatter plot with x-axis: p-value from analyzed gse30870 by MEAL and y-axis: p-value from analyzed gse40279 by linear regression model.
plot(log10(filterd_data$p_30870), log10(filterd_data$p_40279), 
     xlab = "log10(P.value) from analysis of gse30870", ylab = "log10(P.value) from analysis of gse40279", 
     main = "Scatter plot for gse30870(p.value<1e-10) and gse40279(p.value<1e-10)", 
     pch = 20, col = "#8bc34a", cex = 1)

# Assuming merged_30870_40279_p_3$p_30870 and merged_30870_40279_p_3$p_40279 are your data vectors
x <- log10(filterd_data$p_30870)
y <- log10(filterd_data$p_40279)

# Fit a linear model
model <- lm(y ~ x)

# Calculate the R-squared and correlation values
ar_squared <- summary.lm(model)$adj.r.squared
correlation <- cor(x, y)

# Add R-squared and correlation values to the plot
text(-22, -180, sprintf("Adjusted r-squared: %.4f", ar_squared), adj = 0)
text(-22, -190, sprintf("Correlation: %.4f", correlation), adj = 0)
