library(limma)
library(MEAL)
age_40279 <- data.frame()

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


# Select data for the current probe
probe_data <- subset(gse40279_matrix, fData(gse40279_matrix)$ID == 'cg00000029')

# Extract the age and beta value as the predictor and response variables
X <- probe_data$age
y <- assayData(probe_data)$exprs[1, ]

# Create and fit the linear regression model
model <- lm(y ~ X)



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



# Replace 'cg16867657' with the specific probe ID you want to plot
probe_id_of_interest <- 'cg16867657'
model_of_interest <- probe_regression_results[[probe_id_of_interest]]$model

# Assuming you have the 'X' and 'y' data for the specific probe
# You can also create a data frame if you have the original data
data_of_interest <- data.frame(
  X = model_of_interest$x,
  y = model_of_interest$y
)

# Calculate the correlation coefficient between 'X' and 'y'
cor_coefficient <- cor(data_of_interest$X, data_of_interest$y)

# Calculate the total number of points
total_points <- nrow(data_of_interest)

library(ggplot2)
# Plot the data points and the linear regression line
plot_lm <- ggplot(data_of_interest, aes(x = X, y = y)) +
  geom_point(color = "blue", size = 1) +  # Scatter plot of data points
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Linear regression line
  geom_text(x = min(data_of_interest$X), y = max(data_of_interest$y),
            label = paste("Correlation:", round(cor_coefficient, 3)),
            hjust = 0, vjust = 1, color = "black", size = 4) +  # Add correlation coefficient
  geom_text(x = min(data_of_interest$X), y = max(data_of_interest$y) - 0.05,  # Adjust y value to place it below the correlation
            label = paste("Total numbers of points:", total_points),
            hjust = 0, vjust = 1, color = "black", size = 4) +  # Add total number of points
  labs(title = paste("Linear Model for Probe:", probe_id_of_interest),
       x = "Age",
       y = "Beta Value")

# Show the plot
print(plot_lm)



# Coefficients
coefficients(model_probe_cg16867657)

# R-squared value
summary(model_probe_cg16867657)
