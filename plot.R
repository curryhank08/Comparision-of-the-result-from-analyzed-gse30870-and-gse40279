### author: Yao


## manhattan plot
library(qqman)
result_30870_sub$CHR <- as.numeric(result_30870_sub$CHR)
# function from qqman to plot manhattan 
manhattan(result_30870_sub, 
          main = "Manhattan Plot for gse30870",
          cex = 0.3,
          ylim = c(0, 50),
          chr="CHR", 
          bp="RANGE_START", 
          snp= "ID", 
          p="P.Value",
          genomewideline = FALSE,
          suggestiveline = -log10(1e-05))

result_limma_x1$CHR <- as.numeric(result_limma_x1$CHR)
manhattan(result_limma_x1, 
          main = "Manhattan Plot for gse40279-p.value_x1",
          cex = 0.3,
          ylim = c(0, 120),
          chr="CHR", 
          bp="RANGE_START", 
          snp= "ID", 
          p="P.Value",
          genomewideline = FALSE,
          suggestiveline = -log10(1e-05))

result_limma_x1x2$CHR <- as.numeric(result_limma_x1x2$CHR)
manhattan(result_limma_x1x2, 
          main = "Manhattan Plot for gse40279-p.value_x1x2",
          cex = 0.3,
          ylim = c(0, 10),
          chr="CHR", 
          bp="RANGE_START", 
          snp= "ID", 
          p="P.Value",
          genomewideline = FALSE,
          suggestiveline = -log10(1e-05))
## histogram
library(ggplot2)
ggplot(result_30870_sub, aes(x = P.Value)) +
  geom_histogram(binwidth = 0.02) +
  ggtitle("Histogram plot for P.Value (gse30870)")

ggplot(result_30870_sub, aes(x = P.Value)) +
  geom_histogram() +
  ggtitle("Distribution of P.Value (gse30870)") +
  scale_y_continuous(
    breaks = seq(0, 50000, by = 10000),     # Specify the desired breaks for y-axis
    labels = seq(0, 50000, by = 10000),     # Specify the corresponding labels
    limits = c(0, 50000),
    expand = c(0, 0)                  # Adjust the expand parameter for y-axis
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = seq(0, 1, by = 0.1),
    limits = c(0, 1),  # Set limits to match xlim
    expand = c(0, 0)      # Adjust the expand parameter to remove extra space
  ) 

ggplot(result_30870_sub, aes(x = P.Value)) +
  geom_histogram() +
  ggtitle("Distribution of P.Value (gse30870)") +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = seq(0, 1, by = 0.1),
    limits = c(0, 1),  # Set limits to match xlim
    expand = c(0, 0)      # Adjust the expand parameter to remove extra space
  ) 

# Create a histogram using base R plotting functions
hist(result_30870_sub$P.Value, 
     breaks = seq(0, 1, by = 0.01),  # Specify the bin breaks
     col = "skyblue",                # Set color of bars
     main = "Histogram plot for gse30870-P.Value",  # Set main title
     xlab = "P.Value",               # Label for x-axis
     ylab = "Count")             # Label for y-axis

# Create a histogram using base R plotting functions
hist(result_limma_x1$P.Value, 
     breaks = seq(0, 1, by = 0.01),  # Specify the bin breaks
     col = "skyblue",                # Set color of bars
     main = "Histogram plot for gse40279-P.Value_x1",  # Set main title
     xlab = "P.Value",               # Label for x-axis
     ylab = "Count",
     ylim = c(0, 140000))             # Label for y-axis

# Create a histogram using base R plotting functions
hist(result_limma_x1x2$P.Value, 
     breaks = seq(0, 1, by = 0.01),  # Specify the bin breaks
     col = "skyblue",                # Set color of bars
     main = "Histogram plot for gse40279-P.Value_x1x2",  # Set main title
     xlab = "P.Value",               # Label for x-axis
     ylab = "Count") 
