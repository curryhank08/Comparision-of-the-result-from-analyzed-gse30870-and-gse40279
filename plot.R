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


## Venndiagram
# Install ggVennDiagram to plot venndiagram
install.packages("ggVennDiagram") # Optional

" Wrong cuz ggVennDiagram() needs data with type of list as x input
library(ggVennDiagram)
# Assuming you have the adjusted p-values for each probe in each dataset
adjusted_p_values_gse30870sub <- result_30870_sub$adj.P.Val
adjusted_p_values_gse40279_x1 <- result_limma_x1$adj.P.Val
adjusted_p_values_gse40279_x1x2 <- result_limma_x1x2$adj.P.Val

# Specify the condition
condition_p <- 0.05

# Create logical vectors indicating if the adjusted p-values meet the condition
logical_vector_gse30870 <- as.integer(adjusted_p_values_gse30870sub < condition_p)
logical_vector_gse40279_x1 <- as.integer(adjusted_p_values_gse40279_x1 < condition_p)
logical_vector_gse40279_x1x2 <- as.integer(adjusted_p_values_gse40279_x1x2 < condition_p)

# data frame of the values for the probes' adj.pvalue under condition
data <- data.frame(
  `gse30870` = logical_vector_gse30870,
  `gse40279_x1` = logical_vector_gse40279_x1,
  `gse40279_x1x2` = logical_vector_gse40279_x1x2
)
row.names(data) <- probe_ids

# list of the values for the probes' adj.pvalue under condition
data_for_venn <- list(gse30870_adjp = logical_vector_gse30870,
                      gse40279_adjp_x1 = logical_vector_gse40279_x1,
                      gse40279_adjp_x1x2 = logical_vector_gse40279_x1x2)

ggVennDiagram(data_for_venn)
"

# Load the necessary library
library(VennDiagram)
library(ggplot2)

# Create your gene sets based on the condition (adj.pvalue < 0.01(Optional))
set1 <- (result_30870_sub[result_30870_sub$adj.P.Val < 1e-8,])$ID
set2 <- (result_limma_x1[result_limma_x1$adj.P.Val < 1e-8,])$ID
set3 <- (result_limma_x1x2[result_limma_x1x2$adj.P.Val < 1e-8,])$ID

# Create a list of sets
sets_list <- list('30870' = set1, '40279-x1' = set2, '40279-x1x2' = set3)

# Plot the Venn diagram
venn1 <- ggVennDiagram(sets_list)
venn1 + 
  labs(title = "Venn Diagram",
       subtitle = "Note: adj.P.value < 1e-8",
       caption = paste(Sys.Date(), "Yao")) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  theme(
    plot.title = element_text(hjust = 0.1, vjust = 10, size = 20, face = "bold"),  # Adjust title positioning
    plot.subtitle = element_text(hjust = 0.1, vjust = 15),  # Adjust subtitle positioning
    plot.margin = margin(t = 0, r = 10, b = 0, l = 0, unit = "pt")  # Adjust the bottom margin
  )

theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
