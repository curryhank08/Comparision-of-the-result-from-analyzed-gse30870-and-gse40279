## author: haku yao
library(limma)

gset <- gse40279_matrix

'
if (length(gset) > 1) idx <- grep("GPL13534", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
'

'
# Create age categories
age <- pData(gset)$characteristics_ch1
# Remove "age (y):" and convert to numeric
age <- sub("^\\s*age \\(y\\): ", "", age)
age <- as.numeric(age)
# Assign age values to a new column in pData of gset
pData(gset)$age <- age

# Define age categories based on specific age ranges
age_categories <- cut(age,
                      breaks = c(0, 30, 65, Inf),
                      labels = c("Young", "Middle", "Old"),
                      include.lowest = TRUE)

# Assign age categories to the pData of gset
pData(gset)$age_category <- age_categories
'

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex)}

## By limma (Not yet solved)
x1 <- gse40279_matrix$age
x2 <- gse40279_matrix$`gender:ch1`
#f <- factor(conditions)
design <- model.matrix(~ x1+x2+x1:x2)
colnames(design)[colnames(design) == "x2M"] <- "x2"
colnames(design)[colnames(design) == "x1:x2M"] <- "x1:x2"
fit <- lmFit(gse40279_matrix, design)
fit <- eBayes(fit)
result_limma_x1 <- topTable(fit, coef="x1", number = Inf, adjust.method = "BH", sort.by = "P")
result_limma_x1x2 <- topTable(fit, coef="x1:x2", number = Inf, adjust.method = "BH", sort.by = "P")

# result_limma_2 <- result_limma[, c("UCSC_RefGene_Name", "P.Value", "logFC")]
# Outcome of each hypothesis test
results <- decideTests(fit, method = "separate", adjust.method = "BH", p.value = 1e-5)

# Showing numbers of genes significant in each comparison
vennDiagram(results, show.include = FALSE) + title('Numbers of genes significant (p.value<1e-5) in each comparison')

# Install ggVennDiagram to plot venndiagram
install.packages("ggVennDiagram") # Optional

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
