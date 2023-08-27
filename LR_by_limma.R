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
