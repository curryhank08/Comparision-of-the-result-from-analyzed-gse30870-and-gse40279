library(devtools)
install_github('GEOquery','seandavi')
library(remotes)
install_github("seandavi/GEOquery", force = TRUE)


# Load modified GEOquery
library(GEOquery)
# Setting the max timeout_seconds
options(timeout=1000000)
# Check the input timeout_seconds
getOption("timeout")

# Download GSE40279 by the fuction getGEO() from modified GEOquery package.
gse56581 <- getGEO("GSE56581", GSEMatrix = FALSE, AnnotGPL = TRUE, destdir = '/Users/yaohank/Library/CloudStorage/OneDrive-國立成功大學NationalChengKungUniversity/R/Comparision_7_18')
gse56581_matrix <- gse56581[[1]]
gse56581_data <- exprs(gse56581_matrix)

# Download GSE40279 by the fuction getGEO() from modified GEOquery package.
gse42861 <- getGEO("GSE42861", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = '/Users/yaohank/Library/CloudStorage/OneDrive-國立成功大學NationalChengKungUniversity/R/Comparision_7_18')
gse42861_matrix <- gse42861[[1]]
gse42861_data <- exprs(gse42861_matrix)
