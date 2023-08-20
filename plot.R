### author: Yao

library(ggplot2)
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
ggplot(result_30870_sub)



