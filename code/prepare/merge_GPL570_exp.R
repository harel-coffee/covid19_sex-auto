setwd("~/Documents/stanford/wars/gender//data")

library(GEOquery)

gpl <- getGEO("GPL570", destdir=".")

mapping = (Table(gpl))

symbols = c("ACE2", "ESR1", "DPP4",  "GPER1", "ESR2", "AR", "PGR", "TMPRSS2")
 
load("Mergedv2.Rdata")

GSMs = sapply(colnames(Mergedv2)[-1], function(x) unlist(strsplit(x, "_"))[1])

expr = NULL
for (symbol in symbols){
     ids = mapping[mapping$`Gene Symbol` == symbol, "ID"]    
     exp_subset = apply(Mergedv2[Mergedv2$X %in% ids, -1], 2, median)
     expr = rbind(expr, exp_subset)
}

rownames(expr) = symbols
colnames(expr) = GSMs
expr_subset = t(expr)[, "TMPRSS2"]
#write.csv(t(expr), "~/Documents/stanford/wars/gender/data/GPL570_exp.csv")

expr_old = read.csv("GPL570_exp_v0.csv", row.names = 1)
expr_new = merge(expr_old, expr_subset, by.x = 0, by.y = 0)
rownames(expr_new) = expr_new$X
expr_new$TMPRSS2 = expr_new$y
expr_new = expr_new[, -c(which(colnames(expr_new) %in% c("Row.names", "y")))]
write.csv(expr_new, "GPL570_exp.csv")
