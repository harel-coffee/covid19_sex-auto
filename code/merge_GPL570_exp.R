setwd("~/Documents/stanford/wars/cmap/data")

library(GEOquery)

gpl <- getGEO("GPL570", destdir=".")

mapping = (Table(gpl))

symbols = c("ACE2", "ESR1", "DPP4",  "GPER1", "ESR2", "AR", "PGR")
 
load("~/Downloads/Merged42.Rdata")

GSMs = sapply(colnames(Merged42)[-1], function(x) unlist(strsplit(x, "_"))[1])

expr = NULL
for (symbol in symbols){
     ids = mapping[mapping$`Gene Symbol` == symbol, "ID"]    
     exp_subset = apply(Merged42[Merged42$X %in% ids, -1], 2, median)
     expr = rbind(expr, exp_subset)
}

rownames(expr) = symbols
colnames(expr) = GSMs

write.csv(t(expr), "~/Documents/stanford/wars/gender/data/GPL570_exp.csv")
