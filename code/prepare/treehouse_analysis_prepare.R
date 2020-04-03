#analyze ACE2 from the treehouse project
setwd("~/Documents/stanford/wars/gender/data")

library(ggplot2)
library(psych)
library(beanplot)
library(cowplot)
library(stringr)

treehouse_meta = read.csv("treehouse_all_meta.csv", stringsAsFactors = F)
treehouse_meta = subset(treehouse_meta, !is.na(age_in_years) & gender %in% c("female", "male"))
treehouse_meta$sample_id = str_replace_all(treehouse_meta$sample_id, "-", ".")
treehouse_meta$age = NA
treehouse_meta$age[treehouse_meta$age_in_years %in% c("0-1", "1-9", "10-19")] = "0-19"
treehouse_meta$age[treehouse_meta$age_in_years %in% c("20-29", "30-39", "40-49","50-59")] = "20-59"
treehouse_meta$age[treehouse_meta$age_in_years %in% c( "60-69", "70-79", "80-100")] = "60-100"

load("output/octad//id_tissue_prediction_final.RData")
id_tissue_prediction$tissue = as.numeric(as.character(id_tissue_prediction$tissue))
tissues = names(id_tissue_prediction)[-c(1,2,3)]
id_tissue_prediction$tissue = sapply(1:nrow(id_tissue_prediction), function(i){
  #keep original
  t = NA
  if (is.na(id_tissue_prediction$tissue[i])){
    #only if prob greater than a threshold
    t_temp = tissues[id_tissue_prediction$pred[i] + 1]
    if (id_tissue_prediction[i, t_temp] > 0.7) t = t_temp
  }else{
    t = tissues[id_tissue_prediction$tissue[i] + 1]
  }
  t
})
id_tissue_prediction = id_tissue_prediction[, c("id", "tissue")]
id_tissue_prediction$id = str_replace_all(id_tissue_prediction$id, "-", ".")

treehouse_meta = merge(treehouse_meta, id_tissue_prediction, by.x = "sample_id", by.y = "id")
treehouse_meta$tissue = treehouse_meta$tissue.y

load("raw/TcgaTargetGtexMet_rsem_gene_tpm.RData")
#colnames(all_exprs) = str_replace_all(colnames(all_exprs), "\\.", "-")
'
gene1 = "ENSG00000130234" #ACE2 59272
gene3 = "ENSG00000184012"  #TMPRSS2 #ENSG00000116678  ENSG00000174697 
gene2 = "ENSG00000091831" #ESR1 2099 
gene1 = "ENSG00000197635" ##DPP-4  ENSG00000197635 1803
gene1 = "ENSG00000159640" #ACE1 1636 ENSG00000159640
gene2 = "ENSG00000164850"  #GPER1 2852
gene2 = "ENSG00000140009" #ESR2, 2100 
gene2 = "ENSG00000169083" #AR 367
gene2 = "ENSG00000082175" #5241 PR
gene2 = "ENSG00000137869" #1588 Identifiers Aliases	CYP19A1 estrogen synthetase
'
Y_genes = c("ENSG00000067048", "ENSG00000198692", "ENSG00000012817", "ENSG00000176728", "ENSG00000129824", "ENSG00000114374")
Y_expr = apply(all_exprs[Y_genes, ], 2, median)

genes = c("ENSG00000130234", "ENSG00000091831", "ENSG00000197635", "ENSG00000159640", "ENSG00000164850", "ENSG00000140009", "ENSG00000169083", "ENSG00000082175", "ENSG00000137869", "ENSG00000184012")
symbols = c("ACE2", "ESR1", "DPP4", "ACE1", "GPER1", "ESR2", "AR", "PGR", "CYP19A1", "TMPRSS2")

expr = t(all_exprs[genes, ])
colnames(expr) = symbols

expr_meta = merge(treehouse_meta, cbind(expr,Y_expr) , by.x = "sample_id", by.y = 0)
expr_meta$tissue = expr_meta$tissue.y


paste("# samples", nrow(expr_meta))
paste("% predicted gender", sum(table(expr_meta$gender))/nrow(expr_meta))
paste("% labeled tissue", sum(table(expr_meta$tissue_old))/nrow(expr_meta))
paste("% predicted tissue", sum(table(expr_meta$tissue))/nrow(expr_meta))
paste("% predicted age", sum(table(expr_meta$age_group))/nrow(expr_meta))


#expr_meta = expr_meta[expr_meta$ENSG00000130234 > 0.1, ]

paste("# samples", nrow(expr_meta))

write.csv(expr_meta, paste0("treehouse_expr_meta.csv"))
###############

