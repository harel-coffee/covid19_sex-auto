####
#post processing
setwd("~/Documents/stanford/wars/gender/data")

library("ggplot2")
library("beanplot")
library("stringr")

#load gender
load("output/archs//id_gender_prediction_final.RData")
#keep the labeled gender
id_gender_prediction$gender = sapply(1:nrow(id_gender_prediction), function(i){
  g = id_gender_prediction$gender[i]
  if (is.na(g)){
    if (id_gender_prediction$male[i] > 0.8){
      g = "male"
    }
    if (id_gender_prediction$female[i] > 0.8){
      g = "female"
    }
  }else{
    if (g == 1){
      g = "male"
    }else if (g == 0){
      g = "female"
    }
  }
  g
})
id_gender_prediction = id_gender_prediction[, c("id", "gender")]

load("output/archs//id_tissue_prediction_final.RData")
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

#age
load("output/archs//id_age_prediction_final.RData")
data_archs$age_group = sapply(1:nrow(data_archs), function(i){
  #keep original
  a = NA
  if (is.na(data_archs$age[i])){
    #only if prob greater than a threshold
    if (max(data_archs[i, c("group1", "group2", "group3")]) > 0.5){
      a = data_archs$pred[i]
    }
  }else{
    a = data_archs$age[i]
  }
  
  if (is.na(a)){
    NA
  }else if (a == 0){
    "0-19"
  }else if (a == 1){
    "20-59"
  }else if (a == 2){
    "60-100"
  }else{
    NA
  }
})

data_archs = (data_archs[, c("id", "age_group")])
names(data_archs) = c("id", "age")
data_archs = unique(data.frame(data_archs))
#using counts
#expr_meta1 = read.csv("~/Downloads/human_matrix_expr.csv", stringsAsFactors = F, row.names = 1)

#using tpm data
#expr_meta2 = read.csv("~/Downloads/sample_gender_expr_sum.csv", stringsAsFactors = F, row.names = 1)
#colnames(expr_meta) = str_replace_all(colnames(expr_meta), "_expr", "")
#rownames(expr_meta) = expr_meta$samples.sample_ids.

#using tpm data from KE
load("ARCHS4.gene.tpm.matrix.RData")
expr_meta3 = t(ARCHS4.log2.gene.tpm.matrix)
other_meta = unique(read.csv("GEO_gsm_all.csv", stringsAsFactors = F)[, c("gsm", "age", "gender", "tissue", "source_name", "characteristics", "title")])

expr_meta = merge(expr_meta3, other_meta, by.x = 0, by.y = "gsm")
expr_meta = merge(expr_meta, id_gender_prediction, by.x = "Row.names", by.y = "id", all.x=T)
expr_meta = merge(expr_meta, id_tissue_prediction, by.x = "Row.names", by.y = "id", all.x = T)
expr_meta = merge(expr_meta, data_archs, by.x = "Row.names", by.y = "id", all.x = T)
expr_meta$gender = expr_meta$gender.y
expr_meta$tissue = expr_meta$tissue.y #using predicted one
expr_meta$age = expr_meta$age.y #using predicted one

paste("# samples", nrow(expr_meta))
paste("% labeled gender", sum(table(expr_meta$gender.x))/nrow(expr_meta))
paste("% predicted gender", sum(table(expr_meta$gender))/nrow(expr_meta))
paste("% labeled tissue", sum(table(expr_meta$tissue.x))/nrow(expr_meta))
paste("% predicted tissue", sum(table(expr_meta$tissue))/nrow(expr_meta))
paste("% labeled age", sum(table(expr_meta$age.x))/nrow(expr_meta))
paste("% predicted age", sum(table(expr_meta$age))/nrow(expr_meta))


#which we should removed those samples with ACE2 == 0? should be OK, as long as we use the same standard for female and male

#expr_meta = expr_meta[!is.na(expr_meta$ACE2) & expr_meta$ACE2 > 0, ] #& expr_meta$ACE2 > 0
paste("# samples", nrow(expr_meta))

hist(log(expr_meta$ACE2 + 0.01))

t.test(ACE2 ~ gender, data = expr_meta)

symbols = c("ACE2", "ESR1", "DPP4",  "GPER1", "ESR2", "AR", "PGR", "TMPRSS2")

write.csv(expr_meta, "ARCHS_expr_meta.csv")
###############
