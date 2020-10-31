####
#post processing
setwd("~/Documents/stanford/wars/gender/data")

library("ggplot2")
library("beanplot")

#load gender
load("output/geo/id_gender_prediction_final.RData")
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

#tissue
load("output/geo/id_tissue_prediction_final.RData")
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
load("output/geo/id_age_prediction_final.RData")
data_GEO$age_group = sapply(1:nrow(data_GEO), function(i){
  #keep original
  a = NA
  if (is.na(data_GEO$age[i])){
    #only if prob greater than a threshold
    if (max(data_GEO[i, c("group1", "group2", "group3")]) > 0.5){
      a = data_GEO$pred[i]
    }
  }else{
    a = data_GEO$age[i]
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

data_GEO = data_GEO[, c("id", "age_group")]
names(data_GEO) = c("id", "age")

meta = read.csv("GPL570_meta_all.csv", stringsAsFactors = F)
expr = read.csv("GPL570_exp.csv", stringsAsFactors = F)
#duplicated GSM
expr = expr[!duplicated(expr$X),]
expr$gsm = sapply(expr$X, function(x) unlist(strsplit(x, "\\."))[1])


expr_meta = merge( meta,expr, by.x= "gsm", by.y = "gsm")
expr_meta = merge(expr_meta, id_gender_prediction, by.x= "gsm", by.y = "id", all.x=T)
expr_meta = merge(expr_meta, id_tissue_prediction, by.x= "gsm", by.y = "id", all.x= T)
expr_meta = merge(expr_meta, data_GEO, by.x= "gsm", by.y = "id", all.x= T)
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


hist((expr_meta$ACE2 ))

t.test(ACE2 ~ gender, data = expr_meta)
t.test(ACE2 ~ smoke, data = expr_meta)
anova(lm(ACE2 ~ smoke + gender, data = expr_meta))
anova(lm(ACE2 ~  gender, data = expr_meta))
anova(lm(ACE2 ~  age, data = expr_meta))

symbols = c("ACE2", "ESR1", "DPP4",  "GPER1", "ESR2", "AR", "PGR", "TMPRSS2")

write.csv(expr_meta, "GPL570_expr_meta.csv")
###############
