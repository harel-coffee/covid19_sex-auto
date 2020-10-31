# setwd("~/Documents/stanford/wars/gender/data")
rm(list=ls())
setwd('/Users/msun/Documents/covid-19/main')

####
#under expression of ACE2/TMPRSS2 between gender groups
#we collected data from three datasets and grouped them into high and normal group
####

#top * expression (test 80%, 85%, 90%, 95%)
p = 0.9

###########
#compile data
ARCHS_expr_meta = read.csv("raw/ARCHS_expr_meta_final_v8.csv")
ARCHS_expr_meta = ARCHS_expr_meta[!is.na(ARCHS_expr_meta$ACE2) & !is.na(ARCHS_expr_meta$TMPRSS2), ]
ARCHS_expr_meta = ARCHS_expr_meta[ARCHS_expr_meta$gpl %in% c("GPL11154"), ] #seems GPL11154 give quite robust signal
ACE_cutoff =  sort(ARCHS_expr_meta$ACE2)[length((ARCHS_expr_meta$ACE2)) * p] #quantile(ARCHS_expr_meta$ACE2)[4]
TMPRSS2_cutoff = sort(ARCHS_expr_meta$TMPRSS2)[length((ARCHS_expr_meta$TMPRSS2)) * p] #quantile(ARCHS_expr_meta$TMPRSS2)[4]
ESR1_cutoff = sort(ARCHS_expr_meta$ESR1)[length((ARCHS_expr_meta$ESR1)) * p] #quantile(ARCHS_expr_meta$ESR1)[4]
ESR2_cutoff = sort(ARCHS_expr_meta$ESR2)[length((ARCHS_expr_meta$ESR2)) * p] #quantile(ARCHS_expr_meta$ESR2)[4]
AR_cutoff = sort(ARCHS_expr_meta$AR)[length((ARCHS_expr_meta$AR)) * p] #quantile(ARCHS_expr_meta$AR)[4]
PGR_cutoff = sort(ARCHS_expr_meta$PGR)[length((ARCHS_expr_meta$PGR)) * p] #quantile(ARCHS_expr_meta$PGR)[4]
# DPP4_cutoff = sort(ARCHS_expr_meta$DPP4)[length((ARCHS_expr_meta$DPP4)) * p] #quantile(ARCHS_expr_meta$DPP4)[4]

ARCHS_expr_meta$ACE2_TMPRSS2 = sapply(1:nrow(ARCHS_expr_meta), function(i) {
  if (ARCHS_expr_meta$ACE2[i] >  ACE_cutoff & ARCHS_expr_meta$TMPRSS2[i] > TMPRSS2_cutoff){
    1
  }else{
    0
  }
})
ARCHS_expr_meta$ACE2_bin = sapply(1:nrow(ARCHS_expr_meta), function(i) {
  if (ARCHS_expr_meta$ACE2[i] >  ACE_cutoff){
    1
  }else{
    0
  }
})
ARCHS_expr_meta$TMPRSS2_bin = sapply(1:nrow(ARCHS_expr_meta), function(i) {
  if (ARCHS_expr_meta$TMPRSS2[i] > TMPRSS2_cutoff){
    1
  }else{
    0
  }
})

ARCHS_expr_meta$ESR1_bin = sapply(1:nrow(ARCHS_expr_meta), function(i) {
  if (ARCHS_expr_meta$ESR1[i] > ESR1_cutoff){
    1
  }else{
    0
  }
})

ARCHS_expr_meta$ESR2_bin = sapply(1:nrow(ARCHS_expr_meta), function(i) {
  if (ARCHS_expr_meta$ESR2[i] > ESR2_cutoff){
    1
  }else{
    0
  }
})

ARCHS_expr_meta$PGR_bin = sapply(1:nrow(ARCHS_expr_meta), function(i) {
  if (ARCHS_expr_meta$PGR[i] > PGR_cutoff){
    1
  }else{
    0
  }
})

ARCHS_expr_meta$AR_bin = sapply(1:nrow(ARCHS_expr_meta), function(i) {
  if (ARCHS_expr_meta$AR[i] > AR_cutoff){
    1
  }else{
    0
  }
})

# ARCHS_expr_meta$DPP4_bin = sapply(1:nrow(ARCHS_expr_meta), function(i) {
#   if (ARCHS_expr_meta$DPP4[i] > DPP4_cutoff){
#     1
#   }else{
#     0
#   }
# })

ARCHS_expr_meta$sample_id = ARCHS_expr_meta$Row.names
ARCHS_expr_meta$source = "ARCHS"
ARCHS_expr_meta$type = NA
table(ARCHS_expr_meta$ACE2_TMPRSS2)/sum(table(ARCHS_expr_meta$ACE2_TMPRSS2) )

GPL570_expr_meta = read.csv("raw/GPL570_expr_meta_final_v8.csv")
GPL570_expr_meta = GPL570_expr_meta[!is.na(GPL570_expr_meta$ACE2) & !is.na(GPL570_expr_meta$TMPRSS2), ]
ACE_cutoff =  sort(GPL570_expr_meta$ACE2)[length((GPL570_expr_meta$ACE2)) * p] #quantile(GPL570_expr_meta$ACE2)[4]
TMPRSS2_cutoff = sort(GPL570_expr_meta$TMPRSS2)[length((GPL570_expr_meta$TMPRSS2)) * p] #quantile(GPL570_expr_meta$TMPRSS2)[4]
ESR1_cutoff = sort(GPL570_expr_meta$ESR1)[length((GPL570_expr_meta$ESR1)) * p] #quantile(GPL570_expr_meta$ESR1)[4]
ESR2_cutoff = sort(GPL570_expr_meta$ESR2)[length((GPL570_expr_meta$ESR2)) * p] #quantile(GPL570_expr_meta$ESR2)[4]
PGR_cutoff = sort(GPL570_expr_meta$PGR)[length((GPL570_expr_meta$PGR)) * p] #quantile(GPL570_expr_meta$PGR)[4]
AR_cutoff = sort(GPL570_expr_meta$AR)[length((GPL570_expr_meta$AR)) * p] #quantile(GPL570_expr_meta$AR)[4]
# DPP4_cutoff = sort(GPL570_expr_meta$DPP4)[length((GPL570_expr_meta$DPP4)) * p] #quantile(GPL570_expr_meta$DPP4)[4]

GPL570_expr_meta$ACE2_TMPRSS2 = sapply(1:nrow(GPL570_expr_meta), function(i) {
  if (GPL570_expr_meta$ACE2[i] > ACE_cutoff & GPL570_expr_meta$TMPRSS2[i] > TMPRSS2_cutoff){
    1
  }else{
    0
  }
})
GPL570_expr_meta$ACE2_bin = sapply(1:nrow(GPL570_expr_meta), function(i) {
  if (GPL570_expr_meta$ACE2[i] > ACE_cutoff){
    1
  }else{
    0
  }
})
GPL570_expr_meta$TMPRSS2_bin = sapply(1:nrow(GPL570_expr_meta), function(i) {
  if (GPL570_expr_meta$TMPRSS2[i] > TMPRSS2_cutoff){
    1
  }else{
    0
  }
})

GPL570_expr_meta$ESR1_bin = sapply(1:nrow(GPL570_expr_meta), function(i) {
  if (GPL570_expr_meta$ESR1[i] > ESR1_cutoff){
    1
  }else{
    0
  }
})

GPL570_expr_meta$ESR2_bin = sapply(1:nrow(GPL570_expr_meta), function(i) {
  if (GPL570_expr_meta$ESR2[i] > ESR2_cutoff){
    1
  }else{
    0
  }
})

GPL570_expr_meta$PGR_bin = sapply(1:nrow(GPL570_expr_meta), function(i) {
  if (GPL570_expr_meta$PGR[i] > PGR_cutoff){
    1
  }else{
    0
  }
})

GPL570_expr_meta$AR_bin = sapply(1:nrow(GPL570_expr_meta), function(i) {
  if (GPL570_expr_meta$AR[i] > AR_cutoff){
    1
  }else{
    0
  }
})

# GPL570_expr_meta$DPP4_bin = sapply(1:nrow(GPL570_expr_meta), function(i) {
#   if (GPL570_expr_meta$DPP4[i] > DPP4_cutoff){
#     1
#   }else{
#     0
#   }
# })

GPL570_expr_meta$sample_id = GPL570_expr_meta$gsm
GPL570_expr_meta$source = "GPL570"
GPL570_expr_meta$type = NA
table(GPL570_expr_meta$ACE2_TMPRSS2)/sum(table(GPL570_expr_meta$ACE2_TMPRSS2) )

treehouse_expr_meta = read.csv("raw/treehouse_expr_meta_final_v8.csv")
treehouse_expr_meta$smoke = NA
#treehouse_expr_meta$tissue = treehouse_expr_meta$tissue.x
#treehouse_expr_meta = treehouse_expr_meta[treehouse_expr_meta$source == "gtex",]
treehouse_expr_meta = treehouse_expr_meta[!is.na(treehouse_expr_meta$ACE2) & !is.na(treehouse_expr_meta$TMPRSS2), ]
ACE_cutoff = sort(treehouse_expr_meta$ACE2)[length((treehouse_expr_meta$ACE2)) * p] # quantile(treehouse_expr_meta$ACE2)[4]
TMPRSS2_cutoff =  sort(treehouse_expr_meta$TMPRSS2)[length((treehouse_expr_meta$TMPRSS2)) * p] # quantile(treehouse_expr_meta$TMPRSS2)[4]
ESR1_cutoff =  sort(treehouse_expr_meta$ESR1)[length((treehouse_expr_meta$ESR1)) * p] # quantile(treehouse_expr_meta$ESR1)[4]
ESR2_cutoff =  sort(treehouse_expr_meta$ESR2)[length((treehouse_expr_meta$ESR2)) * p] # quantile(treehouse_expr_meta$ESR2)[4]
PGR_cutoff =  sort(treehouse_expr_meta$PGR)[length((treehouse_expr_meta$PGR)) * p] # quantile(treehouse_expr_meta$PGR)[4]
AR_cutoff =  sort(treehouse_expr_meta$AR)[length((treehouse_expr_meta$AR)) * p] # quantile(treehouse_expr_meta$AR)[4]
#DPP4_cutoff =  sort(treehouse_expr_meta$DPP4)[length((treehouse_expr_meta$DPP4)) * p] # quantile(treehouse_expr_meta$DPP4)[4]

treehouse_expr_meta$ACE2_TMPRSS2 = sapply(1:nrow(treehouse_expr_meta), function(i) {
  if (treehouse_expr_meta$ACE2[i] > ACE_cutoff & treehouse_expr_meta$TMPRSS2[i] > TMPRSS2_cutoff){
    1
  }else{
    0
  }
})
treehouse_expr_meta$ACE2_bin = sapply(1:nrow(treehouse_expr_meta), function(i) {
  if (treehouse_expr_meta$ACE2[i] > ACE_cutoff ){
    1
  }else{
    0
  }
})
treehouse_expr_meta$TMPRSS2_bin = sapply(1:nrow(treehouse_expr_meta), function(i) {
  if (treehouse_expr_meta$TMPRSS2[i] > TMPRSS2_cutoff){
    1
  }else{
    0
  }
})

treehouse_expr_meta$ESR1_bin = sapply(1:nrow(treehouse_expr_meta), function(i) {
  if (treehouse_expr_meta$ESR1[i] > ESR1_cutoff){
    1
  }else{
    0
  }
})

treehouse_expr_meta$ESR2_bin = sapply(1:nrow(treehouse_expr_meta), function(i) {
  if (treehouse_expr_meta$ESR2[i] > ESR2_cutoff){
    1
  }else{
    0
  }
})

treehouse_expr_meta$PGR_bin = sapply(1:nrow(treehouse_expr_meta), function(i) {
  if (treehouse_expr_meta$PGR[i] > PGR_cutoff){
    1
  }else{
    0
  }
})

treehouse_expr_meta$AR_bin = sapply(1:nrow(treehouse_expr_meta), function(i) {
  if (treehouse_expr_meta$AR[i] > AR_cutoff){
    1
  }else{
    0
  }
})

# treehouse_expr_meta$DPP4_bin = sapply(1:nrow(treehouse_expr_meta), function(i) {
#   if (treehouse_expr_meta$DPP4[i] > DPP4_cutoff){
#     1
#   }else{
#     0
#   }
# })

table(treehouse_expr_meta$ACE2_TMPRSS2)/sum(table(treehouse_expr_meta$ACE2_TMPRSS2) )
treehouse_expr_meta$source = "treehouse"

expr_meta = rbind(ARCHS_expr_meta[, c("source", "sample_id", "age", "tissue", "gender", "smoke", "type", "ACE2_TMPRSS2", "ACE2_bin", "TMPRSS2_bin", "ESR1_bin", "ESR2_bin", "PGR_bin", "AR_bin")], #, "DPP4_bin"
                  GPL570_expr_meta[, c("source","sample_id", "age", "tissue", "gender","smoke",  "type","ACE2_TMPRSS2", "ACE2_bin", "TMPRSS2_bin", "ESR1_bin", "ESR2_bin", "PGR_bin", "AR_bin")], # , "DPP4_bin"
                  treehouse_expr_meta[, c("source","sample_id", "age", "tissue", "gender", "smoke", "type","ACE2_TMPRSS2", "ACE2_bin", "TMPRSS2_bin", "ESR1_bin", "ESR2_bin", "PGR_bin", "AR_bin")]) # , "DPP4_bin"
expr_meta = expr_meta[expr_meta$source %in% c("GPL570", "treehouse", "ARCHS"),]
expr_meta$tissue = as.character(expr_meta$tissue)
expr_meta[is.na(expr_meta$tissue), "tissue"] = "unknown"

# remove gender specific tissue
#tmp = which(expr_meta$tissue=='breast' | expr_meta$tissue=='prostate' | expr_meta$tissue=='testis')
#expr_meta = expr_meta[-tmp,]

expr_meta$tissue = as.factor(expr_meta$tissue)

# only keep with gender 
# expr_meta = expr_meta[which(!is.na(expr_meta$gender)), ]

#####################
#basic analysis
gpl_dist = as.matrix(table(expr_meta$tissue[expr_meta$source == "GPL570"], expr_meta$gender[expr_meta$source == "GPL570"]))
gpl_dist[,1]/gpl_dist[,2]
A_dist = as.matrix(table(expr_meta$tissue[expr_meta$source == "ARCHS"], expr_meta$gender[expr_meta$source == "ARCHS"]))
A_dist[,1]/A_dist[,2]
T_dist = table(expr_meta$tissue[expr_meta$source == "treehouse"], expr_meta$gender[expr_meta$source == "treehouse"])
T_dist[,1]/T_dist[,2]

table(expr_meta$tissue[expr_meta$source == "GPL570"])/sum(table(expr_meta$tissue[expr_meta$source == "GPL570"]))
table(expr_meta$tissue[expr_meta$source == "ARCHS"])/sum(table(expr_meta$tissue[expr_meta$source == "ARCHS"]))
table(expr_meta$tissue[expr_meta$source == "treehouse"])/sum(table(expr_meta$tissue[expr_meta$source == "treehouse"]))


table(expr_meta$age)/sum(table(expr_meta$age))
table(expr_meta$tissue)/sum(table(expr_meta$tissue))
table(expr_meta$gender)/sum(table(expr_meta$gender))

sum(is.na(expr_meta$tissue))/nrow(expr_meta) #missing tissue


#####
#age difference
#old group has much higher ACE2/TMPRSS2 expression
####
############################
#compare young and old without/with  adjustment
expr_meta_subset = subset(expr_meta, age %in% c("0-19", "60-100"))
expr_meta_subset$age = as.factor(as.character(expr_meta_subset$age))
(fisher.test(table(expr_meta_subset$ACE2_TMPRSS2, expr_meta_subset$age)))
(fisher.test(table(expr_meta_subset$ACE2_bin, expr_meta_subset$age)))
(fisher.test(table(expr_meta_subset$TMPRSS2_bin, expr_meta_subset$age)))
mylogit = (glm(ACE2_TMPRSS2 ~ gender + tissue + age + source, expr_meta_subset, family = binomial()))
summary(mylogit)
exp(cbind(OR = coef(mylogit), confint(mylogit)))


#############
#table 1
#analyze normal samples from treehouse
#all age groups, adjust  age and tissue
#did not observe expression difference between female and male
result = NULL
expr_meta_subset = expr_meta[expr_meta$source %in% c("treehouse") & expr_meta$type %in% c("normal", "adjacent"), ]
dim(expr_meta_subset)

mylogit = (glm(ACE2_bin ~ gender + tissue + age, expr_meta_subset, family = binomial()))
summary(mylogit)
odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
print(odds)
result = rbind(result, data.frame(dataset = "healthy", gene = "ACE2", age_group = "all", OR = round(odds["gendermale", 1], 2), 
                                                CI_25 =  round(odds["gendermale", 2], 2) , 
                                                CI_75 =  round(odds["gendermale", 3], 2) , 
                                                p =  summary(mylogit)$coefficients["gendermale", 4]))

mylogit = (glm(TMPRSS2_bin ~ gender + tissue + age, expr_meta_subset, family = binomial()))
summary(mylogit)
odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
print(odds)
result = rbind(result, data.frame(dataset = "healthy", gene = "TMPRSS2", age_group = "all", OR = round(odds["gendermale", 1], 2), 
                                                CI_25 =  round(odds["gendermale", 2], 2) , 
                                                CI_75 =  round(odds["gendermale", 3], 2) , 
                                                p =  summary(mylogit)$coefficients["gendermale", 4]))

mylogit = (glm(ACE2_TMPRSS2 ~ gender + tissue + age, expr_meta_subset, family = binomial()))
summary(mylogit)
odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
print(odds)
result = rbind(result, data.frame(dataset = "healthy", gene = "ACE2 & TMPRSS2", age_group = "all", OR = round(odds["gendermale", 1], 2), 
                                                CI_25 =  round(odds["gendermale", 2], 2) , 
                                                CI_75 =  round(odds["gendermale", 3], 2) , 
                                                p =  summary(mylogit)$coefficients["gendermale", 4]))

#old group, adjust  tissue
expr_meta_subset = expr_meta[expr_meta$source %in% c("treehouse") & expr_meta$type %in% c("normal", "adjacent") & expr_meta$age %in% c("60-100"), ]
dim(expr_meta_subset)
mylogit = (glm(ACE2_bin ~ gender+ tissue, expr_meta_subset, family = binomial()))
summary(mylogit)
odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
print(odds)
result = rbind(result, data.frame(dataset = "healthy", gene = "ACE2", age_group = "old", OR = round(odds["gendermale", 1], 2), 
                                                CI_25 =  round(odds["gendermale", 2], 2) , 
                                                CI_75 =  round(odds["gendermale", 3], 2) , 
                                                p =  summary(mylogit)$coefficients["gendermale", 4]))

mylogit = (glm(TMPRSS2_bin ~ gender+ tissue, expr_meta_subset, family = binomial()))
summary(mylogit)
odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
print(odds)
result = rbind(result, data.frame(dataset = "healthy", gene = "TMPRSS2", age_group = "old", OR = round(odds["gendermale", 1], 2), 
                                                CI_25 =  round(odds["gendermale", 2], 2) , 
                                                CI_75 =  round(odds["gendermale", 3], 2) , 
                                                p =  summary(mylogit)$coefficients["gendermale", 4]))

mylogit = (glm(ACE2_TMPRSS2 ~ gender+ tissue, expr_meta_subset, family = binomial()))
summary(mylogit)
odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
print(odds)
result = rbind(result, data.frame(dataset = "healthy", gene = "ACE2 & TMPRSS2", age_group = "old", OR = round(odds["gendermale", 1], 2), 
                                                CI_25 =  round(odds["gendermale", 2], 2) , 
                                                CI_75 =  round(odds["gendermale", 3], 2) , 
                                                p =  summary(mylogit)$coefficients["gendermale", 4]))

################
#analysis all samples from merged dataset
###############
expr_meta_subset = expr_meta
dim(expr_meta_subset)

mylogit = (glm(ACE2_bin ~ gender + source + tissue + age, expr_meta_subset, family = binomial()))
summary(mylogit)
odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
print(odds)
result = rbind(result, data.frame(dataset = "All", gene = "ACE2", age_group = "all", OR = round(odds["gendermale", 1], 2), 
                                                CI_25 =  round(odds["gendermale", 2], 2) , 
                                                CI_75 =  round(odds["gendermale", 3], 2) , 
                                                p =  summary(mylogit)$coefficients["gendermale", 4]))

mylogit = (glm(TMPRSS2_bin ~ gender+ source+ tissue + age, expr_meta_subset, family = binomial()))
summary(mylogit)
odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
print(odds)
result = rbind(result, data.frame(dataset = "All", gene = "TMPRSS2", age_group = "all", OR = round(odds["gendermale", 1], 2), 
                                                CI_25 =  round(odds["gendermale", 2], 2) , 
                                                CI_75 =  round(odds["gendermale", 3], 2) , 
                                                p =  summary(mylogit)$coefficients["gendermale", 4]))

mylogit = (glm(ACE2_TMPRSS2 ~ gender+ source+ tissue + age, expr_meta_subset, family = binomial()))
summary(mylogit)
odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
print(odds)
result = rbind(result, data.frame(dataset = "All", gene = "ACE2 & TMPRSS2", age_group = "all", OR = round(odds["gendermale", 1], 2), 
                                                CI_25 =  round(odds["gendermale", 2], 2) , 
                                                CI_75 =  round(odds["gendermale", 3], 2) , 
                                                p =  summary(mylogit)$coefficients["gendermale", 4]))

#analyze old group
expr_meta_subset = expr_meta[expr_meta$age %in% c("60-100"), ]
dim(expr_meta_subset)
mylogit = (glm(ACE2_bin ~ gender+ source+ tissue, expr_meta_subset, family = binomial()))
summary(mylogit)
odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
print(odds)
result = rbind(result, data.frame(dataset = "All", gene = "ACE2", age_group = "old", OR = round(odds["gendermale", 1], 2), 
                                                CI_25 =  round(odds["gendermale", 2], 2) , 
                                                CI_75 =  round(odds["gendermale", 3], 2) , 
                                                p =  summary(mylogit)$coefficients["gendermale", 4]))

mylogit = (glm(TMPRSS2_bin ~ gender+ source+ tissue, expr_meta_subset, family = binomial()))
summary(mylogit)
odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
print(odds)
result = rbind(result, data.frame(dataset = "All", gene = "TMPRSS2", age_group = "old", OR = round(odds["gendermale", 1], 2), 
                                                CI_25 =  round(odds["gendermale", 2], 2) , 
                                                CI_75 =  round(odds["gendermale", 3], 2) , 
                                                p =  summary(mylogit)$coefficients["gendermale", 4]))

mylogit = (glm(ACE2_TMPRSS2 ~ gender+ source+ tissue, expr_meta_subset, family = binomial()))
summary(mylogit)
odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
print(odds)
result = rbind(result, data.frame(dataset = "All", gene = "ACE2 & TMPRSS2", age_group = "old", OR = round(odds["gendermale", 1], 2), 
                                                CI_25 =  round(odds["gendermale", 2], 2) , 
                                                CI_75 =  round(odds["gendermale", 3], 2) , 
                                                p =  summary(mylogit)$coefficients["gendermale", 4]))

fn = paste0('output/labeled/table1_p_', p, '.csv')
write.csv(result, fn)

######################
#table 2
#####################
tissues = c("blood","bone","brain","colon","heart","kidney","liver","lung","pancreas","skin","small intestine", "unknown")

results_all = NULL
for (t in tissues){
  print(t)
  expr_meta_subset = expr_meta[expr_meta$tissue %in% t,]
  mylogit = (glm(TMPRSS2_bin ~ gender+ source + age, expr_meta_subset, family = binomial()))
  
  odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
  
  results_all = rbind(results_all, data.frame(group = "all", tissue = t, OR = round(odds["gendermale", 1], 2), CI_25 =  round(odds["gendermale", 2],2) , 
                                      CI_75 =  round(odds["gendermale", 3],2) , p =  summary(mylogit)$coefficients["gendermale", 4]))
}
results_all$p_adj  = p.adjust(results_all$p)

#old group
results_old = NULL
for (t in tissues){
  print(t)
  expr_meta_subset = expr_meta[expr_meta$tissue %in% t & expr_meta$age %in% c("60-100"),]
  
  if (sum((table(expr_meta_subset[, "TMPRSS2_bin"], expr_meta_subset[, "gender"]) > 3)) != 4) next
    
  mylogit = (glm(TMPRSS2_bin ~ gender+ source , expr_meta_subset, family = binomial()))
  
  odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
  
  results_old = rbind(results_old, data.frame(group = "60-100", tissue = t, OR = round(odds["gendermale", 1], 2), CI_25 =  round(odds["gendermale", 2], 2) , 
                                      CI_75 =  round(odds["gendermale", 3],2) , p =  summary(mylogit)$coefficients["gendermale", 4]))
}
results_old$p_adj  = p.adjust(results_old$p)

results = rbind(results_all, results_old)

write.csv(results, "table2_TMPRSS2.csv")

#################
#table 3: ACE2 vs Hormone
#old group has negative correlation between hormone and ACE2, while the whole group has positive correlation between hormone and ACE2
#################

results = NULL
  for (hormone in c("ESR1_bin", "ESR2_bin", "PGR_bin", "AR_bin")){
    for (virus_gene in c("ACE2_bin", "TMPRSS2_bin", "ACE2_TMPRSS2")){
        expr_meta$hormone = expr_meta[, hormone]
        expr_meta$virus = expr_meta[, virus_gene]
        
          #all age 
          mylogit = (glm(virus ~ hormone +  gender + source + age + tissue, expr_meta, family = binomial()))
          odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
          results = rbind(results, data.frame(age_group = "all", hormone, virus_gene, gender = "all",  OR = round(odds["hormone", 1], 2), CI_25 =  round(odds["hormone", 2], 2) , 
                                                      CI_75 =  round(odds["hormone", 3],2) , p =  summary(mylogit)$coefficients["hormone", 4]))

          mylogit = (glm(virus ~ hormone + source + age + tissue, expr_meta[expr_meta$gender %in% c("male"),], family = binomial()))
          odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
          results = rbind(results, data.frame(age_group = "all", hormone, virus_gene,gender = "male",  OR = round(odds["hormone", 1], 2), CI_25 =  round(odds["hormone", 2], 2) , 
                                              CI_75 =  round(odds["hormone", 3],2) , p =  summary(mylogit)$coefficients["hormone", 4]))
          
          mylogit = (glm(virus ~ hormone +   source + age + tissue, expr_meta[expr_meta$gender %in% c("female"),], family = binomial()))
          odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
          results = rbind(results, data.frame(age_group = "all", hormone, virus_gene,gender = "female",  OR = round(odds["hormone", 1], 2), CI_25 =  round(odds["hormone", 2], 2) , 
                                              CI_75 =  round(odds["hormone", 3],2) , p =  summary(mylogit)$coefficients["hormone", 4]))
          
          #old
          mylogit = (glm(virus ~ hormone +  gender + source +   tissue, expr_meta[expr_meta$age %in% c("60-100") ,], family = binomial()))
          odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
          results = rbind(results, data.frame(age_group = "old", hormone, virus_gene,gender = "all",  OR = round(odds["hormone", 1], 2), CI_25 =  round(odds["hormone", 2], 2) , 
                                              CI_75 =  round(odds["hormone", 3],2) , p =  summary(mylogit)$coefficients["hormone", 4]))
          
          mylogit = (glm(virus ~ hormone + source +   tissue, expr_meta[expr_meta$age %in% c("60-100") & expr_meta$gender %in% c("male"),], family = binomial()))
          odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
          results = rbind(results, data.frame(age_group = "old",hormone, virus_gene, gender = "male",  OR = round(odds["hormone", 1], 2), CI_25 =  round(odds["hormone", 2], 2) , 
                                              CI_75 =  round(odds["hormone", 3],2) , p =  summary(mylogit)$coefficients["hormone", 4]))
          
          mylogit = (glm(virus ~ hormone +   source + tissue , expr_meta[expr_meta$age %in% c("60-100") & expr_meta$gender %in% c("female"),], family = binomial()))
          odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
          results = rbind(results, data.frame(age_group = "old", hormone, virus_gene,gender = "female",  OR = round(odds["hormone", 1], 2), CI_25 =  round(odds["hormone", 2], 2) , 
                                              CI_75 =  round(odds["hormone", 3],2) , p =  summary(mylogit)$coefficients["hormone", 4]))
    }
  }

fn = paste0('output/labeled/table3_p_', p, '.csv')
write.csv(results, fn)

##################
# table 4, different organs
#################
tissues = names(table(expr_meta$tissue))
tissues = tissues[!tissues %in% c('breast', "prostate", "testis")]
results = NULL
for (hormone in c("ESR1_bin", "ESR2_bin", "PGR_bin", "AR_bin")){
  for (virus_gene in c("ACE2_bin", "TMPRSS2_bin", "ACE2_TMPRSS2")){
    for (t in tissues){
      print(c(hormone, virus_gene, t))
      expr_meta$hormone = expr_meta[, hormone]
      expr_meta$virus = expr_meta[, virus_gene]
      
      #all age 
      expr_meta_subset = expr_meta[expr_meta$tissue %in% t,]
      if (length(unique(expr_meta_subset$virus))!=1 & length(unique(expr_meta_subset$hormone))!=1){
        mylogit = (glm(virus ~ hormone +  gender + source + age , expr_meta_subset, family = binomial()))
        odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
        results = rbind(results, data.frame(age_group = "all", hormone, virus_gene, gender = "all", tissue=t, OR = round(odds["hormone", 1], 2), CI_25 =  round(odds["hormone", 2], 2) , 
                                            CI_75 =  round(odds["hormone", 3],2) , p =  summary(mylogit)$coefficients["hormone", 4]))
      }
      
      expr_meta_subset =  expr_meta[expr_meta$gender %in% c("male") & expr_meta$tissue %in% t,] 
        if (length(unique(expr_meta_subset$virus))!=1 & length(unique(expr_meta_subset$hormone))!=1){
        mylogit = (glm(virus ~ hormone + source + age ,expr_meta_subset, family = binomial()))
        odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
        results = rbind(results, data.frame(age_group = "all", hormone, virus_gene,gender = "male",  tissue=t, OR = round(odds["hormone", 1], 2), CI_25 =  round(odds["hormone", 2], 2) , 
                                            CI_75 =  round(odds["hormone", 3],2) , p =  summary(mylogit)$coefficients["hormone", 4]))
      }
      
      expr_meta_subset =  expr_meta[expr_meta$gender %in% c("female") & expr_meta$tissue %in% t,]
        if (length(unique(expr_meta_subset$virus))!=1 & length(unique(expr_meta_subset$hormone))!=1){
        mylogit = (glm(virus ~ hormone +  source + age , expr_meta_subset, family = binomial()))
        odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
        results = rbind(results, data.frame(age_group = "all", hormone, virus_gene,gender = "female", tissue=t, OR = round(odds["hormone", 1], 2), CI_25 =  round(odds["hormone", 2], 2) , 
                                            CI_75 =  round(odds["hormone", 3],2) , p =  summary(mylogit)$coefficients["hormone", 4]))
      }
      
      #old
      expr_meta_subset =   expr_meta[expr_meta$age %in% c("60-100") & expr_meta$tissue %in% t ,]
      if (nrow(expr_meta_subset) > 30 & !t %in% c( "prostate", "testis", "small intestine", "heart") & length(unique(expr_meta_subset$virus))!=1 & length(unique(expr_meta_subset$hormone))!=1){
        mylogit = (glm(virus ~ hormone +  gender + source , expr_meta_subset, family = binomial()))
        odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
        results = rbind(results, data.frame(age_group = "old", hormone, virus_gene,gender = "all",  tissue=t, OR = round(odds["hormone", 1], 2), CI_25 =  round(odds["hormone", 2], 2) , 
                                            CI_75 =  round(odds["hormone", 3],2) , p =  summary(mylogit)$coefficients["hormone", 4]))
      }
      
      expr_meta_subset =   expr_meta[expr_meta$age %in% c("60-100") & expr_meta$gender %in% c("male")  & expr_meta$tissue %in% t,] 
      if (nrow(expr_meta_subset) > 30 & !t %in% c("bone", "brain", "breast", "small intestine", "heart") & length(unique(expr_meta_subset$virus))!=1 & length(unique(expr_meta_subset$hormone))!=1){
        mylogit = (glm(virus ~ hormone + source , expr_meta_subset, family = binomial()))
        odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
        results = rbind(results, data.frame(age_group = "old",hormone, virus_gene, gender = "male", tissue=t, OR = round(odds["hormone", 1], 2), CI_25 =  round(odds["hormone", 2], 2) , 
                                            CI_75 =  round(odds["hormone", 3],2) , p =  summary(mylogit)$coefficients["hormone", 4]))
      }
      
      expr_meta_subset =   expr_meta[expr_meta$age %in% c("60-100") & expr_meta$gender %in% c("female")  & expr_meta$tissue %in% t,]
      if (nrow(expr_meta_subset) > 30  & !t %in% c("bone", "brain", "heart", "prostate", "testis", "small intestine", "skin") & length(unique(expr_meta_subset$virus))!=1 & length(unique(expr_meta_subset$hormone))!=1){
        mylogit = (glm(virus ~ hormone +  source  ,expr_meta_subset, family = binomial()))
        odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
        results = rbind(results, data.frame(age_group = "old", hormone, virus_gene,gender = "female", tissue=t, OR = round(odds["hormone", 1], 2), CI_25 =  round(odds["hormone", 2], 2) , 
                                            CI_75 =  round(odds["hormone", 3],2) , p =  summary(mylogit)$coefficients["hormone", 4]))
      }
    }
  }
}

write.csv(results, "output/table4.csv")






###################
mylogit = (glm(DPP4_bin ~ ESR1_bin  + gender + source + age + tissue, expr_meta, family = binomial()))
summary(mylogit)
exp(cbind(OR = coef(mylogit), confint(mylogit)))

mylogit = (glm(DPP4_bin ~ ESR1_bin  + gender + source + tissue, expr_meta[expr_meta$age %in% c("60-100"),], family = binomial()))
summary(mylogit)
exp(cbind(OR = coef(mylogit), confint(mylogit)))

###
#hormone receptor vs age
expr_meta$age <- relevel(expr_meta$age, ref = "20-59")
mylogit = (glm(AR_bin ~  tissue + age + source, expr_meta[expr_meta$gender %in% c("male"), ], family = binomial()))
summary(mylogit)
exp(cbind(OR = coef(mylogit), confint(mylogit)))

######
#is smoking a factor?
mylogit = (glm(ACE2_bin ~   age + source + smoke + gender, expr_meta, family = binomial()))
summary(mylogit)
exp(cbind(OR = coef(mylogit), confint(mylogit)))

mylogit = (glm(ACE2_bin ~   age + source + smoke + smoke * gender, expr_meta, family = binomial()))
summary(mylogit)
exp(cbind(OR = coef(mylogit), confint(mylogit)))

mylogit = (glm(ACE2_bin ~   age + source + smoke , expr_meta[expr_meta$gender %in% c("male"),], family = binomial()))
summary(mylogit)
exp(cbind(OR = coef(mylogit), confint(mylogit)))


mylogit = (glm(TMPRSS2_bin ~   age + source + smoke, expr_meta, family = binomial()))
summary(mylogit)
exp(cbind(OR = coef(mylogit), confint(mylogit)))

mylogit = (glm(ACE2_TMPRSS2 ~   age + source + smoke, expr_meta, family = binomial()))
summary(mylogit)
exp(cbind(OR = coef(mylogit), confint(mylogit)))

