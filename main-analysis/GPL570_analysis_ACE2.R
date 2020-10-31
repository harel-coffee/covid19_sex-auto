####
#post processing
#setwd("~/Documents/stanford/wars/gender/data")

library("ggplot2")
library("beanplot")
library(OneR)
library(plyr)

####################
## default mean
#####################

rm(list=ls())
root = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root)
fn = paste0(root, '/raw/GPL570_expr_meta_final_v9.csv')
expr_meta = read.csv(fn)
expr_meta = expr_meta[!is.na(expr_meta$gender), ]

expr_meta$tissue = as.character(expr_meta$tissue)
expr_meta[is.na(expr_meta$tissue), "tissue"] = "unknown"
expr_meta$tissue = as.factor(expr_meta$tissue)
expr_meta$gender = factor(expr_meta$gender, levels = c("male", "female"))


#####################
nrow(expr_meta)
length(which(!is.na(expr_meta$gender.x)))
length(which(!is.na(expr_meta$gender.x)))/nrow(expr_meta)
length(which(!is.na(expr_meta$age.x)))
length(which(!is.na(expr_meta$age.x)))/nrow(expr_meta)
length(which(!is.na(expr_meta$tissue.x)))
length(which(!is.na(expr_meta$tissue.x)))/nrow(expr_meta)


length(which(!is.na(expr_meta$gender)))
length(which(!is.na(expr_meta$gender)))/nrow(expr_meta)
length(which(!is.na(expr_meta$age)))
length(which(!is.na(expr_meta$age)))/nrow(expr_meta)
length(which(!is.na(expr_meta$tissue)))
length(which(!is.na(expr_meta$tissue)))/nrow(expr_meta)

#######################
paste("# samples", nrow(expr_meta))
paste("% labeled gender", sum(table(expr_meta$gender.x))/nrow(expr_meta))
paste("% predicted gender", sum(table(expr_meta$gender))/nrow(expr_meta))
paste("% labeled tissue", sum(table(expr_meta$tissue.x))/nrow(expr_meta))
paste("% predicted tissue", sum(table(expr_meta$tissue))/nrow(expr_meta))
paste("% labeled age", sum(table(expr_meta$age.x))/nrow(expr_meta))
paste("% predicted age", sum(table(expr_meta$age))/nrow(expr_meta))
#######################

#####################
expr_meta_subset = expr_meta
t.test(ACE2 ~ gender, data = expr_meta_subset)
wilcox.test(ACE2 ~  gender,  conf.int = TRUE, data = expr_meta_subset)

anova(lm(ACE2 ~ smoke + gender, data = expr_meta_subset))
anova(lm(ACE2 ~  gender, data = expr_meta_subset))
summary(lm(ACE2 ~  gender, data = expr_meta_subset))
summary(lm(ACE2 ~  gender + age, data = expr_meta_subset))
summary(lm(ACE2 ~  gender + tissue, data = expr_meta_subset))

summary(lm(ACE2 ~    age, data = expr_meta_subset))
summary(lm(TMPRSS2 ~    age, data = expr_meta_subset))

########################################
## Part 1 ## overall gender difference
########################################
expr_meta_subset = expr_meta[!is.na(expr_meta$gender), ] 
expr_meta_subset$two_group = paste("all", as.numeric(as.factor(expr_meta_subset$gender)))
group_patients = table(expr_meta_subset$two_group)
expr_meta_subset = subset(expr_meta_subset, two_group %in% names(group_patients[group_patients > 30]))

#find group statistics
#ignore age
ages = "all"
test = sapply(ages, function(age){
  expr_meta_subset_age = expr_meta_subset
  
  # bootstrap for ratio and CI
  n_bs = 1000
  set.seed(1)
  ratio_list = rep(NA, n_bs)
  for (i in 1:n_bs){
    n = nrow(expr_meta_subset_age)
    bootstrap_idx = sample(seq(1, n), size=n, replace=TRUE)
    expr_meta_subset_age_bootstrap = expr_meta_subset_age[bootstrap_idx, ]
    df = aggregate(ACE2~two_group, data = expr_meta_subset_age_bootstrap, FUN=mean)
    ratio_list[i] = df[1, 2]/df[2, 2]
    if (i %% 100==0){print(i)}
  }
  ratio = mean(ratio_list)
  CI = quantile(ratio_list, probs=c(0.025, 0.975), na.rm=TRUE, names=FALSE)
  
  Ttest = t.test(ACE2 ~ two_group, expr_meta_subset_age)
  Wtest = wilcox.test(ACE2 ~  two_group,  conf.int = TRUE, data = expr_meta_subset_age)
  paste0(# "ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
        
        "M/F: ", format(ratio, digit=3), " ", "[", format(CI[1], digits=3)," ", format(CI[2], digits=3), "]", "\n",
        "M-F: ", format(Wtest$estimate[1], digit=2), " ", "[", format(Wtest$conf.int[1], digits=2), " ", format(Wtest$conf.int[2], digits=2), "]", "\n",
         
         # "effect size: ", format(Wtest$estimate[1], digit=2), "\n",
         # "CI: [", format(Wtest$conf.int[1], digits=2)," ", format(Wtest$conf.int[2], digits=2), "]", "\n", 

         "p: ", format(Wtest$p.value, digits = 2, scientific=T), "\n",
         "n: ", nrow(expr_meta_subset_age)
  )
  
  })


pdf(paste("figure/GPL570_ACE2_gender.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((ACE2 )  ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(0,15),what = c(T, T, T, F),
         border = NA, innerborder = "gray", col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "mean", overalllin="mean")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("male", "female"))
text(x=c(1:1), y = 14, labels = test)
dev.off()

##########################################
## Part 2 ## age specific gender difference
##########################################
expr_meta_subset = expr_meta[!is.na(expr_meta$age) & !is.na(expr_meta$gender), ]
expr_meta_subset$two_group = paste(expr_meta_subset$age, as.numeric(as.factor(expr_meta_subset$gender)))
group_patients = table(expr_meta_subset$two_group)
expr_meta_subset = subset(expr_meta_subset, two_group %in% names(group_patients[group_patients > 30]))

anova(lm(ACE2 ~ age, data = expr_meta_subset))
by(expr_meta_subset$ACE2, expr_meta_subset$age, mean)

#find group statistics
ages = sort(unique(expr_meta_subset$age))
test = sapply(ages, function(a){
  expr_meta_subset_age = subset(expr_meta_subset, age == a)
  
  # bootstrap for ratio and CI
  n_bs = 1000
  set.seed(1)
  ratio_list = rep(NA, n_bs)
  for (i in 1:n_bs){
    n = nrow(expr_meta_subset_age)
    bootstrap_idx = sample(seq(1, n), size=n, replace=TRUE)
    expr_meta_subset_age_bootstrap = expr_meta_subset_age[bootstrap_idx, ]
    df = aggregate(ACE2~two_group, data = expr_meta_subset_age_bootstrap, FUN=mean)
    ratio_list[i] = df[1, 2]/df[2, 2]
    if (i %% 100==0){print(i)}
  }
  ratio = mean(ratio_list)
  CI = quantile(ratio_list, probs=c(0.025, 0.975), na.rm=TRUE, names=FALSE)
  
  Ttest = t.test(ACE2 ~ two_group, expr_meta_subset_age)
  Wtest = wilcox.test(ACE2 ~  two_group,  conf.int = TRUE, data = expr_meta_subset_age)
  paste0(# "ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
    
    "M/F: ", format(ratio, digit=3), " ", "[", format(CI[1], digits=3)," ", format(CI[2], digits=3), "]", "\n",
    "M-F: ", format(Wtest$estimate[1], digit=2), " ", "[", format(Wtest$conf.int[1], digits=2), " ", format(Wtest$conf.int[2], digits=2), "]", "\n",
    
    # "effect size: ", format(Wtest$estimate[1], digit=2), "\n",
    # "CI: [", format(Wtest$conf.int[1], digits=2)," ", format(Wtest$conf.int[2], digits=2), "]", "\n", 
    
    "p: ", format(Wtest$p.value, digits = 2, scientific=T), "\n",
    "n: ", nrow(expr_meta_subset_age)
  )
  
  })


pdf(paste("figure/GPL570_ACE2_ages.pdf", sep="_"), width=10)
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((ACE2)  ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(0,15),what = c(T, T, T, F),
         border = NA, col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "mean", overalllin="mean")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("male", "female"))
text(x=c(1:length(unique(expr_meta_subset$age))), y = 14, labels = test)
dev.off()
#}
##################
##########################################
## Part 4 ## Receptor Gene vs Homone gene
##########################################
symbols = c("ESR1", "ESR2", "PGR", "AR")
for (a in c("60-100", "all"))
for (gender in c("female", "male")){
  print(gender)
for (i in 1:length(symbols)){
gene = symbols[i]
  print(gene)
if (gender == "all"){
  expr_meta_subset =  expr_meta
}else{
  expr_meta_subset =  expr_meta[expr_meta$gender %in% gender, ]
}
  
if (a == "60-100"){
    expr_meta_subset =  expr_meta_subset[expr_meta_subset$age %in% a, ]
}
  
  
cutoff1 = sort(expr_meta[,gene], decreasing = T)[round(nrow(expr_meta) * 0.05)] #qnorm(0.01, mean(expr_meta_subset[,gene]), sd(expr_meta_subset[, gene]), lower.tail=F)  #quantile(expr_meta[, "ACE2"])[4]
cutoff2 = sort(expr_meta[,"ACE2"], decreasing = T)[round(nrow(expr_meta) * 0.05)] #qnorm(0.01, mean(expr_meta_subset[,gene]), sd(expr_meta_subset[, gene]), lower.tail=F)  #quantile(expr_meta[, "ACE2"])[4]

expr_meta_subset$gene = expr_meta_subset[, gene]

expr_meta_subset$expr = NA
expr_meta_subset$expr[expr_meta_subset[, gene] > cutoff1 & expr_meta_subset[, "ACE2"] > cutoff2] = "high high"
expr_meta_subset$expr[expr_meta_subset[, gene] < cutoff1 & expr_meta_subset[, "ACE2"] < cutoff2] = "low low"
expr_meta_subset$expr[expr_meta_subset[, gene] < cutoff1 & expr_meta_subset[, "ACE2"] > cutoff2] = "low high"
expr_meta_subset$expr[expr_meta_subset[, gene] > cutoff1 & expr_meta_subset[, "ACE2"] < cutoff2] = "high low"

Ftest = fisher.test(matrix(table(expr_meta_subset$expr), nrow=2))
text0 = paste0("OR: ", round(Ftest$estimate, 2), "[", round(Ftest$conf.int[1], 2)," ", round(Ftest$conf.int[2],2), "]\n",
               "p:", format(Ftest$p.value, digits = 2, scientific=T)
)

# Ttest = t.test(expr_meta_subset[, "ACE2"] ~ expr_meta_subset$expr)
# Wtest = wilcox.test(expr_meta_subset[, "ACE2"] ~ expr_meta_subset$expr,  conf.int = TRUE)
# df = aggregate(expr_meta_subset[, "ACE2"] ~ expr_meta_subset$expr, FUN = median)
# 
# text = paste0(" ", nrow(expr_meta_subset), ", ",
#               format(df[1, 2] / df[2, 2], digit=3), ", ",
#               format(Wtest$p.value, digits = 2, scientific=T), ", "
# )
# print(text)
# 
# # print(c(df[1, 2], df[2, 2], df[1, 2] / df[2, 2]))
# # df = aggregate(expr_meta_subset[, "ACE2"] ~ expr_meta_subset$expr, FUN = mean)
# # print(c(df[1, 2], df[2, 2], df[1, 2] / df[2, 2]))
# 
# # linear model
# model = lm(expr_meta_subset[, "ACE2"]~expr_meta_subset[,gene])
# model_summary = summary(model)
# print(model_summary)
# 
# text = paste0(" ", nrow(expr_meta_subset), ", ",
#               format(model_summary$coefficients[2, 1], digit=3), ", ",
#               format(model_summary$coefficients[2, 4], digits = 2, scientific=T), ", "
# )
# print(text)

 pdf(paste("figure/GPL570_ACE2_", symbols[i], gender, a, ".pdf", sep="_"))
print(ggplot(expr_meta_subset, aes(x= (gene), y = (ACE2), colour = expr )) +  theme_bw()  +
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +
  geom_point(size=0.5) +
  annotate("text", label = text0,  x = 10, y = 14, size = 6) + 
  xlab(gene) + guides(shape=FALSE, size=FALSE) +
  ylab("ACE2") + coord_cartesian(xlim = c(0, 15), ylim=c(0, 15))
)
dev.off()
}
}
########

##########################################
## Part 3 ## tissue specific gender difference
##########################################
expr_meta_subset = expr_meta
main_tissues = c("lung", "liver", "heart", "pancreas", "blood",  "brain", "colon",   "skin", "bone", "kidney", "intestine", "breast", "testis", "prostate", "unknown") #,"breast",
expr_meta_subset$tissue = sapply(expr_meta_subset$tissue, function(x){
  main_tissue = NA
  for (t in main_tissues){
    if (length(grep(t, x, ignore.case = T)) > 0){
      main_tissue = t
      break
    }
  }
  main_tissue
})
top_tissues = names(table(expr_meta_subset$tissue))
expr_meta_subset = subset(expr_meta_subset, tissue %in% top_tissues & !is.na(gender))
expr_meta_subset$two_group = paste(expr_meta_subset$tissue, as.numeric(as.factor(expr_meta_subset$gender)))
table(expr_meta_subset$two_group)


expr_meta_subset$two_group = factor(expr_meta_subset$two_group, levels=c(
  "blood 1", "blood 2",  "bone 1", "bone 2", "brain 1", "brain 2", "breast 1", "breast 2", "colon 1", "colon 2", "heart 1", "heart 2",  "intestine 1", 
  "intestine 2", "kidney 1", "kidney 2", "liver 1", "liver 2", "lung 1", "lung 2",
  "pancreas 1", "pancreas 2", "prostate 1", "prostate 2",  "skin 1", "skin 2", "testis 1", "testis 2", "unknown 1", "unknown 2"))


group_patients = table(expr_meta_subset$two_group)
expr_meta_subset = subset(expr_meta_subset, two_group %in% names(group_patients[group_patients > 20]))

#find group statistics
tissues = as.character(sort(unique(expr_meta_subset$tissue)))
test = sapply(tissues, function(t){
  expr_meta_subset_age = subset(expr_meta_subset, tissue == t)
  if (sum(table(expr_meta_subset_age$two_group) > 10) > 1){
    
    # bootstrap for ratio and CI
    n_bs = 1000
    set.seed(1)
    ratio_list = rep(NA, n_bs)
    for (i in 1:n_bs){
      n = nrow(expr_meta_subset_age)
      bootstrap_idx = sample(seq(1, n), size=n, replace=TRUE)
      expr_meta_subset_age_bootstrap = expr_meta_subset_age[bootstrap_idx, ]
      df = aggregate(ACE2~two_group, data = expr_meta_subset_age_bootstrap, FUN=mean)
      ratio_list[i] = df[1, 2]/df[2, 2]
      if (i %% 100==0){print(i)}
    }
    ratio = mean(ratio_list)
    CI = quantile(ratio_list, probs=c(0.025, 0.975), na.rm=TRUE, names=FALSE)
    
    Ttest = t.test(ACE2 ~ two_group, expr_meta_subset_age)
    Wtest = wilcox.test(ACE2 ~  two_group,  conf.int = TRUE, data = expr_meta_subset_age)
    paste0(# "ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
      
      "M/F: ", format(ratio, digit=3), " ", "[", format(CI[1], digits=3)," ", format(CI[2], digits=3), "]", "\n",
      "M-F: ", format(Wtest$estimate[1], digit=2), " ", "[", format(Wtest$conf.int[1], digits=2), " ", format(Wtest$conf.int[2], digits=2), "]", "\n",
      
      # "effect size: ", format(Wtest$estimate[1], digit=2), "\n",
      # "CI: [", format(Wtest$conf.int[1], digits=2)," ", format(Wtest$conf.int[2], digits=2), "]", "\n", 
      
      "p: ", format(Wtest$p.value, digits = 2, scientific=T), "\n",
      "n: ", nrow(expr_meta_subset_age)
    )
  }else{
    NA
  }
})


pdf(paste("figure/GPL570_ACE2_tissue_20.pdf", sep="_"), width = 35)
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot(ACE2   ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(0,15), what = c(T, T, T, F),
         border = NA, col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "mean", overalllin="mean")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("male", "female"))
text(x=c(1:length(unique(expr_meta_subset$tissue))), y = 14, labels = test)
dev.off()
#}

##########################################
## Part 3 ## tissue specific gender difference (old group)
##########################################
# expr_meta_subset = expr_meta[expr_meta$age %in% c("60-100"),]
main_tissues = c("lung", "liver", "heart", "pancreas", "blood",  "brain", "colon",   "skin", "bone", "kidney", "intestine", "breast", "testis", "prostate", "unknown") #,"breast",
expr_meta_subset$tissue.x = sapply(expr_meta_subset$tissue, function(x){
  main_tissue = NA
  for (t in main_tissues){
    if (length(grep(t, x, ignore.case = T)) > 0){
      main_tissue = t
      break
    }
  }
  main_tissue
})
top_tissues = names(table(expr_meta_subset$tissue.x))
expr_meta_subset = subset(expr_meta_subset, tissue.x %in% top_tissues & !is.na(gender.x))
expr_meta_subset$two_group = paste(expr_meta_subset$tissue.x, as.numeric(as.factor(expr_meta_subset$gender.x)))
table(expr_meta_subset$two_group)


expr_meta_subset$two_group = factor(expr_meta_subset$two_group, levels=c(
  "blood 1", "blood 2",  "bone 1", "bone 2", "brain 1", "brain 2", "breast 1", "breast 2", "colon 1", "colon 2", "heart 1", "heart 2", "kidney 1", "kidney 2", "liver 1", "liver 2", "lung 1", "lung 2",
  "pancreas 1", "pancreas 2", "prostate 1", "prostate 2",  "skin 1", "skin 2",  "unknown 1", "unknown 2"))


group_patients = table(expr_meta_subset$two_group)
expr_meta_subset = subset(expr_meta_subset, two_group %in% names(group_patients[group_patients > 20]))

#find group statistics
tissues = as.character(sort(unique(expr_meta_subset$tissue)))
test = sapply(tissues, function(t){
  expr_meta_subset_age = subset(expr_meta_subset, tissue == t)
  if (sum(table(expr_meta_subset_age$two_group) > 10) > 1){
    
    # bootstrap for ratio and CI
    n_bs = 1000
    set.seed(1)
    ratio_list = rep(NA, n_bs)
    for (i in 1:n_bs){
      n = nrow(expr_meta_subset_age)
      bootstrap_idx = sample(seq(1, n), size=n, replace=TRUE)
      expr_meta_subset_age_bootstrap = expr_meta_subset_age[bootstrap_idx, ]
      df = aggregate(ACE2~two_group, data = expr_meta_subset_age_bootstrap, FUN=mean)
      ratio_list[i] = df[1, 2]/df[2, 2]
      if (i %% 100==0){print(i)}
    }
    ratio = mean(ratio_list)
    CI = quantile(ratio_list, probs=c(0.025, 0.975), na.rm=TRUE, names=FALSE)
    
    Ttest = t.test(ACE2 ~ two_group, expr_meta_subset_age)
    Wtest = wilcox.test(ACE2 ~  two_group,  conf.int = TRUE, data = expr_meta_subset_age)
    paste0(# "ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
      
      "M/F: ", format(ratio, digit=3), " ", "[", format(CI[1], digits=3)," ", format(CI[2], digits=3), "]", "\n",
      "M-F: ", format(Wtest$estimate[1], digit=2), " ", "[", format(Wtest$conf.int[1], digits=2), " ", format(Wtest$conf.int[2], digits=2), "]", "\n",
      
      # "effect size: ", format(Wtest$estimate[1], digit=2), "\n",
      # "CI: [", format(Wtest$conf.int[1], digits=2)," ", format(Wtest$conf.int[2], digits=2), "]", "\n", 
      
      "p: ", format(Wtest$p.value, digits = 2, scientific=T), "\n",
      "n: ", nrow(expr_meta_subset_age)
    )
  }else{
    NA
  }
})


pdf(paste("figure/GPL570_ACE2_tissue_20_old.pdf", sep="_"), width = 35)
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot(ACE2   ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(0,15), what = c(T, T, T, F),
         border = NA, col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "mean", overalllin="mean")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("male", "female"))
text(x=c(1:length(unique(expr_meta_subset$tissue))), y = 14, labels = test)
dev.off()
#}

#####################
###
#check airway
airway = sapply(1:nrow(expr_meta), function(x){
  if (length(grep("airway", paste(expr_meta$source_name_ch1[x], expr_meta$characteristics_ch1[x]), ignore.case = T)) > 0){
    if (length(grep("epithe", paste(expr_meta$source_name_ch1[x], expr_meta$characteristics_ch1[x]), ignore.case = T)) > 0){
      T
    }else{
      F
    }
  }else{
    F
  }
})

expr_meta_subset = expr_meta[airway , ]
Ttest = t.test(ACE2 ~ gender, expr_meta_subset, alternative = "two.sided")
text =     paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
                  "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
                  "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
                  "n: ", nrow(expr_meta_subset)
)

pdf(paste("figure/GPL570_ACE2_airway.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((ACE2)  ~ gender, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(0,15),what = c(T, T, T, F),
         border = NA, col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1), y = 14, labels = text)
dev.off()

##
#maternal blood
maternal = sapply(1:nrow(expr_meta), function(x){
  if (length(grep("blood", paste(expr_meta$source_name[x], expr_meta$characteristics[x]), ignore.case = T)) > 0){
    if (length(grep("maternal", paste(expr_meta$source_name[x], expr_meta$characteristics[x]), ignore.case = T)) > 0){
      T
    }else{
      F
    }
  }else{
    F
  }
})


############
#smoker 
expr_meta_subset = expr_meta[!is.na(expr_meta$smoke) & !is.na(expr_meta$gender), ]
expr_meta_subset$two_group = paste(expr_meta_subset$smoke, as.numeric(as.factor(expr_meta_subset$gender)))
group_patients = table(expr_meta_subset$two_group)
expr_meta_subset = subset(expr_meta_subset, two_group %in% names(group_patients[group_patients > 30]))

anova(lm(ACE2 ~ smoke, data = expr_meta_subset))
by(expr_meta_subset$ACE2, expr_meta_subset$smoke, mean)

#find group statistics
smokes = sort(unique(expr_meta_subset$smoke))
test = sapply(smokes, function(a){
  expr_meta_subset_smoke = subset(expr_meta_subset, smoke == a)
  Ttest = t.test(ACE2 ~ two_group, expr_meta_subset_smoke, alternative = "two.sided")
  paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
         "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
         "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
         "n: ", nrow(expr_meta_subset_smoke)
  )
})


pdf(paste("figure/GPL570_ACE2_smokes.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((ACE2)  ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(0,15),what = c(T, T, T, F),
         border = NA, col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1:length(unique(expr_meta_subset$smoke))), y = 14, labels = test)
dev.off()



############# generate table ######################
expr_meta_subset = expr_meta[!is.na(expr_meta$age) & !is.na(expr_meta$gender), ]
expr_meta_subset$two_group = paste(expr_meta_subset$age, as.numeric(as.factor(expr_meta_subset$gender)))
group_patients = table(expr_meta_subset$two_group)
expr_meta_subset = subset(expr_meta_subset, two_group %in% names(group_patients[group_patients > 30]))

#find group statistics
ages = sort(unique(expr_meta_subset$age))
test = sapply(ages, function(a){
  expr_meta_subset_age = subset(expr_meta_subset, age == a)
  
  Wtest = wilcox.test(ACE2 ~ two_group, data = expr_meta_subset_age, conf.int = TRUE)
  df = aggregate(ACE2 ~ two_group, data = expr_meta_subset_age, FUN=median)
  text = paste0(" ", nrow(expr_meta_subset_age), ", ",
         format(df[1, 2] / df[2, 2], digit=3), ", ", 
         format(Wtest$p.value, digits = 2, scientific=T), ", "
         )
  print(text)
  
  tissues = c('liver', 'lung', 'kidney', 'colon')
  test_tissue = sapply(tissues, function(t){
    expr_meta_subset_age_tissue = subset(expr_meta_subset_age, tissue==t)
    
    Wtest = wilcox.test(ACE2 ~ gender, data = expr_meta_subset_age_tissue, conf.int = TRUE)
    df = aggregate(ACE2 ~ two_group, data = expr_meta_subset_age_tissue, FUN=median)
    text = paste0(" ", nrow(expr_meta_subset_age_tissue), ", ",
                  format(df[1, 2] / df[2, 2], digit=3), ", ", 
                  format(Wtest$p.value, digits = 2, scientific=T), ", ")
    print(text)
    
  })
  # print("")
  
})
print(ages)


########## for reference code ########
#compare tissue
expr_meta_subset = expr_meta
main_tissues = c("lung", "liver",  "pancrea", "blood",  "brain", "colon",   "skin", "bone", "kidney") #,"breast",
expr_meta_subset$tissue = sapply(expr_meta_subset$tissue, function(x){
  main_tissue = NA
  for (t in main_tissues){
    if (length(grep(t, x, ignore.case = T)) > 0){
      main_tissue = t
      break
    }
  }
  main_tissue
})
top_tissues = names(table(expr_meta_subset$tissue))
expr_meta_subset = subset(expr_meta_subset, tissue %in% top_tissues & !is.na(gender))
expr_meta_subset$two_group = paste(expr_meta_subset$tissue, as.numeric(as.factor(expr_meta_subset$gender)))
group_patients = table(expr_meta_subset$two_group)
expr_meta_subset = subset(expr_meta_subset, two_group %in% names(group_patients[group_patients > 20]))

#find group statistics
tissues = as.character(sort(unique(expr_meta_subset$tissue)))
test = sapply(tissues, function(t){
  expr_meta_subset_age = subset(expr_meta_subset, tissue == t)
  if (sum(table(expr_meta_subset_age$two_group) > 10) > 1){
    Ttest = t.test(ACE2 ~ two_group, expr_meta_subset_age, alternative = "two.sided")
    Wtest = wilcox.test(ACE2 ~ two_group, data = expr_meta_subset_age, conf.int = TRUE)
    paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
           "effect size: ", format(Wtest$estimate[1], digits=2), "\n",
           "CI: [", format(Wtest$conf.int[1], digits=2)," ", format(Wtest$conf.int[2], digits=2), "]", "\n", 
           "p: ", format(Wtest$p.value, digits = 2, scientific=T), "\n",
           "n: ", nrow(expr_meta_subset_age)
    )  }else{
      NA
    }
})

################ receptor vs hormone #############################
ages = c("0-19", "20-59", "60-100", "all")
ages = c("all")
for (age_group in ages){
  for (target_gene in c("ACE2")){# , "TMPRSS2"
    cdata_all = NULL
    
    #target_gene = "TMPRSS2" #TMPRSS2
    expr_meta = expr_meta[!is.na(expr_meta[, target_gene]),]
    expr_meta$bin = cut(expr_meta[, target_gene], breaks=c(quantile(expr_meta[, target_gene], probs = c(0.0, 0.7, 0.9, 1.0))), 
        labels=c("Normal","High", "Very High"), include.lowest=TRUE) # bin(expr_meta[, target_gene], 3)
    
    symbols = c("ESR1", "ESR2", "AR", "PGR")
    for (gender in c("female", "male")){
      for (i in 1:length(symbols)){
        gene = symbols[i]
        
        if (gender == "all"){
          expr_meta_subset =  expr_meta #[expr_meta$gender %in% c("female", "male"), ]
        }else{
          expr_meta_subset =  expr_meta[expr_meta$gender %in% gender, ]
        }
        if (age_group != "all"){
          expr_meta_subset= expr_meta_subset[!is.na(expr_meta_subset$age) & expr_meta_subset$age == age_group ,]
        }
        expr_meta_subset$gene = expr_meta_subset[, gene]
        cdata <- ddply(expr_meta_subset, .(bin), summarise, 
                       N    = length(gene),
                       mean = mean(gene),
                       se   = sd(gene, na.rm =T) / sqrt(length(gene))
        )
        
        # anova_test = anova(lm(expr_meta_subset[, gene] ~ expr_meta_subset$bin + expr_meta_subset$tissue))
        # a1 <- aov(expr_meta_subset[, gene] ~ expr_meta_subset$bin + expr_meta_subset$tissue)
        # posthoc <- TukeyHSD(x=a1, 'expr_meta_subset$bin', conf.level=0.95)
        
        anova_test = anova(lm(expr_meta_subset[, gene] ~ expr_meta_subset$bin))
        a1 <- aov(expr_meta_subset[, gene] ~ expr_meta_subset$bin)
        posthoc <- TukeyHSD(x=a1, 'expr_meta_subset$bin', conf.level=0.95)
        
        # for (t in unique(expr_meta_subset$tissue)){
        #   expr_meta_subset_tissue = expr_meta_subset[expr_meta_subset$tissue %in% t, ]
        #   a1 <- aov(expr_meta_subset_tissue[, gene] ~ expr_meta_subset_tissue$bin )
        #   posthoc <- TukeyHSD(x=a1, 'expr_meta_subset_tissue$bin', conf.level=0.95)
        #   print(t)
        #   print(posthoc)
        # }
        
        p = format(anova_test$`Pr(>F)`[1], digits = 2, scientific=T) 
        cdata_all = rbind(cdata_all, data.frame(gene, gender, p, cdata))
      }
    }
    
    cdata_all  = cdata_all[!is.na(cdata_all$bin), ]
    
    # write.csv(cdata_all, paste("figure/bin_plot/GPL570", target_gene, age_group, "_data.csv", sep="_"))
    cdata_all = cdata_all[cdata_all$gender %in% c("male", "female"), ]
    cdata_all$gender_gene = paste(cdata_all$gender, cdata_all$gene)
    #pdf(paste("figure/bin_plot/GPL570", target_gene, age_group, "bin.pdf", sep="_"))
    print(
      ggplot(cdata_all, aes(x=bin, y=mean, group = gender_gene, color = gene, shape = gender)) + coord_cartesian(ylim=c(3, 6))  +  theme_bw() +
        geom_point(position=position_dodge(), stat="identity", size = 4) + geom_line(aes(linetype=gender)) + 
        geom_errorbar(aes(ymin=mean, ymax=mean+se),
                      width=.2) +
        xlab(target_gene) +
        ylab("Expression") + coord_cartesian(ylim=c(2, 10)) 
    ) 
    #dev.off()
    
  }
}

##################### version 2 : transpose   #########################
library(binr)
ages = c("all")#"0-19", "20-59", "60-100",
# expr_meta = expr_meta[which(expr_meta$tissue=='kidney'), ]  # only look at kidney
for (age_group in ages){
  for (target_gene in c("ESR1")){  # "ACE2", "TMPRSS2", , "ESR2", "AR", "PGR"
    cdata_all = NULL
    
    #target_gene = "TMPRSS2" #TMPRSS2
    expr_meta = expr_meta[!is.na(expr_meta[, target_gene]),]
    expr_meta$bin = cut(expr_meta[, target_gene], breaks=c(quantile(expr_meta[, target_gene], probs = c(0.0, 0.7, 0.98, 1.0))), 
                        labels=c("Normal","High", "Very High"), include.lowest=TRUE) # bin(expr_meta[, target_gene], 3)
    
    symbols = c("ACE2", "TMPRSS2", "DPP4") # "ESR1", "ESR2", "AR", "PGR"
    for (gender in c("female", "male")){  #"female", "male"
      for (i in 1:length(symbols)){
        gene = symbols[i]
        
        if (gender == "all"){
          expr_meta_subset =  expr_meta #[expr_meta$gender %in% c("female", "male"), ]
        }else{
          expr_meta_subset =  expr_meta[expr_meta$gender %in% gender, ]
        }
        if (age_group != "all"){
          expr_meta_subset= expr_meta_subset[!is.na(expr_meta_subset$age) & expr_meta_subset$age == age_group ,]
        }
        expr_meta_subset$gene = expr_meta_subset[, gene]
        cdata <- ddply(expr_meta_subset, .(bin), summarise, 
                       N    = length(gene),
                       mean = mean(gene),
                       se   = sd(gene, na.rm =T) / sqrt(length(gene))
        )
        
        # anova_test = anova(lm(expr_meta_subset[, gene] ~ expr_meta_subset$bin + expr_meta_subset$tissue))
        # a1 <- aov(expr_meta_subset[, gene] ~ expr_meta_subset$bin + expr_meta_subset$tissue)
        # posthoc <- TukeyHSD(x=a1, 'expr_meta_subset$bin', conf.level=0.95)
        # print(posthoc)

        anova_test = anova(lm(expr_meta_subset[, gene] ~ expr_meta_subset$bin))
        a1 <- aov(expr_meta_subset[, gene] ~ expr_meta_subset$bin)
        posthoc <- TukeyHSD(x=a1, 'expr_meta_subset$bin', conf.level=0.95)
        
        # for (t in unique(expr_meta_subset$tissue)){
        #   expr_meta_subset_tissue = expr_meta_subset[expr_meta_subset$tissue %in% t, ]
        #   a1 <- aov(expr_meta_subset_tissue[, gene] ~ expr_meta_subset_tissue$bin )
        #   posthoc <- TukeyHSD(x=a1, 'expr_meta_subset_tissue$bin', conf.level=0.95)
        #   print(t)
        #   print(posthoc)
        # }
        
        p = format(anova_test$`Pr(>F)`[1], digits = 2, scientific=T) 
        cdata_all = rbind(cdata_all, data.frame(gene, gender, p, cdata))
      }
    }
    
    cdata_all  = cdata_all[!is.na(cdata_all$bin), ]
    
    # write.csv(cdata_all, paste("figure/bin_plot/GPL570", target_gene, age_group, "_data.csv", sep="_"))
    cdata_all = cdata_all[cdata_all$gender %in% c("male", "female"), ] # "male", "female"
    cdata_all$gender_gene = paste(cdata_all$gender, cdata_all$gene)
    #pdf(paste("figure/bin_plot/GPL570", target_gene, age_group, "bin.pdf", sep="_"))
    print(
      ggplot(cdata_all[which(cdata_all$bin != 'group2'), ], aes(x=bin, y=mean, group = gender_gene, color = gene, shape = gender)) + coord_cartesian(ylim=c(3, 7))  +  theme_bw() +
        geom_point(position=position_dodge(), stat="identity", size = 4) + geom_line(aes(linetype=gender)) + 
        geom_errorbar(aes(ymin=mean, ymax=mean+se),
                      width=.2) +
        xlab(target_gene) +
        ylab("Expression") + coord_cartesian(ylim=c(3, 7)) 
    ) 
    #dev.off()
    
  }
}

####################### binary RG vs Categorical HG #######################
expr_meta_subset = expr_meta[!is.na(expr_meta$age) & !is.na(expr_meta$tissue) & !is.na(expr_meta$gender), ]

receptor_gene = 'ACE2'
expr_meta_subset$bin_RG = cut(expr_meta_subset[, receptor_gene], breaks=c(quantile(expr_meta[, receptor_gene], probs = c(0.0, 0.7, 0.9, 1.0))), 
    labels=c("low","high","very high"), include.lowest=TRUE) 
hormone_gene = 'ESR1'
expr_meta_subset$bin_HG = cut(expr_meta_subset[, hormone_gene], breaks=c(quantile(expr_meta[, receptor_gene], probs = c(0.0, 0.7, 0.9, 1.0))), 
                              labels=c("low","high", "very high"), include.lowest=TRUE) 
mod1 = glm(bin_RG ~ gender, data=expr_meta_subset, family=binomial)
summary(mod1)
mod2 = glm(bin_RG ~ gender + age + gender*age, data=expr_meta_subset, family=binomial)
summary(mod2)
mod3 = glm(bin_RG ~ gender + tissue + gender*tissue, data=expr_meta_subset, family=binomial)
summary(mod3)
mod4 = glm(bin_RG ~ gender + tissue + gender*tissue + gender*age + gender*age*tissue, data=expr_meta_subset, family=binomial)
summary(mod4)

mod20 = glm(bin_HG ~ gender, data=expr_meta_subset, family=binomial)
summary(mod20)
mod21 = glm(bin_RG ~ bin_HG, data=expr_meta_subset, family=binomial)
summary(mod21)
mod22 = glm(bin_RG ~ bin_HG + gender, data=expr_meta_subset, family=binomial)
summary(mod22)
mod23 = glm(bin_RG ~ gender + bin_HG + gender*bin_HG, data=expr_meta_subset, family=binomial)
summary(mod23)

mod = glm(gender ~ bin_HG + bin_RG, data=expr_meta_subset, family=binomial)
summary(mod)

age_list=c("0-19", "20-59", "60-100")
for (age in age_list){
  print(age)
  d = expr_meta_subset[expr_meta_subset$age==age,]
  mod = glm(gender ~ bin_HG + bin_RG, data=d, family=binomial)
  print(summary(mod))
}

