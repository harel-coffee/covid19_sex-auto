####
#post processing
#setwd("~/Documents/stanford/wars/gender/data")

library("ggplot2")
library("beanplot")

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


###############
paste("# samples", nrow(expr_meta))
paste("% labeled gender", sum(table(expr_meta$gender.x))/nrow(expr_meta))
paste("% predicted gender", sum(table(expr_meta$gender))/nrow(expr_meta))
paste("% labeled tissue", sum(table(expr_meta$tissue.x))/nrow(expr_meta))
paste("% predicted tissue", sum(table(expr_meta$tissue))/nrow(expr_meta))
paste("% labeled age", sum(table(expr_meta$age.x))/nrow(expr_meta))
paste("% predicted age", sum(table(expr_meta$age))/nrow(expr_meta))
###############

expr_meta = expr_meta[!is.na(expr_meta$gender), ]   #GEO GPL570 no need to remove
symbols = c( "ESR1", "DPP4",  "GPER1", "ESR2", "AR", "PGR", "ACE2")

#####################
expr_meta_subset = expr_meta
t.test(TMPRSS2 ~ gender, data = expr_meta_subset)
anova(lm(TMPRSS2 ~ smoke + gender, data = expr_meta_subset))
anova(lm(TMPRSS2 ~  gender, data = expr_meta_subset))
summary(lm(TMPRSS2 ~  gender, data = expr_meta_subset))
summary(lm(TMPRSS2 ~  gender + age, data = expr_meta_subset))
summary(lm(TMPRSS2 ~  gender + tissue, data = expr_meta_subset))

summary(lm(TMPRSS2 ~    age, data = expr_meta_subset))
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
    df = aggregate(TMPRSS2~two_group, data = expr_meta_subset_age_bootstrap, FUN=mean)
    ratio_list[i] = df[1, 2]/df[2, 2]
    if (i %% 100==0){print(i)}
  }
  ratio = mean(ratio_list)
  CI = quantile(ratio_list, probs=c(0.025, 0.975), na.rm=TRUE, names=FALSE)
  
  Ttest = t.test(TMPRSS2 ~ two_group, expr_meta_subset_age)
  Wtest = wilcox.test(TMPRSS2 ~  two_group,  conf.int = TRUE, data = expr_meta_subset_age)
  paste0(# "ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
    
    "M/F: ", format(ratio, digit=3), " ", "[", format(CI[1], digits=3)," ", format(CI[2], digits=3), "]", "\n",
    "M-F: ", format(Wtest$estimate[1], digit=2), " ", "[", format(Wtest$conf.int[1], digits=2), " ", format(Wtest$conf.int[2], digits=2), "]", "\n",
    
    # "effect size: ", format(Wtest$estimate[1], digit=2), "\n",
    # "CI: [", format(Wtest$conf.int[1], digits=2)," ", format(Wtest$conf.int[2], digits=2), "]", "\n", 
    
    "p: ", format(Wtest$p.value, digits = 2, scientific=T), "\n",
    "n: ", nrow(expr_meta_subset_age)
  )
  
  })


pdf(paste("figure/GPL570_TMPRSS2_gender.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((TMPRSS2 )  ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "TMPRSS2 Expression", side = "both", ylim=c(1,15),what = c(T, T, T, F),
         border = NA, innerborder = "gray", col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "mean", overalllin="mean")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("male", "female"))
text(x=c(1:1), y = 14, labels = test)
dev.off()

########################################
## Part 2 ## age specific gender difference
########################################
expr_meta_subset = expr_meta[!is.na(expr_meta$age) & !is.na(expr_meta$gender), ]
expr_meta_subset$two_group = paste(expr_meta_subset$age, as.numeric(as.factor(expr_meta_subset$gender)))
group_patients = table(expr_meta_subset$two_group)
expr_meta_subset = subset(expr_meta_subset, two_group %in% names(group_patients[group_patients > 30]))

anova(lm(TMPRSS2 ~ age, data = expr_meta_subset))
by(expr_meta_subset$TMPRSS2, expr_meta_subset$age, mean)

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
    df = aggregate(TMPRSS2~two_group, data = expr_meta_subset_age_bootstrap, FUN=mean)
    ratio_list[i] = df[1, 2]/df[2, 2]
    if (i %% 100==0){print(i)}
  }
  ratio = mean(ratio_list)
  CI = quantile(ratio_list, probs=c(0.025, 0.975), na.rm=TRUE, names=FALSE)
  
  Ttest = t.test(TMPRSS2 ~ two_group, expr_meta_subset_age)
  Wtest = wilcox.test(TMPRSS2 ~  two_group,  conf.int = TRUE, data = expr_meta_subset_age)
  paste0(# "ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
    
    "M/F: ", format(ratio, digit=3), " ", "[", format(CI[1], digits=3)," ", format(CI[2], digits=3), "]", "\n",
    "M-F: ", format(Wtest$estimate[1], digit=2), " ", "[", format(Wtest$conf.int[1], digits=2), " ", format(Wtest$conf.int[2], digits=2), "]", "\n",
    
    # "effect size: ", format(Wtest$estimate[1], digit=2), "\n",
    # "CI: [", format(Wtest$conf.int[1], digits=2)," ", format(Wtest$conf.int[2], digits=2), "]", "\n", 
    
    "p: ", format(Wtest$p.value, digits = 2, scientific=T), "\n",
    "n: ", nrow(expr_meta_subset_age)
  )
  
  })


pdf(paste("figure/GPL570_TMPRSS2_ages.pdf", sep="_"), width=10)
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((TMPRSS2)  ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "TMPRSS2 Expression", side = "both", ylim=c(0,15),what = c(T, T, T, F),
         border = NA, col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "mean", overalllin="mean")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("male", "female"))
text(x=c(1:length(unique(expr_meta_subset$age))), y = 14, labels = test)
dev.off()
#}
##################
########################################
## Part 4 ## receptor gene vs homorne gene
########################################
symbols = c("ESR1", "ESR2", "PGR", "AR", "GPER1")
for (gender in c("female", "male")){
for (i in 1:length(symbols)){
gene = symbols[i]

if (gender == "all"){
  expr_meta_subset =  expr_meta #[expr_meta$gender %in% c("female", "male"), ]
}else{
  expr_meta_subset =  expr_meta[expr_meta$gender %in% gender, ]
}
cutoff = sort(expr_meta_subset[,gene], decreasing = T)[round(nrow(expr_meta_subset) * 0.1)] #qnorm(0.01, mean(expr_meta_subset[,gene]), sd(expr_meta_subset[, gene]), lower.tail=F)  #quantile(expr_meta[, "TMPRSS2"])[4]

expr_meta_subset$gene = expr_meta_subset[, gene]
expr_meta_subset$expr = "high"
expr_meta_subset$expr[expr_meta_subset[, gene] < cutoff] = "low"

Ttest = t.test(expr_meta_subset[, "TMPRSS2"] ~ expr_meta_subset$expr)
Wtest = wilcox.test(expr_meta_subset[, "TMPRSS2"] ~ expr_meta_subset$expr,  conf.int = TRUE)
df = aggregate(expr_meta_subset[, "TMPRSS2"] ~ expr_meta_subset$expr, FUN = median)

text = paste0(" ", nrow(expr_meta_subset), ", ",
              format(df[1, 2] / df[2, 2], digit=3), ", ", 
              format(Wtest$p.value, digits = 2, scientific=T), ", "
)
print(text)


pdf(paste("figure/GPL570_TMPRSS2_", symbols[i], gender, ".pdf", sep="_"))
print(ggplot(expr_meta_subset, aes(x= (gene), y = (TMPRSS2 ), colour = expr )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
  geom_point(size=0.5) + 
  annotate("text", label = paste( "p = ", format(Wtest$p.value, digit=2), ", ratio = ",  format(df[1, 2] / df[2, 2], digit=3), ":1", sep=""),  x = 10, y = 14, size = 6) + 
  xlab(gene) + guides(shape=FALSE, size=FALSE) +
  ylab("TMPRSS2") + coord_cartesian(xlim = c(0, 15), ylim=c(0, 15)) 
)
dev.off()

}
}
########

###############
########################################
## Part 3 ## tissue specific gender difference
########################################
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
      df = aggregate(TMPRSS2~two_group, data = expr_meta_subset_age_bootstrap, FUN=mean)
      ratio_list[i] = df[1, 2]/df[2, 2]
      if (i %% 100==0){print(i)}
    }
    ratio = mean(ratio_list)
    CI = quantile(ratio_list, probs=c(0.025, 0.975), na.rm=TRUE, names=FALSE)
    
    Ttest = t.test(TMPRSS2 ~ two_group, expr_meta_subset_age)
    Wtest = wilcox.test(TMPRSS2 ~  two_group,  conf.int = TRUE, data = expr_meta_subset_age)
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


pdf(paste("figure/GPL570_TMPRSS2_tissue_20.pdf", sep="_"), width = 35)
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot(TMPRSS2   ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "TMPRSS2 Expression", side = "both", ylim=c(0,15), what = c(T, T, T, F),
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
Ttest = t.test(TMPRSS2 ~ gender, expr_meta_subset, alternative = "two.sided")
text =     paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
                  "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
                  "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
                  "n: ", nrow(expr_meta_subset)
)

pdf(paste("figure/GPL570_TMPRSS2_airway.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((TMPRSS2)  ~ gender, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "TMPRSS2 Expression", side = "both", ylim=c(0,15),what = c(T, T, T, F),
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

anova(lm(TMPRSS2 ~ smoke, data = expr_meta_subset))
by(expr_meta_subset$TMPRSS2, expr_meta_subset$smoke, mean)

#find group statistics
smokes = sort(unique(expr_meta_subset$smoke))
test = sapply(smokes, function(a){
  expr_meta_subset_smoke = subset(expr_meta_subset, smoke == a)
  Ttest = t.test(TMPRSS2 ~ two_group, expr_meta_subset_smoke, alternative = "two.sided")
  paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
         "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
         "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
         "n: ", nrow(expr_meta_subset_smoke)
  )
})


pdf(paste("figure/GPL570_TMPRSS2_smokes.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((TMPRSS2)  ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "TMPRSS2 Expression", side = "both", ylim=c(0,15),what = c(T, T, T, F),
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

  Wtest = wilcox.test(TMPRSS2 ~ two_group, data = expr_meta_subset_age, conf.int = TRUE)
  df = aggregate(TMPRSS2 ~ two_group, data = expr_meta_subset_age, FUN=median)
  text = paste0(" ", nrow(expr_meta_subset_age), ", ",
                format(df[1, 2] / df[2, 2], digit=3), ", ", 
                format(Wtest$p.value, digits = 2, scientific=T), ", "
  )
  print(text)
  
  tissues = c('liver', 'lung', 'kidney', 'colon')
  test_tissue = sapply(tissues, function(t){
    expr_meta_subset_age_tissue = subset(expr_meta_subset_age, tissue==t)
    
    Wtest = wilcox.test(TMPRSS2 ~ gender, data = expr_meta_subset_age_tissue, conf.int = TRUE)
    df = aggregate(TMPRSS2 ~ two_group, data = expr_meta_subset_age_tissue, FUN=median)
    text = paste0(" ", nrow(expr_meta_subset_age_tissue), ", ",
                  format(df[1, 2] / df[2, 2], digit=3), ", ", 
                  format(Wtest$p.value, digits = 2, scientific=T), ", ")
    print(text)
    
  })
  # print("")
  
})
print(ages)
