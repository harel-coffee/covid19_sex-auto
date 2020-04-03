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


hist((expr_meta$TMPRSS2 ))

t.test(TMPRSS2 ~ gender, data = expr_meta)
t.test(TMPRSS2 ~ smoke, data = expr_meta)
anova(lm(TMPRSS2 ~ smoke + gender, data = expr_meta))
anova(lm(TMPRSS2 ~  gender, data = expr_meta))
anova(lm(TMPRSS2 ~  age, data = expr_meta))

symbols = c("TMPRSS2", "ESR1", "DPP4",  "GPER1", "ESR2", "AR", "PGR", "TMPRSS2")

#write.csv(expr_meta, "GPL570_expr_meta.csv")
###############
##
expr_meta_subset = expr_meta[!is.na(expr_meta$gender), ] 
expr_meta_subset$two_group = paste("all", as.numeric(as.factor(expr_meta_subset$gender)))
group_patients = table(expr_meta_subset$two_group)
expr_meta_subset = subset(expr_meta_subset, two_group %in% names(group_patients[group_patients > 30]))

#find group statistics
#ignore age
ages = "all"
test = sapply(ages, function(age){
  expr_meta_subset_age = expr_meta_subset
  Ttest = t.test(TMPRSS2 ~ two_group, expr_meta_subset_age)
  paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
         "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
         "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
         "n: ", nrow(expr_meta_subset_age)
  )
  })


pdf(paste("figure/GPL570_TMPRSS2_gender.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((TMPRSS2 )  ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "TMPRSS2 Expression", side = "both", ylim=c(1,15),what = c(T, T, T, F),
         border = NA, innerborder = "gray", col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1:1), y = 14, labels = test)
dev.off()

ages = "all"
test = sapply(ages, function(age){
  expr_meta_subset_age = expr_meta_subset
  Ttest = t.test(DPP4 ~ two_group, expr_meta_subset_age)
  paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
         "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
         "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
         "n: ", nrow(expr_meta_subset_age)
  )})


pdf(paste("figure/GPL570_DPP4_gender.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((DPP4 )  ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "DPP4 Expression", side = "both", ylim=c(0,15),what = c(T, T, T, F),
         border = NA, innerborder = "gray", col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1:1), y = 14, labels = test)
dev.off()
####################
#####################
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
  Ttest = t.test(TMPRSS2 ~ two_group, expr_meta_subset_age, alternative = "two.sided")
  paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
         "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
         "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
         "n: ", nrow(expr_meta_subset_age)
  )
  })


pdf(paste("figure/GPL570_TMPRSS2_ages.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((TMPRSS2)  ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "TMPRSS2 Expression", side = "both", ylim=c(0,15),what = c(T, T, T, F),
         border = NA, col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1:length(unique(expr_meta_subset$age))), y = 14, labels = test)
dev.off()
#}
##################
##############
for (gender in c("female", "male", "all")){
for (i in 1:length(symbols)){
gene = symbols[i]

if (gender == "all"){
  expr_meta_subset =  expr_meta #[expr_meta$gender %in% c("female", "male"), ]
}else{
  expr_meta_subset =  expr_meta[expr_meta$gender %in% gender, ]
}
cutoff = qnorm(0.01, mean(expr_meta_subset[,gene]), sd(expr_meta_subset[, gene]), lower.tail=F) # quantile(expr_meta[, gene])[4]

expr_meta_subset$gene = expr_meta_subset[, gene]
expr_meta_subset$expr = "high"
expr_meta_subset$expr[expr_meta_subset[, gene] < cutoff] = "low"

Ttest = t.test(expr_meta_subset[, "TMPRSS2"] ~ expr_meta_subset$expr)

pdf(paste("figure/GPL570_TMPRSS2_", symbols[i], gender, ".pdf", sep="_"))
print(ggplot(expr_meta_subset, aes(x= (gene), y = (TMPRSS2 ), colour = expr )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
  geom_point(size=0.5) + 
  annotate("text", label = paste( "p = ", format(Ttest$p.value, digit=2), ", ratio = ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=3), ":1", sep=""),  x = 10, y = 14, size = 6) + 
  xlab(gene) + guides(shape=FALSE, size=FALSE) +
  ylab("TMPRSS2") + coord_cartesian(xlim = c(0, 15), ylim=c(0, 15)) 
)
dev.off()
}
}
########
#DPP4
for (gender in c("female", "male", "all")){
for (i in 1:length(symbols)){
  if (gender == "all"){
    expr_meta_subset =  expr_meta #[expr_meta$gender %in% c("female", "male"), ]
  }else{
    expr_meta_subset =  expr_meta[expr_meta$gender %in% gender, ]
  }
  
  gene = symbols[i]
  cutoff = qnorm(0.01, mean(expr_meta_subset[,gene]), sd(expr_meta_subset[, gene]), lower.tail=F) # quantile(expr_meta[, gene])[4]
  
  expr_meta_subset =  expr_meta[expr_meta$gender %in% c("female", "male"), ]
  expr_meta_subset$gene = expr_meta_subset[, gene]
  expr_meta_subset$expr = "high"
  expr_meta_subset$expr[expr_meta_subset[, gene] < cutoff] = "low"
  
  Ttest = t.test(expr_meta_subset[, "DPP4"] ~ expr_meta_subset$expr)
  
  pdf(paste("figure/GPL570_DPP4_", symbols[i], gender, ".pdf", sep="_"))
  print(ggplot(expr_meta_subset, aes(x= (gene ), y = (DPP4 ), colour = expr )) +  theme_bw()  + 
          theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
          geom_point(size=0.5) + 
          annotate("text", label = paste( "p = ", format(Ttest$p.value, digit=2), ", ratio = ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=3), ":1", sep=""),  x = 10, y = 14, size = 6) + 
          xlab(gene) + guides(shape=FALSE, size=FALSE) +
          ylab("DPP4") + coord_cartesian(xlim = c(0, 15), ylim=c(0, 15)) 
  )
  dev.off()
}
}

###############
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
  Ttest = t.test(TMPRSS2 ~ two_group, expr_meta_subset_age, alternative = "two.sided")
  paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
         "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
         "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
         "n: ", nrow(expr_meta_subset_age)
  )  }else{
    NA
  }
})



pdf(paste("figure/GPL570_TMPRSS2_tissue_30.pdf", sep="_"), width = 20)
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot(TMPRSS2   ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "TMPRSS2 Expression", side = "both", ylim=c(0,15), what = c(T, T, T, F),
         border = NA, col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("female", "male"))
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

