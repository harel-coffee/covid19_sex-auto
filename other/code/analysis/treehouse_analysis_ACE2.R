#analyze ACE2 from the treehouse project
setwd("~/Documents/stanford/wars/gender/data")

library(ggplot2)
library(psych)
library(beanplot)
library(cowplot)
library(stringr)


expr_meta = read.csv("treehouse_expr_meta.csv")
###############
paste("# samples", nrow(expr_meta))
paste("% predicted gender", sum(table(expr_meta$gender))/nrow(expr_meta))
paste("% labeled tissue", sum(table(expr_meta$tissue_old))/nrow(expr_meta))
paste("% predicted tissue", sum(table(expr_meta$tissue))/nrow(expr_meta))
paste("% predicted age", sum(table(expr_meta$age))/nrow(expr_meta))

expr_meta = expr_meta[!is.na(expr_meta$gender) & expr_meta$ACE2 > 0.5, ] 
symbols = c("ACE2", "ESR1", "DPP4",  "GPER1", "ESR2", "AR", "PGR", "TMPRSS2")

###################
expr_meta_subset = expr_meta

t.test(ACE2 ~ gender, data = expr_meta_subset)
anova(lm(ACE2 ~  gender, data = expr_meta_subset))
summary(lm(ACE2 ~  gender, data = expr_meta_subset))
summary(lm(ACE2 ~  gender + age, data = expr_meta_subset))
summary(lm(ACE2 ~  gender + tissue, data = expr_meta_subset))

summary(lm(ACE2 ~    age, data = expr_meta_subset))
summary(lm(TMPRSS2 ~    age, data = expr_meta_subset))

###################
table(expr_meta$age)
table(expr_meta$type)

#visualize gtex normal
expr_meta_subset = expr_meta[expr_meta$source == "gtex", ]
#describe.by(expr_meta_subset[, c("ACE2", "gender")], group = "gender")
t.test(ACE2 ~ gender, data = expr_meta)
tissues = table(tolower(expr_meta_subset$tissue))
tissues = names(tissues[tissues > 30])
p = sapply(tissues, function(x){
  expr_meta_subset = expr_meta[tolower(expr_meta$tissue) == x,]
  if (sum(table(expr_meta_subset$gender) > 10) > 1){
    t.test(ACE2 ~ gender, data = expr_meta_subset)$p.value
  }else{
    NA
  }
})

#+ geom_jitter()

##########################
#plot ACE2 expression between male and female using all data
expr_meta_subset = expr_meta[!is.na(expr_meta$gender), ] 
expr_meta_subset$two_group = paste("all", as.numeric(as.factor(expr_meta_subset$gender)))
group_patients = table(expr_meta_subset$two_group)
expr_meta_subset = subset(expr_meta_subset, two_group %in% names(group_patients[group_patients > 30]))

#find group statistics
#ignore age
ages = "all"
test = sapply(ages, function(age){
  expr_meta_subset_age = expr_meta_subset
  Ttest = t.test(ACE2 ~ two_group, expr_meta_subset_age)
  paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
         "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
         "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
         "n: ", nrow(expr_meta_subset_age)
  )
  })


pdf(paste("figure/treehouse_ACE2_gender.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot(ACE2  ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(-1,6),what = c(T, T, T, F),
         border = NA, innerborder = "gray", col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1:1), y = 5, labels = test)
dev.off()

########################
#########################
#plot ACE2 express across age
expr_meta_subset = expr_meta[expr_meta$type %in% c("normal", "adjacent", "cancer"), ] #, "cancer"
#expr_meta_subset = expr_meta_subset[unique(c(grep("leu", expr_meta_subset$disease, ignore.case = T),
#                                     grep("blood", expr_meta_subset$tissue, ignore.case = T))), ]
expr_meta_subset$two_group = paste(expr_meta_subset$age, as.numeric(as.factor(expr_meta_subset$gender)))
group_patients = table(expr_meta_subset$two_group)
expr_meta_subset = subset(expr_meta_subset, two_group %in% names(group_patients[group_patients > 20]))

#find group statistics
ages = sort(unique(expr_meta_subset$age))
test = sapply(ages, function(a){
  expr_meta_subset_age = subset(expr_meta_subset, age %in% a)
  Ttest = t.test(ACE2 ~ two_group, expr_meta_subset_age, alternative = "two.sided")
  paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
         "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
         "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
         "n: ", nrow(expr_meta_subset_age)
  )
  })

  
pdf(paste("figure/treehouse_ACE2_ages.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((ACE2)  ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(-1,6),what = c(T, T, T, F),
         border = NA, col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1:length(unique(expr_meta_subset$age))), y = 5, labels = test)
dev.off()
#}

#############################

##
#plot ESR1 vs ACE2
for (gender in c("female", "male", "all")){
  
for (i in 1:length(symbols)){
gene = symbols[i]
expr_meta_subset = expr_meta[expr_meta_subset$type %in% c("normal", "adjacent", "cancer"), ]

if (gender == "all"){
  expr_meta_subset =  expr_meta_subset #[expr_meta$gender %in% c("female", "male"), ]
}else{
  expr_meta_subset =  expr_meta_subset[expr_meta_subset$gender %in% gender, ]
}

#build on normal distribution, may not be appropriate
cutoff = sort(expr_meta_subset[,gene], decreasing = T)[round(nrow(expr_meta_subset) * 0.1)] #qnorm(0.01, mean(expr_meta_subset[,gene]), sd(expr_meta_subset[, gene]), lower.tail=F)  #quantile(expr_meta[, "ACE2"])[4]

expr_meta_subset$gene = expr_meta_subset[, gene]
expr_meta_subset$expr = "high"
expr_meta_subset$expr[expr_meta_subset[, gene] < cutoff] = "low"

if (sum(table(expr_meta_subset$expr) > 10) == 1) next

Ttest = t.test(expr_meta_subset[, "ACE2"] ~ expr_meta_subset$expr)

pdf(paste("figure/treehouse_ACE2_", symbols[i], gender, ".pdf", sep="_"))
print(ggplot(expr_meta_subset, aes(x= gene, y = ACE2, colour = expr )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
   geom_point(size=0.5) + 
  annotate("text", label = paste( "p = ", format(Ttest$p.value, digit=2), ", ratio = ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=3), ":1", sep=""),  x = 3, y = 4, size = 6) + 
  xlab(symbols[i]) + guides(shape=FALSE, size=FALSE) +
  ylab("ACE2") + coord_cartesian(xlim = c(0, 4), ylim=c(0, 4)) 
)
dev.off()
}
}

##############
#tissue
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
expr_meta_subset = subset(expr_meta_subset, two_group %in% names(group_patients[group_patients > 30]))

#find group statistics
tissues = as.character(sort(unique(expr_meta_subset$tissue)))
test = sapply(tissues, function(t){
  expr_meta_subset_age = subset(expr_meta_subset, tissue == t)
  if (sum(table(expr_meta_subset_age$two_group) > 10) > 1){
    Ttest = t.test(ACE2 ~ two_group, expr_meta_subset_age, alternative = "two.sided")
    paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
           "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
           "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
           "n: ", nrow(expr_meta_subset_age)
    )
    }else{
    NA
  }
})


pdf(paste("figure/treehouse_ACE2_tissue_30.pdf", sep="_"), width = 20)
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((ACE2 )  ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(-1,6), what = c(T, T, T, F),
         border = NA, col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1:length(unique(expr_meta_subset$tissue))), y = 5, labels = test)
dev.off()

