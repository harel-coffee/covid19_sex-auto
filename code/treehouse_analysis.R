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
treehouse_meta$age_group = NA
treehouse_meta$age_group[treehouse_meta$age_in_years %in% c("0-1", "1-9", "10-19")] = "0-19"
treehouse_meta$age_group[treehouse_meta$age_in_years %in% c("20-29", "30-39", "40-49")] = "20-49"
treehouse_meta$age_group[treehouse_meta$age_in_years %in% c("50-59", "60-69", "70-79", "80-100")] = "50-100"

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

genes = c("ENSG00000130234", "ENSG00000091831", "ENSG00000197635", "ENSG00000159640", "ENSG00000164850", "ENSG00000140009", "ENSG00000169083", "ENSG00000082175", "ENSG00000137869")
symbols = c("ACE2", "ESR1", "DPP4", "ACE1", "GPER1", "ESR2", "AR", "PGR", "CYP19A1")

expr = t(all_exprs[genes, ])
expr_meta = merge(treehouse_meta, cbind(expr,Y_expr) , by.x = "sample_id", by.y = 0)

write.csv(expr_meta, paste0("treehouse_expr_meta.csv"))

##################
#check gender genes
boxplot(Y_expr ~ gender, data = expr_meta)
nrow(expr_meta[expr_meta$Y_expr > 0.5 & expr_meta$gender == "female", ])/nrow(expr_meta)
nrow(expr_meta[expr_meta$Y_expr < 0.5 & expr_meta$gender == "male", ])/nrow(expr_meta)

#View(expr_meta[expr_meta$Y_expr > 0.5 & expr_meta$gender == "female", ])
#View(expr_meta[expr_meta$Y_expr < 0.5 & expr_meta$gender == "male", ])
###################
table(expr_meta$age_group)
table(expr_meta$type)

#visualize gtex normal
expr_meta_subset = expr_meta[expr_meta$source == "gtex", ]
describe.by(expr_meta_subset[, c("ENSG00000130234", "gender")], group = "gender")
t.test(ENSG00000130234 ~ gender, data = expr_meta)
tissues = table(tolower(expr_meta_subset$tissue))
tissues = names(tissues[tissues > 30])
p = sapply(tissues, function(x){
  expr_meta_subset = expr_meta[tolower(expr_meta$tissue) == x,]
  if (sum(table(expr_meta_subset$gender) > 10) > 1){
    t.test(ENSG00000130234 ~ gender, data = expr_meta_subset)$p.value
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
  Ttest = t.test(ENSG00000130234 ~ two_group, expr_meta_subset_age)
  paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
         "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
         "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
         "n: ", nrow(expr_meta_subset_age)
  )
  })


pdf(paste("figure/treehouse_ACE2_gender.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot(ENSG00000130234  ~ two_group, data = expr_meta_subset, ll = 0.15,
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(-1,6),what = c(T, T, T, F),
         border = NA, innerborder = "gray", col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1:1), y = 5, labels = test)
dev.off()

ages = "all"
test = sapply(ages, function(age){
  expr_meta_subset_age = expr_meta_subset
  Ttest = t.test(ENSG00000197635 ~ two_group, expr_meta_subset_age)
  paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
         "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
         "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
         "n: ", nrow(expr_meta_subset_age)
  )
})


pdf(paste("figure/treehouse_DPP4_gender.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot(ENSG00000197635  ~ two_group, data = expr_meta_subset, ll = 0.15,
         main = "", ylab = "DPP4 Expression", side = "both", ylim=c(-1,6),what = c(T, T, T, F),
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
expr_meta_subset$two_group = paste(expr_meta_subset$age_group, as.numeric(as.factor(expr_meta_subset$gender)))
group_patients = table(expr_meta_subset$two_group)
expr_meta_subset = subset(expr_meta_subset, two_group %in% names(group_patients[group_patients > 30]))

#find group statistics
ages = sort(unique(expr_meta_subset$age_group))
test = sapply(ages, function(age){
  expr_meta_subset_age = subset(expr_meta_subset, age_group == age)
  Ttest = t.test(ENSG00000130234 ~ two_group, expr_meta_subset_age, alternative = "two.sided")
  paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
         "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
         "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
         "n: ", nrow(expr_meta_subset_age)
  )
  })

  
pdf(paste("figure/treehouse_ACE2_ages.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((ENSG00000130234)  ~ two_group, data = expr_meta_subset, ll = 0.15,
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(-1,6),what = c(T, T, T, F),
         border = NA, col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1:length(unique(expr_meta_subset$age_in_years))), y = 5, labels = test)
dev.off()
#}

#############################

##
#plot ESR1 vs ACE2
for (i in 1:length(genes)){
gene = genes[i]
expr_meta_subset = expr_meta[expr_meta_subset$type %in% c("normal", "adjacent", "cancer"), ]

cutoff = qnorm(0.01, mean(expr_meta_subset[,gene]), sd(expr_meta_subset[, gene]), lower.tail=F)  #quantile(expr_meta[, "ENSG00000130234"])[4]

expr_meta_subset$gene = expr_meta_subset[, gene]
expr_meta_subset$expr = "high"
expr_meta_subset$expr[expr_meta_subset[, gene] < cutoff] = "low"

Ttest = t.test(expr_meta_subset[, "ENSG00000130234"] ~ expr_meta_subset$expr)

pdf(paste("figure/treehouse_ACE2_", symbols[i], ".pdf", sep="_"))
print(ggplot(expr_meta_subset, aes(x= gene, y = ENSG00000130234, colour = expr )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
   geom_point(size=0.5) + 
  annotate("text", label = paste( "p = ", format(Ttest$p.value, digit=2), ", ratio = ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=3), ":1", sep=""),  x = 3, y = 4, size = 6) + 
  xlab(symbols[i]) + guides(shape=FALSE, size=FALSE) +
  ylab("ACE2") + coord_cartesian(xlim = c(0, 4), ylim=c(0, 4)) 
)
dev.off()
}

#########
#plot ESR1 vs DPP4
for (i in 1:length(genes)){
gene = gene = genes[i] 
expr_meta_subset = expr_meta[expr_meta$type %in% c("normal", "adjacent", "cancer"), ]

cutoff = qnorm(0.01, mean(expr_meta_subset[,gene]), sd(expr_meta_subset[, gene]), lower.tail=F)  #quantile(expr_meta[, "ENSG00000130234"])[4]

expr_meta_subset$gene = expr_meta_subset[, gene]
expr_meta_subset$expr = "high"
expr_meta_subset$expr[expr_meta_subset[, gene] < cutoff] = "low"

Ttest = t.test(expr_meta_subset[, "ENSG00000197635"] ~ expr_meta_subset$expr)

pdf(paste("figure/treehouse_DPP4_", symbols[i], ".pdf", sep="_"))
print(ggplot(expr_meta_subset, aes(x= gene, y = ENSG00000197635, colour = expr )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
  geom_point(size=0.5) + 
  annotate("text", label = paste( "p = ", format(Ttest$p.value, digit=2), ", ratio = ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=3), sep=""),  x = 3, y = 4, size = 6) + 
  xlab(symbols[i]) + guides(shape=FALSE, size=FALSE) +
  ylab("DPP4") + coord_cartesian(xlim = c(0, 4), ylim=c(0, 4)) 
)
dev.off()
}


##############################
#ACE2 expression with other genes across ages
ACE2_quantile = quantile(expr_meta$ENSG00000130234)[3]
for (age in unique(expr_meta$age_group)){
  for (sex in unique(expr_meta$gender)){
    expr_meta_subset = expr_meta[expr_meta$age_group == age & expr_meta$gender == sex, ]
    print(paste(age, sex))
    print(anova(lm(expr_meta_subset$ENSG00000130234 ~ expr_meta_subset$ENSG00000137869)))
    plot(expr_meta_subset$ENSG00000130234 ~ expr_meta_subset$ENSG00000137869, 
         main = paste(round(cor(expr_meta_subset$ENSG00000130234, expr_meta_subset$ENSG00000137869), 2),
                      age, sex))
    
    Ttest = t.test(expr_meta_subset$ENSG00000137869[expr_meta_subset$ENSG00000130234 > ACE2_quantile],
                   expr_meta_subset$ENSG00000137869[expr_meta_subset$ENSG00000130234 < ACE2_quantile])
    print(paste(age, sex, Ttest$p.value))
  }
}

#gene

#######################
gene = "ENSG00000091831"
cutoff = quantile(expr_meta[, gene])[3]

expr_meta_subset = expr_meta[expr_meta$type %in% c("normal", "adjacent", "cancer"), ]
expr_meta_subset$expr = 1
expr_meta_subset$expr[expr_meta_subset[, gene] < cutoff] = 2

expr_meta_subset$two_group = paste(paste(expr_meta_subset$gender, expr_meta_subset$age_group, sep = "_"), expr_meta_subset$expr)
group_patients = table(expr_meta_subset$two_group)
expr_meta_subset = subset(expr_meta_subset, two_group %in% names(group_patients[group_patients > 20]))

#find group statistics
two_groups = sort(unique(paste(expr_meta_subset$gender, expr_meta_subset$age_group)))
test = sapply(two_groups, function(g){
  expr_meta_subset_age = subset(expr_meta_subset, age_group == unlist(strsplit(g, " "))[2] &
                                  gender == unlist(strsplit(g, " "))[1])
  Ttest = t.test(expr_meta_subset_age$ENSG00000130234[expr_meta_subset_age$expr == 1],
                 expr_meta_subset_age$ENSG00000130234[expr_meta_subset_age$expr == 2])
  paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
         "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
         "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
         "n: ", nrow(expr_meta_subset_age)
  )
  print(Ttest)
})

pdf(paste0("age_gender_", gene, "_ACE2.pdf"), width = 10)
beanplot((ENSG00000130234)  ~ two_group, data = expr_meta_subset, ll = 0.15,
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(-1,6), xlab = gene,what = c(T, T, T, F),
         border = NA, col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "mean", overalllin="mean")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("high", "low"))
text(x=c(1:length(unique(paste(expr_meta_subset$gender, expr_meta_subset$age_group)))), y = 5, labels = test)
dev.off()


##############
#tissue
expr_meta_subset = expr_meta
main_tissues = c("lung", "liver", "pancreas", "blood", "intestine", "brain", "skin", "colon", "kidney")
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
    Ttest = t.test(ENSG00000130234 ~ two_group, expr_meta_subset_age, alternative = "two.sided")
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
beanplot((ENSG00000130234 )  ~ two_group, data = expr_meta_subset, ll = 0.15,
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(-1,6), what = c(T, T, T, F),
         border = NA, col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1:length(unique(expr_meta_subset$tissue))), y = 5, labels = test)
dev.off()

