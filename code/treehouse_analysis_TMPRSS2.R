#analyze TMPRSS2 from the treehouse project
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
treehouse_meta$age_group[treehouse_meta$age_in_years %in% c("20-29", "30-39", "40-49","50-59")] = "20-59"
treehouse_meta$age_group[treehouse_meta$age_in_years %in% c( "60-69", "70-79", "80-100")] = "60-100"

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

genes = c("ENSG00000184012", "ENSG00000091831", "ENSG00000197635", "ENSG00000159640", "ENSG00000164850", "ENSG00000140009", "ENSG00000169083", "ENSG00000082175", "ENSG00000137869", "ENSG00000184012")
symbols = c("TMPRSS2", "ESR1", "DPP4", "ACE1", "GPER1", "ESR2", "AR", "PGR", "CYP19A1", "TMPRSS2")

expr = t(all_exprs[genes, ])
expr_meta = merge(treehouse_meta, cbind(expr,Y_expr) , by.x = "sample_id", by.y = 0)
expr_meta$tissue = expr_meta$tissue.y


paste("# samples", nrow(expr_meta))
paste("% predicted gender", sum(table(expr_meta$gender))/nrow(expr_meta))
paste("% labeled tissue", sum(table(expr_meta$tissue_old))/nrow(expr_meta))
paste("% predicted tissue", sum(table(expr_meta$tissue))/nrow(expr_meta))
paste("% predicted age", sum(table(expr_meta$age_group))/nrow(expr_meta))


expr_meta = expr_meta[expr_meta$ENSG00000184012 > 0.1, ]

paste("# samples", nrow(expr_meta))

#write.csv(expr_meta, paste0("treehouse_expr_meta.csv"))
###############

###################
table(expr_meta$age_group)
table(expr_meta$type)

#visualize gtex normal
expr_meta_subset = expr_meta[expr_meta$source == "gtex", ]
#describe.by(expr_meta_subset[, c("ENSG00000184012", "gender")], group = "gender")
t.test(ENSG00000184012 ~ gender, data = expr_meta)
tissues = table(tolower(expr_meta_subset$tissue))
tissues = names(tissues[tissues > 30])
p = sapply(tissues, function(x){
  expr_meta_subset = expr_meta[tolower(expr_meta$tissue) == x,]
  if (sum(table(expr_meta_subset$gender) > 10) > 1){
    t.test(ENSG00000184012 ~ gender, data = expr_meta_subset)$p.value
  }else{
    NA
  }
})

#+ geom_jitter()

##########################
#plot TMPRSS2 expression between male and female using all data
expr_meta_subset = expr_meta[!is.na(expr_meta$gender), ] 
expr_meta_subset$two_group = paste("all", as.numeric(as.factor(expr_meta_subset$gender)))
group_patients = table(expr_meta_subset$two_group)
expr_meta_subset = subset(expr_meta_subset, two_group %in% names(group_patients[group_patients > 30]))

#find group statistics
#ignore age
ages = "all"
test = sapply(ages, function(age){
  expr_meta_subset_age = expr_meta_subset
  Ttest = t.test(ENSG00000184012 ~ two_group, expr_meta_subset_age)
  paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
         "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
         "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
         "n: ", nrow(expr_meta_subset_age)
  )
  })


pdf(paste("figure/treehouse_TMPRSS2_gender.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot(ENSG00000184012  ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "TMPRSS2 Expression", side = "both", ylim=c(-1,6),what = c(T, T, T, F),
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
beanplot(ENSG00000197635  ~ two_group, data = expr_meta_subset, ll = 0.15,  log = "",
         main = "", ylab = "DPP4 Expression", side = "both", ylim=c(-1,6),what = c(T, T, T, F),
         border = NA, innerborder = "gray", col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1:1), y = 5, labels = test)
dev.off()
########################
#########################
#plot TMPRSS2 express across age
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
  Ttest = t.test(ENSG00000184012 ~ two_group, expr_meta_subset_age, alternative = "two.sided")
  paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
         "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
         "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
         "n: ", nrow(expr_meta_subset_age)
  )
  })

  
pdf(paste("figure/treehouse_TMPRSS2_ages.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((ENSG00000184012)  ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "TMPRSS2 Expression", side = "both", ylim=c(-1,6),what = c(T, T, T, F),
         border = NA, col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1:length(unique(expr_meta_subset$age_in_years))), y = 5, labels = test)
dev.off()
#}

#############################

##
#plot ESR1 vs TMPRSS2
for (gender in c("female", "male", "all")){
  
for (i in 1:length(genes)){
gene = genes[i]
expr_meta_subset = expr_meta[expr_meta_subset$type %in% c("normal", "adjacent", "cancer"), ]

if (gender == "all"){
  expr_meta_subset =  expr_meta_subset #[expr_meta$gender %in% c("female", "male"), ]
}else{
  expr_meta_subset =  expr_meta_subset[expr_meta_subset$gender %in% gender, ]
}

#build on normal distribution, may not be appropriate
cutoff = qnorm(0.01, mean(expr_meta_subset[,gene]), sd(expr_meta_subset[, gene]), lower.tail=F)  #quantile(expr_meta[, "ENSG00000184012"])[4]

expr_meta_subset$gene = expr_meta_subset[, gene]
expr_meta_subset$expr = "high"
expr_meta_subset$expr[expr_meta_subset[, gene] < cutoff] = "low"

if (sum(table(expr_meta_subset$expr) > 10) == 1) next

Ttest = t.test(expr_meta_subset[, "ENSG00000184012"] ~ expr_meta_subset$expr)

pdf(paste("figure/treehouse_TMPRSS2_", symbols[i], gender, ".pdf", sep="_"))
print(ggplot(expr_meta_subset, aes(x= gene, y = ENSG00000184012, colour = expr )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
   geom_point(size=0.5) + 
  annotate("text", label = paste( "p = ", format(Ttest$p.value, digit=2), ", ratio = ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=3), ":1", sep=""),  x = 3, y = 4, size = 6) + 
  xlab(symbols[i]) + guides(shape=FALSE, size=FALSE) +
  ylab("TMPRSS2") + coord_cartesian(xlim = c(0, 4), ylim=c(0, 4)) 
)
dev.off()
}
}

#########
#plot ESR1 vs DPP4
for (gender in c("female", "male", "all")){
  
for (i in 1:length(genes)){
gene = genes[i] 
expr_meta_subset = expr_meta[expr_meta$type %in% c("normal", "adjacent", "cancer"), ]

if (gender == "all"){
  expr_meta_subset =  expr_meta_subset #[expr_meta$gender %in% c("female", "male"), ]
}else{
  expr_meta_subset =  expr_meta_subset[expr_meta_subset$gender %in% gender, ]
}

cutoff = qnorm(0.01, mean(expr_meta_subset[,gene]), sd(expr_meta_subset[, gene]), lower.tail=F)  #quantile(expr_meta[, "ENSG00000184012"])[4]

expr_meta_subset$gene = expr_meta_subset[, gene]
expr_meta_subset$expr = "high"
expr_meta_subset$expr[expr_meta_subset[, gene] < cutoff] = "low"

if (sum(table(expr_meta_subset$expr) > 10) == 1) next

Ttest = t.test(expr_meta_subset[, "ENSG00000197635"] ~ expr_meta_subset$expr)

pdf(paste("figure/treehouse_DPP4_", symbols[i], gender, ".pdf", sep="_"))
print(ggplot(expr_meta_subset, aes(x= gene, y = ENSG00000197635, colour = expr )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
  geom_point(size=0.5) + 
  annotate("text", label = paste( "p = ", format(Ttest$p.value, digit=2), ", ratio = ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=3), sep=""),  x = 3, y = 4, size = 6) + 
  xlab(symbols[i]) + guides(shape=FALSE, size=FALSE) +
  ylab("DPP4") + coord_cartesian(xlim = c(0, 4), ylim=c(0, 4)) 
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
    Ttest = t.test(ENSG00000184012 ~ two_group, expr_meta_subset_age, alternative = "two.sided")
    paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
           "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
           "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
           "n: ", nrow(expr_meta_subset_age)
    )
    }else{
    NA
  }
})


pdf(paste("figure/treehouse_TMPRSS2_tissue_30.pdf", sep="_"), width = 20)
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((ENSG00000184012 )  ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "TMPRSS2 Expression", side = "both", ylim=c(-1,6), what = c(T, T, T, F),
         border = NA, col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1:length(unique(expr_meta_subset$tissue))), y = 5, labels = test)
dev.off()

