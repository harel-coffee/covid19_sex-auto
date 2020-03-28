####
#post processing
setwd("~/Documents/stanford/wars/gender/data")

library("ggplot2")
library("beanplot")

meta = read.csv("GPL570_meta_all.csv", stringsAsFactors = F)
expr = read.csv("GPL570_exp.csv", stringsAsFactors = F)
#duplicated GSM
expr = expr[!duplicated(expr$X),]
expr$gsm = sapply(expr$X, function(x) unlist(strsplit(x, "\\."))[1])

expr_meta = merge( meta,expr, by.x= "gsm", by.y = "gsm")

#expr_meta = expr_meta[expr_meta$ACE2 > 1, ]
hist((expr_meta$ACE2 ))

t.test(ACE2 ~ gender, data = expr_meta)

symbols = c("ACE2", "ESR1", "DPP4",  "GPER1", "ESR2", "AR", "PGR")

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
  Ttest = t.test(ACE2 ~ two_group, expr_meta_subset_age)
  paste0("p= ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
         format(Ttest$estimate[1] / Ttest$estimate[2], digit=3), ":1" )
})


pdf(paste("figure/GPL570_ACE2_gender.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((ACE2 )  ~ two_group, data = expr_meta_subset, ll = 0.15,
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(1,15),what = c(T, T, T, F),
         border = NA, innerborder = "gray", col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("topleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1:1), y = 14, labels = test)
dev.off()

ages = "all"
test = sapply(ages, function(age){
  expr_meta_subset_age = expr_meta_subset
  Ttest = t.test(DPP4 ~ two_group, expr_meta_subset_age)
  paste0("p= ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
         format(Ttest$estimate[1] / Ttest$estimate[2], digit=3), ":1" )
})


pdf(paste("figure/GPL570_DPP4_gender.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((DPP4 )  ~ two_group, data = expr_meta_subset, ll = 0.15,
         main = "", ylab = "DPP4 Expression", side = "both", ylim=c(1,15),what = c(T, T, T, F),
         border = NA, innerborder = "gray", col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("topleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1:1), y = 14, labels = test)
dev.off()
####################
#####################
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
  Ttest = t.test(ACE2 ~ two_group, expr_meta_subset_age, alternative = "less")
  paste0("p= ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
         format(Ttest$estimate[1] / Ttest$estimate[2], digit=3), ":1" )
})


pdf(paste("figure/GPL570_ACE2_ages.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((ACE2)  ~ two_group, data = expr_meta_subset, ll = 0.15,
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(1,15),what = c(T, T, T, F),
         border = NA, col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("topleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1:length(unique(expr_meta_subset$age))), y = 2, labels = test)
dev.off()
#}
##################
##############
for (i in 1:length(symbols)){
gene = symbols[i]

expr_meta_subset =  expr_meta #[expr_meta$gender %in% c("female", "male"), ]
cutoff = qnorm(0.01, mean(expr_meta[,gene]), sd(expr_meta[, gene]), lower.tail=F) # quantile(expr_meta[, gene])[4]

expr_meta_subset$gene = expr_meta_subset[, gene]
expr_meta_subset$expr = "high"
expr_meta_subset$expr[expr_meta_subset[, gene] < cutoff] = "low"

Ttest = t.test(expr_meta_subset[, "ACE2"] ~ expr_meta_subset$expr)

pdf(paste("figure/GPL570_ACE2_", symbols[i], ".pdf", sep="_"))
print(ggplot(expr_meta_subset, aes(x= (gene), y = (ACE2 ), colour = expr )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
  geom_point(size=0.5) + 
  annotate("text", label = paste( "p = ", format(Ttest$p.value, digit=2), ", ratio = ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=3), ":1", sep=""),  x = 10, y = 14, size = 6) + 
  xlab(gene) + guides(shape=FALSE, size=FALSE) +
  ylab("ACE2") + coord_cartesian(xlim = c(0, 15), ylim=c(0, 15)) 
)
dev.off()
}

########
#DPP4

for (i in 1:length(symbols)){
  gene = symbols[i]
  cutoff = qnorm(0.01, mean(expr_meta[,gene]), sd(expr_meta[, gene]), lower.tail=F) # quantile(expr_meta[, gene])[4]
  
  expr_meta_subset =  expr_meta[expr_meta$gender %in% c("female", "male"), ]
  expr_meta_subset$gene = expr_meta_subset[, gene]
  expr_meta_subset$expr = "high"
  expr_meta_subset$expr[expr_meta_subset[, gene] < cutoff] = "low"
  
  Ttest = t.test(expr_meta_subset[, "DPP4"] ~ expr_meta_subset$expr)
  
  pdf(paste("figure/GPL570_DPP4_", symbols[i], ".pdf", sep="_"))
  print(ggplot(expr_meta_subset, aes(x= (gene ), y = (DPP4 ), colour = expr )) +  theme_bw()  + 
          theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
          geom_point(size=0.5) + 
          annotate("text", label = paste( "p = ", format(Ttest$p.value, digit=2), ", ratio = ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=3), ":1", sep=""),  x = 10, y = 14, size = 6) + 
          xlab(gene) + guides(shape=FALSE, size=FALSE) +
          ylab("DPP4") + coord_cartesian(xlim = c(0, 15), ylim=c(0, 15)) 
  )
  dev.off()
}


###############
#compare tissue
expr_meta_subset = expr_meta
main_tissues = c("lung", "liver",  "pancrea", "blood", "intestine", "brain", "colon",   "skin", "bone") #,"breast",
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
  paste0("p= ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
         format(Ttest$estimate[1] / Ttest$estimate[2], digit=3), ":1" )
  }else{
    NA
  }
})



pdf(paste("figure/GPL570_ACE2_tissue_30.pdf", sep="_"), width = 20)
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot(ACE2   ~ two_group, data = expr_meta_subset, ll = 0.15,
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(1,15), what = c(T, T, T, F),
         border = NA, col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("topleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1:length(unique(expr_meta_subset$tissue))), y = 2, labels = test)
dev.off()
#}

#####################
#investigate highly expressed ACE2
tissue_info = read.csv("GEO_gsm_all.csv")
expr_meta_subset = merge(expr_meta, tissue_info[, c("gsm", "gse", "tissue")], by.x= "Row.names", by.y = "gsm")
tail(sort(by(expr_meta_subset$ACE2, expr_meta_subset$gse, median)), 10)
expr_meta_subset = expr_meta_subset[expr_meta_subset$gender == "female" & expr_meta_subset$tissue %in%  c("tissue: blood","tissue: cord blood"),]
t.test(ACE2 ~ tissue, data = expr_meta_subset)

expr_meta_subset = expr_meta_subset[ expr_meta_subset$tissue %in%  c("tissue: bronchial brushing"),]
t.test(ACE2 ~ gender, data = expr_meta_subset)



