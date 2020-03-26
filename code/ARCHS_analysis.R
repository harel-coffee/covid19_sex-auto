####
#post processing
setwd("~/Documents/stanford/wars/gender/data")

library("ggplot2")
library("beanplot")

load("~/Downloads/y_pred0.RData")
load("~/Downloads/y_pred.RData")
y_pred0 = y_pred0[, c("score", "pred_gender")]
names(y_pred0) = c("score", "gender")
y = rbind(y_pred[, c("score", "gender")], y_pred0)

expr_meta = read.csv("~/Downloads/human_matrix_expr.csv", stringsAsFactors = F, row.names = 1)
expr_meta = merge(expr_meta, y, by = 0)
expr_meta = expr_meta[expr_meta$ACE2 > 1, ]
hist(log(expr_meta$ACE2 + 0.01))

t.test(ACE2 ~ gender, data = expr_meta)

symbols = c("ACE2", "ESR1", "DPP4",  "GPER1", "ESR2", "AR", "PGR")

###############
##
#
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


pdf(paste("figure/ARCHS_ACE2_gender.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot(log(ACE2 + 0.01)  ~ two_group, data = expr_meta_subset, ll = 0.15,
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(1,15),what = c(T, T, T, F),
         border = NA, innerborder = "gray", col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("topleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1:1), y = 14, labels = test)
dev.off()

##############
for (i in 1:length(symbols)){
gene = symbols[i]
cutoff = quantile(expr_meta[, "ACE2"])[4]

expr_meta_subset =  expr_meta[expr_meta$gender %in% c("female", "male"), ]
expr_meta_subset$gene = expr_meta_subset[, gene]
expr_meta_subset$expr = "high"
expr_meta_subset$expr[expr_meta_subset[, "ACE2"] < cutoff] = "low"

Ttest = t.test(expr_meta_subset[, gene] ~ expr_meta_subset$expr)

pdf(paste("figure/ARCHS_ACE2_", symbols[i], ".pdf", sep="_"))
print(ggplot(expr_meta_subset, aes(x= log(gene + 0.01), y = log(ACE2 + 0.01), colour = expr )) +  theme_bw()  + 
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
  cutoff = quantile(expr_meta[, "DPP4"])[4]
  
  expr_meta_subset =  expr_meta[expr_meta$gender %in% c("female", "male"), ]
  expr_meta_subset$gene = expr_meta_subset[, gene]
  expr_meta_subset$expr = "high"
  expr_meta_subset$expr[expr_meta_subset[, "DPP4"] < cutoff] = "low"
  
  Ttest = t.test(expr_meta_subset[, gene] ~ expr_meta_subset$expr)
  
  pdf(paste("figure/ARCHS_DPP4_", symbols[i], ".pdf", sep="_"))
  print(ggplot(expr_meta_subset, aes(x= log(gene + 0.01), y = log(DPP4 + 0.01), colour = expr )) +  theme_bw()  + 
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
tissue_info = read.csv("GEO_gsm_all.csv")
expr_meta_subset = merge(expr_meta, tissue_info[!is.na(tissue_info$tissue), c("gsm", "tissue")], by.x= "Row.names", by.y = "gsm")
top_tissues = names(tail(sort(table(expr_meta_subset$tissue)), 30))
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
  paste0("p= ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
         format(Ttest$estimate[1] / Ttest$estimate[2], digit=3), ":1" )
  }else{
    NA
  }
})



pdf(paste("figure/ARCHS_ACE2_tissue_30.pdf", sep="_"), width = 20)
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot(log(ACE2 + 0.01)  ~ two_group, data = expr_meta_subset, ll = 0.15,
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(-1,10), what = c(T, T, T, F),
         border = NA, col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "median", overalllin="median")
legend("topleft", fill = c("blue", "red"),
       legend = c("female", "male"))
text(x=c(1:length(unique(expr_meta_subset$tissue))), y = -0.5, labels = test)
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



