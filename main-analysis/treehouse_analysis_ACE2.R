#analyze ACE2 from the treehouse project
library(ggplot2)
library(psych)
library(beanplot)
library(cowplot)
library(stringr)

rm(list=ls())
root = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root)
fn = paste0(root, '/raw/treehouse_expr_meta_final_v8.csv')
expr_meta = read.csv(fn)
expr_meta = expr_meta[!is.na(expr_meta$gender) & expr_meta$ACE2 > 0.5, ]
symbols = c("ACE2", "ESR1", "DPP4",  "GPER1", "ESR2", "AR", "PGR", "TMPRSS2")

expr_meta$tissue = as.character(expr_meta$tissue)
expr_meta[is.na(expr_meta$tissue), "tissue"] = "unknown"
expr_meta$tissue = as.factor(expr_meta$tissue)
expr_meta$gender = factor(expr_meta$gender, levels = c("male", "female"))

expr_meta = expr_meta[expr_meta$type %in% c("normal", "adjacent"), ]

nrow(expr_meta)
length(which(!is.na(expr_meta$gender)))
length(which(!is.na(expr_meta$gender)))/nrow(expr_meta)
length(which(!is.na(expr_meta$age)))
length(which(!is.na(expr_meta$age)))/nrow(expr_meta)
length(which(!is.na(expr_meta$tissue)))
length(which(!is.na(expr_meta$tissue)))/nrow(expr_meta)




###############
paste("# samples", nrow(expr_meta))
paste("% predicted gender", sum(table(expr_meta$gender))/nrow(expr_meta))
paste("% labeled tissue", sum(table(expr_meta$tissue_old))/nrow(expr_meta))
paste("% predicted tissue", sum(table(expr_meta$tissue))/nrow(expr_meta))
paste("% predicted age", sum(table(expr_meta$age))/nrow(expr_meta))
################


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
    df = aggregate(ACE2~two_group, data = expr_meta_subset_age_bootstrap, FUN=mean) #mean
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


pdf(paste("figure/treehouse_ACE2_gender.pdf", sep="_"))
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot(ACE2  ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(-1,6),what = c(T, T, T, F),
         border = NA, innerborder = "gray", col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "mean", overalllin="mean")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("male", "female"))
text(x=c(1:1), y = 5, labels = test)
dev.off()

########################################
## Part 2 ## age specific gender difference
########################################
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

  
pdf(paste("figure/treehouse_ACE2_ages.pdf", sep="_"), width=10)
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

########################################
## Part 4 ## hormone gene
########################################
#plot ESR1 vs ACE2
symbols = c("ESR1", "ESR2", "PGR", "AR", "GPER1")
for (gender in c("female", "male")){
  print(gender)
for (i in 1:length(symbols)){
gene = symbols[i]
print(gene)
expr_meta_subset = expr_meta[expr_meta$type %in% c("normal", "adjacent", "cancer"), ]

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

Wtest = wilcox.test(expr_meta_subset[, "ACE2"] ~ expr_meta_subset$expr,  conf.int = TRUE)
df = aggregate(expr_meta_subset[, "ACE2"] ~ expr_meta_subset$expr, FUN = median)

# text = paste0(" ", nrow(expr_meta_subset), ", ",
#               format(df[1, 2] / df[2, 2], digit=3), ", ", 
#               format(Wtest$p.value, digits = 2, scientific=T), ", "
# )
# print(text)

print(c(df[1, 2], df[2, 2], df[1, 2] / df[2, 2]))
df = aggregate(expr_meta_subset[, "ACE2"] ~ expr_meta_subset$expr, FUN = mean)
print(c(df[1, 2], df[2, 2], df[1, 2] / df[2, 2]))


#pdf(paste("figure/treehouse_ACE2_", symbols[i], gender, ".pdf", sep="_"))
print(ggplot(expr_meta_subset, aes(x= (gene), y = (ACE2 ), colour = expr )) +  theme_bw()  +
        theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +
        geom_point(size=0.5) +
        annotate("text", label = paste( "p = ", format(Wtest$p.value, digit=2), ", ratio = ", format(df[1, 2] / df[2, 2], digit=3), ":1", sep=""),  x = 10, y = 14, size = 6) +
        xlab(gene) + guides(shape=FALSE, size=FALSE) +
        ylab("ACE2") + coord_cartesian(xlim = c(0, 6), ylim=c(0, 6))
)
#dev.off()

}
}

########################################
## Part 3 ## tissue specific gender diff
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
  "brain 1", "brain 2", "breast 1", "breast 2", "colon 1", "colon 2", "heart 1", "heart 2",  "intestine 1", 
  "intestine 2", "kidney 1", "kidney 2", "liver 1", "liver 2", "lung 1", "lung 2",
  "pancreas 1", "pancreas 2", "prostate 1", "prostate 2",  "skin 1", "skin 2", "testis 1", "testis 2", "unknown 1", "unknown 2"))


group_patients = table(expr_meta_subset$two_group)
expr_meta_subset = subset(expr_meta_subset, two_group %in% names(group_patients[group_patients > 20])) # intestine only has 23 for female

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
      df = aggregate(ACE2~two_group, data = expr_meta_subset_age_bootstrap, FUN=mean) #mean
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


pdf(paste("figure/treehouse_ACE2_tissue_20.pdf", sep="_"), width = 30)
# par(mar = c(14, 4, 4, 2) + 0.1, las=2)
beanplot((ACE2 )  ~ two_group, data = expr_meta_subset, ll = 0.15, log = "",
         main = "", ylab = "ACE2 Expression", side = "both", ylim=c(-1,15), what = c(T, T, T, F),
         border = NA, col = list(c("blue","blue"), c("red", "red")), 
         beanlines = "mean", overalllin="mean")
legend("bottomleft", fill = c("blue", "red"),
       legend = c("male", "female"))
text(x=c(1:length(unique(expr_meta_subset$tissue))), y = 14, labels = test)
dev.off()


############# generate table ######################
expr_meta_subset = expr_meta[!is.na(expr_meta$age) & !is.na(expr_meta$gender), ]
expr_meta_subset$two_group = paste(expr_meta_subset$age, as.numeric(as.factor(expr_meta_subset$gender)))
group_patients = table(expr_meta_subset$two_group)
expr_meta_subset = subset(expr_meta_subset, two_group %in% names(group_patients[group_patients > 30]))

#find group statistics
ages = sort(unique(expr_meta_subset$age))
ages = c('20-59', '60-100')
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


##################### version 2 : transpose   #########################
library(binr) # ACE2, age 0-19 all male, excluded
ages = c("all") # "20-59", "60-100",
for (age_group in ages){
  for (target_gene in c("ESR1")){  # "ACE2", "TMPRSS2", "ESR2", "AR", "PGR"
    cdata_all = NULL
    
    #target_gene = "TMPRSS2" #TMPRSS2
    expr_meta = expr_meta[!is.na(expr_meta[, target_gene]),]
    expr_meta$bin = cut(expr_meta[, target_gene], breaks=c(quantile(expr_meta[, target_gene], probs = c(0.0, 0.7, 0.98, 1.0))), 
                        labels=c("Normal","High", "Very High"), include.lowest=TRUE) # bin(expr_meta[, target_gene], 3)
    
    symbols = c("ACE2", "TMPRSS2", "DPP4") # "ESR1", "ESR2", "AR", "PGR"
    for (gender in c("female", "male")){#"female", "male"
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
    cdata_all = cdata_all[cdata_all$gender %in% c("female", "male"), ]# "male", "female" # "all"
    cdata_all$gender_gene = paste(cdata_all$gender, cdata_all$gene)
    print(cdata_all)
    #pdf(paste("figure/bin_plot/GPL570", target_gene, age_group, "bin.pdf", sep="_"))
    print(
      ggplot(cdata_all[which(cdata_all$bin != 'High'), ], aes(x=bin, y=mean, group = gender_gene, color = gene, shape = gender)) + coord_cartesian(ylim=c(3, 6))  +  theme_bw() +
        geom_point(position=position_dodge(), stat="identity", size = 4) + geom_line(aes(linetype=gender)) + 
        geom_errorbar(aes(ymin=mean, ymax=mean+se),
                      width=.2) +
        xlab(target_gene) +
        ylab("Expression") + coord_cartesian(ylim=c(0, 4)) 
    ) 
    #dev.off()
    
  }
}

