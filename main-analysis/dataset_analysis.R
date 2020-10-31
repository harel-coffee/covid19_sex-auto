expr_meta = read.csv("GPL570_expr_meta_final_v6.csv")
expr_meta = expr_meta[!is.na(expr_meta$gender), ]
expr_meta$gender = factor(expr_meta$gender, levels = c("female", "male"))
ACE_cutoff =  sort(expr_meta$ACE2)[length((expr_meta$ACE2)) * 0.9] #quantile(GPL570_expr_meta$ACE2)[4]

expr_meta$ACE2_bin = sapply(1:nrow(expr_meta), function(i) {
  if (expr_meta$ACE2[i] > ACE_cutoff){
    1
  }else{
    0
  }
})

gender_sets = data.frame(table(expr_meta$series_id, expr_meta$gender))

gender_sets = gender_sets[gender_sets$Freq > 5, ]

gse = table(gender_sets$Var1)
gse = names(gse[gse> 1])

gender_sets = gender_sets[gender_sets$Var1 %in% gse, ]

results = NULL
for (g in gse){
  expr_meta_subset = expr_meta[expr_meta$series_id %in% g, ]
  
  if (sum((table(expr_meta_subset[, "ACE2_bin"], expr_meta_subset[, "gender"]) > 1)) != 4) next
  
  mylogit = (glm(ACE2_bin ~ gender, expr_meta_subset, family = binomial()))
  odds = exp(cbind(OR = coef(mylogit), confint(mylogit)))
  results = rbind(results, data.frame(g, OR = round(odds["gendermale", 1], 2), CI_25 =  round(odds["gendermale", 2], 2) , 
                                      CI_75 =  round(odds["gendermale", 3],2) , p =  summary(mylogit)$coefficients["gendermale", 4]))
  #Wtest = wilcox.test(ACE2 ~ gender, expr_meta_subset, conf.int = T) #(M-F)
  #Ttest = t.test(ACE2 ~ gender, expr_meta_subset) #(M-F)

  #results = rbind(results, data.frame(g, diff = Wtest$estimate, p = Wtest$p.value, diff_t = Ttest$estimate[1] - Ttest$estimate[2], diff_p = Ttest$p.value))
}

results$q = p.adjust(results$p)
sum(results$diff[results$p < 0.05] > 0)
sum(results$diff[results$p < 0.05]  < 0)

sum(results$diff_t[results$diff_p < 0.05] > 0)
sum(results$diff_t[results$diff_p < 0.05]  < 0)


results = results[order(results$diff_t),]

dim(results[results$OR > 1 & results$p < 0.05,])
dim(results[results$OR < 1 & results$p < 0.05,])

###ARCHS

#

#fix gender age tissue