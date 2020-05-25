#visualize odds ratio
setwd('/Users/msun/Documents/covid-19/main')
forestplot <- function(d, xlab="OR", ylab="Study"){
  require(ggplot2)
  p <- ggplot(d, aes(x=x, y=y, ymin=ylo, ymax=yhi)) +  theme_bw()  + 
    theme(axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
    geom_pointrange() + 
    coord_flip() +
    geom_hline(aes(x=0, yintercept= 0), lty=2) + coord_cartesian( ylim=c(-10, 10))  +
    ylab(xlab) +
    xlab(ylab) #switch because of the coord_flip() above
  return(p)
}

odds = read.csv("output/table3_v2.csv", stringsAsFactors = F)
odds$OR = log(odds$OR)
odds$CI_25 = log(odds$CI_25)
odds$CI_75 = log(odds$CI_75)

d = odds[odds$age_group == "all" & odds$gender == "male" & odds$virus_gene == "ACE2",  c("hormone", "OR", "CI_25", "CI_75")] #need to change
colnames(d) = c("x", "y", "ylo", "yhi")
forestplot(d)

# table 1 visualization
rm(list=ls())
setwd('/Users/msun/Documents/covid-19/main')
thresholds = c(0.75, 0.8, 0.85, 0.9, 0.95)
odds = NULL
for (t in thresholds){
  fn = paste0('output/table1_p_', t, '.csv')
  d = read.csv(fn)
  d$threshold = rep(t, nrow(d))
  odds = rbind(odds, d)
  print(t)
}

# p-value bool
p_boo = rep(NA, nrow(odds))
for (i in 1:nrow(odds)){
  if (odds$p[i] < 0.001){
    p_boo[i]=2
  }else{
    p_boo[i]=1
  }
}
odds$p_boo = p_boo

genes = c('ACE2', 'TMPRSS2', 'ACE2 & TMPRSS2')
# g='ACE2'
# ages = c('all', 'old')
color = c(2, 4)
shape = c(16, 17)
for (g in genes){
  df = odds[which(odds$gene==g & odds$age_group=='all'), ]
  plot(df$threshold, df$OR, col=color[df$dataset], ylim=c(0.8, 1.7), pch=shape[df$p_boo], xlab='threshold', ylab='OR', main=g)
  lines(x=seq(0.75, 0.95, 0.01), y=rep(1,length(seq(0.75, 0.95, 0.01))), lty=2)
  legend('topleft', c('All', 'Healthy', 'All p<0.001', 'Healthy p<0.001'), col=c(2, 4, 2, 4), pch=c(16, 16, 17, 17), bty='n')
  print(g)
}
for (g in genes){
  df = odds[which(odds$gene==g & odds$age_group=='old'), ]
  plot(df$threshold, df$OR, col=color[df$dataset], ylim=c(0.3, 1.7), pch=shape[df$p_boo], xlab='threshold', ylab='OR', main=g)
  lines(x=seq(0.75, 0.95, 0.01), y=rep(1,length(seq(0.75, 0.95, 0.01))), lty=2)
  legend('topleft', c('All', 'Healthy', 'All p<0.001', 'Healthy p<0.001'), col=c(2, 4, 2, 4), pch=c(16, 16, 17, 17), bty='n')
  print(g)
}

# table 3 visualization
rm(list=ls())
setwd('/Users/msun/Documents/covid-19/main')
thresholds = c(0.75, 0.8, 0.85, 0.9, 0.95)
odds = NULL
for (t in thresholds){
  fn = paste0('output/table3_p_', t, '.csv')
  d = read.csv(fn)
  d$threshold = rep(t, nrow(d))
  odds = rbind(odds, d)
  print(t)
}

# p-value bool
p_boo = rep(NA, nrow(odds))
for (i in 1:nrow(odds)){
  if (odds$p[i] < 0.001){
    p_boo[i]=2
  }else{
    p_boo[i]=1
  }
}
odds$p_boo = p_boo

genes = c('ESR1_bin', 'ESR2_bin', 'PGR_bin', 'AR_bin')
color = c(1, 2, 4)
shape = c(16, 17)
for (g in genes){
  df = odds[which(odds$virus_gene=='ACE2_bin' & odds$gender != 'all' & odds$hormone==g & odds$age_group=='all'), ]
  plot(df$threshold, df$OR, col=color[df$gender], ylim=c(0, 4), pch=shape[df$p_boo], xlab='threshold', ylab='OR', main=g)
  lines(x=seq(0.75, 0.95, 0.01), y=rep(1,length(seq(0.75, 0.95, 0.01))), lty=2)
  legend('topleft', c('female', 'male', 'female p<0.001', 'male p<0.001'), col=c(2, 4, 2, 4), pch=c(16, 16, 17, 17), bty='n')
  print(g)
}
for (g in genes){
  df = odds[which(odds$virus_gene=='ACE2_bin' & odds$gender != 'all' & odds$hormone==g & odds$age_group=='old'), ]
  plot(df$threshold, df$OR, col=color[df$gender], ylim=c(0, 4), pch=shape[df$p_boo], xlab='threshold', ylab='OR', main=g)
  lines(x=seq(0.75, 0.95, 0.01), y=rep(1,length(seq(0.75, 0.95, 0.01))), lty=2)
  legend('topleft', c('female', 'male', 'female p<0.001', 'male p<0.001'), col=c(2, 4, 2, 4), pch=c(16, 16, 17, 17), bty='n')
  print(g)
}




