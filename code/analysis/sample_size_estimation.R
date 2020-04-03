
setwd("~/Documents/stanford/wars/gender/data")

treehouse_meta = read.csv("treehouse_expr_meta.csv")

sig_ratio = NULL
for (i in seq(1000, 10000, 100)){
  ps = sapply(1:100, function(j){
  t.test(ENSG00000130234 ~ gender, treehouse_meta[sample(1:nrow(treehouse_meta), i), ])$p.value
  })
  sig_ratio = c(sig_ratio, sum(ps < 0.01)/100)
}

pdf("figure//treehouse_samples.pdf")
plot(seq(1000, 10000, 100), sig_ratio, xlab = "# samples", ylab = "pass rate")
dev.off()

##########
GPL570_meta = read.csv("GPL570_expr_meta.csv")

sig_ratio = NULL
for (i in seq(1000, 100000, 100)){
  ps = sapply(1:100, function(j){
    t.test(ACE2 ~ gender, GPL570_meta[sample(1:nrow(GPL570_meta), i), ])$p.value
  })
  sig_ratio = c(sig_ratio, sum(ps < 0.01)/100)
}
pdf("figure/GPL570_samples.pdf")
plot(seq(1000, 100000, 100), sig_ratio, xlab = "# samples", ylab = "pass rate")
dev.off()


######
ARCHS_meta = read.csv("ARCHS_expr_meta.csv")

sig_ratio = NULL
for (i in seq(1000, 10000, 100)){
  ps = sapply(1:100, function(j){
    t.test(ACE2 ~ gender, ARCHS_meta[sample(1:nrow(ARCHS_meta), i), ])$p.value
  })
  sig_ratio = c(sig_ratio, sum(ps < 0.01)/100)
}

pdf("figure//treehouse_samples.pdf")
plot(seq(1000, 10000, 100), sig_ratio, xlab = "# samples", ylab = "pass rate")
dev.off()