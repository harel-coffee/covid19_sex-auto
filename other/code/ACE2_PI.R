TCGA_PI = read.csv("~/Downloads/TCGA_PI.csv", stringsAsFactors = F)
TCGA_PI$type = sapply(TCGA_PI$Sample.ID, function(x) paste(unlist(strsplit(x, "-"))[4], collapse = "."))
TCGA_PI$Sample.ID = sapply(TCGA_PI$Sample.ID, function(x) paste(unlist(strsplit(x, "-"))[1:3], collapse = "."))

TCGA_PI = TCGA_PI[TCGA_PI$type %in% c("01A", "01B"), ]

expr_meta$Sample.ID = sapply(expr_meta$sample_id, function(x) paste(unlist(strsplit(x, "\\."))[1:3], collapse = "."))
TCGA_PI_ACE2 = merge(TCGA_PI, expr_meta, by.x= "Sample.ID", by.y = "Sample.ID")

TCGA_PI_ACE2_subset = TCGA_PI_ACE2[TCGA_PI_ACE2$ENSG00000130234 > 2, ]

sapply(unique(TCGA_PI_ACE2$Cancer.type), function(x){
  TCGA_PI_ACE2_subset = TCGA_PI_ACE2[TCGA_PI_ACE2$Cancer.type %in% x, ]
  cor = cor.test(TCGA_PI_ACE2_subset$ENSG00000130234, TCGA_PI_ACE2_subset$PI.score )
  
  plot(TCGA_PI_ACE2_subset$ENSG00000130234, TCGA_PI_ACE2_subset$PI.score, main = paste(x , cor$p.value, cor$estimate))
})

