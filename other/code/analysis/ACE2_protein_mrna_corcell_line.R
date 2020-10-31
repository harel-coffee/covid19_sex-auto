load("~/Downloads/octad_cell_line/OCTAD_cell_line.RData")
table(octad_cell_line_features$type)

as.character(octad_cell_line_features[octad_cell_line_features$type == "protein","id"])

plot(octad_cell_line_matrix[,"expression_ACE2"], 
     octad_cell_line_matrix[,"expression_ESR1"])

mRNA_protein = format(cor(octad_cell_line_matrix[,"expression_ESR1"], 
           octad_cell_line_matrix[,"protein_ER-alpha"], use = "pairwise.complete.obs"), digit=2) 

pdf("figure/ESR1_mRNA_protein.pdf")
plot(octad_cell_line_matrix[,"expression_ESR1"], 
     octad_cell_line_matrix[,"protein_ER-alpha"], xlab = "ESR1 mRNA expression", ylab = "ER-alpha protein",
     main = paste("cor:", mRNA_protein))
dev.off()

pdf("figure/AR_mRNA_protein.pdf")
mRNA_protein = format(cor(octad_cell_line_matrix[,"expression_AR"], 
                          octad_cell_line_matrix[,"protein_AR"], use = "pairwise.complete.obs"), digit=2) 

plot(octad_cell_line_matrix[,"expression_AR"], 
     octad_cell_line_matrix[,"protein_AR"], xlab = "AR mRNA expression", ylab = "AR protein",
     main = paste("cor:", mRNA_protein))
dev.off()



plot(octad_cell_line_matrix[,"expression_ESR1"], 
     octad_cell_line_matrix[,"protein_ER-alpha_pS118"])
cor(octad_cell_line_matrix[,"expression_ESR1"], 
    octad_cell_line_matrix[,"protein_ER-alpha_pS118"], use = "pairwise.complete.obs")


(cor(octad_cell_line_matrix[,"expression_ACE2"], 
    octad_cell_line_matrix[,as.character(octad_cell_line_features[octad_cell_line_features$type == "sensitivity","id"])
], use = "pairwise.complete.obs"))
head(sort(cor(octad_cell_line_matrix[,"gene_effect_ACE2"], 
     octad_cell_line_matrix[,as.character(octad_cell_line_features[octad_cell_line_features$type == "sensitivity","id"])
                            ], use = "pairwise.complete.obs")[1,]))
tail(sort(cor(octad_cell_line_matrix[,"gene_effect_ACE2"], 
              octad_cell_line_matrix[,as.character(octad_cell_line_features[octad_cell_line_features$type == "sensitivity","id"])
                                     ], use = "pairwise.complete.obs")[1,]))
plot(octad_cell_line_matrix[,"expression_ACE2"], 
     octad_cell_line_matrix[,"protein_E-Cadherin"])
