library(data.table)

data = fread("~/Downloads/Adult-Trachea2_dge.txt", sep = "\t")
data_subset = data.frame(data[data$GENE %in% c("ESR1", "ACE2", "ESR2", "AR", "PGR", "DPP4", "GPER1"), ])
rownames(data_subset) = data_subset$GENE
data_subset = data.frame(t(data_subset[, -1]))

data_subset = data_subset[data_subset$ESR1 > 0 & data_subset$ACE2 > 0, ]


data = fread("~/Downloads/Adult-Lung1_dge.txt", sep = "\t")
data_subset = data.frame(data[data$GENE %in% c("ESR1", "ACE2", "ESR2", "AR", "PGR", "DPP4", "GPER1"), ])
rownames(data_subset) = data_subset$GENE
data_subset = data.frame(t(data_subset[, -1]))
sum(data_subset$ACE2 > 0) / nrow(data_subset)

data = fread("~/Downloads/Adult-Lung2_dge.txt", sep = "\t")
data_subset = data.frame(data[data$GENE %in% c("ESR1", "ACE2", "ESR2", "AR", "PGR", "DPP4", "GPER1"), ])
rownames(data_subset) = data_subset$GENE
data_subset = data.frame(t(data_subset[, -1]))
sum(data_subset$ACE2 > 0) / nrow(data_subset)

data = fread("~/Downloads/Adult-JeJunum2_dge.txt", sep = "\t")
data_subset = data.frame(data) #[data$GENE %in% c("ESR1", "ACE2", "ESR2", "AR", "PGR", "DPP4", "GPER1"), ])
rownames(data_subset) = data_subset$GENE
data_subset = data.frame(t(data_subset[, -1]))
sum(data_subset$ACE2 > 0) / nrow(data_subset)
data_subset = data_subset[data_subset$ACE2 > 0, ]
cors = (sort(cor(data_subset)["ACE2",]))
tail(cors, 30)

write.csv(data.frame(cors), "~/Downloads/ACE2_Adult-JeJunum2.csv")
