GSE = "GSE59491"
data <- (getGEO(GSE, destdir=".", GSEMatrix=TRUE)[[1]])
expr_matrix <- exprs(data)
nn <- phenoData(data)@data

case <- as.character(nn$geo_accession[grep("_2", nn$title)])
control <- as.character(nn$geo_accession[grep("_1", nn$title)])

#

platform <- annotation(data)
gpl <- getGEO(platform, destdir=".")

mapping = Table(gpl)

gene1_ids = mapping$ID[mapping$ENTREZ_GENE_ID %in% c(59272)]
gene2_ids = mapping$ID[mapping$ENTREZ_GENE_ID %in% c(2099)]

plot(expr_matrix[gene1_ids, ], expr_matrix[gene2_ids,])

cor.test(expr_matrix[gene1_ids, ], expr_matrix[gene2_ids,])

t.test(expr_matrix[gene2_ids, control], expr_matrix[gene2_ids,case])

mean(expr_matrix[gene1_ids, ])/ mean(expr_matrix[gene2_ids,])




#########
GSE = "GSE46510"
data <- (getGEO(GSE, destdir=".", GSEMatrix=TRUE)[[1]])
expr_matrix <- exprs(data)
nn <- phenoData(data)@data


platform <- annotation(data)
gpl <- getGEO(platform, destdir=".")

mapping = Table(gpl)

gene1_ids = mapping$ID[mapping$SPOT_ID  %in% c(59272)]
gene2_ids = mapping$ID[mapping$SPOT_ID %in% c(2099)]

plot(expr_matrix[gene1_ids, ], expr_matrix[gene2_ids,])

cor.test(expr_matrix[gene1_ids, ], expr_matrix[gene2_ids,])

t.test(expr_matrix[gene2_ids, control], expr_matrix[gene2_ids,case])

mean(expr_matrix[gene1_ids, ])/ mean(expr_matrix[gene2_ids,])
