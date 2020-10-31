#given one matrix
#compute ACE2 co-expressed gene with ESR1

phyper.test <- function(gene_list1, gene_list2, background_genes){
  #gene list1: balls drawn, gene_list2: white balls in background
  
  #first step should remove the genes not shown in background, the step is moved to the parent step to increase efficiency
  #gene_list1 <- gene_list1[gene_list1 %in% background_genes]
  #gene_list2 <- gene_list2[gene_list2 %in% background_genes]
  
  q = sum(gene_list1 %in% gene_list2)
  m = length(gene_list2)
  n = length(background_genes) - length(gene_list2)
  k = length(gene_list1)
  # p = 1 - phyper(q,m,n,k)
  p = phyper(q,m,n,k,lower.tail=F,log.p=T)
  
  return(exp(p))
}

require(httr)
require(dplyr)
enrichGeneList <- function (gene.list, databases=c("ChEA_2016", "ENCODE_TF_ChIP-seq_2015"), fdr.cutoff=NULL) {
  ######Step 1: Post gene list to EnrichR
  req.body <- list(list=paste(gene.list, collapse="\n"))
  post.req <- httr::POST("http://amp.pharm.mssm.edu/Enrichr/addList", 
                         encode="multipart", 
                         body=I(req.body))
  ids = content(post.req,type='application/json')
  results = data.frame(matrix(ncol = 6,nrow=0))
  colnames(results) = c('database',
                        'category',
                        'pval_adj',
                        'Odds_Ratio',
                        'Combined Score',
                        'Genes')
  for(db in databases){
    getEnrich = GET('http://amp.pharm.mssm.edu/Enrichr/export',
                    query=list(userListId = ids$userListId,
                               filename='example',
                               backgroundType= db))
    query.results <- readr::read_tsv(content(getEnrich,type='text/tsv'))
    query.results$database = db
    query.results = query.results %>% select(database,
                                             category=Term,
                                             pval_adj=`Adjusted P-value`,
                                             Odds_Ratio = `Odds Ratio`,
                                             `Combined Score`,
                                             Genes)
    results = rbind(results,query.results)
  }
  results
}

setwd("~/Documents/stanford/wars/gender/data")

load("raw/TcgaTargetGtexMet_rsem_gene_tpm.RData")

drug_gene_set_up = read.delim("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.txt",  sep ="\t", header=F, stringsAsFactors = F)
drug_gene_set = read.delim("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.txt",  sep ="$", header=F, stringsAsFactors = F)
all_genes = table(unlist(strsplit(drug_gene_set$V1, "\t")))
all_genes = names(all_genes[all_genes > 2])
  
ESR1_gene = drug_gene_set_up[drug_gene_set_up$V1 ==  "ESR1 CHEA", ][1,]
ESR1_gene = ESR1_gene[ESR1_gene != ""][-1]

AR_gene = drug_gene_set_up[drug_gene_set_up$V1 ==  "AR CHEA", ][1,]
AR_gene = AR_gene[AR_gene != ""][-1]


mapping = read.csv("gencode.v23.annotation.gene.probeMap.csv", stringsAsFactors = F)[, c("ensembl", "gene")]
ESR1_ensembl = mapping$ensembl[mapping$gene %in% ESR1_gene]
AR_ensembl = mapping$ensembl[mapping$gene %in% AR_gene]
all_genes_ensembl = mapping$ensembl[mapping$gene %in% all_genes]

res = enrichGeneList(ESR1_gene)
###
#gtext normal
ESR_target_genes = ESR1_ensembl
expr_meta = read.csv("treehouse_expr_meta_final_v6.csv", stringsAsFactors = F)
expr_meta$tissue = as.character(expr_meta$tissue)
expr_meta[is.na(expr_meta$tissue), "tissue"] = "unknown"
expr_meta$tissue = as.factor(expr_meta$tissue)

expr_meta$tissue_old = tolower(expr_meta$tissue_old)
tissues = sort(table(expr_meta$tissue_old[expr_meta$gender %in% c("female","male")])) #[expr_meta$ == "gtex" ])
tissues = names(tissues[tissues > 30])

tissue_ESR1_AR = NULL

for (tissue in tissues){
background =  all_genes_ensembl #rownames(expr_subet) #names(head(sort(cors[1, ]), 20000))
  
select_samples = expr_meta$sample_id[expr_meta$tissue_old %in% c(tissue) & expr_meta$gender %in% c('female', "male")]
expr_subet = all_exprs[, colnames(all_exprs) %in% select_samples]
cors = cor(t(expr_subet[c("ENSG00000130234", "ENSG00000184012",  "ENSG00000197635", "ENSG00000159640", "ENSG00000091831"), ]), t(expr_subet[rownames(expr_subet) %in% all_genes_ensembl, ]))
top_negative = names(head(sort(cors[2, ]), 200))
top_positive = names(tail(sort(cors[2, ]), 200))

ESR1_ACE2_neg = phyper.test(top_negative, ESR1_ensembl, background )
ESR1_ACE2_pos = phyper.test(top_positive, ESR1_ensembl, background )

AR_ACE2_neg = phyper.test(top_negative, AR_ensembl, background )
AR_ACE2_pos = phyper.test(top_positive, AR_ensembl, background )

tissue_ESR1_AR = rbind(tissue_ESR1_AR, data.frame( ESR1_ACE2_neg, ESR1_ACE2_pos, AR_ACE2_neg, AR_ACE2_pos))
#res = enrichGeneList(mapping$gene[mapping$ensembl %in% top_negative])
#print(head(res))
}

rownames(tissue_ESR1_AR) = tolower(tissues)
colnames(tissue_ESR1_AR) = c("ESR1 TF vs neg co-expressed ACE2", "ESR1 TF vs pos co-expressed ACE2", "AR TF vs neg co-expressed ACE2", "AR TF vs pos co-expressed ACE2")
library(pheatmap)
pheatmap(t(-log(tissue_ESR1_AR, 10)), cellheight = 10, cellwidth = 10, show_colnames=T, legend=T, show_rownames=T, cluster_cols=F, cluster_rows = F, filename = "figure/ACE2_ESR1_AR_female.pdf") #, filename = "ritonavir_dz_cor_dose.pdf, filename = "ritonavir_selected_dz_cor.pdf"

