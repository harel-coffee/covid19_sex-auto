setwd("~/Documents/stanford/wars/cmap/data//")

drug_set_name = "Disease_Perturbations_from_GEO"
drug_gene_set_up = read.delim(paste0( drug_set_name, "_up.txt"),  sep ="$", header=F, stringsAsFactors = F)
drug_gene_set_down = read.delim(paste0(drug_set_name, "_down.txt"),  sep ="$", header=F, stringsAsFactors = F)


ACE2_up = NULL
ACE2_down = NULL

for (drug_id in 1:nrow(drug_gene_set_up)){  #nrow(drug_gene_set) nrow(drug_gene_set_up
  drugs_up = unlist(strsplit(drug_gene_set_up[drug_id,1], "\t\t"))
  drug_name_up = drugs_up[1]
  
  drugs_down = unlist(strsplit(drug_gene_set_down[drug_id,1], "\t\t"))
  drug_name_down = drugs_down[1]
  
  drug_genes_up = as.character(sapply(unlist(strsplit(drugs_up[2], "\t")), function(x) {unlist(strsplit(x, ","))[1]}))
  drug_genes_up_value = as.character(sapply(unlist(strsplit(drugs_up[2], "\t")), function(x) {unlist(strsplit(x, ","))[2]}))
  
  drug_genes_down = as.character(sapply(unlist(strsplit(drugs_down[2], "\t")), function(x) {unlist(strsplit(x, ","))[1]}))
  drug_genes_down_value = as.character(sapply(unlist(strsplit(drugs_down[2], "\t")), function(x) {unlist(strsplit(x, ","))[2]}))
  
  if (drug_name_up != drug_name_down){
    stop(paste(drug_id, " name inconsitent"))
  }else{
    drug = drug_name_up
  }
  dataset = str_replace_all(drug, " ", "_")
  drug = str_trim(str_replace_all(unlist(strsplit(dataset, "_"))[1], ",", ""))
 
  if (sum(drug_genes_down %in% c("ACE2"))>0)  ACE2_down = c(ACE2_down, dataset)
  if (sum(drug_genes_up %in% c("ACE2"))>0)  ACE2_up = c(ACE2_up, dataset) #TMPRSS2
  
}


