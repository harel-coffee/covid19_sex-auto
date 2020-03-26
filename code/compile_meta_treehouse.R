#compare meta info
#age distribution: 0-1, 1-10, 20-20, 20-30, etc
setwd("~/Documents/stanford/wars/gender/data/")
treehouse_meta = read.csv("raw/treehouse_public_samples_clinical_metadata.2017-09-11.tsv", sep ="\t", stringsAsFactors = F)
treehouse_meta$tissue = NA
treehouse_meta$source = "treehouse"
treehouse_meta$type = "cancer"
treehouse_meta$age_in_years = sapply(treehouse_meta$age_in_years, function(x){
  if (is.na(x)){
    NA
  }else if (x <= 1){
    "0-1"
  }else if (x < 10){
    "1-9"
  }else if (x < 20){
    "10-19"
  }else if (x < 30){
    "20-29"
  }else if (x < 40){
    "30-39"
  }else if (x < 50){
    "40-49"
  }else if (x < 60){
    "50-59"
  }else if (x < 70){
    "60-69"
  }else if (x < 80){
    "70-79"
  }else{
    "80-100"
  }
})

gtex_meta = read.csv("raw/GTEX_phenotype", stringsAsFactors = F, sep = "\t")
gtex_meta_more = read.csv("raw/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", stringsAsFactors = F, sep = "\t")
gtex_meta_all = merge(gtex_meta, gtex_meta_more, by.x = "X_patient", by.y = "SUBJID")
gtex_meta_all = gtex_meta_all[, c("Sample", "AGE", "X_gender", "DTHHRDY", "X_primary_site")]
names(gtex_meta_all) = c("sample_id", "age_in_years", "gender", "disease", "tissue")
gtex_meta_all$source = "gtex"
gtex_meta_all$type = "normal"

CBTTC = read.csv("raw/Participants.txt", stringsAsFactors = F, sep = "\t")
CBTTC_diagnose = read.csv("raw/Histological_Diagnoses.txt", stringsAsFactors = F, sep = "\t")
CBTTC_meta = merge(CBTTC, CBTTC_diagnose, by.x = "xena_sample", by.y = "xena.id")
CBTTC_meta$age = CBTTC_meta$Age.at.Diagnosis..Days./365
CBTTC_meta$age_in_years = sapply(CBTTC_meta$age, function(x){
  if (is.na(x)){
    NA
  }else if (x <= 1){
    "0-1"
  }else if (x < 10){
    "1-9"
  }else if (x < 20){
    "10-19"
  }else if (x < 30){
    "20-29"
  }else if (x < 40){
    "30-39"
  }else if (x < 50){
    "40-49"
  }else if (x < 60){
    "50-59"
  }else if (x < 70){
    "60-69"
  }else if (x < 80){
    "70-79"
  }else{
    "80-100"
  }
})
CBTTC_meta = CBTTC_meta[, c("xena_sample", "age_in_years", "Gender", "Histological.Diagnosis..NCIT.", "Histological.Tumor.Location" )]
names(CBTTC_meta) = c("sample_id", "age_in_years", "gender", "disease", "tissue")
CBTTC_meta$source = "CBTTC"
CBTTC_meta$type = "cancer"

treehouse_all_meta = rbind(treehouse_meta, gtex_meta_all, CBTTC_meta)

#merge with processed ones
#load from other sources
load("TcgaTargetGtexMet_rsem_gene_pheno.RData")
library("stringr")
pheno$sample = str_replace_all(pheno$sample, "\\.", "-")
pheno = pheno[!pheno$sample %in% treehouse_all_meta$sample_id,]
#retrieve age and gender info from TCGA
pheno$patient = sapply(pheno$sample, function(x) {
  if (length(unlist(strsplit(x, "-"))) > 3){
  paste(unlist(strsplit(x, "-"))[1:3], collapse ="-")
  }else{
    x
  }
})
treehouse_all_meta_subset = treehouse_all_meta
treehouse_all_meta_subset$patient =  sapply(treehouse_all_meta_subset$sample_id, function(x) {
  if (length(unlist(strsplit(x, "-"))) > 3){
    paste(unlist(strsplit(x, "-"))[1:3], collapse ="-")
  }else{
    x
  }
})
pheno = merge(pheno, treehouse_all_meta_subset, by = "patient", all.x = T)

pheno = data.frame(sample_id = pheno$sample, age_in_years = pheno$age_in_years, gender= pheno$gender,
                   disease = pheno$cancer, tissue = pheno$biopsy, source = "others", type = pheno$type.x )

treehouse_all_meta = rbind(treehouse_all_meta, pheno )

treehouse_all_meta$gender = tolower(treehouse_all_meta$gender)
treehouse_all_meta$gender[!treehouse_all_meta$gender %in% c("female", "male")] = "others"

table(treehouse_all_meta$gender)
table(treehouse_all_meta$age_in_years)
table(treehouse_all_meta$type)
write.csv(treehouse_all_meta, "treehouse_all_meta.csv")
