#
setwd("/Users/ChenB1/Documents/stanford/wars/gender/data/")
library(stringr)
load("raw/human_gsm_meta.rda")

gsm_gender = NULL


output = paste("gsm", "gender", "tissue", "age", "source_name",
               "gse" , "title" , "characteristics" , sep = "\t")
write(output, "GEO_gsm_all.txt", append = T)
for (i in 1:length(gsmMeta)){
  print(i)
  gsm = gsmMeta[[i]]$Sample_geo_accession
  
  gender = NA
  if (length(grep("male", gsmMeta[[i]]$Sample_characteristics_ch1, ignore.case = T)) > 0){
     if (length(grep("female", gsmMeta[[i]]$Sample_characteristics_ch1, ignore.case = T)) > 0){
       gender = "female"
     }else{
       gender = "male"
     }
  }
  
  tissue = NA
  if (length(grep("tissue:", gsmMeta[[i]]$Sample_characteristics_ch1, ignore.case = T)) > 0){
      tissue = gsmMeta[[i]]$Sample_characteristics_ch1[grep("tissue:", gsmMeta[[i]]$Sample_characteristics_ch1, ignore.case = T)[1]]
  }

  age = NA
  if (length(grep("age:", gsmMeta[[i]]$Sample_characteristics_ch1, ignore.case = T)) > 0){
    age = gsmMeta[[i]]$Sample_characteristics_ch1[grep("age:", gsmMeta[[i]]$Sample_characteristics_ch1, ignore.case = T)[1]]
  }
  
  output = paste(gsm, gender, tissue, age, source_name = gsmMeta[[i]]$Sample_source_name_ch1, 
                 gse =  gsmMeta[[i]]$Sample_series_id, title = gsmMeta[[i]]$Sample_title, characteristics = paste(gsmMeta[[i]]$Sample_characteristics_ch1, collapse = "; "), sep = "\t")
  write(output, "GEO_gsm_all.txt", append = T)
}

#write.csv(gsm_gender, "GEO_gsm_all.csv")

#####
#rs <- dbGetQuery(con,'select gsm, title, characteristics_ch1 from gsm where gpl == "GPL570" and LOWER(characteristics_ch1) LIKE "%male%" ')

geo_gsm_all = read.csv("GEO_gsm_all.txt", sep = "\t")
geo_gsm_all$tissue = tolower(geo_gsm_all$tissue)
write.csv(geo_gsm_all, "geo_gsm_all.csv")
