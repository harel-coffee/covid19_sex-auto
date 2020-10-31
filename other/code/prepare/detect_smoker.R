#reformat tissues
setwd("~/Documents/stanford/wars/gender/data")


geo_rna_seq = read.csv("GEO_gsm_all.csv")

smoke = sapply(1:nrow(geo_rna_seq), function(i){
   smoker = NA
    if (length(grep("smoke", paste(geo_rna_seq$characteristics[i], geo_rna_seq$source_name[i]), ignore.case = T)) > 0){
      if (length(grep("non| no ", paste(geo_rna_seq$characteristics[i], geo_rna_seq$source_name[i]), ignore.case = T)) > 0){
        smoker = F
      }else{
        smoker = T
      }
   }
  smoker
})

geo_rna_seq$smoke = smoke
write.csv(geo_rna_seq, "GEO_gsm_all.csv")

####################
#####################
gpl570 = read.csv("GPL570_meta_all.csv")
smoke = sapply(1:nrow(gpl570), function(i){
  smoker = NA
  if (length(grep("smoke", paste(gpl570$characteristics_ch1[i], gpl570$source_name_ch1[i]), ignore.case = T)) > 0){
    if (length(grep("non| no ", paste(gpl570$characteristics_ch1[i], gpl570$source_name_ch1[i]), ignore.case = T)) > 0){
      smoker = F
    }else{
      smoker = T
    }
  }
  smoker
})

gpl570$smoke = smoke
write.csv(gpl570, "GPL570_meta_all.csv")

