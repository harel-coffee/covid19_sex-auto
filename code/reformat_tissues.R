#reformat tissues
setwd("~/Documents/stanford/wars/gender/data")

main_tissues = c("airway", "lung", "liver",  "pancrea", "blood", "intestine", "brain", "colon", "breast",  "skin", "bone", "kidney", "prostate", "gland") #,

geo_rna_seq = read.csv("GEO_gsm_all.csv")

tissue_new = sapply(1:nrow(geo_rna_seq), function(i){
  main_tissue = NA
  for (t in main_tissues){
    if (length(grep(t, paste(geo_rna_seq$tissue[i], geo_rna_seq$source_name[i]), ignore.case = T)) > 0){
      main_tissue = t
      break
    }
  }
  main_tissue
})

geo_rna_seq$tissue_old = geo_rna_seq$tissue
geo_rna_seq$tissue = tissue_new
write.csv(geo_rna_seq, "GEO_gsm_all.csv")

####################
#####################
gpl570 = read.csv("GPL570_meta_all.csv")
tissue_new = sapply(1:nrow(gpl570), function(i){
  main_tissue = NA
  for (t in main_tissues){
    if (length(grep(t, paste(gpl570$tissue[i], gpl570$source_name_ch1[i]), ignore.case = T)) > 0){
      main_tissue = t
      break
    }
  }
  main_tissue
})

gpl570$tissue_old = gpl570$tissue
gpl570$tissue = tissue_new
write.csv(gpl570, "GPL570_meta_all.csv")

####
treehouse = read.csv("treehouse_all_meta.csv")

tissue_new = sapply(1:nrow(treehouse), function(i){
  main_tissue = NA
  for (t in main_tissues){
    if (length(grep(t, paste(treehouse$tissue[i], treehouse$disease[i]), ignore.case = T)) > 0){
      main_tissue = t
      break
    }
  }
  if (is.na(main_tissue)){
    if (length(grep("glio|neuroblastoma|medulloblastoma", treehouse$disease[i])) > 0){
      main_tissue = "brain"
    }
    if (length(grep("leukemia", treehouse$disease[i])) > 0){
      main_tissue = "blood"
    }    
    if (length(grep("hepa|cholan", treehouse$disease[i])) > 0){
      main_tissue = "liver"
    }    
    if (length(grep("renal", treehouse$disease[i])) > 0){
      main_tissue = "kidney"
    }    
  }
  main_tissue
})
treehouse$tissue_old = treehouse$tissue
treehouse$tissue = tissue_new
write.csv(treehouse, "treehouse_all_meta.csv")
