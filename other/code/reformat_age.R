#reformat tissues
setwd("~/Documents/stanford/wars/gender/data")

library("stringr")

geo_rna_seq = read.csv("GEO_gsm_all.csv")

ages = sapply(geo_rna_seq$age, function(age){
    new_age = NA
    #print(age)
    if (!is.na(age)){
      if (str_starts(age, "age:")){
        age = str_trim(str_replace(age, "age:", ""))
        age = str_trim(str_replace(age, "year", ""))
        if (age == "adult"){
          new_age = "20-49"
        }
        age = as.numeric(age)
        if (is.na(age)){
          new_age = NA
        }else if ( age < 20){
          new_age = "0-19"
        }else if (age >= 20 & age < 50){
          new_age = "20-49"
        }else if (age >= 49 & age < 100){
          new_age = "50-100"
        }else{
          new_age = NA
        }
      }
    }
    new_age
})

geo_rna_seq$age_old = geo_rna_seq$age
geo_rna_seq$age = ages
write.csv(geo_rna_seq, "GEO_gsm_all.csv")

############
gpl570 = read.csv("GPL570_meta_all.csv")
ages = sapply(gpl570$age, function(age){
  new_age = NA
  #print(age)
  if (!is.na(age)){
    if (str_starts(age, "age:")){
      age = str_trim(str_replace(age, "age:", ""))
      age = str_trim(str_replace(age, "year", ""))
      if (age == "adult"){
        new_age = "20-49"
      }
      age = as.numeric(age)
      if (is.na(age)){
        new_age = NA
      }else if ( age < 20){
        new_age = "0-19"
      }else if (age >= 20 & age < 50){
        new_age = "20-49"
      }else if (age >= 49 & age < 100){
        new_age = "50-100"
      }else{
        new_age = NA
      }
    }
  }
  new_age
})

gpl570$age_old = gpl570$age
gpl570$age = ages
write.csv(gpl570, "GPL570_meta_all.csv")
