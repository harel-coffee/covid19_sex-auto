#
setwd("/Users/ChenB1/Documents/stanford/wars/gender/data/")

#
'library(GEOmetadb)
con <- dbConnect(SQLite(),"GEOmetadb.sqlite")
dbListFields(con,"gsm")
rs <- dbGetQuery(con,‘select gsm,series_id, title, source_name_ch1, characteristics_ch1 from gsm where gpl == “GPL570” ’)
'

library(stringr)

gsmMeta = read.csv("GPL570_meta.csv", stringsAsFactors = F)

#output = paste("gsm", "gender", "tissue", "age", "source_name",
#               "gse" , "title" , "characteristics" , sep = "\t")
#write(output, "GPL570_gsm_all.txt", append = T)
gsmMeta$gender = NA
for (i in 1:nrow(gsmMeta)){
  print(i)
  gsm = gsmMeta$gsm[i]
  
  gender = NA
  if (length(grep("male", gsmMeta$characteristics_ch1[i], ignore.case = T)) > 0){
     if (length(grep("female", gsmMeta$characteristics_ch1[i], ignore.case = T)) > 0){
       gender = "female"
     }else{
       gender = "male"
     }
  }
  
  gsmMeta$gender[i] = gender
  tissue = NA

  age = NA
 # if (length(grep("age:", gsmMeta[[i]]$Sample_characteristics_ch1, ignore.case = T)) > 0){
 #   age = gsmMeta[[i]]$Sample_characteristics_ch1[grep("age:", gsmMeta[[i]]$Sample_characteristics_ch1, ignore.case = T)[1]]
 # }
  
 # output = paste(gsm, gender, tissue, age, source_name = gsmMeta[[i]]$Sample_source_name_ch1, 
 #                gse =  gsmMeta[[i]]$Sample_series_id, title = gsmMeta[[i]]$Sample_title, characteristics = paste(gsmMeta[[i]]$Sample_characteristics_ch1, collapse = "; "), sep = "\t")
 # write(output, "GPL570_gsm_all.txt", append = T)
}

#write.csv(gsm_gender, "GEO_gsm_all.csv")

#####
#rs <- dbGetQuery(con,'select gsm, title, characteristics_ch1 from gsm where gpl == "GPL570" and LOWER(characteristics_ch1) LIKE "%male%" ')

tail(sort(table(gsmMeta$source_name_ch1)))

