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
gsmMeta$age = NA
gsmMeta$tissue = NA

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
  if (is.na(gender)){
    if (length(grep("Sex F", gsmMeta$characteristics_ch1[i], ignore.case = T)) > 0){
      gender = "female"
    }
    if (length(grep("Sex M", gsmMeta$characteristics_ch1[i], ignore.case = T)) > 0){
      gender = "male"
    }    
  }
  
  gsmMeta$gender[i] = gender
  
  tissue = NA
  if (length(grep("tissue:", gsmMeta$characteristics_ch1[i], ignore.case = T)) > 0){
      items = unlist(strsplit(gsmMeta$characteristics_ch1[i], ";|,"))
      tissue = str_trim(tolower(items[grep("tissue:", items, ignore.case = T)]))
  }
  gsmMeta$tissue[i] = tissue
  
  age = NA
  if (length(grep("age:", gsmMeta$characteristics_ch1[i], ignore.case = T)) > 0){
    items = unlist(strsplit(gsmMeta$characteristics_ch1[i], ";|,"))
    age = str_trim(tolower(items[grep("age:", items, ignore.case = T)]))
  }
  gsmMeta$age[i] = age
  
# if (length(grep("age:", gsmMeta[[i]]$Sample_characteristics_ch1, ignore.case = T)) > 0){
 #   age = gsmMeta[[i]]$Sample_characteristics_ch1[grep("age:", gsmMeta[[i]]$Sample_characteristics_ch1, ignore.case = T)[1]]
 # }
  
 # output = paste(gsm, gender, tissue, age, source_name = gsmMeta[[i]]$Sample_source_name_ch1, 
 #                gse =  gsmMeta[[i]]$Sample_series_id, title = gsmMeta[[i]]$Sample_title, characteristics = paste(gsmMeta[[i]]$Sample_characteristics_ch1, collapse = "; "), sep = "\t")
 # write(output, "GPL570_gsm_all.txt", append = T)
}

write.csv(gsmMeta, "GPL570_meta_all.csv")

#####
#rs <- dbGetQuery(con,'select gsm, title, characteristics_ch1 from gsm where gpl == "GPL570" and LOWER(characteristics_ch1) LIKE "%male%" ')

tail(sort(table(gsmMeta$source_name_ch1)))
tail(sort(table(gsmMeta$gender)))
tail(sort(table(gsmMeta$tissue)), 30)
tail(sort(table(gsmMeta$age)))

#########
#expression
load('~/Downloads/Merged42.Rdata')


