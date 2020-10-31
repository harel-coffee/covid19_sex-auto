
#download data here https://github.com/beoutbreakprepared/nCoV2019/tree/master/latest_data
#as of 04/02/2020
data1 = read.csv("latestdata.csv", stringsAsFactors = F)
data1 = data1[!is.na(data1$age) & !is.na(data1$sex) & data1$age != "" & data1$sex %in% c("male", "female"), ]
data1$age = sapply(data1$age, function(x){
  a = unlist(strsplit(x, ""))
  if (length(a) > 1){
    paste0(a[1], "0s")
  }else{
    "0s"
  }
})

data1$age_group = NA
data1$age_group[data1$age %in% c("00s", "0s", "10s")] = "0-19"
data1$age_group[data1$age %in% c("20s", "30s", "40s", "50s")] = "20-59"
data1$age_group[data1$age %in% c("60s", "70s", "80s", "90s")] = "60-100"

cases = data.frame(t(table(data1$sex,  data1$age_group)))
cases = cases[cases$Var2 %in% c("female", "male") & cases$Var1 !="", ]
cases$number_per_M = (cases$Freq/sum(cases$Freq)) * 1000000

qplot(cases$Var1, cases$Freq, color = cases$Var2, xlab = "age", ylab = "# cases", main = "all")

pdf("figure/gender_case.pdf")
ggplot(cases, aes(x= Var1 , y = number_per_M, fill = Var2 )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
  geom_bar(stat="identity") + scale_fill_manual(values=c('blue','red')) + 
  ylab("# cases in one million patients") + xlab ("") 
dev.off()

