data = read.csv("~/Downloads/coronavirusdataset/PatientInfo.csv")

unlist(by(data$sex, data$age, count))

data = data[data$age != "", ]
cases = data.frame(t(table(data$sex,  data$age)))
cases = cases[cases$Var2 %in% c("female", "male") & cases$Var1 !="", ]

library(ggplot2)

qplot(cases$Var1, cases$Freq, color = cases$Var2, xlab = "age", ylab = "# cases", main = "korea")
data1$group = "hormone"
data1$group[data1$age %in% c("00s", "0s", "10s",  "70s","80s", "90s")] = "no hormone"
table(data1$group, data1$sex)
chisq.test(table(data1$group, data1$sex))


####
#global

data1 = read.csv("~/Downloads/Copy of COVID19_2020_open_line_list - outside_Hubei.csv", stringsAsFactors = F)
data1 = data1[!is.na(data1$age) & !is.na(data1$sex) & data1$age != "" & data1$sex %in% c("male", "female"), ]
data1$age = sapply(data1$age, function(x){
  a = unlist(strsplit(x, ""))
  if (length(a) > 1){
    if (a[1] == 0){
      paste0( "0s")
    }else{
      paste0(a[1], "0s")
    }
  }else{
    "0s"
  }
})

cases = data.frame(t(table(data1$sex,  data1$age)))
cases = cases[cases$Var2 %in% c("female", "male") & cases$Var1 !="", ]

qplot(cases$Var1, cases$Freq, color = cases$Var2, xlab = "age", ylab = "# cases", main = "all")

pdf("global_male_female_cases.pdf")
cases$gender = cases$Var2
ggplot(cases, aes(x= Var1, y = Freq, colour = gender )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
  geom_point(size=2) + 
  xlab("age") + guides(shape=FALSE, size=FALSE) +
  ylab("# cases") + coord_cartesian( ylim=c(0, 300))
dev.off()

data1$group = "hormone"
data1$group[data1$age %in% c("00", "0s", "10s",  "80s", "90s")] = "no hormone"
table(data1$group, data1$sex)
chisq.test(table(data1$group, data1$sex))

fisher.test(data1$age, data1$sex)

ratio = cases$Freq[cases$Var2 == "male"]/cases$Freq[cases$Var2 == "female"]
hist(data1_male$age_new)