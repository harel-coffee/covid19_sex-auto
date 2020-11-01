library(xCell)
library(GEOquery)
library(dplyr)
library(stringi)
library(stringr)
setwd("~/Desktop/coronavirus_dz_signature/covid_gender/Shreya_new_data/Cell_composition_work/")
mapping = read.csv("~/Desktop/ECMO_Vs_MODS_data/transcriptome_analysis/gencode.v23.annotation.gene.probeMap.csv")[, c("gene", "ensembl")]

##################
GSE_id="GSE157103"
gse <- getGEO(GSE_id)
Metadata=pData(gse[[1]])
names(Metadata)[11]= "age"
for (i in 1:nrow(Metadata)){
  Metadata$age_new[i]= as.numeric(str_replace(as.character(Metadata$age[i]), "age \\(years\\): ","")) #correct >89
}

SRP279280_COVID_Metadata=SRP279280_COVID_Metadata%>% left_join(Metadata %>% select(geo_accession, age_new), by=c("Sample_Name"="geo_accession"))

exp = merge(mapping, SRP279280_COVID_log2.tpm.matrix, by.x= "ensembl", by.y = 0)
exp = exp[, -1]
exp = aggregate(. ~ gene, exp, median)
rownames(exp) = exp$gene
exp = exp[,-1]
cells  = xCellAnalysis(exp)

cell_meta = merge(SRP279280_COVID_Metadata, t(cells), by.x = "Sample_Name", by.y = 0)
#### cell enrichment for male
cell_meta = cell_meta[cell_meta$Sex %in% c("male", "female"),]
cell_meta$icu = as.factor(cell_meta$icu)
#cell_meta$gender=as.factor(ifelse(cell_meta$Sex=="female", "yes", "no"))
cell_meta=cell_meta[-6, -4]
OR_pvalue_g <- data.frame("names"=colnames(cell_meta[,c(18:85)]), OR=rep(NA, ncol(cell_meta[,c(18:85)])), lower=rep(NA, ncol(cell_meta[,c(18:85)])),
                          upper=rep(NA, ncol(cell_meta[,c(18:85)])), pvalue=rep(NA, ncol(cell_meta[,c(18:85)])), Adjp=rep(NA, ncol(cell_meta[,c(18:85)])))
mylogit= list()
factors=colnames(cell_meta[,c(18:85)])
for (i in 2:(nrow(cells))){
  data= cell_meta[cell_meta$disease_state == "COVID-19" & cell_meta$Sex=="male", c("icu", factors[i], "age_new")]
  data[,2]=scale(as.numeric(data[,2]))
  tests <- coef(summary(glm(icu ~  . , data = data, family = "binomial")))
  mylogit[[i]] <- glm(icu ~  . , data = data, family = "binomial")
  test <- exp(cbind(OR = coef(mylogit[[i]]), confint(mylogit[[i]])))
  OR_pvalue_g$OR[i]=test[2,1]
  OR_pvalue_g$lower[i]=test[2,2]
  OR_pvalue_g$upper[i]=test[2,3]
  OR_pvalue_g$pvalue[i] <- tests[,4][2]

}
OR_pvalue_g$Adjp=p.adjust(OR_pvalue_g$pvalue, method = "fdr")

df=OR_pvalue_g[!is.na(OR_pvalue_g$Adjp) & OR_pvalue_g$Adjp<=0.05   ,]
df = df[!df$names %in% c("ImmuneScore"),]

df = df[order(df$OR), ]
OR_value_male_icu_nonicu=df
write.csv(df, "male_OR_with_age_adjusted.csv")
df$yAxis = nrow(df):1
yAxis = nrow(df):1
names = df$names
p <- ggplot(df, aes(x = OR, y = yAxis))
p + geom_vline(aes(xintercept = 1), size = 0.5, linetype = "dashed") +
  geom_errorbarh(aes(xmax = upper, xmin = lower), size = 0.8, height = .3, color = "gray50") +
  geom_point(size = 5.5, color = "orange") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.border = element_rect(fill=NA, colour = "black", size=0.8),
        axis.text.x = element_text(size=16, face="bold"),
        axis.text.y = element_text(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        plot.title = element_text(size=16, face="bold", hjust = 0.7, vjust = -15)) +
  scale_y_continuous(breaks = yAxis, labels = names) +
  scale_x_continuous(breaks = c(0.1, 0.5, 1, 1.5,  2.5,  4,  6, 9)) +
  coord_trans(x = "log10") +
  ylab("") +
  xlab("Odds ratio (log scale)")+ggtitle("Male")

ggsave("Male_OR_adjusted_with_age.pdf", height = 5, width = 8)


######## Cell enrichment for Female patients
#cell_meta=cell_meta[-6, -4]
OR_pvalue_g <- data.frame("names"=colnames(cell_meta[,c(18:85)]), OR=rep(NA, ncol(cell_meta[,c(18:85)])), lower=rep(NA, ncol(cell_meta[,c(18:85)])),
                          upper=rep(NA, ncol(cell_meta[,c(18:85)])), pvalue=rep(NA, ncol(cell_meta[,c(18:85)])), Adjp=rep(NA, ncol(cell_meta[,c(18:85)])))
mylogit= list()
factors=colnames(cell_meta[,c(18:85)])
for (i in 2:(nrow(cells))){
  data= cell_meta[cell_meta$disease_state == "COVID-19" & cell_meta$Sex=="female", c("icu", factors[i], "age_new")]
  data[,2]=scale(as.numeric(data[,2]))
  tests <- coef(summary(glm(icu ~  . , data = data, family = "binomial")))
  mylogit[[i]] <- glm(icu ~  . , data = data, family = "binomial")
  test <- exp(cbind(OR = coef(mylogit[[i]]), confint(mylogit[[i]])))
  OR_pvalue_g$OR[i]=test[2,1]
  OR_pvalue_g$lower[i]=test[2,2]
  OR_pvalue_g$upper[i]=test[2,3]
  OR_pvalue_g$pvalue[i] <- tests[,4][2]

}
OR_pvalue_g$Adjp=p.adjust(OR_pvalue_g$pvalue, method = "fdr")

df=OR_pvalue_g[!is.na(OR_pvalue_g$Adjp) & OR_pvalue_g$Adjp<=0.05   ,]
df = df[!df$names %in% c("ImmuneScore", "StromaScore", "Endothelial cells"),]##, "Endothelial cells"

df = df[order(df$OR), ]
OR_value_female_icu_nonicu=df
write.csv(df, "Female_OR_with_age_adjusted.csv")
df$yAxis = nrow(df):1
yAxis = nrow(df):1
names = df$names
p <- ggplot(df, aes(x = OR, y = yAxis))
p + geom_vline(aes(xintercept = 1), size = 0.5, linetype = "dashed") +
  geom_errorbarh(aes(xmax = upper, xmin = lower), size = 0.8, height = .3, color = "gray50") +
  geom_point(size = 5.5, color = "orange") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.border = element_rect(fill=NA, colour = "black", size=0.8),
        axis.text.x = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        plot.title = element_text(size=16, face="bold", hjust = 0.9, vjust = -15)) +
  scale_y_continuous(breaks = yAxis, labels = names) +
  scale_x_continuous(breaks = c(0.01, 0.05, 1, 10, 60)) +
  coord_trans(x = "log10") +
  ylab("") +
  xlab("Odds ratio (log scale)")+ggtitle("Female")
ggsave("Female_OR_adjusted_with_age.pdf", height = 4, width = 7)


##Heatmaps for cells in male
library(ComplexHeatmap)
# For male
exp_data=cbind(cell_meta[,c(1:10)],cell_meta[,colnames(cell_meta)%in% OR_value_male_icu_nonicu$names])
exp_data=cbind(cell_meta[,c(1:10)], apply(cell_meta[,colnames(cell_meta)%in% OR_value_male_icu_nonicu$names], 2, function(x){scale(x)}))
exp_data1= exp_data[exp_data$Sex=="male" &exp_data$disease_state == "COVID-19",]
exp_data1= exp_data1[order(exp_data1$icu),]
pdf("Heatmaps_male_OR_1_cells.pdf")
Heatmap(as.matrix(exp_data1[,c(11:18)]), show_row_names = F, cluster_rows = F)+
  Heatmap(exp_data1$icu, name="ICU",c("cyan4", "coral1"))
dev.off()


## Heatmaps for cells in female
## For female

exp_data=cbind(cell_meta[,c(1:10)],cell_meta[,colnames(cell_meta)%in% OR_value_female_icu_nonicu$names])
exp_data=cbind(cell_meta[,c(1:10)],apply(cell_meta[,colnames(cell_meta)%in% OR_value_female_icu_nonicu$names], 2, function(x){scale(x)}))
exp_data1= exp_data[exp_data$Sex=="female" & exp_data$disease_state=="COVID-19",]
exp_data1= exp_data1[order(exp_data1$icu),]
pdf("Heatmaps_female_OR_1_cells.pdf")
Heatmap(as.matrix(exp_data1[,c(11:15)]), show_row_names = F, cluster_rows = F)+
  Heatmap(exp_data1$icu, name="ICU",c("cyan4", "coral1"))
dev.off()

