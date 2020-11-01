##########
# for GSE156063
##########
setwd("../GSE156063_age//")
study="GSE156063"
Male_CT_SARS2 = read.csv(paste(study,"_age_above60_male_CT_SARS-Cov-2_DE_gene.csv", sep = ""))
Female_CT_SARS2= read.csv(paste(study,"_age_above60_female_CT_SARS-Cov-2_DE_gene.csv",sep = ""))
Female_male_SARS2 = read.csv(paste(study,"_age_above60_female_male_SARS-Cov-2_DE_gene.csv", sep = ""))
Female_male_CT = read.csv(paste(study,"_age_above60_female_male_CT_DE_gene.csv",sep = ""))

### prepare the data list for comparisons
data_list = list("Male_CT_vs_SARS2"= Male_CT_SARS2$gene[!is.na(Male_CT_SARS2$gene)], "Female_CT_vs_SARS2"= Female_CT_SARS2$gene[!is.na(Female_CT_SARS2$gene)],
                 "Female_vs_male_CT"=Female_male_CT$gene[!is.na(Female_male_CT$gene)],"Female_vs_male_SARS2"=Female_male_SARS2$gene[!is.na(Female_male_SARS2$gene)])

## create the venn diagramm file and plot them
if(length(data_list) >3){
  venn.plot=venn.diagram(data_list, NULL, fill=c("darkmagenta", "grey","darkblue","orange"), draw.quad.venn(area1 = length(Male_CT_vs_SARS2), area2 = length(Female_CT_vs_SARS2), area3 = length(Female_vs_male_CT),
                                                                                                            area4 = length(Female_vs_male_SARS2)), alpha=0.5,cex = 1.1, cat.fontface=2, cat.dist = c(0.25,0.25,0.1, 0.1), cat.pos = c(-7,6,0,-15),cat.cex=0.9,
                         main =paste(study, "age_above_60", sep = "_"), main.fontface = 2)
}else{
  venn.plot=venn.diagram(data_list, NULL, fill=c("darkmagenta", "grey","darkblue"), draw.triple.venn(area1 = length(Male_CT_vs_SARS2), area2 = length(Female_CT_vs_SARS2), area3 = length(Female_vs_male_SARS2))
                         , alpha=0.5,cex = 1.1, cat.fontface=2, cat.dist = c(0.05,0.05, -0.47), cat.pos = c(-10,10, 0),cat.cex=0.9,
                         main = paste(study, "age_above_60", sep = "_"), main.fontface = 2)
}

dev.off()
grid.draw(venn.plot)


## save the image
pdf(paste(study,"Age_above_60_venn_diagram.pdf", sep = ""), width = 7, height = 4)
grid.draw(venn.plot)
dev.off()

## Find the genes in each group
library(gplots)
ItemsList <- venn(data_list, show.plot = FALSE)
GSE156063_above60_Gene_list= attributes(ItemsList)$intersections
#GSE152075_Gene_list= data.frame(matrix(unlist(Gene_list), nrow=length(Gene_list), byrow=T))
#write.csv(GSE152075_Gene_list)
## Below 60
study="GSE156063"
Male_CT_SARS2 = read.csv(paste(study,"_age_age_below60_male_CT_SARS-Cov-2_DE_gene.csv", sep = ""))
Female_CT_SARS2= read.csv(paste(study,"_age_age_below60_female_CT_SARS-Cov-2_DE_gene.csv",sep = ""))
Female_male_SARS2 = read.csv(paste(study,"_age_age_below60_female_male_SARS-Cov-2_DE_gene.csv", sep = ""))
Female_male_CT = read.csv(paste(study,"_age_age_below60_female_male_CT_DE_gene.csv",sep = ""))

### prepare the data list for comparisons
data_list = list("Male_CT_vs_SARS2"= Male_CT_SARS2$gene[!is.na(Male_CT_SARS2$gene)], "Female_CT_vs_SARS2"= Female_CT_SARS2$gene[!is.na(Female_CT_SARS2$gene)],
                 "Female_vs_male_CT"=Female_male_CT$gene[!is.na(Female_male_CT$gene)],"Female_vs_male_SARS2"=Female_male_SARS2$gene[!is.na(Female_male_SARS2$gene)])

## create the venn diagramm file and plot them
if(length(data_list) >3){
  venn.plot=venn.diagram(data_list, NULL, fill=c("darkmagenta", "grey","darkblue","orange"), draw.quad.venn(area1 = length(Male_CT_vs_SARS2), area2 = length(Female_CT_vs_SARS2), area3 = length(Female_vs_male_CT),
                                                                                                            area4 = length(Female_vs_male_SARS2)), alpha=0.5,cex = 1.1, cat.fontface=2, cat.dist = c(0.25,0.25,0.1, 0.1), cat.pos = c(-7,6,0,-15),cat.cex=0.9,
                         main =paste(study, "age_below_60", sep = "_"), main.fontface = 2)
}else{
  venn.plot=venn.diagram(data_list, NULL, fill=c("darkmagenta", "grey","darkblue"), draw.triple.venn(area1 = length(Male_CT_vs_SARS2), area2 = length(Female_CT_vs_SARS2), area3 = length(Female_vs_male_SARS2))
                         , alpha=0.5,cex = 1.1, cat.fontface=2, cat.dist = c(0.05,0.05, -0.47), cat.pos = c(-10,10, 0),cat.cex=0.9,
                         main = paste(study, "age_below_60", sep = "_"), main.fontface = 2)
}

dev.off()
grid.draw(venn.plot)


## save the image
pdf(paste(study,"Age_below_60_venn_diagram.pdf", sep = ""), width = 7, height = 4)
grid.draw(venn.plot)
dev.off()

## Find the genes in each group
library(gplots)
ItemsList <- venn(data_list, show.plot = FALSE)
GSE156063_below60_Gene_list= attributes(ItemsList)$intersections


DE_file= list()
#DE_file=list(Male_CT_SARS2,Female_CT_SARS2)
DE_file= list("GSE156063_Male_CT_vs_SARS"=Male_CT_SARS2, "GSE156063_female_CT_vs_SARS"=Female_CT_SARS2)

# for (GSE in GSE_ids){
#   DE_file[[paste(GSE, "_Male_CT_vs_SARS",sep = "")]]= read.csv(paste("../",GSE,"/", GSE,"_male_CT_SARS-Cov2_DE_gene.csv", sep = ""))
#   DE_file[[paste(GSE, "_female_CT_vs_SARS",sep = "")]]= read.csv(paste("../",GSE,"/", GSE,"_female_CT_SARS-Cov2_DE_gene.csv", sep = ""))
# }

# DE_file[["GSE152075_female_CT_vs_SARS"]]$gene=DE_file[["GSE152075_female_CT_vs_SARS"]]$identifier
#
# comparison_list= list("GSE152075"=GSE152075_Gene_list, "GSE156063"=GSE156063_Gene_list, "SRP267176"=SRP267176_Gene_list)

## Extract the DE genes for the specific genes
spc_DE_list = list()
for (GSE in study){
  spc_DE_list[[paste(GSE,"_Male_CT_vs_SARS", sep = "")]]= DE_file[[paste(GSE,"_Male_CT_vs_SARS", sep = "")]][DE_file[[paste(GSE,"_Male_CT_vs_SARS", sep = "")]]$gene %in% GSE156063_above60_Gene_list[["Male_CT_vs_SARS2"]],]
  spc_DE_list[[paste(GSE,"_female_CT_vs_SARS", sep = "")]]= DE_file[[paste(GSE,"_female_CT_vs_SARS", sep = "")]][DE_file[[paste(GSE,"_female_CT_vs_SARS", sep = "")]]$gene %in% GSE156063_above60_Gene_list[["Female_CT_vs_SARS2"]],]
}





#Perform enrichment analysis of the samples
library(ggplot2)
library(ggpubr)
library(ggthemes)
#if you want to plot up-regulated and downregulted both togther, Assign upregulated pathways in regulation column= "1" and for downregulation= "2

library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018")

GSE_ids= c("GSE156063", "GSE152075", "SRP267176")
GSE_ids="SRP279280"
## for female
for (GSE in study){{
  up_enriched <- enrichr(spc_DE_list[[paste(GSE,"_female_CT_vs_SARS", sep = "")]]$identifier[spc_DE_list[[paste(GSE,"_female_CT_vs_SARS", sep = "")]]$log2FoldChange >1], dbs)
  #printEnrich(up_enriched, paste(GSE,"up_female_CT_vs_SARS_GO.txt",sep="_") , sep = "\t", columns = c(1:9))
  up_bp <- up_enriched[["GO_Biological_Process_2018"]]
  up_bp_sign = up_bp[up_bp$Adjusted.P.value <= 0.2,]
  up_bp_sign=up_bp_sign %>% separate(Overlap, c("A", "B"))
  up_bp_sign$GeneRatio= as.numeric(up_bp_sign$A)/as.numeric(up_bp_sign$B)
  if(nrow(up_bp_sign) >1){
    up_bp_sign$log10_adjp <- -log10(up_bp_sign$Adjusted.P.value)
    ggplot(up_bp_sign, aes(x = log10_adjp, y= reorder(Term,log10_adjp))) + geom_point(aes(color= as.factor(log10_adjp), size=GeneRatio))+ylab(label = NULL)+ xlab(label = "-log10(Adjusted P-value)")+
      theme_few()+theme(axis.text.x = element_text(size=10, face="bold"),
                        axis.text.y = element_text(size=10, face="bold"),
                        axis.title =  element_text(size=10, face="bold"),
                        legend.position = "none",
                        panel.border = element_rect(fill=NA, colour = "black", size=0.8),
                        axis.ticks = element_line(size = 0.5))
    ggsave(paste(GSE,"Above60_up_female_CT_vs_SARS_GO.pdf", sep="_"), dpi = 350, width = 7, height = 6)
    write.csv(up_bp_sign, paste(GSE,"Above60_up_female_CT_vs_SARS_GO.csv", sep="_"))
  }else{
    print("No significant pathways were found")
  }
}

  for (GSE in study){
    down_enriched <- enrichr(spc_DE_list[[paste(GSE,"_female_CT_vs_SARS", sep = "")]]$identifier[spc_DE_list[[paste(GSE,"_female_CT_vs_SARS", sep = "")]]$log2FoldChange <1], dbs)
    # printEnrich(down_enriched, paste(GSE,"down_female_CT_vs_SARS_GO.txt",sep="_") , sep = "\t", columns = c(1:9))
    down_bp <- down_enriched[["GO_Biological_Process_2018"]]
    down_bp_sign = down_bp[down_bp$Adjusted.P.value <= 0.2,]
    down_bp_sign=down_bp_sign %>% separate(Overlap, c("A", "B"))
    down_bp_sign$GeneRatio= as.numeric(down_bp_sign$A)/as.numeric(down_bp_sign$B)
    if(nrow(down_bp_sign) >1){
      down_bp_sign$log10_adjp <- -log10(down_bp_sign$Adjusted.P.value)
      ggplot(down_bp_sign[c(1:20),], aes(x = log10_adjp, y= reorder(Term,log10_adjp))) + geom_point(aes(color= as.factor(log10_adjp), size=GeneRatio))+ylab(label = NULL)+ xlab(label = "-log10(Adjusted P-value)")+
        theme_few()+theme(axis.text.x = element_text(size=10, face="bold"),
                          axis.text.y = element_text(size=10, face="bold"),
                          axis.title =  element_text(size=10, face="bold"),
                          legend.position = "none",
                          panel.border = element_rect(fill=NA, colour = "black", size=0.8),
                          axis.ticks = element_line(size = 0.5))
      ggsave(paste(GSE,"Above60_down_female_CT_vs_SARS_GO.pdf", sep="_"), dpi = 350, width = 8, height = 6)
      write.csv(down_bp_sign, paste(GSE,"Above60_down_female_CT_vs_SARS_GO.csv", sep="_"))
    }else{
      print("No significant pathways were found")
    }
  }
}


## for male

for (GSE in study){{
  up_enriched <- enrichr(spc_DE_list[[paste(GSE,"_Male_CT_vs_SARS", sep = "")]]$identifier[spc_DE_list[[paste(GSE,"_Male_CT_vs_SARS", sep = "")]]$log2FoldChange >1], dbs)
  #printEnrich(up_enriched, paste(GSE,"up_Male_CT_vs_SARS_GO.txt",sep="_") , sep = "\t", columns = c(1:9))
  up_bp <- up_enriched[["GO_Biological_Process_2018"]]
  up_bp_sign = up_bp[up_bp$Adjusted.P.value <= 0.2,]
  up_bp_sign=up_bp_sign %>% separate(Overlap, c("A", "B"))
  up_bp_sign$GeneRatio= as.numeric(up_bp_sign$A)/as.numeric(up_bp_sign$B)
  if(nrow(up_bp_sign) >1){
    up_bp_sign$log10_adjp <- -log10(up_bp_sign$Adjusted.P.value)
    ggplot(up_bp_sign, aes(x = log10_adjp, y= reorder(Term,log10_adjp))) + geom_point(aes(color= as.factor(log10_adjp), size=GeneRatio))+ylab(label = NULL)+ xlab(label = "-log10(Adjusted P-value)")+
      theme_few()+theme(axis.text.x = element_text(size=10, face="bold"),
                        axis.text.y = element_text(size=10, face="bold"),
                        axis.title =  element_text(size=10, face="bold"),
                        legend.position = "none",
                        panel.border = element_rect(fill=NA, colour = "black", size=0.8),
                        axis.ticks = element_line(size = 0.5))
    ggsave(paste(GSE,"Above60_up_Male_CT_vs_SARS_GO.pdf", sep="_"), dpi = 350, width = 7, height = 6)
    write.csv(up_bp_sign, paste(GSE,"Above60_up_Male_CT_vs_SARS_GO.csv", sep="_"))
  }else{
    print("No significant pathways were found")
  }
}

  for (GSE in study){
    down_enriched <- enrichr(spc_DE_list[[paste(GSE,"_Male_CT_vs_SARS", sep = "")]]$identifier[spc_DE_list[[paste(GSE,"_Male_CT_vs_SARS", sep = "")]]$log2FoldChange <1], dbs)
    #printEnrich(down_enriched, paste(GSE,"down_Male_CT_vs_SARS_GO.txt",sep="_") , sep = "\t", columns = c(1:9))
    down_bp <- down_enriched[["GO_Biological_Process_2018"]]
    down_bp_sign = down_bp[down_bp$Adjusted.P.value <= 0.2,]
    down_bp_sign=down_bp_sign %>% separate(Overlap, c("A", "B"))
    down_bp_sign$GeneRatio= as.numeric(down_bp_sign$A)/as.numeric(down_bp_sign$B)
    if(nrow(down_bp_sign) >1){
      down_bp_sign$log10_adjp <- -log10(down_bp_sign$Adjusted.P.value)
      ggplot(down_bp_sign, aes(x = log10_adjp, y= reorder(Term,log10_adjp))) + geom_point(aes(color= as.factor(log10_adjp), size=GeneRatio))+ylab(label = NULL)+ xlab(label = "-log10(Adjusted P-value)")+
        theme_few()+theme(axis.text.x = element_text(size=10, face="bold"),
                          axis.text.y = element_text(size=10, face="bold"),
                          axis.title =  element_text(size=10, face="bold"),
                          legend.position = "none",
                          panel.border = element_rect(fill=NA, colour = "black", size=0.8),
                          axis.ticks = element_line(size = 0.5))
      ggsave(paste(GSE,"Above60_down_Male_CT_vs_SARS_GO.pdf", sep="_"), dpi = 350, width = 7, height = 6)
      write.csv(down_bp_sign, paste(GSE,"Above60_down_Male_CT_vs_SARS_GO.csv", sep="_"))
    }else{
      print("No significant pathways were found")
    }
  }
}
