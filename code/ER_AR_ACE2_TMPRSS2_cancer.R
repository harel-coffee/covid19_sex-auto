#examine ACE2/TMPRSS2 expression in ER positive breast cancer AR positive prostate
library(cgdsr)
library(survival)
library(survminer)

mycgds = CGDS("https://www.cbioportal.org/")

test(mycgds)

# Get list of cancer studies at server
a = getCancerStudies(mycgds)
#retrieve gene symbol
symbols = c("ACE2", "TMPRSS2", "DPP4")


###########l
#breast cancer
# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = "brca_tcga_pub" #prad_tcga_pub 

# Get data slices for a specified list of genes, genetic profile and case list
symbol_set = symbols #split(symbols, ceiling(seq_along(symbols)/90))

mrna_matrix = data.frame()
  mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]
  print(mycancerstudy)
  # Get available genetic profiles
  mygeneticprofiles = getGeneticProfiles(mycgds,mycancerstudy)[,1]
  
  #other two types missed ACE2 expression
  mygeneticprofile =  mygeneticprofiles[grep("mrna", mygeneticprofiles)][1]

  print(symbol_set)
  mrna = getProfileData(mycgds,symbol_set ,mygeneticprofile,mycaselist)

  myclinicaldata = getClinicalData(mycgds,mycaselist)

gene1 = "ACE2"
a = rbind(data.frame(score = mrna[rownames(myclinicaldata[myclinicaldata$ER_STATUS == "Positive",]), gene1], group="ER positive"),
            data.frame(score = mrna[rownames(myclinicaldata[myclinicaldata$ER_STATUS == "Negative",]), gene1], group="ER negative"))
a = a[!is.na(a$score), ]
Ttest = t.test(score ~ group, a)
paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
       "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
       "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
       "n: ", nrow(a)
)
pdf(paste0("figure/", gene1, "ER_status.pdf"))
createGGPlot(a, measureVar="score", groupVars = "group", main = "Primary Breast Cancer", ylab = paste0(gene1, " expression"), xlab = "", "two.sided")
dev.off()

#######
##prostate cancer
# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = "prad_tcga_pub" #prad_tcga_pub 

# Get data slices for a specified list of genes, genetic profile and case list
symbol_set = symbols #split(symbols, ceiling(seq_along(symbols)/90))

mrna_matrix = data.frame()
mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]
print(mycancerstudy)
# Get available genetic profiles
mygeneticprofiles = getGeneticProfiles(mycgds,mycancerstudy)[,1]

#match term can be changed.
mygeneticprofile = "prad_tcga_pub_rna_seq_v2_mrna_median_Zscores" # mygeneticprofiles[grep("mrna", mygeneticprofiles)][2]

print(symbol_set)
mrna = getProfileData(mycgds,symbol_set ,mygeneticprofile,mycaselist)

myclinicaldata = getClinicalData(mycgds,mycaselist)

gene1 = "TMPRSS2"
a = rbind(data.frame(score = mrna[rownames(myclinicaldata[myclinicaldata$AR_SCORE > 0,]), gene1], group="AR positive"),
          data.frame(score = mrna[rownames(myclinicaldata[myclinicaldata$AR_SCORE < 0,]), gene1], group="AR negative"))
a = a[!is.na(a$score), ]
Ttest = t.test(score ~ group, a)
paste0("ratio: ", format(Ttest$estimate[1] / Ttest$estimate[2], digit=2), ":1", "\n",
       "diff: [", format(Ttest$conf.int[1], digit=2)," ", format(Ttest$conf.int[2], digit=2), "]", "\n", 
       "p: ", format(Ttest$p.value, digits = 2, scientific=T), "\n",
       "n: ", nrow(a)
)
pdf(paste0("figure/", gene1, "AR_protein.pdf"))
createGGPlot(a, measureVar="score", groupVars = "group", main = "Primary Prostate Cancer", ylab = paste0(gene1, " expression"), xlab = "", "two.sided")
dev.off()

