#drug induced gene expression
source("../code/createGGPlot.R") #used for plotting

drugs = read.csv("drug_consensus_rank_BC.csv", stringsAsFactors = F)
agonist = drugs$name[drugs$moa == "estrogen receptor agonist"]
antagonist = drugs$name[drugs$moa == "estrogen receptor antagonist"]

load("cmap_signatures.RData")
#big rank numerbs means down regulation

meta = read.csv("cmap_drug_experiments.csv")
meta = meta[meta$cell_line == "MCF7", ]
#59272 ACE2
#7113 TMPRSS2
#7153 TOP2A
#2099 ESR1
#DPP-4  ENSG00000197635 1803
#ACE1 1636 ENSG00000159640
meta[tolower(meta$name) %in% agonist,]
meta[tolower(meta$name) %in% antagonist,]

cmap_signatures = t(cmap_signatures[cmap_signatures$V1 %in% c("59272", "2099"),])

agonist_ids = paste0("V", meta$id[tolower(meta$name) %in% agonist] + 1)
antagonist_ids = paste0("V", meta$id[tolower(meta$name) %in% antagonist] + 1)

wilcox.test(cmap_signatures[agonist_ids, 2], cmap_signatures[antagonist_ids, 2]) #ACE2
wilcox.test(cmap_signatures[agonist_ids, 1], cmap_signatures[antagonist_ids, 1]) # ESR1

a = rbind(data.frame(score = cmap_signatures[agonist_ids, 2], group="ESR agonist"),
          data.frame(score = cmap_signatures[antagonist_ids, 2], group="ESR antagonist"))
createGGPlot(a, measureVar="score", groupVars = "group", main = "", ylab = "ACE expression", xlab = "", "two.sided")

meta$id1 = paste0("V", meta$id + 1)
View(merge(meta, data.frame(cmap_signatures[antagonist_ids, 2]), by.x="id1", by.y= 0))
View(merge(meta, data.frame(cmap_signatures[agonist_ids, 2]), by.x="id1", by.y= 0))

