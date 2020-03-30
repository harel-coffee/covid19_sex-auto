setwd("~/Documents/stanford/wars/cmap/data//")

meta_files = list.files(pattern =  "GSE")

dz_all = NULL
meta_files_valid = NULL
for (meta_file  in meta_files){
  if (file.exists(paste0(getwd(), "/", meta_file, "/dz_signature_all.txt"))) {
    dz = read.csv(paste0(getwd(), "/", meta_file, "/dz_signature_all.txt"), sep = "\t")
    if (nrow(dz) > 0 & ncol(dz) == 6){
    dz_all = rbind(dz_all, data.frame(meta_file, dz))
    }
  }
}

dz_all_subset = dz_all[grep("SAR", dz_all$meta_file), ]
#ACE2 expressed higher after infection
View(dz_all_subset[dz_all_subset$Symbol == "ACE2", ])
View(dz_all_subset[dz_all_subset$Symbol == "DPP4", ]) #but not for MERS
View(dz_all_subset[dz_all_subset$Symbol == "ESR1", ]) #but not for MERS
View(dz_all_subset[dz_all_subset$Symbol == "ESR2", ]) #but not for MERS
View(dz_all_subset[dz_all_subset$Symbol == "PGR", ]) #but not for MERS
View(dz_all_subset[dz_all_subset$Symbol == "AR", ]) #but not for MERS

sum(dz_all_subset$Symbol == "ACE2" & dz_all_subset$up_down == "up")/length(meta_files)
sum(dz_all_subset$Symbol == "ACE2" & dz_all_subset$up_down == "down")/length(meta_files)
sum(dz_all_subset$Symbol == "DPP4" & dz_all_subset$up_down == "up")/length(meta_files)
sum(dz_all_subset$Symbol == "DPP4" & dz_all_subset$up_down == "down")/length(meta_files)
sum(dz_all_subset$Symbol == "ESR1" & dz_all_subset$up_down == "up")/length(meta_files)
sum(dz_all_subset$Symbol == "ESR2" & dz_all_subset$up_down == "down")/length(meta_files)
sum(dz_all_subset$Symbol == "PGR" & dz_all_subset$up_down == "up")/length(meta_files)
sum(dz_all_subset$Symbol == "PGR" & dz_all_subset$up_down == "down")/length(meta_files)
sum(dz_all_subset$Symbol == "AR" & dz_all_subset$up_down == "up")/length(meta_files)
sum(dz_all_subset$Symbol == "AR" & dz_all_subset$up_down == "down")/length(meta_files)

