library(tidyr)

args = commandArgs(trailingOnly=TRUE)

setwd("~/Documents/PhD/haplotype-project")
source("Scripts/Haplotype-projekt-ver2.R")

readData()

# runMutation(mut = "c.7617+1G>A", gene = "BRCA2", country = NULL, F)
runMutation(mut = "c.5266dupC", gene = "BRCA1", country = NULL, F)
# runMutation(mut = "c.1016dupA", gene = "BRCA1", country = NULL, F)
# runMutation(mut = "c.1016delA", gene = "BRCA1", country = NULL, F)

#prepare_dataframe(mut = "c.1016delA", gene = "BRCA1")
prepare_dataframe(mut = "c.5266dupC", gene = "BRCA1")

famHaplotypes <- findFamHaplotypes(mut = "c.1016delA", gene = "BRCA1")
famHaplotypes <- findFamHaplotypes(mut = "c.1016dupA", gene = "BRCA1")
famHaplotypes <- findFamHaplotypes(mut = "c.5266dupC", gene = "BRCA1")
famHaplotypes <- findFamHaplotypes(mut = "c.7617+1G>A", gene = "BRCA2")

dim(famHaplotypes)

# if (args[1] == 1) {
#     gene = "BRCA1"
# 
#     muts = c("c.-200-?_80+?del", "c.1687C>T", "c.181T>G", "c.211A>G", "c.2475delC",
#              "c.2681_2682delAA", "c.3319G>T", "c.3331_3334delCAAG", "c.3481_3491del11", 
#              "c.3700_3704del5", "c.3756_3759delGTCT") 
#     
#     runMutations(muts, gene, T)
# } else if (args[1] == 2){
#     gene = "BRCA1"
#     
#     muts = c("c.4035delA", "c.4065_4068delTCAA", 
#      "c.4186-?_4357+?dup", "c.427G>T", "c.4327C>T", "c.5123C>A", "c.5266dupC", 
#      "c.5333-36_5406+400del510", "c.5503C>T", "c.68_69delAG")
#     
#     runMutations(muts, gene, T)
# } else if (args[1] == 3){
#     gene = "BRCA2"
#     muts = c("c.1310_1313delAAGA", "c.1813dupA", "c.2808_2811delACAA", "c.3847_3848delGT",
#              "c.4478_4481delAAAG", "c.5645C>A", "c.5682C>G", "c.5722_5723delCT")
#     
#     runMutations(muts, gene, T)
# } else if (args[1] == 4){
#     gene = "BRCA2"
#     muts = c("c.5946delT", "c.6275_6276delTT", "c.7069_7070delCT", "c.7617+1G>A", 
#              "c.771_775del5", "c.7934delG", "c.8537_8538delAG")
#     
#     runMutations(muts, gene, T)
# }



# Reformat genotype data
#brca1_geno_merged_twoRows = read.table("input/139_ouh_june_2017/139_ouh_brca1_onco_geno_twoRows.txt", header = T, stringsAsFactors = F) 
#save(brca1_geno_merged_twoRows, file = "cache/brca1_geno_merged_twoRows.RData")