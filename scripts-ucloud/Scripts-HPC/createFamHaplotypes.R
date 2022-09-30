source("Haplotype-projekt-ver2.R")

readData()

# BRCA2
l=findAllFamHaplotypes2(pheno_data=brca2_pheno_merged, geno_data=brca2_geno_merged, gene="BRCA2")
# BRCA1
#l=findAllFamHaplotypes2(pheno_data=brca1_pheno_merged, geno_data=brca1_geno_merged, gene="BRCA1")
