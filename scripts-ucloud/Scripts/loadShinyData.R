initShinyData <- function(){
    print("loading BRCA2 toy dataset...")
    
    loadedGenes <<- c("BRCA2_toySample")
    
    brca2_toy_pheno <<- readRDS("input/139_ouh_june_2017/BRCA2_toySample/BRCA2_toySample_pheno.Rdata")
    brca2_toy_geno <<- readRDS("input/139_ouh_june_2017/BRCA2_toySample/BRCA2_toySample_geno.Rdata")
    brca2_toy_fam_geno <<- readRDS("input/139_ouh_june_2017/BRCA2_toySample/BRCA2_toySample_fam_geno.Rdata")
    
    # Read SNP genomic coordinates
    chr_coords_all <<- read.table("input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed-with_centiMorgan_coords_from_linear_interpolation.txt", header = T, sep = "\t")
    
    # Read Genetic coordinates (centimorgan)
    morgan.brca2 <<- read.table2("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr13_combined_b37.txt", header = T) # phase3 map
    
    # BRCA mutation coordinates
    brca2_mutation_coords <<- read.table("input/139_ouh_june_2017/BRCA2_mut_position.txt", header = T, sep = "\t")
    
    # Initiate project lookup table
    open_projects <<- list()
    open_projects[["BRCA2_toySample"]] <<- list(project_path = paste0("cache/BRCA2_toySample/"), project = F)#paste0("cache/BRCA2_toySample/"), project = F)
    
    # create projects' cache folder
    dir.create("cache_projects", showWarnings = F)
    
    # Load genotype frequencies
    popFreqs_brca1 <<- read.table2("cache/pop_freqs_genotypes_brca1_ordered.txt", header = T)
    popFreqs_brca2 <<- read.table2("cache/pop_freqs_genotypes_brca2_ordered.txt", header = T)
    
    print("Done loading BRCA2 toy dataset")
}

loadBRCA1 <- function(){
    print("loading BRCA1 dataset...")
    brca1_pheno <- read.csv("input/139_ouh_june_2017/B1_Onco_phenotype_distribution_311215.csv")
    brca1_pheno_hisp <- read.csv("input/139_ouh_june_2017/B1_Hispanic_OncoArray_phenotypes_230817.csv")
    brca1_pheno_hisp <<- brca1_pheno_hisp %>% mutate(Country="Hispanic")
    
    # I think these are not important and therefore removed in the merge step.
    brca1_pheno_merged <<- rbind(brca1_pheno[, 1:47], brca1_pheno_hisp) %>% arrange(FamCode)
    
    if (!file.exists("cache/brca1_geno_plink.RData")){
        brca1_geno_eu_plink <- read.table2("input/139_ouh_june_2017/139_ouh_brca1_onco_geno_plink_format.txt", header = T)
        brca1_geno_hisp_plink <- read.table2("input/139_ouh_june_2017/139_ouh_female_hisp_brca1_onco_geno_plink_format.txt", header = T)
        
        brca1_geno_plink <- rbind(brca1_geno_eu_plink, brca1_geno_hisp_plink)
        
        save(brca1_geno_plink, file = "cache/brca1_geno_plink.RData")
    }
    load("cache/brca1_geno_plink.RData", envir = .GlobalEnv)
    
    # Load family geno data for brca1
    if (file.exists("cache/famHaplotypes-BRCA1-geno.RData") || file.exists("cache/famHaplotypes-BRCA1-geno_plink_format.RData")){
        load("cache/famHaplotypes-BRCA1-geno_plink_format.RData", envir = .GlobalEnv)
    }else{
        famHaplotypes_brca1_geno <<- read.table(file = "input/139_ouh_june_2017/famHaplotypes-BRCA1-geno.txt", sep = "\t", header = T)
        save(famHaplotypes_brca1_geno, file = "cache/famHaplotypes-BRCA1-geno.RData")
        famHaplotypes_brca1_geno_plink <<- read.table(file = "input/139_ouh_june_2017/famHaplotypes-BRCA1-geno_fam_plink_format.txt", sep = "\t", header = T)
        save(famHaplotypes_brca1_geno_plink, file = "cache/famHaplotypes-BRCA1-geno_plink_format.RData")
    }
    
    # Read SNP genomic coordinates
    chr_coords_all <<- read.table("input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed-with_centiMorgan_coords_from_linear_interpolation.txt", header = T, sep = "\t")
    
    # Read Genetic coordinates (centimorgan)
    morgan.brca1 <<- read.table2("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr17_combined_b37.txt", header = T) # phase3 map
    
    # BRCA mutation coordinates
    brca1_mutation_coords <<- read.table("input/139_ouh_june_2017/BRCA1_mut_position.txt", header = T, sep = "\t")
    
    gene <<- "BRCA1"
    
    ### Add project to open projects
    open_projects[[gene]] <<- list(project_path = paste0("cache/", gene, "/"), project = F)
    
    print("Done loading BRCA1 dataset")
}

loadBRCA2 <- function(){
    print("loading BRCA2 dataset...")
    # Read phenotype data
    brca2_pheno <- read.csv("input/139_ouh_june_2017/B2_Onco_phenotype_distribution_180816.csv")
    brca2_pheno_hisp <- read.csv("input/139_ouh_june_2017/B2_Hispanic_OncoArray_phenotypes_230817.csv")
    
    brca2_pheno_hisp <<- brca2_pheno_hisp %>% mutate(Country="Hispanic")
    
    # Suppress warnings of NAs introduced - non-important data
    brca2_pheno_merged <<- suppressWarnings(rbind(brca2_pheno[, 1:47], brca2_pheno_hisp) %>% arrange(FamCode))
    
    # Read genotype data in plink format
    if (!file.exists("cache/brca2_geno_plink.RData")){
        brca2_geno_eu_plink <- read.table2("input/139_ouh_june_2017/139_ouh_brca2_onco_geno_plink_format.txt", header = T)
        brca2_geno_hisp_plink <- read.table2("input/139_ouh_june_2017/139_ouh_female_hisp_brca2_onco_geno_plink_format.txt", header = T)
        
        brca2_geno_plink <- rbind(brca2_geno_eu_plink, brca2_geno_hisp_plink)
        
        save(brca2_geno_plink, file = "cache/brca2_geno_plink.RData")
    }
    load("cache/brca2_geno_plink.RData", envir = .GlobalEnv)
    
    # Load family geno data for brca2
    if (file.exists("cache/famHaplotypes-BRCA2-geno.RData") || file.exists("cache/famHaplotypes-BRCA2-geno_plink_format.RData")){
        load("cache/famHaplotypes-BRCA2-geno_plink_format.RData", envir = .GlobalEnv)
    }else{
        famHaplotypes_brca2_geno <<- read.table(file = "input/139_ouh_june_2017/famHaplotypes-BRCA2-geno.txt", sep = "\t", header = T)
        save(famHaplotypes_brca2_geno, file = "cache/famHaplotypes-BRCA2-geno.RData")
        famHaplotypes_brca2_geno_plink <<- read.table(file = "input/139_ouh_june_2017/famHaplotypes-BRCA2-geno_fam_plink_format.txt", sep = "\t", header = T)
        save(famHaplotypes_brca2_geno_plink, file = "cache/famHaplotypes-BRCA2-geno_plink_format.RData")
    }
    
    # Read SNP genomic coordinates
    chr_coords_all <<- read.table("input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed-with_centiMorgan_coords_from_linear_interpolation.txt", header = T, sep = "\t")
    
    # Read Genetic coordinates (centimorgan)
    morgan.brca2 <<- read.table2("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr13_combined_b37.txt", header = T) # phase3 map
    
    # BRCA mutation coordinates
    brca2_mutation_coords <<- read.table("input/139_ouh_june_2017/BRCA2_mut_position.txt", header = T, sep = "\t")
    
    gene <<- "BRCA2"
    
    ### Add project to open projects
    open_projects[[gene]] <<- list(project_path = paste0("cache/", gene, "/"), project = F)
    
    print("Done loading BRCA2 dataset")
}

loadBRCA1Sim <- function(){
    print("loading BRCA1_simulated dataset...")
    #filename = paste0("classify-simulated-population-age/BRCA1_simulated-starGenealogy-100_samples-1_simulations-generations_100_100_step10-seed_9012/simulated_population/BRCA1_simulated-starGenealogy-100_samples-1_simulations-generations_100-seed_9012")
    #filename = paste0("classify-simulated-population-age/BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_10_500_step10-seed_42/simulated_population/BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_10-seed_42")
    #brca1_simulated_pheno_merged <<- read.table(file = paste0(filename, "-pheno.txt"), header = T)
    #brca1_simulated_geno_plink <<- read.table(file = paste0(filename, "-geno.txt"), header = T)
    
    # Read SNP genomic coordinates
    chr_coords_all <<- read.table("input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed-with_centiMorgan_coords_from_linear_interpolation.txt", header = T, sep = "\t")
    #chr_coords_all <<- read.table("input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed-with_centiMorgan_coords_from_linear_interpolation_hapmap-phase2.txt", header = T, sep = "\t")
    
    # Read Genetic coordinates (centimorgan)
    morgan.brca1 <<- read.table2("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr17_combined_b37.txt", header = T) # phase3 map
    
    # BRCA mutation coordinates
    brca1_mutation_coords <<- read.table("input/139_ouh_june_2017/BRCA1_mut_position.txt", header = T, sep = "\t")
    
    gene <<- "BRCA1_simulated"
    
    ### Add project to open projects
    open_projects[[gene]] <<- list(project_path = paste0("cache/", gene, "/"), project = F)
    
    print("Done loading BRCA1_simulated dataset")
}

loadBRCA2Sim <- function(){
    print("loading BRCA2_simulated dataset...")
    
    gene <<- "BRCA2_simulated"
    
    #filename = paste0("classify-simulated-population-age/BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_10_500_step10-seed_42/simulated_population/BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_100-seed_42")
    #brca2_simulated_pheno_merged <<- read.table(file = paste0(filename, "-pheno.txt"), header = T)
    #brca2_simulated_geno_plink <<- read.table(file = paste0(filename, "-geno.txt"), header = T)
    
    # Read SNP genomic coordinates
    chr_coords_all <<- read.table("input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed-with_centiMorgan_coords_from_linear_interpolation.txt", header = T, sep = "\t")
    
    # Read Genetic coordinates (centimorgan)
    morgan.brca2 <<- read.table2("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr13_combined_b37.txt", header = T) # phase3 map
    
    # BRCA mutation coordinates
    brca2_mutation_coords <<- read.table("input/139_ouh_june_2017/BRCA2_mut_position.txt", header = T, sep = "\t")
    
    ### Add project to open projects
    open_projects[[gene]] <<- list(project_path = paste0("cache/", gene, "/"), project = F)
    
    print("Not yet available")
    
    #print("Done loading BRCA2_simulated dataset")
}

loadProject <- function(project_id){
    if (project_id == ""){
        return(c(1,"Project ID is missing"))
    }
    if (toupper(project_id) == "BRCA1"){
        loadBRCA1()
        return(c(0, paste0("Project loaded successfully")))
    }
    if (toupper(project_id) == "BRCA2"){
        loadBRCA2()
        return(c(0, paste0("Project loaded successfully")))
    } 
    project_path = paste0("cache_projects/", project_id, "/")
    if (!dir.exists(project_path)){
        print(paste0("Project does not exists. You can create a new project using option 3."))
        return(c(1, paste0("Project does not exists. You can create a new project using option 3.")))
    } else {
        ### Read data
        project_pheno <<- readRDS(file = paste0(project_path, project_id, "_pheno.RData"))
        project_geno <<- readRDS(file = paste0(project_path, project_id, "_geno.RData"))
        project_coords <<- readRDS(file = paste0(project_path, project_id, "_coords.RData"))
        project_mut_info <<- readRDS(file = paste0(project_path, project_id, "_mut_info.RData"))
        
        gene <- project_mut_info$gene[1]
        gene <<- gene
        
        ### Add project to open projects
        open_projects[[gene]] <<- list(project_path = project_path,
                                      project_pheno = project_pheno,
                                      project_geno = project_geno,
                                      project_coords = project_coords,
                                      project_mut_info = project_mut_info,
                                      project = T)
        
        return(c(0, paste0("Project loaded successfully")))
    }
}

createProject <- function(project_id, pheno_file, geno_file, coords_file, mut_info_file){
    print(paste(project_id, pheno_file$datapath, geno_file$datapath, coords_file$datapath, mut_info_file$datapath))
    if (project_id == ""){
        return(c(1,"Project ID is missing"))
    }
    # if (gene == ""){
    #     return(c(1,"Gene name is missing"))
    # }
    if (length(pheno_file) == 0){
        return(c(1,"Phenotype file is missing"))
    }
    if (length(geno_file) == 0){
        return(c(1,"Genotype file is missing"))
    }
    if (length(coords_file) == 0){
        return(c(1,"Coordinates file is missing"))
    }
    if (length(mut_info_file) == 0){
        return(c(1,"Coordinates file is missing"))
    }
    
    # Deal with starting with number
    project_path = paste0("cache_projects/", project_id, "/")
    if (dir.exists(project_path) || toupper(project_id) %in% c("BRCA1","BRCA2")){
        print("Project already exists. Try loading it using the 'Load Project' tool.")
        return(c(1,"Project already exists. Try loading it using the 'Load Project' tool."))
    } else {
        
        ### Load data - supported formats: txt (tabs- or space-seperated), tsv, csv
        pheno_ex <- strsplit(basename(pheno_file$datapath), split="\\.")[[1]][2]
        if (pheno_ex %in% c("txt", "tsv")){
            project_pheno <<- read.table(file = pheno_file$datapath, header = T, fill = T, sep = "\t", quote = "")
        } else if (pheno_ex == "csv"){
            project_pheno <<- read.csv(file = pheno_file$datapath, header = T, fill = T)
        } else if (pheno_ex %in% c("RData", "rdata", "Rdata")){
            project_pheno <<- readRDS(file = pheno_file$datapath)
        } else {
            return(c(1, "Failed to read phenotype file. File format not supported."))
        }
        geno_ex <- strsplit(basename(geno_file$datapath), split="\\.")[[1]][2]
        if (geno_ex %in% c("txt", "tsv")){
            project_geno <<- read.table(file = geno_file$datapath, header = T, sep = "\t", quote = "")
        } else if (geno_ex == "csv"){
            project_geno <<- read.csv(file = geno_file$datapath, header = T)
        } else if (geno_ex %in% c("RData", "rdata", "Rdata")){
            project_geno <<- readRDS(file = geno_file$datapath)
        } else {
            return(c(1, "Failed to read genotype file. File format not supported."))
        }
        coords_ex <- strsplit(basename(coords_file$datapath), split="\\.")[[1]][2]
        if (coords_ex %in% c("txt", "tsv")){
            project_coords <<- read.table(file = coords_file$datapath, header = T, sep = "\t", quote = "")
        } else if (coords_ex == "csv"){
            project_coords <<- read.csv(file = coords_file$datapath, header = T)
        } else if (coords_ex %in% c("RData", "rdata", "Rdata")){
            project_coords <<- readRDS(file = coords_file$datapath)
        } else {
            return(c(1, "Failed to read coordinates file. File format not supported."))
        }
        mut_info_ex <- strsplit(basename(mut_info_file$datapath), split="\\.")[[1]][2]
        if (mut_info_ex %in% c("txt", "tsv")){
            project_mut_info <<- read.table(file = mut_info_file$datapath, header = T, stringsAsFactors = F)
        } else if (mut_info_ex == "csv"){
            project_mut_info <<- read.csv(file = mut_info_file$datapath, header = T, stringsAsFactors = F)
        } else if (mut_info_ex %in% c("RData", "rdata", "Rdata")){
            project_mut_info <<- readRDS(file = mut_info_file$datapath)
        } else {
            return(c(1, "Failed to read coordinates file. File format not supported."))
        }
        
        ### Check files
        if (sum(c("Onc_ID", "FamCode", "Country", "Mut1HGVS") %in% names(project_pheno)) != 4){
            if (sum(c("Onc_ID", "Mut1HGVS") %in% names(project_pheno)) != 2){
                return(c(1, "Phenotype file does not contain the headers: Onc_ID, Mut1HGVS."))
            }
            if (!("FamCode" %in% names(project_pheno))){
                project_pheno$FamCode = project_pheno$Onc_ID
            }
            if (!("Country" %in% names(project_pheno))){
                project_pheno$Country = "No country info"
            }
        }
        if (names(project_geno)[1] != "SNP"){
            return(c(1, "Genotype file does not have 'SNP' header as first column"))
        }
        if (sum(c("SNP", "Chr_numeric", "position_b37")  %in% names(project_coords)) != 3){
            return(c(1, "Coordinates file does not contain the headers: SNP, Chr_numeric, position_b37"))
        }
        
        if (sum(c("Mut1HGVS", "gene", "chr", "start", "stop")  %in% names(project_mut_info)) != 5){
            return(c(1, "Mutation info file does not contain the headers: Mut1HGVS, gene, chr, start, stop"))
        }
        gene <- unique(project_mut_info$gene)
        if (length(gene) != 1 || length(unique(project_mut_info$chr)) != 1) {
            return(c(1, "Mutation info file may not contain information for more than one gene (or chromosome)."))
        }
        if (sum(project_mut_info$Mut1HGVS %in% project_pheno$Mut1HGVS) == 0){
            # return(c(1, "No overlap of mutation names in phenotype file and mutation info file!"))
            project_pheno$Mut1HGVS = project_mut_info$Mut1HGVS
        }
        if (sum(project_geno$SNP %in% project_pheno$Onc_ID) == 0){
            return(c(1, "No overlap of samples names in phenotype file (Onc_ID) and genotype file (SNP)!"))
        }
        # Add cM info if missing
        if (sum(c("cM") %in% names(project_coords)) != 1){
            chr = project_mut_info$chr
            project_coords <- filter(project_coords, Chr_numeric == chr)
            morgan.coords <<- read.table2(paste0("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr",chr,"_combined_b37.txt"), header = T) # phase3 map
            project_coords$cM <- approx(x = morgan.coords[,1], y = morgan.coords[,3], xout = project_coords$position_b37)$y
            project_coords <<- project_coords
        }
        
        # Make gene global
        gene <<- gene
        
        
        ### Create cache structure
        createCacheStructure(project_path)
        
        ### Save data
        saveRDS(project_pheno, file = paste0(project_path, project_id, "_pheno.RData"))
        saveRDS(project_geno, file = paste0(project_path, project_id, "_geno.RData"))
        saveRDS(project_coords, file = paste0(project_path, project_id, "_coords.RData"))
        saveRDS(project_mut_info, file = paste0(project_path, project_id, "_mut_info.RData"))
        
        ### Add project to open projects
        open_projects[[gene]] <<- list(project_path = project_path,
                                      project_pheno = project_pheno,
                                      project_geno = project_geno,
                                      project_coords = project_coords,
                                      project_mut_info = project_mut_info,
                                      project = T)

        
        return(c(0, "Project created successfully"))
    }
}

createCacheStructure <- function(project_path){
    ### Create cache structure
    dir.create(project_path, showWarnings = F)
    dir.create(paste0(project_path, "breakpointPlots"), showWarnings = F)
    dir.create(paste0(project_path, "breakpointPlotsCountry"), showWarnings = F)
    dir.create(paste0(project_path, "ClusterValidationPlots"), showWarnings = F)
    dir.create(paste0(project_path, "ClusterValidationValues"), showWarnings = F)
    dir.create(paste0(project_path, "distMatrices"), showWarnings = F)
    dir.create(paste0(project_path, "haplotypePlots"), showWarnings = F)
    dir.create(paste0(project_path, "haplotypePlotsCountry"), showWarnings = F)
    dir.create(paste0(project_path, "mutationAge_country"), showWarnings = F)
    dir.create(paste0(project_path, "mutationAge_groups"), showWarnings = F)
    dir.create(paste0(project_path, "nearestBreakpointPlots"), showWarnings = F)
    dir.create(paste0(project_path, "nearestBreakpointPlotsCountry"), showWarnings = F)
    dir.create(paste0(project_path, "nearestBreaksStatistics"), showWarnings = F)
}
