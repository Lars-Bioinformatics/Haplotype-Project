library(readxl)
library(dplyr)
library(forcats)
library(fastmatch)
library(ggplot2)

setwd("~/Documents/PhD/haplotype-project")

### Input: mutation name, or vector of mutation names ###
run <- function(mut_name, brca_name){
    "
    Run for BRCA1 or BRCA2
    Input: 
        - mut_name:  mutation name, or vector of mutation names
        - brca_name: the name of BRCA gene containing the mut_name mutation

    note: '<<-' makes global variables, convinient for testing purposes.
    "
    # Read data
    if (brca_name == "BRCA1"){
        # Load onco information
        brca_data <- read_excel("input/139_ouh_june_2017/Onco_pheno_BC1.xlsx")
        
        # Load SNP data
        if(file.exists("cache/haplotypes_brca1.RData")){
            load("cache/haplotypes_brca1.RData")
            haplotypes_brca <<- haplotypes_brca
        } else {
            haplotypes_brca <<- read.table("input/139_ouh_june_2017/139_ouh_brca1_onco_geno.txt", header = T)
            save(haplotypes_brca, file = "cache/haplotypes_brca1.RData")
        }
        
    } else { # BRCA2
        # Load onco information
        brca_data <- read_excel("input/139_ouh_june_2017/Onco_pheno_BC2.xlsx")
        
        # Load SNP data
        if(file.exists("cache/haplotypes_brca2.RData")){
            load("cache/haplotypes_brca2.RData")
            haplotypes_brca <<- haplotypes_brca
        } else {
            haplotypes_brca <<- read.table("input/139_ouh_june_2017/139_ouh_brca2_onco_geno.txt", header = T)
            save(haplotypes_brca, file = "cache/haplotypes_brca2.RData")
        }
    }

    chr_coords <<- read.table("input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position.txt", header = T, sep = ",")
    
    # BRCA gene information
    if (brca_name == "BRCA1"){
        brca_name <<- "BRCA1"
        brca_start <<- 41197695; brca_stop <<- 41276113
        brca_middle <<- (brca_stop - brca_start)/2 + brca_start
        brca_chr <<- 17
    } else { # BRCA2
        brca_name <<- "BRCA2"
        brca_start <<- 32890598; brca_stop <<- 32972907
        brca_middle <<- (brca_stop - brca_start)/2 + brca_start
        brca_chr <<- 13
    }

    # all_mut_names_brca <<- unique(brca_data$Mut1HGVS)
    # mut_name <- sample(all_mut_names_brca, 20)
    # mut_name <- all_mut_names_brca

    # Required number of families
    num.fam = 1
    
    # Compute haplotype breaks for each mutation
    for (mut in mut_name){
        brca_data_mut <- filter(brca_data, Mut1HGVS %in% mut)
        if (length(unique(brca_data_mut$FamCode)) >= num.fam){
            print(length(unique(brca_data_mut$FamCode)))
            haplotype_mutation(brca_data_mut, haplotypes_brca, chr_coords, brca_start,
                               brca_stop, brca_middle, brca_chr, brca_name, mut)
        }
    }
}

## ---- haplotype_mutation
haplotype_mutation <- function(brca_data, haplotypes_brca, chr_coords,
	brca_start, brca_stop, brca_middle, brca_chr, brca_name, mut,
	country = NULL, distinct = F){
	"
	Haplotype samples with a given brca mutation
	"

    # Prints for checking output
    print(dim(brca_data)) # 8x67
    print(dim(haplotypes_brca)) # 15679x12468
    print(dim(chr_coords)) # 12467x3
    print(mut) # c.1016delA

    # Make brca_data a global variable (for testing and debugging)
    brca_data_temp <<- brca_data
    brca_data <<- brca_data_temp

	# Filter out the SNPs needed i.e. on current chromosome
	chr_coords = filter(chr_coords, chr_coords$Chr_numeric == brca_chr)

    # Extract samples with given mutation,
	# and remove SNPs not related to actual BRCA gene
    matched <<- extractSamples(brca_data, haplotypes_brca, chr_coords, distinct)
    #   matched <<- extractSamples(brca_data, haplotypes_brca, chr_coords, F)
    # matched <- matched[, c(1, na.omit(fmatch(chr_coords$SNP, colnames(matched))))]

    # If the haplotype should be defined across all samples
    if (is.null(country)){
    	haplotypes <<- findConsensus(matched)
    # If haplotype should be based on samples from specific country
    } else {
    	haplotypes <<- findConsensus(filter(matched,
    		SNP %in% subset(brca_data, Country == "SPAIN")$Onc_ID))
    }
    #haplotypes <<- haplotypes

    # Find SNPs breaking the haplotype
    breaks <<- findHaplotypeBreaks(matched, haplotypes)

	# If no samples breaking haplotype, then stop analysis of current mutation
    if (length(breaks) == 0){
        print(paste("No haplotype breaks found in families for mut:", mut))
        return()
    }

	# Map found breaks (SNPs) to their genomic position
    chr_pos <<- mapSNPs(breaks, chr_coords, brca_middle,
						brca_start, brca_stop, brca_chr, brca_name)
    
    # Save information required to plot and find nearest breaks
    save(chr_pos, file = paste0("cache/chr_pos-", brca_name, "-", paste(mut, collapse = ","), ".RData"))
    # Load information, if already computed (SHOULD BE MOVED BEFORE COMPUTATIONS)
    # load("cache/chr_pos-BRCA1-c.1016delA.RData")

    # Find the nearest break on each side of the genome
    # and compute the mean length of the breaks for each country
    #### haplotypeDistance <- findNearestBreaks(chr_pos, mut)
    
	# Visualize the breaks around the BRCA gene
    p <- plot_haplotype(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut)
    
    # ggplot(haplotypeDistance, aes(x=sample_id, y=positive_pos)) + geom_bar()
    
    return(p)
        
}

## ---- extractSamples
extractSamples <- function(brca_data, haplotypes_brca, chr_coords, distinct = F){
    "
    Extract relevant samples from table of all genotypes
    "

    # Keep only one member of a family and remove the others, if chosen
    if (distinct){
        brca_data <- brca_data[!duplicated(brca_data$FamCode), ]
        # Equivalent to above
        # brca_data <- subset(brca_data, !duplicated(FamCode))
        # brca_data <- distinct(brca_data, FamCode, .keep_all = T)
    }

    # Extract samples (with given brca mutation) from the big haplotypes dataset
    # Faster to use fmatch than match
    index <- fmatch(brca_data$Onc_ID, haplotypes_brca$SNP)
    matched <- haplotypes_brca[index, ]

    # as above, but slower
    # matched <- filter(haplotypes_brca, SNP %in% brca_data$Onc_ID)
    print(dim(matched))
	# Only keep SNPs (columns) that are related to current chromosome
    matched <- matched[, c(1, na.omit(fmatch(chr_coords$SNP, colnames(matched))))]

	print(dim(matched))
    return(matched)
}

## ---- findConsensus
findConsensus <- function(matched){
    "
    Find haplotype consensus for each SNP
    Haplotype consensus is defined as the most
    frequent letter (A,C,G,T,I,D) in all samples
    "
    print("Find haplotype consensus for each SNP")
    haplotypes <- c("SNP")
    for (pos in 2:ncol(matched)){
        if (pos %% 500 == 0) print(paste("Process:", pos, "of", ncol(matched), "SNPs"))

        # Combine a SNP column to one string
        seq <- paste0(matched[, pos], collapse = "")
        # Split string into one-character vector
        seq2 <- strsplit(seq, split = "")[[1]]
        # Counts the number of occurences of each element
        counts <- table(seq2)
        # print(counts)

        # Skip the '-' character - as it's unknown value
        start_index = 1
        if (names(counts[1]) == "-") {
            start_index = 2
        }

        # List of candidates for most frequent haplotype
        candidates = names(counts[start_index])
        # Highest number of occurences
        best = counts[start_index]

        # if only '-' character in string, return no consensus
        if (length(counts) == 1 && names(counts[1]) == "-"){
            new = "FAIL"
            names(new) <- names(matched[pos])
            haplotypes <- c(haplotypes, new)
            next
        }
        # If only one element (except -), return that as consensus
        if (length(counts) == start_index){
            new = candidates
            names(new) <- names(matched[pos])
            haplotypes <- c(haplotypes, new)
            next
        }

        # Check if more than one element has the highest frequency
        for (i in (start_index+1):length(counts)){
            if (counts[i] == best){
                candidates = c(candidates, names(counts[i]))
            }
            if (counts[i] > best){
                candidates = names(counts[i])
                best = counts[i]
            }
            # print(candidates)
        }

        # If only one candidate, return as consensus
        if (length(candidates) == 1){
            new = candidates
            names(new) <- names(matched[pos])
            haplotypes <- c(haplotypes, new)
        # When multiple possible haplotypes, return FAIL
        } else {
            new = "FAIL"
            names(new) <- names(matched[pos])
            haplotypes <- c(haplotypes, new)
        }
    }
    # Return list representing consensus haplotype
    return(haplotypes)
}

## ---- findHaplotypeBreaks
findHaplotypeBreaks <- function(matched, haplotypes){
    breaks <- list()
    for (pos in 2:ncol(matched)){

        if (pos %% 500 == 0) print(paste("Process:", pos, "of", ncol(matched), "SNPs"))

        haplotype = haplotypes[pos]
        # print(haplotype)
        if (haplotype == "FAIL"){
            next
        }

        for (i in 1:nrow(matched)){
            # i = 1
            #row <- 1
            row <- strsplit(as.character(matched[i,pos]), split = "")[[1]]
            #row <- c("A", "A")
            # print(row)
            if (row[1] == "-" || row[2] == "-"){
                next
            }
            if (row[1] != haplotype && row[2] != haplotype){
                #print(row)
                #print(haplotype)
                index = as.character(matched[i,1])
                breaks[[index]] = c(breaks[[index]], colnames(matched[pos]))
            }
        }
    }
    print(paste("Num samples:", length(breaks)))
    print(paste("Num breaks:", length(breaks[[1]])))
    return(breaks)
}

## ---- mapSNPs
### Map found SNPs to their genomic position ###
mapSNPs <- function(breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name){
#mapSNPs <- function(){
    chr_pos = NULL
    # For each sample (patient) with given mutation
    for (i in 1:length(breaks)){
        # Extract rows with SNPs in brca1 chromosome (17) and where SNP is homozygote for that sample
        #chr_pos2 = filter(chr_coords, chr_coords$Chr_numeric == brca_chr, chr_coords$SNP %in% breaks[[i]])
		chr_pos2 = filter(chr_coords, chr_coords$SNP %in% breaks[[i]])
		print(dim(chr_pos2))
        # Removes SNPs in the brca1 gene
        chr_pos2 = filter(chr_pos2, as.numeric(as.character(position_b37)) < brca_start | as.numeric(as.character(position_b37)) > brca_stop)
        # Add sample_id column
        chr_pos2["sample_id"] = names(breaks[i])
        # Add position_b37_adjusted column, where brca1 = 0
        chr_pos2 = mutate(chr_pos2, position_b37_adjusted = as.numeric(as.character(position_b37)) - brca_middle)
        # Combine all samples into big data frame for plotting (likely not most efficient way)
        chr_pos = rbind(chr_pos, chr_pos2)
    }
    # Add columns Country and FamCode to data frame
    index <- fmatch(chr_pos$sample_id, brca_data_temp$Onc_ID)
    chr_pos[c("Country", "FamCode")] <- brca_data_temp[index, ] %>% select(Country, FamCode)

    # Print some data information - here counting samples by Country
    print(chr_pos %>% distinct(sample_id, FamCode, Country) %>% count(Country))
    print(dim(chr_pos))
    # unique(chr_pos$sample_id)
    # print(chr_pos)

    return(chr_pos)
}

## ---- plot
### Visual output of haplotypes and breaks ###
plot_haplotype <- function(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut){

    # Define y axis step distance
    tick_steps = 3000000
    # Set the adjusted beginning and ending position of brca gene 
    brca_max = (brca_stop-brca_start)/2; brca_min = -brca_max
    # Number of samples
    num_samples = length(unique(chr_pos$sample_id))
    
    # Plot countries together that has less than 5% samples of total datasets
    p = chr_pos %>% mutate(Country = fct_lump(Country, prop = 0.05)) %>% # group_by(Country) %>%
        
        # Set up the ggplot, with x sorted by most frequent country, and sample_id as tie breaker
        # The genomic position defines the y axis
        ggplot(aes(x=interaction(sample_id, fct_infreq(Country)), y=position_b37_adjusted)) +
        
        # Plots the breakpoints around the brca gene
        geom_point(aes(color = Country), size = 0.8) +
        #geom_point(size = 0.8, color = "darkblue") +
        
        # Defines the color of the dots
        scale_color_brewer(palette="Set1") +
        
        # Plots the area of the brca gene
        #geom_hline(yintercept = 0, color = "red", size = 1) +
        geom_hline(aes(yintercept = brca_max), color = "darkblue", size = 0.4) +
        geom_hline(aes(yintercept = brca_min), color = "darkblue", size = 0.4) +
        geom_rect(aes(xmin = 0, xmax = num_samples+1, ymin=brca_min, ymax=brca_max), fill="darkblue") +
        
        # Inserts the plot title and centers it
        ggtitle(paste("Mutation:", paste(mut, collapse = ","))) + # Adds plot title
        theme(plot.title = element_text(hjust = 0.5)) + # Centering of plot title
        
        # Change the angle of the x-labels
        # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # 45 degree x-labels
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  # vertical x-labels
        
        # Define x-ticks and their label names
        scale_x_discrete(breaks=interaction(chr_pos$sample_id, chr_pos$Country), labels = chr_pos$sample_id) +
        #scale_x_discrete(breaks=interaction(chr_pos$sample_id, chr_pos$Country)) +

        # Define y-ticks and their label names (sets BRCA name in y axis)
        scale_y_continuous(breaks = c(seq(from = tick_steps, to = 5*tick_steps, by = tick_steps),
                                      seq(from = -tick_steps, to = -5*tick_steps, by = -tick_steps), 0),
                           labels = c(seq(from = tick_steps, to = 5*tick_steps, by = tick_steps),
                                      seq(from = -tick_steps, to = -5*tick_steps, by = -tick_steps), brca_name)) +
        
        
        # Define visible y-axis
        #ylim(-1000000,1000000) + 
        
        # Set x and y axis label names
        ylab("Genomic position") + xlab("Sample")
    
    print(p)
    # print(p + ylim(-1000000,1000000))

    # Save plots to SVG and PNG files
    # ggsave(paste0("plots/", brca_name, "-", paste(mut, collapse = ","), "-", length(num_samples), "_samples.svg"))
    ggsave(paste0("plots/", brca_name, "-", paste(mut, collapse = ","), "-", num_samples, "_samples.png"))
    
    return(p)
}

## ---- findNearestBreaks
findNearestBreaks <- function(chr_pos, mut){
    # Extract (positive position) breakpoint closest to brca gene
    positive <- chr_pos %>%
        group_by(sample_id) %>%
        filter(position_b37_adjusted > 0) %>%
        slice(which.min(position_b37_adjusted))
    # Extract (negative position) breakpoint closest to brca gene
    negative <- chr_pos %>%
        group_by(sample_id) %>%
        filter(position_b37_adjusted < 0) %>%
        slice(which.max(position_b37_adjusted))

    # Merge breakpoints into one data frame
    haplotypeDistance <- merge(select(positive, sample_id, Country, positive_pos = position_b37_adjusted),
                               select(negative, sample_id, negative_pos = position_b37_adjusted))
    # Compute haplotype length for each sample
    haplotypeDistance <- haplotypeDistance %>% mutate(haplo_distance = positive_pos + abs(negative_pos))

    # Arrange samples according to plot i.e. most frequent country samples first
    ordered_haploDist <- haplotypeDistance %>% arrange(fct_infreq(Country), sample_id)

    # Compute mean haplotype length for all countries
    ordered_mean_table <- haplotypeDistance %>%
        group_by(Country) %>%
        #summarize_each(funs(mean(., na.rm = TRUE)), haplo_distance) %>%
        #summarise_at(.vars=c("haplo_distance"), .funs=c(Mean="mean")) %>%
        summarise(mean_haplo_distance=mean(haplo_distance, na.rm=TRUE)) %>%
        mutate(mean_haplo_distance = round(mean_haplo_distance)) %>%
        arrange(desc(mean_haplo_distance))

    # Write tables to file
    post_filename <- paste0(brca_name, "-", paste(mut, collapse = ","), "-", length(breaks), "_samples.txt")
    write.table(ordered_haploDist, file = paste0("distance_statistics/haplotype_dist-", post_filename),
                sep = "\t", col.names = T, quote = F, row.names = F)
    write.table(ordered_mean_table, file = paste0("distance_statistics/haplotype_dist_mean-", post_filename),
                sep = "\t", col.names = T, quote = F, row.names = F)

    print("DONE")
    
    return(haplotypeDistance)
}


countHetAndHomo <- function(matched){
    het = 0
    homo = 0
    #row=1
    #pos = 2
    for (row in 1:nrow(matched)){
        het = 0
        homo = 0
        for (pos in 2:ncol(matched)){
            #if (pos %% 500 == 0) print(paste("Process:", pos, "of", ncol(matched), "SNPs"))
        
            # Split string into one-character vector
            seq <- strsplit(paste0(matched[row, pos]), split = "")[[1]]
            
            if (seq[1] == "-" || seq[2] == "-"){
                next
            }
            
            if (seq[1] == seq[2]) {
                homo = homo + 1
            } else {
                het = het + 1
            }
            
            
        }
        print(rownames(matched[row,]))
        print(het)
        print(homo)
        print("")
    }
}
countHetAndHomo(matched)

### Run BRCA1 samples ###
#brca_name = "BRCA1"
# mut_name = "c.68_69delAG"
# mut_name = "c.1016delA"
# mut_name = "c.1016dupA"
# mut_name = "c.5266dupC"
# # mut_name = c("c.1016delA", "c.1016dupA")
#mut <- mut_name
mut = "c.1016delA"
#run_brca1(mut_name
run("c.1016delA", "BRCA1")
#run_brca1(mut_name)


### Run BRCA2 samples ###
#brca_name = "BRCA2"
#mut_name <- "c.1128delT"
# mut_name <- "c.7617+1G>A"
# mut <- mut_name
# run_brca2(mut_name)
run("c.1128delT", "BRCA2")
