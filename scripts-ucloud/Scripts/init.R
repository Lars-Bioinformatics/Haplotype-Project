setwd("~/Documents/PhD/haplotype-project")

library(shiny)

source("Scripts/Haplotype-projekt-ver2.R")
source("Scripts/Clustering.R")
source("Scripts/Mutation-Age.R")
source("Scripts/loadShinyData.R")
source("Scripts/Gandolfo_Speed_Mutation_Age_Estimation.R")
source("Scripts/Mutation-age-tripleA.R")
source("Scripts/maxLikelihood.R")

## Read input data into global variables
system.time(readData())
system.time(readFamData())
#readSimulatedData()

runApp("Scripts/", host = "0.0.0.0")

# init mutation
mut = "c.7617+1G>A"; gene = "BRCA2"
readIndvAndFamData(gene, mut)
