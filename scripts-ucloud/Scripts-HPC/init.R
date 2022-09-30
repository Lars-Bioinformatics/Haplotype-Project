setwd("~/Documents/PhD/haplotype-project")

library(shiny)

source("Scripts/Haplotype-projekt-ver2.R")
source("ScriptsClustering.R")

## Read input data into global variables
readData()
readFamData()

runApp()

# init mutation
mut = "c.7617+1G>A"; gene = "BRCA2"
readIndvAndFamData(gene, mut)
