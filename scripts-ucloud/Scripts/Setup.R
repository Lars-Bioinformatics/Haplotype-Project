# Setup for haplotype project

# Installs packages needed - Note: these commanda will overwrite already installed packages with the same name
packages <- c("readxl", "dplyr", "forcats", "fastmatch", "ggplot2", "amap", "ggdendro",
              "dendextend", "RColorBrewer", "shiny", "shinycssloaders", "shinyjs", "ggpubr",
              "proxy", "pvclust", "factoextra", "NbClust", "rmarkdown", "shinythemes",
              "gridExtra", "gtable", "dbscan", "DT", "data.table", "e1071", "randomForest", 
              "Metrics", "caret", "tsne", "doParallel", "statmod", "tweedie", "doSnow", "boot")

# Install for local users
install.packages(packages)


### NOT FINISHED ###
# Install for all users (Ubuntu)
# rcmd = paste0("R -e \"install.packages('", packages,"', repos='https://cran.rstudio.com/')\"")
# paste0("sudo su - -c ", rcmd)
# system(paste0("sudo su - -c ", rcmd))
