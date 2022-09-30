library(shiny)
library(shinycssloaders)
library(shinythemes)
library(shinyjs)
library(rmarkdown)
library(DT)
#library(svglite)
#library(foreach)
#library(doParallel)


useSubCluster = F # should be removed sometime

if (Sys.info()["sysname"] == "Linux"){
    ### Deployment setup ###

    source("Haplotype-projekt-ver2.R")
    source("Clustering.R")
    source("Mutation-Age.R")
    source("loadShinyData.R")
    source("Gandolfo_Speed_Mutation_Age_Estimation.R")
    source("Mutation-age-tripleA.R")
    source("maxLikelihood.R")

    # Read input data into global variables
    # readData()
    # readFamData()

    ### END Deployment setup ###
} else { # Darwin
    ### Development setup ###
    setwd("~/Documents/PhD/haplotype-project")

    source("Scripts/Haplotype-projekt-ver2.R")
    source("Scripts/Clustering.R")
    source("Scripts/Mutation-Age.R")
    source("Scripts/loadShinyData.R")
    source("Scripts/Gandolfo_Speed_Mutation_Age_Estimation.R")
    source("Scripts/Mutation-age-tripleA.R")
    source("Scripts/maxLikelihood.R")
    
    #loadBRCA1()
    #loadedGenes <<- c("BRCA1")

    #readSimulatedData()

    ### END Development setup ###
}

# Load toy data
initShinyData()

ui <- navbarPage("", inverse = F, windowTitle = "MutFounder", theme = shinytheme("flatly"), 
                 tabPanel("MutFounder: Haplotype Analysis of Mutations"),
    #fluidPage(theme = shinytheme("flatly"),#("lumen"),#("yeti"),#("united"),#("readable"),#("paper"),#("lumen"),#("cosmo"),#("flatly"), ("cerulean"),
    #shinythemes::themeSelector(),
    
    tags$head(
        # h1{color: #28699C;}
        # h2{color: #28699C;}
        # h3{color: #28699C;}
        # h4{color: #28699C;}
        # #sidebar {
        # background-color: #02274e;
        #     color: white;
        # }
        # header {
        #     background-color: #02274e;
        #     color: white;
        # }
        # #tabsetPanel {
        # color: #28699C;
        #     }
        # .navbar-default .navbar-nav > .active{
        #     background: #2f3e4e;
        # }
        # green progress bar
        # .progress-bar{background-color:#18BC9C;}
        # .shiny-progress-notification .progres-bar{
        #     color: blue;
        #     background-color: red;
        # }   
        
        # default flatly theme colors for navigation panel
        # .navbar-default .navbar-nav > .active > a,
        # .navbar-default .navbar-nav > .active > a:hover,
        # .navbar-default .navbar-nav > .active > a:focus {
        #     background: #2f3e4e;
        # }
        tags$style(HTML("

            h1{color: #28699C;}
            h2{color: #28699C;}
            h3{color: #28699C;}
            h4{color: #28699C;}
    
            .navbar-brand, .navbar-nav li a {
                line-height: 80px;
                height: 80px;
                padding-top: 0px;
                font-size: 46px;
                padding-left: 15px;
            }

            .navbar-default .navbar-nav > .active > a,
            .navbar-default .navbar-nav > .active > a:hover,
            .navbar-default .navbar-nav > .active > a:focus {
                background: #003851;
            }
            
            .navbar-default {
                background-color: #003851;
            }

            .btn-default {
                color: white;
                background-color: #4c8af5;
                border-color: #4c8af5;
            }
            .btn-default:hover, .btn-default:focus, .btn-default:active{
                color: white;
                background-color: #0c53cf;
                border-color: #0c53cf;
            }

            .progress-bar{background-color:#0c53cf;}

            .tabbable > .nav > li > a {
                
            }

            #clust_info {
                font-style: bold;
            }
            

            pre {
                color: #2f3e4e;
            }

            #sidebar {
                padding-top: 0px;
            }

            ul {
                padding: 0px;
                margin: 0px;
                border: 0px;
            }
        "))
    ),
    
    # tags$head(tags$style(type = 'text/css', '.navbar { background-color: #262626;
    #                                            font-family: Arial;
    #                                            font-size: 13px;
    #                                            color: #FF0000; }')),
    
    #titlePanel("Haplotype Analysis of BRCA1/2 Mutations"),
    #headerPanel(title = "Haplotype Analysis of Mutations", windowTitle = "Haplotype Analysis"),
    #navbarPage("Haplotype Analysis of Mutations", inverse = F, id = "navbar"),
    
    sidebarLayout(
        sidebarPanel(id = "sidebar",
            h3("Choose Mutation"),
            # selectInput("mut",
            #     label = "Select mutation:",
            #     #choices = unique(brca1_pheno_merged$Mut1HGVS),
            #     choices = list("BRCA1"=count(brca1_pheno_merged, Mut1HGVS) %>% filter(n>100) %>% arrange(n) %>% .$Mut1HGVS,
            #                    "BRCA2"=count(brca2_pheno_merged, Mut1HGVS) %>% filter(n>100) %>% arrange(n) %>% .$Mut1HGVS),
            #     selected = "c.5123C>A"),
            
            selectInput("gene",
                        label = "Select gene:",
                        # choices = list("BRCA1",
                        #                "BRCA2",
                        #                "BRCA1_simulated"),
                        #                "BRCA2_simulated"),
                        choices = loadedGenes,
                        #selected = "BRCA1"),
                        #selected = "BRCA1"),
                        selected = loadedGenes[1]),
            
            #uiOutput("selectGene"),
            #uiOutput("selectMut"),

            selectInput("mut",
                label = "Select mutation:",
                choices = "c.7617+1G>A",
                selected = "c.7617+1G>A"),
            
            # selectInput("mut",
            #     label = "Select mutation:",
            #     choices = get_muts(),
            #     selected = get_selected_mut()),
            # conditionalPanel(
            #     condition = "input.gene=='BRCA1'",
            #     selectInput("mut1",
            #                 label = "Select mutation:",
            #                 #choices = count(brca1_pheno_merged, Mut1HGVS) %>% filter(n>5) %>% arrange(Mut1HGVS) %>% .$Mut1HGVS,
            #                 choices = count(brca1_pheno_merged, Mut1HGVS) %>% filter(n>100) %>% arrange(n) %>% .$Mut1HGVS,
            #                 #choices = c("c.68_69delAG", "c.5123C>A"),
            #                 selected = "c.68_69delAG")
            #                 #selected = "c.1016delA")
            # ),
            # conditionalPanel(
            #     condition = "input.gene=='BRCA2'",
            #     selectInput("mut2",
            #                 label = "Select mutation:",
            #                 #choices = count(brca2_pheno_merged, Mut1HGVS) %>% filter(n>5) %>% arrange(n) %>% .$Mut1HGVS,
            #                 choices = count(brca2_pheno_merged, Mut1HGVS) %>% filter(n>100) %>% arrange(n) %>% .$Mut1HGVS,
            #                 #choices = c("c.7617+1G>A", "c.1310_1313delAAGA"),
            #                 selected = "c.7617+1G>A")
            #                 #selected = "c.2808_2811delACAA")
            #                 #selected = "c.1310_1313delAAGA")
            # ),
            # conditionalPanel(
            #     condition = "input.gene=='BRCA1_simulated'",
            #     selectInput("mut3",
            #                 label = "Select mutation:",
            #                 choices = paste0("mutSim_BRCA1_",rep(seq(10,500,10), each=10),"_", 1:10),
            #                 selected = paste0("mutSim_BRCA1_100_", 1))
            # ),
            # conditionalPanel(
            #     condition = "input.gene=='BRCA2_simulated'",
            #     selectInput("mut4",
            #                 label = "Select mutation:",
            #                 choices = paste0("mutSim_BRCA2_", 1:20),
            #                 selected = paste0("mutSim_BRCA2_", 1))
            # ),
            h3("Grouping Data"),
            radioButtons("grouping",
                         label = "Group data by:",
                         choices = list('Clustering', 'Country'),
                         selected = 'Clustering'),
            selectInput("fam",
                        label = "Combine family haplotypes or consider individually:",
                        choices = list('Family' = TRUE,
                                       'Individually' = FALSE),
                        #selected = FALSE),
                        selected = TRUE),
            h3("Clustering"),
            selectInput("distMethod",
                        label = "Select distance metric:",
                        choices = list('firstBreakDist',
                                       'firstBreakDist_cM',
                                       "firstUpperBreakDist",
                                       # "firstUpperBreakDist_cM",
                                       "firstUnderBreakDist",
                                       # "firstUnderBreakDist_cM",
                                       #"firstBreakEuclideanDist",
                                       #"firstBreakEuclideanDist_cM",
                                       "matchstates",
                                       "genpofad",
                                       "mrca",
                                       "nei",
                                       "euclidean weighted according to brca dist"="euclideanDist_BRCAdistWeighted",
                                       'euclidean',
                                       "centered pearson"="correlation",
                                       "pearson"
                                       #"canberra",
                                       #"spearman",
                                       #"manhattan",
                                       #"kendall"
                                       ),
                        selected = 'firstBreakDist'),
            selectInput("clustMethod",
                        label = "Clustering method:",
                        choices = list('complete',
                                       # 'ward.D',
                                       'ward'='ward.D2',
                                       'single',
                                       'average',
                                       # 'mcquitty',
                                       # 'median',
                                       'centroid',
                                       'k-means',
                                       'DBSCAN',
                                       'HDBSCAN'),
                        #selected = 'complete'),
                        selected = 'ward.D2'),
            # selectInput("k",
            #             label = "Define number of clusters (k):",
            #             choices = list(1,2,3,4,5,6,7,8,9,10),
            #             selected = 2),
            # sliderInput("k",
            #             label = "Define number of clusters (k):",
            #             min=1, max=10, value = 2),
            # numericInput("k",
            #             label = "Define number of clusters (k):",
            #             min=1, max=10, value = 2),
            selectInput("k",
                        label = "Select number of clusters (k):",
                        # choices = list("Automatic", "Automatic (extended)", 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
                        choices = list("Automatic", 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25),
                        #selected = "Automatic"),
                        selected = 2),
            # h4("Ancestral plots"),
            # numericInput("group_selection",
            #             label = "Select group to show plot:",
            #             min = 1, max = 10,
            #             value = 1)
            h5("(H)DBSCAN parameters"),
            numericInput("minPts",
                         "Set minPts value",
                         min=0, value=5),
            numericInput("eps",
                         "Set eps value (DBSCAN only)",
                         min=0, value=0.5),
            h3("Ancestral Haplotype"),
            selectInput("ancestral",
                        label = "Ancestral reconstruction method:",
                        choices = list('Branch and bound - indep. sides' = "branchBoundIndep",
                                       'Branch and bound' = "branchBound",
                                       'BranchBoundIndepSides-alleleFreqs' = "branchBoundIndep_alleleFreq",
                                       'BranchBound-alleleFreqs' = "branchBound_alleleFreq",
                                       'Most frequent base' = 'mostFreqBase',
                                       'Simulated founder' = 'simulatedFounder'),
                        selected = "branchBoundIndep"),
            selectInput("cutoff",
                        label = "Select frequency cutoff for ancestral haplotype", #(default: 0.5)
                        choices = seq(0.5, 0.8, 0.1),
                        selected = 0.5),
            numericInput("min_anc_samples",
                         "Set min required samples for branch and bound",
                         min=2, value=0.25),
            h3("Miscellaneous"),
            selectInput("genetic_measure",
                        label = "Select genetic distance measure",
                        choices = list("Physical distance (bp) - relative to mutation"="physical_dist",
                                       #"Physical distance (bp)"="absolute_physical_dist",
                                       "Genetic distance (cM) - relative to mutation"="genetic_dist"
                                       #"Genetic distance (cM)"="absolute_genetic_dist"
                                       ),
                        selected = "physical_dist"),
                        #selected = "centimorgan"),
            numericInput("req_fam_members",
                         label = "Minimum required family members", #(higher value produces better haplotypes but decrease data points),
                         min = 2, max = 10,
                         value = 2),
            numericInput("req_group_members",
                        label = "Minimum required members in group for plotting", #(higher value saves more computation time),
                        min = 1, max = 100,
                        value = 3),
            selectInput("zoom",
                        label = "Zoom in on area around mutation (Mb or cM)",
                        choices = c("No zoom"=0, 1:10),
                        selected = "No zoom")
        , width = 3 ),
        mainPanel(
            tabsetPanel(id = "tabsetPanel", 
                tabPanel("Load Data",
                    h2("Load Data"),
                    HTML("Below are three options to load data into the shiny webapp. 
                         You can load multiple datasets and projects into the same session."),
                    hr(),
                    h3("Option 1: Load 'Build-in' Dataset"),
                    #splitLayout(#cellWidths = c("25%", "75%"),
                    fluidRow(
                        column(3, selectInput("dataset",
                                label = "Select dataset:",
                                choices = list("BRCA1",
                                               "BRCA2",
                                               "BRCA1_simulated",
                                               "BRCA2_simulated"),
                                selected = "BRCA1")), 
                        column(2, actionButton("buttonLoad", "Load data", style = "margin-top: 25px")),
                        column(6, hidden(div(id="text_div", uiOutput("loadingDone", style = "margin-top: 30px"))))
                        #withSpinner(uiOutput("loadingDone"))
                    ),
                    hr(),
                    h3("Option 2: Load Project"),
                    fluidRow(
                        column(3, textInput("project_id_exist", label = "Enter name of project to load")), 
                        column(2, actionButton("buttonLoadProject", "Load project", style = "margin-top: 25px")),
                        column(6, hidden(div(id="div_loadProject", uiOutput("loadProjectState", style = "margin-top: 30px"))))
                    ),
                    hr(),
                    h3("Option 3: Create Project"),
                    fluidRow(
                        column(3, textInput("project_id", label = "Enter project name"))
                        #column(3, textInput("project_gene", label = "Enter name of gene to investigate"))
                    ),
                    #br(),
                    fluidRow(
                        column(3, fileInput("pheno_file", "Choose phenotype file (csv, tsv, txt)",
                                  multiple = FALSE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv"))),
                        column(3, fileInput("geno_file", "Choose genotype file (csv, tsv, txt)",
                                  multiple = FALSE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")))
                    ),
                    fluidRow(
                        column(3, fileInput("coords_file", "Choose coordinates file (csv, tsv, txt)", 
                                  multiple = FALSE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv"))),
                        column(3, fileInput("mut_info_file", "Choose mutation info file (csv, tsv, txt)", 
                                            multiple = FALSE,
                                            accept = c("text/csv",
                                                       "text/comma-separated-values,text/plain",
                                                       ".csv")))
                    ),
                    fluidRow(
                        column(2, actionButton("buttonCreateProject", "Create project")),#, style = "margin-top: 25px")),#, class="btn-info")),#style = "color: white; background-color: #4c8af5; border-color: #4c8af5; margin-top: 25px")),
                        column(6, hidden(div(id="div_createProject", uiOutput("createProjectState", style = "margin-top: 10px"))))
                    )
                    #withSpinner(uiOutput("createProject"))
                ),
                #tabPanel("Hierachical Clustering",
                tabPanel("Clustering",
                        #h2("Hierachical clustering"),
                        h2("Clustering"),
                        textOutput("mutInfo"),
                        withSpinner(plotOutput("hc_plot", height = 600), type = 1, color.background = "#FFFFFF"),
                        htmlOutput("clust_description"),
                        h4("Clustering information:"),
                        htmlOutput("clust_selected"),
                        #h5("Groups are read from right to left in the dendrogram above i.e. group 1 is the rightmost branch etc."),
                        verbatimTextOutput("clust_info")),
                tabPanel("Breakpoints",
                        #h2("Ancestral plots"),
                        h2("Visualizing breakpoints"),
                        h5("Plots all homozygous breakpoints found on either side of the mutation"),
                        br(),
                        withSpinner(uiOutput("breakpoints"))
                ),
                tabPanel("Haplotypes",
                         #h2("Ancestral plots"),
                         h2("Visualizing Haplotypes"),
                         h5("Plots all the haplotypes of each sample based on given ancestral haplotype"),
                         HTML(text = "<pre><b>Plot descriptions:</b>\n    Left-top:     Ordered by group\n    Right-top:    Ordered by total haplotype length\n    Left-bottom:  Ordered by haplotype length above mutation\n    Right-bottom: Ordered by haplotype length below mutation</pre>"),
                         br(),
                         withSpinner(uiOutput("haplotypes_plot"))
                ),
                tabPanel("Age Estimation",
                    h2("Estimating the Age of the Mutation"),
                    #tableOutput("nearestBreaksStatistics"),
                    # h4("Method 1"),
                    # withMathJax(),
                    # h5("Assuming \\(\\theta\\) is the recombination fraction and the probability of recombination between the brca gene and median distance to the first breakpoint among all the samples in a given group, then"),
                    # h5("$$(1-\\theta \ )^n = 0.5$$"),
                    # h5("where \\(n\\) is the age (in generations) of the mutation, that we want to find, thus the age of the mutation is given by"),
                    # h5("$$n = \\frac{\\log(0.5)}{\\log(1-\\theta )}$$"),
                    #withSpinner(DT::dataTableOutput("age_estimation")),
                    #withSpinner(tableOutput("age_estimation")),

                    # h4("Method 2"),
                    # h5("Assuming 1 Mb = 1 CM, then"),
                    # h5("$$n = \\frac{100}{\\text{Median of shared haplotype length}}$$"),
                    # HTML("<br><br>"),
                    #withSpinner(uiOutput("age_estimation_old"))
                    tabsetPanel(id = "mutation_age_table",
                                #tabPanel("Method 1", withSpinner(tableOutput("age_estimation1"))),
                                tabPanel("Method 1",
                                         #h4("Method 1"),
                                         withMathJax(),
                                         h5("Assuming \\(\\theta\\) is the recombination fraction and the probability of recombination between the brca gene and median distance to the first breakpoint among all the samples in a given group, then"),
                                         h5("$$(1-\\theta \ )^n = 0.5$$"),
                                         h5("where \\(n\\) is the age (in generations) of the mutation, that we want to find, thus the age of the mutation is given by"),
                                         h5("$$n = \\frac{\\log(0.5)}{\\log(1-\\theta )}$$"),
                                         HTML("<br>"),
                                         h4("Age in generations:"),
                                         withSpinner(DT::dataTableOutput("age_estimation1"))),
                                #tabPanel("Method 2", withSpinner(tableOutput("age_estimation2"))))
                                # tabPanel("Method 2",
                                #          #h4("Method 2"),
                                #          withMathJax(),
                                #          h5("Assuming 1 Mb = 1 CM, then"),
                                #          h5("$$n = \\frac{100}{\\text{Median of shared haplotype length in Mb}}$$"),
                                #          HTML("<br><br><br>"),
                                #          h4("Age in generations:"),
                                #          withSpinner(DT::dataTableOutput("age_estimation2"))),
                                tabPanel("Method 2",
                                         withMathJax(),
                                         h5(""),
                                         h5("$$n = \\frac{100}{\\text{Mean of shared haplotype length in cM}}$$"),
                                         HTML("<br><br><br>"),
                                         h4("Age in generations:"),
                                         withSpinner(DT::dataTableOutput("age_estimation3"))),
                                # tabPanel("Method 4",
                                #          h5("Same as method 1 but computing age seperate for each side of mutation"),
                                #          HTML("<br>"),
                                #          h4("Age in generations (above mut):"),
                                #          withSpinner(DT::dataTableOutput("age_estimation4_1")),
                                #          HTML("<br>"),
                                #          h4("Age in generations (below mut):"),
                                #          withSpinner(DT::dataTableOutput("age_estimation4_2"))),
                                tabPanel("Method 3",
                                         withMathJax(),
                                         h5("Gandolfo and Speed's Mutation Age Estimation (2014)"),
                                         HTML("<br><br><br>"),
                                         h4("Age in generations:"),
                                         h5("Assuming independent genealogy:"),
                                         withSpinner(DT::dataTableOutput("age_estimation6")),
                                         h5("Assuming correlated genealogy:"),
                                         withSpinner(DT::dataTableOutput("age_estimation6_cor"))),
                                tabPanel("Method 4",
                                         withMathJax(),
                                         h5("Age estimation method based on the compound Poisson-gamma distributions and the gamma distribution defined by the Tweedie class of distributions. Thus,"),
                                         h5("$$Y \\sim Tw_p(\\mu,\\phi)$$"),
                                         h5("with mean and variance"),
                                         h5("$$E[Y] = \\mu \\quad \\text{ and } \\quad Var[Y] = \\phi\\mu^p$$"),
                                         h5("where \\(\\phi\\) is the dispersion paramter and \\(p\\) in the range (1,2]. The optimal parameters are chosen using maximum likelihood."),
                                         h5("The age is then:"),
                                         h5("$$\\frac{100}{\\mu}$$"),
                                         h5("Which is similar to the age estimates from method 2 and 3. However, this method might improve the bounds of the confidence interval."),
                                         HTML("<br><br><br>"),
                                         h4("Age in generations:"),
                                         withSpinner(DT::dataTableOutput("age_estimation8"))),
                                tabPanel("Method 5",
                                         withMathJax(),
                                         h5("Maximum-likelihood estimation"),
                                         h5("$$R_x(n) = L_x(n) = (1-\\theta_{x-1})^n \\cdot (1-(1-(\\theta_x - \\theta_{x-1}))^n)$$"),
                                         h5("$$L = \\prod_{s=1}^{n} L_x(n) \\cdot R_x(n)$$"),
                                         HTML("<br><br><br>"),
                                         h4("Age in generations:"),
                                         withSpinner(DT::dataTableOutput("age_estimation5")))
                    )
                ),
                tabPanel("Migration Pattern",
                         h2("Migration Pattern of Founder Mutations"),
                         h5("Visualizing suggested migration pattern"),
                         br(),
                         withSpinner(uiOutput("migration"))
                ),
                tabPanel("Additional Tools",
                        h2("Additional Tools"),
                        tabsetPanel(id = "additional_tools",
                                    tabPanel("Haplotype data",
                                             h2("Haplotype Length Information"),
                                             br(),
                                             withSpinner(DT::dataTableOutput("mean_haplotype_data")),
                                             withSpinner(uiOutput("haplotype_data"))
                                    ),
                                    tabPanel("Sub-Clustering",
                                             #h2("Hierachical clustering"),
                                             h2("Sub-Clustering"),
                                             h4("Zooming in on clusters from Clustering tab"),
                                             wellPanel(
                                                 #h4("Choose options for second dendrogram"),
                                                 selectInput("group_subClust",
                                                             label = "Select group to zoom in on:",
                                                             choices = c(1),
                                                             selected = 1),
                                                 #uiOutput("group_subCluster"),
                                                 selectInput("k_subClust",
                                                             label = "Select number of clusters (k):",
                                                             choices = c(1:20),
                                                             selected = 1)
                                             ),
                                             textOutput("mutInfo_subCluster"),
                                             withSpinner(plotOutput("hc_plot_subCluster", height = 600)),
                                             htmlOutput("clust_description_subCluster"),
                                             h4("Clustering information:"),
                                             htmlOutput("clust_selected_subCluster"),
                                             verbatimTextOutput("clust_info_subCluster")
                                    ),
                                    tabPanel("Cluster Validity Checking",
                                             h2("Estimate number of groups in clustering"),
                                             h3("TESTING! THIS SECTION IS STILL UNDER DEVELOPMENT. USE WITH CARE."),
                                             tags$em("Note: This tab is only updated/valid when k is set to either 'Automatic' or 'Automatic (extended)'"),
                                             hr(),
                                             h3("Summary of validation scores"),
                                             textOutput("numScores"),
                                             tableOutput("summary_cluster_estimation"),
                                             htmlOutput("optimalClusters"),
                                             # fluidRow(
                                             #     column(width = 6, withSpinner(uiOutput("clust_estimation")),
                                             #     column(width = 6, withSpinner(uiOutput("clust_estimation"))
                                             # ),
                                             withSpinner(uiOutput("clust_estimation"))
                                    ),
                                    tabPanel("Compare Dendrograms",
                                             h2("Visual comparison of two dendrograms"),
                                             wellPanel(
                                                 h4("Choose options for second dendrogram"),
                                                 selectInput("distMethod2",
                                                             label = "Select distance metric:",
                                                             choices = list('firstBreakDist',
                                                                            "firstUpperBreakDist",
                                                                            "firstUnderBreakDist",
                                                                            "matchstates",
                                                                            "genpofad",
                                                                            "mrca",
                                                                            "nei",
                                                                            "euclidean weighted according to brca dist"="euclideanDist_BRCAdistWeighted"
                                                             ),
                                                             selected = 'firstBreakDist'),
                                                 selectInput("clustMethod2",
                                                             label = "Clustering method:",
                                                             choices = list('complete',
                                                                            'ward.D',
                                                                            'ward.D2',
                                                                            'single'),
                                                             selected = 'complete')
                                             ),
                                             withSpinner(plotOutput("dendro_comparison_plot", height = 800))
                                    ),
                                    tabPanel("Nearest Breakpoints",
                                             h2("Visualizing nearest breakpoints"),
                                             h5("Plots only the nearest homozygous breakpoint on either side of the BRCA gene. Similar to Haplotypes plot, just showing nearest break as point."),
                                             #HTML("<pre><b>Plot description:</b>\n    Left:   Ordered by group\n    Middle: Ordered by positive genomic position\n    Right:  Ordered by negative genomic position</pre>"),
                                             HTML("<pre><span><b>Plot descriptions:</b>\n    Left-top:     Ordered by group\n    Right-top:    Ordered by haplotype length\n    Left-bottom:  Ordered by positive genomic position\n    Right-bottom: Ordered by negative genomic position</span></pre>"),
                                             br(),
                                             withSpinner(uiOutput("nearest_breakpoint"))
                                    )
                        )
                ),
                tabPanel("How To Use",
                    includeMarkdown("instructions.rmd")
                )
            )

        , width = 9 )
    )
)

server <- function(input, output, session) {
    # Allow file size up to 300 MB
    options(shiny.maxRequestSize=300*1024^2)
    options(shiny.reactlog=T)

    ###########################################################################
    ### Load data
    ###########################################################################
    load_dataset <- observeEvent(input$buttonLoad, {
        
        # gene should probably de independent from dataset name
        gene <<- input$dataset
    
        switch (gene,
                "BRCA1" = loadBRCA1(),
                "BRCA2" = loadBRCA2(),
                "BRCA1_simulated" = loadBRCA1Sim(),
                "BRCA2_simulated" = loadBRCA2Sim()
        )
        update_geneList(gene)
        
        #"Dataset loaded"
        toggle('text_div')
        output$loadingDone <- renderText({paste("<font color=\"#00e500\"><b>Dataset", input$dataset, "loaded!</b></font>")})
        
    })
    
    update_geneList <- function(gene){
        if (!(gene %in% loadedGenes)) loadedGenes <<- c(loadedGenes, gene)
        updateSelectInput(session, "gene", label = "Select gene:", choices = loadedGenes, selected = gene)
    }
    
    # observeEvent(input$buttonLoad, {
    #     toggle('text_div')
    #     output$loadingDone <- renderText({paste("<b>Dataset", input$dataset, "loaded!</b>")})
    # })
    
    # Register updates to project
    observe({
        input$dataset
        input$project_id
        input$project_id_exist
        #hide('text_div')
        output$loadingDone <- renderText({""})
        output$loadProjectState <- renderText({""})
        output$createProjectState <- renderText({""})
    })
    
    # observe({
    #     input$buttonLoad
    #     updateSelectInput(session, "gene", label = "Select gene:", choices = loadedGenes, selected = input$gene)
    # })
    
    observeEvent(input$buttonLoadProject, {
        status <- loadProject(input$project_id_exist)
        toggle("div_loadProject")
        if (status[1] == 0){
            update_geneList(gene)
            output$loadProjectState <- renderText(paste0("<font color=\"#00e500\"><b>", status[2], "</b></font>"))
        } else {
            output$loadProjectState <- renderText(paste0("<font color=\"#FF0000\"><b>", status[2], "</b></font>"))
        }
    })
    
    observeEvent(input$buttonCreateProject, {
        status = createProject(input$project_id, #input$project_gene, 
                               input$pheno_file, input$geno_file,
                               input$coords_file, input$mut_info_file)
        toggle("div_createProject")
        if (status[1] == 0){
            update_geneList(gene)
            output$createProjectState <- renderText(paste0("<font color=\"#00e500\"><b>", status[2], "</b></font>"))
        } else {
            output$createProjectState <- renderText(paste0("<font color=\"#FF0000\"><b>", status[2], "</b></font>"))
        }
    })
    
    ###########################################################################
    ### Prepare data
    ###########################################################################
        
    #output$selectMut <- renderUI({
    observe({
    #selectMut <- reactive({
        if (is.null(input$gene)){
            print("kagemand")
            choices = "Not yet available"
            selected = "Not yet available"
        } else if (input$gene == "BRCA2_toySample"){
            choices = "c.7617+1G>A"
            selected = "c.7617+1G>A"
        } else if (open_projects[[gene]]$project == T){
            choices = open_projects[[gene]]$project_mut_info$Mut1HGVS
            selected = choices[1]
        } else if (length(input$gene) == 0 || input$gene == "BRCA1"){
            choices = count(brca1_pheno_merged, Mut1HGVS) %>% filter(n>=80) %>% arrange(n) %>% .$Mut1HGVS
            selected = "c.68_69delAG"
        } else if (input$gene == "BRCA2"){
            choices = count(brca2_pheno_merged, Mut1HGVS) %>% filter(n>=80) %>% arrange(n) %>% .$Mut1HGVS
            selected = "c.7617+1G>A"
        } else if (input$gene == "BRCA1_simulated"){
            #choices = unique(brca1_simulated_pheno_merged$Mut1HGVS)
            #choices = paste0("mutSim_BRCA1_",rep(seq(10,500,10), each=10),"_", 1:10)
            choices = paste0("mutSim_BRCA1_",rep(seq(10,2010,50), each=10),"_", 1:10)
            selected = choices[1]
        } else if (input$gene == "BRCA2_simulated"){
                choices = "Not yet available"#unique(brca1_simulated_pheno_merged$Mut1HGVS)
                selected = "Not yet available"#choices[1]
        }
        #selectInput("mut",
        updateSelectInput(session, "mut",
                    label = "Select mutation:",
                    choices = choices,
                    selected = selected)
    })
    
    get_mut <- function(){
        if (input$gene == "BRCA1") {
            input$mut1
        } else if (input$gene == "BRCA2") {
            input$mut2
        } else if (input$gene == "BRCA1_simulated") {
            input$mut3
        } else if (input$gene == "BRCA2_simulated") {
            input$mut4
        }
    }

    loadBRCA1Sim_gen <- reactive({
        gen = strsplit(input$mut, "_")[[1]][3]
        filename = paste0("classify-simulated-population-age/BRCA1_simulated-famHaplotypes-starGenealogy-100_samples-10_simulations-generations_10_2010_step50-seed_42/simulated_population/BRCA1_simulated-starGenealogy-100_samples-10_simulations-generations_",gen,"-seed_42")
        #filename = paste0("classify-simulated-population-age/BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_10_500_step10-seed_42/simulated_population/BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_",gen,"-seed_42")
        #filename = paste0("classify-simulated-population-age/BRCA1_simulated-starGenealogy-100_samples-10_simulations-generations_10_10010_step50-seed_42/simulated_population/BRCA1_simulated-starGenealogy-100_samples-10_simulations-generations_",gen,"-seed_42")
        filename_sim <<- filename
        brca1_simulated_pheno_merged <<- read.table2(file = paste0(filename, "-pheno.txt"), header = T)
        brca1_simulated_geno_plink <<- read.table2(file = paste0(filename, "-geno.txt"), header = T)
    })
    
    updatedField <- function(){
        # Basic fields are updated
        gene <<- input$gene
        mut <<- input$mut
        fam <<- if (input$fam) "fam" else "individual"
        ancestral_method <<- input$ancestral
        
        validate(
            need(!(fam == "fam" && grepl("BRCA1_simulated", gene)),
                 "Family haplotype is not possible for simulated data. Switch to consider individually."),
            need(!(grepl("BRCA2_simulated", gene)),
                 "Not yet available")
        )
        
        # Required fam members updated
        req_fam_members <<- input$req_fam_members
        fam_members <<- if (input$fam) paste0("-famMembers_", req_fam_members) else ""
        
        # Load current project based on gene selected
        project_path <<- open_projects[[gene]]$project_path
        if (open_projects[[gene]]$project==F && grepl("BRCA1_simulated", gene)) loadBRCA1Sim_gen()
    }

    prepare <- reactive({
        updatedField()
        if (open_projects[[gene]]$project){
            prepare_dataframe(mut, gene, gene_info_obj = open_projects[[gene]]$project_mut_info)
        } else {
            prepare_dataframe(mut, gene)
        }
        extractFamPlink(req_fam_members)
    })

    ###########################################################################
    ### Specify number of clusters - either manually or automatically
    ###########################################################################
    getK <- reactive({
        if (input$grouping == "Country"){
            countries = get_countries()
            return(length(countries))
        }
        
        if (input$clustMethod == "DBSCAN" || input$clustMethod == "HDBSCAN"){
            hc()
            return(max(brca_data$cluster_groups))
        }

        if (input$k == "Automatic"){
            dst = dst()
            #fam = if (input$fam) "fam" else "individual"
            index <- getScores()
            k = numCluster_estimation(max.clust = 12,
                                      distMethod = input$distMethod,
                                      clustMethod = input$clustMethod,
                                      fam = fam, dist = dst,
                                      req_members = input$req_group_members,
                                      index = index)
        } else if (input$k == "Automatic (extended)"){
            k = 2
        } else {
            k = as.integer(input$k)
            numCluster <<- t(c(0,0))
        }
        #print(k)
        return(k)
    })

    ###########################################################################
    ### Clustering tab
    ###########################################################################

    output$mutInfo <- renderText({
        prepare()
        paste("Number of samples:", nrow(matched_plink_single), "- Number of families (with 2 members or more):", nrow(matched_plink_fam))
    })

    getDist <- function(distMethod){
        fam = if (input$fam) "fam" else "individual"
        filename = paste0(project_path, "distMatrices/", gene, "-", mut_name, "-distMethod-", distMethod, "-", fam, fam_members, ".rds")

        if (distMethod == "firstBreakDist"){
            dst <- firstBreakDist(fam = input$fam, filename = filename)
        } else if (distMethod %in% c("firstBreakDist_cM", "firstUpperBreakDist", "firstUnderBreakDist", "firstUpperBreakDist_cM", "firstUnderBreakDist_cM", 
                                     "firstBreakEuclideanDist","firstBreakEuclideanDist_cM","matchstates", "genpofad", "mrca", "nei")){
            dst <- computeCustomDist(get(distMethod), input$fam, filename)
        } else if (distMethod == "euclideanDist_BRCAdistWeighted"){
            if (file.exists(filename)){
                dst <- readRDS(filename)
            } else {
                dst <- euclideanDist_BRCAdistWeighted(input$fam)
                saveRDS(dst, file = filename)
            }
        } else {
            if (file.exists(filename)){
                dst <- readRDS(filename)
            } else {
                dst <- amapDist(distMethod, input$fam)
                saveRDS(dst, file = filename)
            }
        }

        print(dim(dst))
        return(dst)
    }

    dst <- reactive({
        # Make sure it trigger, if updated
        updatedField()
        distMethod = input$distMethod
        getDist(distMethod)
    })

    dst2 <- reactive({
        distMethod = input$distMethod2
        getDist(distMethod)
    })

    hc <- reactive({
        #print(2)
        dst <- dst()
        if (input$clustMethod == "DBSCAN"){
            db <<- fpc::dbscan(data = dst, eps = input$eps, MinPts = input$minPts, method = "dist")
            cluster_groups <<- db$cluster
            brca_data.filtered$cluster_groups <<- db$cluster
            brca_data <<- brca_data.filtered
        } else if (input$clustMethod == "HDBSCAN") {
            hdb <<- dbscan::hdbscan(x = dst, minPts = input$minPts)
            cluster_groups <<- hdb$cluster
            brca_data.filtered$cluster_groups <<- hdb$cluster
            brca_data <<- brca_data.filtered
        } else if (input$clustMethod == "k-means"){
            km <<- kmeans(dst, getK())
            cluster_groups <<- km$cluster
            brca_data.filtered$cluster_groups <<- km$cluster
            brca_data <<- brca_data.filtered
        } else {
            hiearchical_clustering(dst, clustMethod = input$clustMethod, k = getK(), doHapPlot = F)
        }

    })

    output$hc_plot <- renderPlot({
        print(1)
        updatedField()
        #prepare()

        # print(input$mut1)
        # print(input$mut2)
        # print(input$gene)
        # print(input$fam)
        # print(input$k)

        k = getK()
        hc <- hc()
        print(3)

        if (input$clustMethod == "DBSCAN"){
            dbscan_plot(db)
        } else if (input$clustMethod == "HDBSCAN"){
            hdbscan_plot(hdb_object = hdb, data = dst())
        } else if (input$clustMethod == "k-means"){
            plot(c(1,2,3),c(1,2,3), main = "Dummy plot")
        } else {
            plot_dendrogram(hc, k)
        }

    }, height = 600)

    output$clust_description <- renderText({
        if (input$clustMethod %in% c("DBSCAN", "HDBSCAN", "k-means")){
            return("")
        } else {
            return("Groups are read from right to left in the dendrogram above i.e. group 1 is the rightmost branch etc.")
        }
    })

    output$clust_selected <- renderText({
        paste("Number of clusters: <b>", getK(), "</b>")
    })

    output$clust_info <- renderPrint({
        # Reactive dependencies
        updatedField()
        distMethod = input$distMethod
        clustMethod <- input$clustMethod
        k = getK()
        hc <- hc()
        # Print group information
        groups <- list()
        for (i in 1:max(cluster_groups)){
            #print(subset(brca_data.filtered, cluster_groups == i) %>% select(FamCode, Onc_ID, Country, cluster_groups))
            cat(paste0("Group ", i, ":\n"))
            #print(subset(brca_data.filtered, cluster_groups == i) %>% count(Country))
            groups[[i]] <- as.data.frame(subset(brca_data, cluster_groups == i) %>% count(Country))# %>% select(Country, number=n))
            print(groups[[i]], row.names = F)
            cat("\n")
            cat(paste("Total:", sum(groups[[i]]$n)))
            if (i < max(cluster_groups)) cat("\n\n\n")
        }

        #print(groups)
    })

    ###########################################################################
    ### Sub-Clustering tab
    ###########################################################################

    use_subCluster <- function(){
        updatedField()
        k = getK()
        hc <- hc()
        brca_data <<- subset(brca_data, cluster_groups == input$group_subClust)
        brca_data.filtered <<- brca_data
        matched_plink <<- extractSamples(brca_data, matched_plink, chr_coords)
        dst <<- firstBreakDist(matched_custom = matched_plink)
        hc <<- hiearchical_clustering(dst, clustMethod = input$clustMethod, k = input$k_subClust)
    }

    get_sub_dendrogram <- reactive({
        updatedField()
        k = getK()
        hc <- hc()
        dend <- as.dendrogram(hc)
        dend_list <- get_subdendrograms(dend = dend, k = k)

        # Remove warning due to slow update of input$group_subClust
        if (is.null(input$group_subClust)){
            sub_dend <<- dend_list[[1]]
        } else {
            sub_dend <<- dend_list[[as.integer(input$group_subClust)]]
        }

        brca_data_subDend <<- brca_data[sort(order.dendrogram(sub_dend)), ]
        #brca_data_subDend <<- brca_data[order.dendrogram(sub_dend), ]

        return(brca_data_subDend)
    })

    #output$group_subCluster <- renderUI({
    observe({
        #prepare()
        k = getK()
        #hc <- hc()
        #selectInput("group_subClust",
        updateSelectInput(session, "group_subClust", 
                    label = "Select group to zoom in on:",
                    choices = c(1:k),
                    selected = input$group_subClust)
    })

    output$mutInfo_subCluster <- renderText({
        prepare()
        brca_data_subDend <- get_sub_dendrogram()
        matched_plink_single2 = filter(matched_plink_single, SNP %in% brca_data_subDend$Onc_ID)
        matched_plink_fam2 = filter(matched_plink_fam, SNP %in% brca_data_subDend$Onc_ID)
        paste("Number of samples:", nrow(matched_plink_single2), "- Number of families (with more than 2 members):", nrow(matched_plink_fam2))
    })

    output$hc_plot_subCluster <- renderPlot({
        prepare()
        #updatedField()
        k = getK()
        hc <- hc()

        # Remove warning due to slow update of input$group_subClust
        if (is.null(input$group_subClust)){
            return()
        }

        plot_dendrogram(hc = hc, k = as.integer(input$k_subClust), leaf = as.integer(input$group_subClust), sub_k = k)

    })

    output$clust_description_subCluster <- renderText({
        if (input$clustMethod %in% c("DBSCAN", "HDBSCAN", "k-means")){
            return("")
        } else {
            return("Groups are read from right to left in the dendrogram above i.e. group 1 is the rightmost branch etc.")
        }
    })

    output$clust_selected_subCluster <- renderText({
        paste("Number of clusters: <b>", input$k_subClust, "</b>")
    })

    output$clust_info_subCluster <- renderPrint({
        # Reactive dependencies
        prepare()
        #updatedField()
        #distMethod = input$distMethod
        #clustMethod <- input$clustMethod
        k = getK()
        hc <- hc()

        brca_data_subDend <- get_sub_dendrogram()

        k_subClust = as.integer(input$k_subClust)
        dend_sublist <- get_subdendrograms(reindex_dend(sub_dend), k = k_subClust)
        #cluster_groups <<- cutree(sub_dend, k=3)

        # Print group information
        groups <- list()
        for (i in 1:k_subClust){
            cat(paste0("Group ", i, ":\n"))
            #groups[[i]] <- as.data.frame(subset(brca_data_subDend, cluster_groups == i) %>% count(Country))
            groups[[i]] <- as.data.frame(brca_data_subDend[order.dendrogram(dend_sublist[[i]]), ] %>% count(Country))
            print(groups[[i]], row.names = F)
            cat("\n")
            cat(paste("Total:", sum(groups[[i]]$n)))
            if (i < max(cluster_groups)) cat("\n\n\n")
        }
        #print(groups)
    })

    ###########################################################################
    ### Breakpoints tab
    ###########################################################################

    get_filename <- function(plot_type = "breakpointPlot", group = NULL, country = NULL, possible_zoom=F, 
                             mutAge_method = NULL, file_ext = "png"){
        k = getK()
        
        db_val = ""
        if (input$clustMethod == "DBSCAN"){
            db_val = paste0("-minPts_", input$minPts, "-eps_", input$eps)
        } else if (input$clustMethod == "HDBSCAN") {
            db_val = paste0("-minPts_", input$minPts)
        }
        
        minSamples=""
        if(ancestral_method=="branchBound" || ancestral_method == "branchBoundIndep"){
            minSamples=paste0("-minSamples_",input$min_anc_samples)
        }
        
        # Mutation age method filename
        if (!is.null(mutAge_method)){
            filename = paste0(project_path, "/mutationAge_groups/", mut_name, "/", gene, "-", mut_name, "-MutationAge-", 
                               "ancestMethod_", ancestral_method, minSamples, "-distMethod-", input$distMethod, "-clustMethod-",
                               input$clustMethod, "-", fam, fam_members, "-k_", k, "-cutoff_", input$cutoff, db_val,
                               "-method_", mutAge_method, ".txt")
            return(filename)
        }
        
        # Plot filename
        zoom = ""
        if (possible_zoom){
            if (as.integer(input$zoom) > 0 && input$genetic_measure == "physical_dist"){
                zoom = paste0("-zoomLevel_", as.integer(input$zoom), "Mb")
            } else if (as.integer(input$zoom) > 0 && input$genetic_measure == "genetic_dist"){
                zoom = paste0("-zoomLevel_", as.integer(input$zoom), "cM")
            }
        }
        
        if (input$grouping == "Clustering"){
            
            dir.create(paste0(project_path, "/", plot_type, "s/", mut_name, "/"), showWarnings = F)
            filename = paste0(project_path, "/", plot_type, "s/", mut_name, "/", gene, "-", mut_name, "-", plot_type,
                              "-ancestMethod_", ancestral_method, minSamples, "-distMethod_", input$distMethod, "-clustMethod_", 
                              input$clustMethod, "-", fam, fam_members, "-k_", k, "-", input$genetic_measure,
                              "-cutoff_", input$cutoff, db_val, zoom,"-group_", group, ".", file_ext)
        } else {
            
            dir.create(paste0(project_path, "/",plot_type, "sCountry/", mut_name, "/"), showWarnings = F)
            filename = paste0(project_path, "/",plot_type, "sCountry/", mut_name, "/", gene, "-", mut_name, "-", plot_type,
                              "-ancestMethod_", ancestral_method, minSamples, "-Country_", country, "-cutoff_",
                              input$cutoff, "-", fam, fam_members, zoom, ".", file_ext)
        }
        return(filename)
    }
    
    get_countries <- reactive({
        updatedField()
        quantile = 0.05
        # Get countries in current dataset
        if (input$fam) {
            brca_data_temp <- brca_data %>% filter(Onc_ID %in% matched_plink_fam$SNP)
            countries <- brca_data_temp %>% group_by(Country) %>% summarise(n=n()/nrow(.)) %>% filter(n > quantile) %>% .$Country
        } else {
            countries <- brca_data %>% group_by(Country) %>% summarise(n=n()/nrow(.)) %>% filter(n > quantile) %>% .$Country
        }
        return(countries)
    })
    
    # Returns chr_pos dataframes for all groups in a list 
    haplotype_data <- reactive({
        prepare()
        k = getK()
        countries = get_countries()
        haplotype_data = list()
        for (i in 1:k){
            if (input$grouping == "Clustering") {
                hc()
                haplotype_data[[i]] = haplotype_mutation(mut, prepare = F, group = i, cutoff = input$cutoff, ancestral_method = ancestral_method, min_samples = input$min_anc_samples)
            } else {
                country = countries[i]
                haplotype_data[[i]] = haplotype_mutation(mut, prepare = F, country = country, cutoff = input$cutoff, ancestral_method = ancestral_method, min_samples = input$min_anc_samples)
            }
        }
        return(haplotype_data)
    })

    plot_breakpoints <- reactive({
        updatedField()
        countries <- get_countries()
        # Test
        if (useSubCluster){
            use_subCluster()
            k = as.integer(input$k_subClust)
        } else {
            k = getK()
        }
        
        withProgress(message = "Computing breakpoints and creating plots", value = 0, {
            for (i in 1:k){
                if (input$grouping == "Clustering" && length(cluster_groups[cluster_groups==i]) < input$req_group_members) {
                    incProgress(amount = 1/k, detail = paste(i, "of", k))
                    next
                }
                local({
                    i <- i
                    
                    zoom_level = 0
                    if (as.integer(input$zoom) > 0 && input$genetic_measure == "physical_dist"){
                        zoom_level = as.integer(input$zoom) * 10^6
                    } else if (as.integer(input$zoom) > 0 && input$genetic_measure == "genetic_dist"){
                        zoom_level = as.integer(input$zoom)
                    }
                    
                    if (input$grouping == "Clustering"){
                        # k2 = k
                        # db_val = ""
                        # if (input$clustMethod == "DBSCAN"){
                        #     db_val = paste0("-minPts_", input$minPts, "-eps_", input$eps)
                        #     k2 = 2 # Fix k to not make many identical plots
                        # } else if (input$clustMethod == "HDBSCAN") {
                        #     db_val = paste0("-minPts_", input$minPts)
                        #     k2 = 2 # Fix k to not make many identical plots
                        # }
                        # 
                        # zoom = ""
                        # if (as.integer(input$zoom) > 0 && input$genetic_measure == "physical_dist"){
                        #     zoom = paste0("-zoomLevel_", as.integer(input$zoom), "Mb")
                        # } else if (as.integer(input$zoom) > 0 && input$genetic_measure == "genetic_dist"){
                        #     zoom = paste0("-zoomLevel_", as.integer(input$zoom), "cM")
                        # }
                        # 
                        # dir.create(paste0(project_path, "/breakpointPlots/", mut_name, "/"), showWarnings = F)
                        # filename = paste0(project_path, "/breakpointPlots/", mut_name, "/", gene, "-", mut_name, "-BreakpointPlot-",
                        #                   "-ancestMethod_", ancestral_method, "-distMethod_", input$distMethod, "-clustMethod_", 
                        #                   input$clustMethod, "-", fam, fam_members, "-k_", k2, "-", input$genetic_measure,
                        #                   "-cutoff_", input$cutoff, db_val, zoom,"-group_", i, ".png")
                        
                        filename = get_filename(plot_type = "breakpointPlot", group = i, possible_zoom = T)
                        plotname = paste0("BreakpointsPlot-Group_",i)
                    } else { #Country
                        country = countries[i]
                        print(country)
                        
                        # dir.create(paste0(project_path, "breakpointPlotsCountry/", mut_name, "/"), showWarnings = F)
                        # filename = paste0(project_path, "/breakpointPlotsCountry/", mut_name, "/", gene, "-", mut_name,
                        #                    "-BreakpointPlot-ancestMethod_", ancestral_method, "-Country_", country, "-cutoff_",
                        #                    input$cutoff, "-", fam, fam_members, ".png")
                        filename = get_filename(country = country, plot_type = "breakpointPlot", possible_zoom = T)
                        plotname = paste0("BreakpointsPlot-", country)
                    }

                    if (!file.exists(filename)){
                        incProgress(amount = 0, detail = paste(i, "of", k))
                        haplotype_data = haplotype_data()
                        
                        #chr_pos = haplotype_mutation(mut, prepare = F, group = i, cutoff = input$cutoff, ancestral_method = ancestral_method)
                        chr_pos = haplotype_data[[i]]
                        
                        incProgress(amount = 1/(k*2), detail = paste(i, "of", k))
                        if (input$grouping == "Clustering"){
                            p = plot_breakpoints_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut,
                                                       group = i, col = input$genetic_measure, zoom = zoom_level)
                            p = add_second_legend(p, chr_pos)
                        } else {
                            p = plot_breakpoints_country(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, 
                                                         otherCountry = T, zoom = zoom_level)
                        }
                        
                        ggsave(filename = filename, plot = p, scale = 0.8, width = 17.5, height = 10, dpi = 300)

                        incProgress(amount = 1/(k*2), detail = paste(i, "of", k))
                        
                        # General ratio calculation for height:
                        #   assuming plots should have form factor width_p x height_p (e.g. 17,5 x height_p)
                        #   and assuming width_w (e.g. 1225) and height_w (e.g. 800) pixels of window available
                        #   Then the height_p is computed as: height_p = width_p / (width_w/height_w)
                        
                        # Ratio conversion: 17,5/(1225/400) = 5,714 - for 3 plots on one row
                        #ggsave(filename3, scale = 1, width = 17.5, height = 5.714)
                        
                        # Ratio conversion: 17,5/(1225/800) = 11,429 - for 2x2 plots on grid
                        # ggsave(filename = filename3, plot = p2, scale = 1, width = 17.5, height = 11.429)
                        
                        # incProgress(amount = 1/(k*2), detail = paste(i, "of", k))
                    }
                    output[[plotname]] <- renderImage({
                        list(src = filename, width = "100%", height = 700)
                    }, deleteFile = F)
                })
            }
        })
    })
    
    plot_nearestBreakpoints <- reactive({
        updatedField()
        countries <- get_countries()
        k = getK()
        
        withProgress(message = "Computing breakpoints and creating plots", value = 0, {
            for (i in 1:k){
                if (input$grouping == "Clustering" && length(cluster_groups[cluster_groups==i]) < input$req_group_members) {
                    incProgress(amount = 1/k, detail = paste(i, "of", k))
                    next
                }
                local({
                    i <- i
                    if (input$grouping == "Clustering"){
                        # k2 = k
                        # db_val = ""
                        # if (input$clustMethod == "DBSCAN"){
                        #     db_val = paste0("-minPts_", input$minPts, "-eps_", input$eps)
                        #     k2 = 2 # Fix k to not make many identical plots
                        # } else if (input$clustMethod == "HDBSCAN") {
                        #     db_val = paste0("-minPts_", input$minPts)
                        #     k2 = 2 # Fix k to not make many identical plots
                        # }
                        # 
                        # dir.create(paste0(project_path, "/nearestBreakpointPlots/", mut_name, "/"), showWarnings = F)
                        # filename = paste0(project_path, "/nearestBreakpointPlots/", mut_name, "/", gene, "-", mut_name, "-NearestBreakpointPlot-",
                        #                   "ancestMethod_", ancestral_method, "-distMethod_", input$distMethod, "-clustMethod_", input$clustMethod, "-",
                        #                   fam, fam_members, "-k_", k2, "-", input$genetic_measure, "-cutoff_", input$cutoff, db_val, "-group_", i, ".png")
                        
                        filename = get_filename(plot_type = "nearestBreakpointPlot", group = i)
                        plotname = paste0("NearestBreakpointsPlot-Group_",i)
                    } else {
                        country = countries[i]
                        print(country)
                        
                        # dir.create(paste0(project_path, "/nearestBreakpointPlotsCountry/", mut_name, "/"), showWarnings = F)
                        # filename = paste0(project_path, "/nearestBreakpointPlotsCountry/", mut_name, "/", gene, "-", mut_name,
                        #                    "-NearestBreakpointPlot-ancestMethod_", ancestral_method, "-Country_", country, "-cutoff_",
                        #                    input$cutoff, "-", fam, fam_members, ".png")
                        
                        filename = get_filename(plot_type = "nearestBreakpointPlot", country = country)
                        plotname = paste0("NearestBreakpointsPlot-",country)
                    }
                    
                    if (!file.exists(filename)){
                        incProgress(amount = 0, detail = paste(i, "of", k))
                        
                        haplotype_data = haplotype_data()
                        
                        #chr_pos <<- haplotype_mutation(mut, prepare = F, group = i, cutoff = input$cutoff, ancestral_method = ancestral_method)
                        chr_pos <- haplotype_data[[i]]
                        incProgress(amount = 1/(k*2), detail = paste(i, "of", k))
                        
                        if (input$grouping == "Clustering"){
                            nearestBreakPlots = plot_nearest_brca_break_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group = i, col = input$genetic_measure)
                            p2 = ggarrange(nearestBreakPlots[[1]], nearestBreakPlots[[4]], nearestBreakPlots[[2]], nearestBreakPlots[[3]], ncol=2, nrow=2)
                            p2 = add_second_legend(p2, chr_pos)
                        } else {
                            nearestBreakPlots = plot_nearest_brca_break_country(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, fixObject = T) #, country = country)
                            p2 = ggarrange(nearestBreakPlots[[1]], nearestBreakPlots[[4]], nearestBreakPlots[[2]], nearestBreakPlots[[3]], ncol=2, nrow=2)
                        }
                        
                        # General ratio calculation for height:
                        #   assuming plots should have form factor width_p x height_p (e.g. 17,5 x height_p)
                        #   and assuming width_w (e.g. 1225) and height_w (e.g. 800) pixels of window available
                        #   Then the height_p is computed as: height_p = width_p / (width_w/height_w)
                        
                        # Ratio conversion: 17,5/(1225/400) = 5,714 - for 3 plots on one row
                        #ggsave(filename3, scale = 1, width = 17.5, height = 5.714)
                        
                        # Ratio conversion: 17,5/(1225/800) = 11,429 - for 2x2 plots on grid
                        ggsave(filename = filename, plot = p2, scale = 1, width = 17.5, height = 11.429)
                        
                        incProgress(amount = 1/(k*2), detail = paste(i, "of", k))
                    }
                    output[[plotname]] <- renderImage({
                        list(src = filename, width = "100%", height = 800)
                    }, deleteFile = F)
                })
            }
        })
        
    })
    
    # Deprecated
    plot_country_breakpoints <- reactive({
        updatedField()
        #prepare() # Needed here as well!
        #k = getK()
        #colors_country()
        #fam = if (input$fam) "fam" else "individual"
        quantile = 0.05
        countries <<- brca_data %>% group_by(Country) %>% summarise(n=n()/nrow(.)) %>% filter(n > quantile) %>% .$Country
        print(countries)

        withProgress(message = "Computing breakpoints and creating plots", value = 0, {
            k = length(countries)
            for (i in 1:length(countries)){
                local({
                    country = countries[i]
                    print(country)
                    
                    dir.create(paste0(project_path, "breakpointPlotsCountry/", mut_name, "/"), showWarnings = F)
                    filename2 = paste0(project_path, "/breakpointPlotsCountry/", mut_name, "/", gene, "-", mut_name,
                                       "-BreakpointPlot-ancestMethod_", ancestral_method, "-Country_", country, "-cutoff_",
                                       input$cutoff, "-", fam, fam_members, ".png")
                    dir.create(paste0(project_path, "/nearestBreakpointPlotsCountry/", mut_name, "/"), showWarnings = F)
                    filename3 = paste0(project_path, "/nearestBreakpointPlotsCountry/", mut_name, "/", gene, "-", mut_name,
                                       "-NearestBreakpointPlot-ancestMethod_", ancestral_method, "-Country_", country, "-cutoff_",
                                       input$cutoff, "-", fam, fam_members, ".png")
                    plotname = paste0("BreakpointsPlot-", country)
                    plotname2 = paste0("NearestBreakpointsPlot-",country)

                    if (!file.exists(filename2)){
                        incProgress(amount = 0, detail = paste(i, "of", k))
                        haplotype_data = haplotype_data()
                        
                        #chr_pos = haplotype_mutation(mut, prepare = F, country = country, cutoff = input$cutoff, ancestral_method = ancestral_method)
                        chr_pos = haplotype_data[[i]]
                        incProgress(amount = 1/(k*3), detail = paste(i, "of", k))

                        p = plot_breakpoints_country(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, otherCountry = T)#, country = country)
                        ggsave(filename2, scale = 1, width = 17.5, height = 10)
                        incProgress(amount = 1/(k*3), detail = paste(i, "of", k))

                        nearestBreakPlots = plot_nearest_brca_break_country(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, fixObject = T) #, country = country)
                        p2 = ggarrange(nearestBreakPlots[[1]], nearestBreakPlots[[4]], nearestBreakPlots[[2]], nearestBreakPlots[[3]], ncol=2, nrow=2)

                        # General ratio calculation for height:
                        #   assuming plots should have form factor width_p x height_p (e.g. 17,5 x height_p)
                        #   and assuming width_w (e.g. 1225) and height_w (e.g. 800) pixels of window available
                        #   Then the height_p is computed as: height_p = width_p / (width_w/height_w)

                        # Ratio conversion: 17,5/(1225/400) = 5,714 - for 3 plots on one row
                        #ggsave(filename3, scale = 1, width = 17.5, height = 5.714)

                        # Ratio conversion: 17,5/(1225/800) = 11,429 - for 2x2 plots on grid
                        ggsave(filename3, scale = 1, width = 17.5, height = 11.429)

                        incProgress(amount = 1/(k*3), detail = paste(i, "of", k))
                    }

                    output[[plotname]] <- renderImage({
                        #list(src = filename2, width = 1225, height = 700)
                        list(src = filename2, width = "100%", height = 700)
                    }, deleteFile = F)
                    output[[plotname2]] <- renderImage({
                        list(src = filename3, width = "100%", height = 800)
                    }, deleteFile = F)
                })
            }
        })
        return(countries)
    })

    # Deprecated
    plot_cluster_breakpoints <- reactive({
        # Make sure we update
        updatedField()

        # Test
        if (useSubCluster){
            use_subCluster()
            k = as.integer(input$k_subClust)
        } else {
            k = getK()
        }
        #k = getK()
        #brca_data <<- brca_data.filtered
        #req_group_members <<- 3

        withProgress(message = "Computing breakpoints and creating plots", value = 0, {
            for (i in 1:k){
                if (length(cluster_groups[cluster_groups==i]) < input$req_group_members) {
                    incProgress(amount = 1/k, detail = paste(i, "of", k))
                    next
                }
                local({
                    i <- i
                    k2 = k
                    db_val = ""
                    if (input$clustMethod == "DBSCAN"){
                        db_val = paste0("-minPts_", input$minPts, "-eps_", input$eps)
                        k2 = 2 # Fix k to not make many identical plots
                    } else if (input$clustMethod == "HDBSCAN") {
                        db_val = paste0("-minPts_", input$minPts)
                        k2 = 2 # Fix k to not make many identical plots
                    }

                    zoom = ""
                    zoom_level = 0
                    if (as.integer(input$zoom) > 0 && input$genetic_measure == "physical_dist"){
                        zoom = paste0("-zoomLevel_", as.integer(input$zoom), "Mb")
                        zoom_level = as.integer(input$zoom) * 10^6
                    } else if (as.integer(input$zoom) > 0 && input$genetic_measure == "genetic_dist"){
                        zoom = paste0("-zoomLevel_", as.integer(input$zoom), "cM")
                        zoom_level = as.integer(input$zoom)
                    }

                    dir.create(paste0(project_path, "/breakpointPlots/", mut_name, "/"), showWarnings = F)
                    filename = paste0(project_path, "/breakpointPlots/", mut_name, "/", gene, "-", mut_name, "-BreakpointPlot-",
                                      "-ancestMethod_", ancestral_method, "distMethod_", input$distMethod, "-clustMethod_", 
                                      input$clustMethod, "-", fam, fam_members, "-k_", k2, "-", input$genetic_measure,
                                      "-cutoff_", input$cutoff, db_val, zoom,"-group_", i, ".pdf")
                    # dir.create(paste0(project_path, "/nearestBreakpointPlots/", mut_name, "/"), showWarnings = F)
                    # filename3 = paste0(project_path, "/nearestBreakpointPlots/", mut_name, "/", gene, "-", mut_name, "-NearestBreakpointPlot-",
                    #                    "distMethod_", input$distMethod, "-clustMethod_", input$clustMethod, "-",
                    #                    fam, "-k_", k2, "-", input$genetic_measure, "-cutoff_", input$cutoff, db_val, "-group_", i, ".png")

                    plotname = paste0("BreakpointsPlot-Group_",i)
                    # plotname2 = paste0("NearestBreakpointsPlot-Group_",i)

                    # if (file.exists(filename)){
                    #     output[[plotname]] <- renderPlot({
                    #         readRDS(filename)
                    #     })
                    # } else {
                    #     output[[plotname]] <- renderPlot({
                    #         p=haplotype_mutation(mut, prepare = F, group = i)
                    #         saveRDS(p, file = filename)
                    #         p
                    #     })
                    # }
                    if (!file.exists(filename)){
                        incProgress(amount = 0, detail = paste(i, "of", k))
                        haplotype_data = haplotype_data()

                        #chr_pos = haplotype_mutation(mut, prepare = F, group = i, cutoff = input$cutoff, ancestral_method = ancestral_method)
                        chr_pos = haplotype_data[[i]]
                        
                        incProgress(amount = 1/(k*2), detail = paste(i, "of", k))

                        p = plot_breakpoints_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut,
                                                   group = i, col = input$genetic_measure, zoom = zoom_level)
                        p = add_second_legend(p, chr_pos)
                        ggsave(filename = filename, plot = p, scale = 0.8, width = 17.5, height = 10, dpi = 300)
                        #ggsave(filename2, scale = 1, width = 35, height = 20, dpi = 200, units = "cm")
                        #ggsave(filename2, scale = 1, width = 17.5, height = 10, dpi = 100)
                        incProgress(amount = 1/(k*2), detail = paste(i, "of", k))

                        # nearestBreakPlots = plot_nearest_brca_break_group(chr_pos, brca_chr, brca_name,
                        #                                                   brca_start, brca_stop, mut, group = i, col = input$genetic_measure)
                        # p2 = ggarrange(nearestBreakPlots[[1]], nearestBreakPlots[[4]], nearestBreakPlots[[2]], nearestBreakPlots[[3]], ncol=2, nrow=2)
                        # p2 = add_second_legend(p2, chr_pos)

                        # General ratio calculation for height:
                        #   assuming plots should have form factor width_p x height_p (e.g. 17,5 x height_p)
                        #   and assuming width_w (e.g. 1225) and height_w (e.g. 800) pixels of window available
                        #   Then the height_p is computed as: height_p = width_p / (width_w/height_w)

                        # Ratio conversion: 17,5/(1225/400) = 5,714 - for 3 plots on one row
                        #ggsave(filename3, scale = 1, width = 17.5, height = 5.714)

                        # Ratio conversion: 17,5/(1225/800) = 11,429 - for 2x2 plots on grid
                        # ggsave(filename = filename3, plot = p2, scale = 1, width = 17.5, height = 11.429)

                        # incProgress(amount = 1/(k*2), detail = paste(i, "of", k))
                    }
                    output[[plotname]] <- renderImage({
                        #list(src = filename, width = 1225, height = 700)
                        list(src = filename, width = "100%", height = 700)
                    }, deleteFile = F)
                    # output[[plotname2]] <- renderImage({
                    #     list(src = filename3, width = "100%", height = 800)
                    # }, deleteFile = F)
                })
            }
        })
    })

    # Deprecated
    plot_cluster_nearestBreakpoints <- reactive({
        # Make sure we update
        updatedField()
        k = getK()
        #brca_data <<- brca_data.filtered
        #colors_clustering()

        withProgress(message = "Computing breakpoints and creating plots", value = 0, {
            for (i in 1:k){
                if (length(cluster_groups[cluster_groups==i]) < input$req_group_members) {
                    incProgress(amount = 1/k, detail = paste(i, "of", k))
                    next
                }
                local({
                    i <- i
                    k2 = k
                    db_val = ""
                    if (input$clustMethod == "DBSCAN"){
                        db_val = paste0("-minPts_", input$minPts, "-eps_", input$eps)
                        k2 = 2 # Fix k to not make many identical plots
                    } else if (input$clustMethod == "HDBSCAN") {
                        db_val = paste0("-minPts_", input$minPts)
                        k2 = 2 # Fix k to not make many identical plots
                    }

                    dir.create(paste0(project_path, "/nearestBreakpointPlots/", mut_name, "/"), showWarnings = F)
                    filename = paste0(project_path, "/nearestBreakpointPlots/", mut_name, "/", gene, "-", mut_name, "-NearestBreakpointPlot-",
                                       "ancestMethod_", ancestral_method, "-distMethod_", input$distMethod, "-clustMethod_", input$clustMethod, "-",
                                       fam, fam_members, "-k_", k2, "-", input$genetic_measure, "-cutoff_", input$cutoff, db_val, "-group_", i, ".png")
                    plotname = paste0("NearestBreakpointsPlot-Group_",i)

                    if (!file.exists(filename)){
                        incProgress(amount = 0, detail = paste(i, "of", k))

                        haplotype_data = haplotype_data()
                        
                        #chr_pos <<- haplotype_mutation(mut, prepare = F, group = i, cutoff = input$cutoff, ancestral_method = ancestral_method)
                        chr_pos <- haplotype_data[[i]]
                        incProgress(amount = 1/(k*2), detail = paste(i, "of", k))

                        nearestBreakPlots = plot_nearest_brca_break_group(chr_pos, brca_chr, brca_name,
                                                                          brca_start, brca_stop, mut, group = i, col = input$genetic_measure)
                        p2 = ggarrange(nearestBreakPlots[[1]], nearestBreakPlots[[4]], nearestBreakPlots[[2]], nearestBreakPlots[[3]], ncol=2, nrow=2)
                        p2 = add_second_legend(p2, chr_pos)

                        # General ratio calculation for height:
                        #   assuming plots should have form factor width_p x height_p (e.g. 17,5 x height_p)
                        #   and assuming width_w (e.g. 1225) and height_w (e.g. 800) pixels of window available
                        #   Then the height_p is computed as: height_p = width_p / (width_w/height_w)

                        # Ratio conversion: 17,5/(1225/400) = 5,714 - for 3 plots on one row
                        #ggsave(filename3, scale = 1, width = 17.5, height = 5.714)

                        # Ratio conversion: 17,5/(1225/800) = 11,429 - for 2x2 plots on grid
                        ggsave(filename = filename, plot = p2, scale = 1, width = 17.5, height = 11.429)

                        incProgress(amount = 1/(k*2), detail = paste(i, "of", k))
                    }
                    output[[plotname]] <- renderImage({
                        list(src = filename, width = "100%", height = 800)
                    }, deleteFile = F)
                })
            }
        })
    })

    plot_haplotypes <- reactive({
        # Make sure we update
        updatedField()
        #brca_data <<- brca_data.filtered
        quantile = 0.05
        if (input$fam) {
            #fam = "fam"
            brca_data_temp <- brca_data %>% filter(Onc_ID %in% matched_plink_fam$SNP)
            countries <<- brca_data_temp %>% group_by(Country) %>% summarise(n=n()/nrow(.)) %>% filter(n > quantile) %>% .$Country
        } else {
            #fam = "individual"
            countries <<- brca_data %>% group_by(Country) %>% summarise(n=n()/nrow(.)) %>% filter(n > quantile) %>% .$Country
        }

        if (useSubCluster){
            k = as.integer(input$k_subClust)
        } else {
            k = getK()
        }

        #registerDoParallel(cores = 4)

        withProgress(message = "Computing haplotypes and creating plots", value = 0, {
            #foreach (i=1:k) %dopar% {
            for (i in 1:k) {
                if ((input$grouping == "Clustering" && length(cluster_groups[cluster_groups==i]) < input$req_group_members)){
                    #|| ((input$ancestral == "branchBoundIndep" || input$ancestral == "branchBound") && length(cluster_groups[cluster_groups==i]) < 3)) {
                    incProgress(amount = 1/k, detail = paste(i, "of", k))
                    next
                }
                local({
                    if (input$grouping == "Clustering") {
                        # k2 = k
                        # db_val = ""
                        # if (input$clustMethod == "DBSCAN"){
                        #     db_val = paste0("-minPts_", input$minPts, "-eps_", input$eps)
                        #     k2 = 2
                        # } else if (input$clustMethod == "HDBSCAN") {
                        #     db_val = paste0("-minPts_", input$minPts)
                        #     k2 = 2 # Fix k to not make many identical plots
                        # }
                        # 
                        # dir.create(paste0(project_path, "/haplotypePlots/", mut_name, "/"), showWarnings = F)
                        # filename = paste0(project_path, "/haplotypePlots/", mut_name, "/", gene, "-", mut_name, "-HaplotypePlot-",
                        #                "ancestMethod_", ancestral_method, "-distMethod-", input$distMethod, "-clustMethod-", input$clustMethod,
                        #                "-distMeasure-", input$genetic_measure, "-", fam, fam_members, "-k_", k2, "-cutoff_", input$cutoff,
                        #                db_val, "-group_", i, ".png")
                        
                        filename = get_filename(plot_type = "haplotypePlot", group = i, file_ext = "RData")
                        plotname = paste0("haplotypePlot-Group_",i)
                    } else { # Country
                        country = countries[i]
                        brca_data <<- filter(brca_data, Mut1HGVS %in% mut)
                        
                        # dir.create(paste0(project_path, "/haplotypePlotsCountry/", mut_name, "/"), showWarnings = F)
                        # filename = paste0(project_path, "/haplotypePlotsCountry/", mut_name, "/", gene, "-", mut_name,
                        #                   "-HaplotypePlot-ancestMethod_", ancestral_method, "-Country_", country,
                        #                   "-distMeasure-", input$genetic_measure, "-cutoff_", input$cutoff, "-", fam, fam_members, ".png")
                        
                        filename = get_filename(plot_type = "haplotypePlot", country = country, file_ext = "pdf")
                        plotname = paste0("haplotypePlot-", country)
                    }

                    if (!file.exists(filename)){
                        incProgress(amount = 1/(k*3), detail = paste(i, "of", k))
                        haplotype_data = haplotype_data()
                        
                        if (input$grouping == "Clustering") {
                            #chr_pos = haplotype_mutation(mut, prepare = F, group = i, cutoff = input$cutoff, ancestral_method = ancestral_method)
                            chr_pos = haplotype_data[[i]]
                            incProgress(amount = 1/(k*3), detail = paste(i, "of", k))
                            haploPlots = plot_haplotype(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group = i, col = input$genetic_measure)
                            p = ggarrange(haploPlots[[1]], haploPlots[[4]], haploPlots[[2]], haploPlots[[3]], ncol=2, nrow=2)
                            p = add_second_legend(p, chr_pos)
                        } else { # Country
                            #chr_pos = haplotype_mutation(mut, prepare = F, country = country, cutoff = input$cutoff, ancestral_method = ancestral_method)
                            chr_pos = haplotype_data[[i]]
                            incProgress(amount = 1/(k*3), detail = paste(i, "of", k))
                            haploPlots = plot_haplotype(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, col = input$genetic_measure)
                            p = ggarrange(haploPlots[[1]], haploPlots[[4]], haploPlots[[2]], haploPlots[[3]], ncol=2, nrow=2)
                            # p = add_second_legend(p, chr_pos, otherCountry = T) # not needed
                        }
                        if (tail(strsplit(basename(filename), split="\\.")[[1]], n = 1) == "RData"){
                            saveRDS(object = p, file = filename)
                            # For easy inspection later
                            filename2 = get_filename(plot_type = "haplotypePlot", group = i)
                            if(!file.exists(filename2)){
                                ggsave(filename = filename2, plot = p, scale = 1, width = 17.5, height = 11.429)
                            }
                        } else {
                            ggsave(filename = filename, plot = p, scale = 1, width = 17.5, height = 11.429)
                        }
                        #cowplot::save_plot(filename = filename, plot = p)

                        incProgress(amount = 1/(k*3), detail = paste(i, "of", k))
                    }

                    if (tail(strsplit(basename(filename), split="\\.")[[1]], n = 1) == "RData"){
                        output[[plotname]] <- renderPlot({
                            p = readRDS(filename)
                            grid.arrange(p)
                        })
                    } else {
                        output[[plotname]] <- renderImage({
                            list(src = filename, width = "100%", height = 800)
                        }, deleteFile = F)
                    }
                })
            }
            # for (i in 1:k) {
            #     local({
            #         if (input$grouping == "Clustering") {
            #             k2 = k
            #             db_val = ""
            #             if (input$clustMethod == "DBSCAN"){
            #                 db_val = paste0("-minPts_", input$minPts, "-eps_", input$eps)
            #                 k2 = 2
            #             } else if (input$clustMethod == "HDBSCAN") {
            #                 db_val = paste0("-minPts_", input$minPts)
            #                 k2 = 2 # Fix k to not make many identical plots
            #             }
            #
            #             dir.create(paste0(project_path, "/haplotypePlots/", mut_name, "/"), showWarnings = F)
            #             filename = paste0(project_path, "/haplotypePlots/", mut_name, "/", gene, "-", mut_name, "-HaplotypePlot-",
            #                               "distMethod-", input$distMethod, "-clustMethod-", input$clustMethod, "-",
            #                               fam, "-k_", k2, "-cutoff_", input$cutoff, db_val, "-group_", i, ".png")
            #             plotname = paste0("haplotypePlot-Group_",i)
            #         } else { # Country
            #             country = countries[i]
            #             brca_data <<- filter(brca2_pheno_merged, Mut1HGVS %in% mut)
            #             print(head(brca_data))
            #             dir.create(paste0(project_path, "/haplotypePlotsCountry/", mut_name, "/"), showWarnings = F)
            #             filename = paste0(project_path, "/haplotypePlotsCountry/", mut_name, "/", gene, "-", mut_name,
            #                               "-HaplotypePlot-Country_", country, "-cutoff_", input$cutoff, "-", fam, ".png")
            #             plotname = paste0("haplotypePlot-", country)
            #         }
            #         output[[plotname]] <- renderImage({
            #             list(src = filename, width = "100%", height = 800)
            #         }, deleteFile = F)
            #     })
            # }
        })
    })

    header <- function(i){
        if (input$grouping == "Clustering"){
            paste("Ancestral haplotype defined by:", "Group", i)
        } else {
            countries <- get_countries()
            paste("Ancestral haplotype defined by:", countries[i])
        }
    }

    output$breakpoints <- renderUI({

        prepare()
        countries <- get_countries()

        plotname <- function(i){
            if (input$grouping == "Clustering"){
                plotname <- paste0("BreakpointsPlot-Group_",i)
            } else {
                plotname <- paste0("BreakpointsPlot-", countries[i])
            }
        }

        if (useSubCluster){
            use_subCluster()
            k = as.integer(input$k_subClust)
        } else {
            k = getK()
        }
        
        plot_cluster_output_list <- lapply(1:k, function(i) {
            hc()
            if (length(cluster_groups[cluster_groups==i]) < input$req_group_members){
                list(tags$h4(header(i)), tags$html(
                    paste("Note: Group consists of less than", input$req_group_members, "members\n\n")),
                    HTML("<br><br><br>"))
            } else {
                plotname = plotname(i)
                list(tags$h4(header(i)), plotOutput(plotname, height = 720))
            }
        })

        plot_country_output_list <- lapply(1:k, function(i){
            if (!(countries[i] %in% names(table(fct_lump(brca_data$Country, prop = 0.05))))){
                list(tags$h4(header(i)), tags$html(
                    paste("Note: Country consists of less than", input$req_group_members, "members\n\n")),
                    HTML("<br><br><br>"))
            } else {
                plotname = plotname(i)
                list(tags$h4(header(i)), plotOutput(plotname, height = 720))
            }
        })

        plot_breakpoints()
        if (input$grouping == "Clustering"){
            #plot_cluster_breakpoints()
            # Convert the list to a tagList
            do.call(tagList, plot_cluster_output_list)
        } else { # Country
            #plot_country_breakpoints()
            # Convert the list to a tagList
            do.call(tagList, plot_country_output_list)
        }
    })

    output$nearest_breakpoint <- renderUI({

        prepare()
        #countries <<- brca_data %>% group_by(Country) %>% summarise(n=n()/nrow(.)) %>% filter(n > 0.05) %>% .$Country
        countries <- get_countries()

        k = getK()
        
        plotname <- function(i){
            if (input$grouping == "Clustering"){
                plotname <- paste0("NearestBreakpointsPlot-Group_",i)
            } else {
                plotname <- paste0("NearestBreakpointsPlot-", countries[i])
            }
        }

        #plot_cluster_output_list <- lapply(1:getK(), function(i) {
        plot_cluster_output_list <- lapply(1:k, function(i) {
            hc()
            if (length(cluster_groups[cluster_groups==i]) < input$req_group_members){
                list(tags$h4(header(i)), tags$html(
                    paste("Note: Group consists of less than", input$req_group_members, "members\n\n")),
                    HTML("<br><br><br>"))
            } else {
                plotname = plotname(i)
                list(tags$h4(header(i)), plotOutput(plotname, height = 820))
            }
        })

        #plot_country_output_list <- lapply(1:length(countries), function(i) {
        plot_country_output_list <- lapply(1:k, function(i) {
            if (!(countries[i] %in% names(table(fct_lump(brca_data$Country, prop = 0.05))))){
                list(tags$h4(header(i)), tags$html(
                    paste("Note: Country consists of less than", input$req_group_members, "members\n\n")),
                    HTML("<br><br><br>"))
            } else {
                plotname = plotname(i)
                list(tags$h4(header(i)), plotOutput(plotname, height = 820))
            }
        })

        # Convert the list to a tagList - this is necessary for the list of items
        # to display properly.
        plot_nearestBreakpoints()
        if (input$grouping == "Clustering"){
            #plot_cluster_nearestBreakpoints()
            do.call(tagList, plot_cluster_output_list)
        } else {
            #plot_country_breakpoints()
            do.call(tagList, plot_country_output_list)
        }
    })

    output$haplotypes_plot <- renderUI({
        prepare()
        countries <<- brca_data %>% group_by(Country) %>% summarise(n=n()/nrow(.)) %>% filter(n > 0.05) %>% .$Country

        plotname <- function(i){
            if (input$grouping == "Clustering"){
                plotname <- paste0("haplotypePlot-Group_",i)
            } else {
                plotname <- paste0("haplotypePlot-", countries[i])
            }
        }


        plot_cluster_output_list <- lapply(1:getK(), function(i) {
            hc()
            if (length(cluster_groups[cluster_groups==i]) < input$req_group_members){
                list(tags$h4(header(i)), tags$html(
                    paste("Note: Group consists of less than", input$req_group_members, "members\n\n")),
                    HTML("<br><br><br>"))
            # } else if ((input$ancestral == "branchBoundIndep" || input$ancestral == "branchBound") && length(cluster_groups[cluster_groups==i]) < 3){
            #     list(tags$h4(header(i)), tags$html(
            #         paste("Note: Group consists of less than 3 members\n\n")),
            #         HTML("<br><br><br>"))
            } else {
                list(tags$h4(header(i)), plotOutput(plotname(i), height = 820))
            }
        })

        plot_country_output_list <- lapply(1:length(countries), function(i){
            if (!(countries[i] %in% names(table(fct_lump(brca_data$Country, prop = 0.05))))){
                list(tags$h4(header(i)), tags$html(
                    paste("Note: Country consists of less than", input$req_group_members, "members\n\n")),
                    HTML("<br><br><br>"))
            } else {
                list(tags$h4(header(i)), plotOutput(plotname(i), height = 820))
            }
        })

        
        
        plot_haplotypes()
        if (input$grouping == "Clustering"){
            # Convert the list to a tagList
            do.call(tagList, plot_cluster_output_list)
        } else {
            # Convert the list to a tagList
            do.call(tagList, plot_country_output_list)
        }
    })

    ###########################################################################
    ### Migration
    ###########################################################################
    output$migration <- renderUI({
        prepare()
        # countries <<- brca_data %>% group_by(Country) %>% summarise(n=n()/nrow(.)) %>% filter(n > 0.05) %>% .$Country
        
        plotname <- function(i){
            plotname <- paste0("haplotypePlot-Group_",i)
        }
        
        
        # plot_cluster_output_list <- lapply(1:getK(), function(i) {
        plot_cluster_output_list <- lapply(1:2, function(i) {
            list(tags$h4(header(i)), div(plotOutput(plotname(i), height = 420), align="center"))
        })
        
        # plot_country_output_list <- lapply(1:length(countries), function(i){
        #     if (!(countries[i] %in% names(table(fct_lump(brca_data$Country, prop = 0.05))))){
        #         list(tags$h4(header(i)), tags$html(
        #             paste("Note: Country consists of less than", input$req_group_members, "members\n\n")),
        #             HTML("<br><br><br>"))
        #     } else {
        #         list(tags$h4(header(i)), plotOutput(plotname(i), height = 820))
        #     }
        # })
        
        output[[plotname(2)]] <- renderImage({
            list(src = "cache/BRCA2/migration/c.7617+1G_A/BRCA2-c.7617+1G_A-Migration-group2.png", height = 400)
        }, deleteFile = F)
        
        output[[plotname(1)]] <- renderImage({
            list(src = "cache/BRCA2/migration/c.7617+1G_A/BRCA2-c.7617+1G_A-Migration-group4.png", height = 400)
        }, deleteFile = F)
        
        do.call(tagList, plot_cluster_output_list)
        
        # plot_haplotypes()
        # if (input$grouping == "Clustering"){
        #     # Convert the list to a tagList
        #     do.call(tagList, plot_cluster_output_list)
        # } else {
        #     # Convert the list to a tagList
        #     do.call(tagList, plot_cluster_output_list)
        #     # do.call(tagList, plot_country_output_list)
        # }
    })
    
    ###########################################################################
    ### Cluster groups estimation
    ###########################################################################
    getScores <- reactive({
        if (input$k == "Automatic"){
            return(c("mcclain", "cindex", "silhouette", "dunn"))
        }
        if (input$k == "Automatic (extended)"){
            return(c("mcclain", "cindex", "silhouette", "dunn", "kl", "ch",
                     "hartigan", "ball", "ptbiserial", "gap", "gamma", "gplus",
                     "tau", "sdindex", "sdbw", "db", "hubert", "dindex"))
        }
    })

    output$numScores <- renderText({
        paste("Using", length(getScores()), "scores to estimate the number of clusters")
    })

    output$optimalClusters <- renderText({
        paste("Suggesting: <b>", getK(), "clusters</b>")
    })

    output$summary_cluster_estimation <- renderTable({
        k = getK()
        df = data.frame(numCluster)
        names(df) = c("# Clusters", "# Scores agreeing")
        df
    }, align = "c")

    output$clust_estimation <- renderUI({
        prepare()
        index <- getScores()
        #index <<- c("silhouette")
        k = getK()
        hc()
        fam = if (input$fam) "fam" else "individual"

        plot_output_list <- lapply(1:(length(index)), function(i) {
            filename = paste0(project_path, "/ClusterValidationPlots/", gene, "-", mut_name, 
                              "-ClusterValidationPlot-ancestMethod_", ancestral_method, "-", index[i],
                              "_score-distMethod-", input$distMethod, "-clustMethod-",
                              input$clustMethod, "-", fam, ".png")
            plotname = paste0("clust_estimate_", i)
            #local({
                # output[[plotname]] <- renderPlot({
                #     numCluster_plots[[i]]
                # })
                output[[plotname]] <- renderImage({
                    list(src = filename, width = "100%", height = 400)
                }, deleteFile = F)
            #})
            list(tags$h4(paste(index[i], "score")), plotOutput(plotname, height = 420))
        })

        # Convert the list to a tagList - this is necessary for the list of items
        # to display properly.
        do.call(tagList, plot_output_list)
    })

    output$nearestBreaksStatistics <- renderTable({
        prepare()
        k = getK()
        hc()
        for (i in 1:k){
            chr_pos = haplotype_mutation(mut, prepare = F, group = i, cutoff = input$cutoff, ancestral_method = ancestral_method, min_samples = input$min_anc_samples)
            nearestBreaksStatistics(chr_pos, mut, country = F, group = i, input = input)
        }
        # i=2
        # chr_pos = haplotype_mutation(mut, prepare = F, group = i, cutoff = input$cutoff)
        # nearestBreaksStatistics(chr_pos, mut, country = F, group = i, input = input)
    }, align = "c")

    ###########################################################################
    ### Visual comparison of two dendrograms
    ###########################################################################
    output$dendro_comparison_plot <- renderPlot({
        prepare()
        dst <- dst()
        dst2 <- dst2()

        k <- getK()
        compare_dendrograms(dst, dst2, input$clustMethod, input$clustMethod2)
    })

    ###########################################################################
    ### 
    ###########################################################################
    get_haplotype_table <- reactive({
        haplotype_data <- haplotype_data()
        k = getK()
        hap_table = list()
        for (i in 1:k){
            chr_pos2 = haplotype_data[[i]]
            
            pos <- chr_pos2 %>%
                group_by(sample_id) %>%
                filter(position_b37_adjusted >= 0) %>%
                slice(1) %>%
                rename(position_pos = position_b37, position_pos_adjusted = position_b37_adjusted,
                       cM_pos_adjusted = cM_adjusted, snp_pos = SNP, cM_pos = cM)
            
            neg <- chr_pos2 %>%
                group_by(sample_id) %>%
                filter(position_b37_adjusted <= 0) %>%
                arrange(sample_id, desc(position_b37_adjusted)) %>%
                slice(1) %>%
                rename(position_neg = position_b37, position_neg_adjusted = position_b37_adjusted,
                       cM_neg_adjusted = cM_adjusted, snp_neg = SNP, cM_neg = cM)
            
            chr_pos_minmax <- merge(pos, neg, by = intersect(names(pos), names(neg))) %>%
                mutate(haplo_length = position_pos-position_neg)
            
            
            if (input$grouping == "Clustering"){
                data = select(chr_pos_minmax, Sample = sample_id, FamCode, Country, 
                              Group = cluster_groups, Pos_right = position_pos_adjusted, 
                              Pos_left = position_neg_adjusted, cM_right=cM_pos_adjusted, 
                              cM_left=cM_neg_adjusted)
            } else {
                data = select(chr_pos_minmax, Sample = sample_id, FamCode, Country, 
                              Pos_right = position_pos_adjusted, Pos_left = position_neg_adjusted, 
                              cM_right=cM_pos_adjusted, cM_left=cM_neg_adjusted)
            }
            
            hap_table[[i]] = data
        }
        return(hap_table)
    })
    
    output$haplotype_data <- renderUI({
        prepare()
        #haplotype_data <- haplotype_data()
        k = getK()
        print(paste("k:", k))
        
        
        pretty_haplotype_table <- reactive({
            hap_table <- get_haplotype_table()
            hap_table2 <<- hap_table
            for (i in 1:k){
                data2 = datatable(hap_table[[i]], rownames = F,
                                  options = list(pageLength = 25,
                                                 #dom = 't',
                                                 columnDefs = list(list(className = 'dt-center', targets = 0:ncol(data)-1)))) %>%
                    formatStyle(1, fontWeight = "bold", borderright = "black") %>%
                    formatRound(columns = c(5:6), digits = 0, interval = 3, mark = " ") %>% 
                    formatRound(columns = c(7:8), digits = 2, interval = 3, mark = " ") #%>% formatStyle(columns = c(1:ncol(data)), 'text-align' = 'dt-center')
                
                # quoted = T needed to keep all the tables
                output[[paste0("haplotype_data_", i)]] = DT::renderDataTable(data2, quoted = T)
            }
        })
        
        haplotype_table_list <- lapply(1:k, function(i){
            list(tags$h4(header(i)), DT::dataTableOutput(paste0("haplotype_data_", i)), br())
            #list(tags$h4(header(i)), tableOutput(paste0("haplotype_data_", i)), br())
        })
        
        
        pretty_haplotype_table()
        do.call(tagList, haplotype_table_list)
    })
    
    output$mean_haplotype_data <- DT::renderDataTable({
        prepare()
        k = getK()
        hap_table <- get_haplotype_table()
        mean_cm_right = unlist(lapply(1:k, function(i) {
            hap=subset(hap_table2[[i]], Group == i)
            mean(hap$cM_right)
        }))
        mean_cm_left = unlist(lapply(1:k, function(i) {
            hap=subset(hap_table2[[i]], Group == i)
            mean(hap$cM_left)
        }))
        data = data.frame(Group=1:k, mean_cm_right=mean_cm_right, mean_cm_left=mean_cm_left)
        
        cols = ncol(data)
        datatable(data, rownames = F,
                  options = list(pageLength = 100,
                                 dom = 't',
                                 columnDefs = list(list(className = 'dt-center', targets = 0:(cols-1))))) %>%
            formatStyle(1, fontWeight = "bold", borderright = "black") %>%
            formatRound(columns = 2:cols, digits = 2, interval = 3, mark = " ") %>%
            formatRound(columns = 1, digits = 0, interval = 3, mark = " ")
    })
    
    ###########################################################################
    ### Estimating Mutation Age
    ###########################################################################
    prettify_datatable <- function(age_table){
        cols = ncol(age_table)
        age_table = datatable(age_table, rownames = F,
                              options = list(pageLength = 100,
                                             dom = 't',
                                             columnDefs = list(list(className = 'dt-center', targets = 0:(cols-1))))) %>%
            formatStyle(1, fontWeight = "bold", borderright = "black") #%>%
            #formatRound(columns = 2:cols, digits = 0, interval = 3, mark = " ")
        
        return(age_table)
    }
    
    compute_age_estimation <- function(method = 1){
        # Make sure we update
        prepare()
        k = getK()
        withProgress(message = "Computing Age Estimations", value = 1/(k^2+1), {
    
            filename1 = get_filename(mutAge_method = 1)
            # filename2 = get_filename(mutAge_method = 2)
            filename2 = ""
            filename3 = get_filename(mutAge_method = 3)
            # filename4_1 = get_filename(mutAge_method = "4_1")
            filename4_1 = ""
            # filename4_2 = get_filename(mutAge_method = "4_2")
            filename4_2 = ""
            filename5 = get_filename(mutAge_method = 5)
            filename6 = get_filename(mutAge_method = "6-gandolfo")
            filename6_cor = get_filename(mutAge_method = "6_cor-gandolfo")
            filename8 = get_filename(mutAge_method = 8)
    
            if (!(file.exists(filename1) || file.exists(filename2) || file.exists(filename6))){
                prepare()
                hc()
                dir.create(paste0(project_path, "/mutationAge_groups/", mut_name, "/"), showWarnings = F)
                computeMutationAge(gene, mut, isFam = input$fam, filename1 = filename1,
                                   filename2 = filename2, filename3 = filename3,
                                   filename4_1 = filename4_1, filename4_2 = filename4_2,
                                   filename5 = filename5, filename6 = filename6,
                                   filename6_cor = filename6_cor, filename8 = filename8, cluster = T, cutoff = input$cutoff,
                                   ancestral_method = ancestral_method, minSamples = input$min_anc_samples,
                                   progress = k^2+1)
            }
    
            if (method == 1){
                age_table <- read.table(filename1, header = T, sep = "\t", check.names = F)
                #print(age_table)
            } else if (method == 2) {
                age_table <- read.table(filename2, header = T, sep = "\t", check.names = F)
            } else if (method == 3) {
                age_table <- read.table(filename3, header = T, sep = "\t", check.names = F)
            } else if (method == "4_1") {
                age_table <- read.table(filename4_1, header = T, sep = "\t", check.names = F)
            } else if (method == "4_2") {
                age_table <- read.table(filename4_2, header = T, sep = "\t", check.names = F)
            } else if (method == 5) {
                age_table <- read.table(filename5, header = T, sep = "\t", check.names = F)
            } else if (method == 6) {
                age_table <- read.table(filename6, header = T, sep = "\t", check.names = F)
            } else if (method == "6_cor") {
                age_table <- read.table(filename6_cor, header = T, sep = "\t", check.names = F)
            } else if (method == 8) {
                age_table <- read.table(filename8, header = T, sep = "\t", check.names = F)
            }
    
            age_table = prettify_datatable(age_table)
            incProgress(amount = 1, detail = "Done")
            
        })
        return(age_table)
    }

    output$age_estimation1 <- DT::renderDataTable({

        #reset("mutation_age_table")
        # Make sure we update
        updatedField()
        #prepare()
        if (useSubCluster){
            k = as.integer(input$k_subClust)
        } else {
            k = getK()
        }

        age_table <- compute_age_estimation(method = 1)
        return(age_table)

    })

    output$age_estimation2 <- DT::renderDataTable({

        # Make sure we update
        updatedField()
        #prepare()
        #fam = if (input$fam) "fam" else "individual"
        k = getK()

        age_table <- compute_age_estimation(method = 2)
        return(age_table)

    })

    output$age_estimation3 <- DT::renderDataTable({

        # Make sure we update
        prepare()
        k = getK()

        age_table <- compute_age_estimation(method = 3)
        return(age_table)

    })

    output$age_estimation4_1 <- DT::renderDataTable({

        # Make sure we update
        prepare()
        k = getK()

        age_table <- compute_age_estimation(method = "4_1")
        return(age_table)

    })

    output$age_estimation4_2 <- DT::renderDataTable({

        # Make sure we update
        prepare()
        k = getK()

        age_table <- compute_age_estimation(method = "4_2")
        return(age_table)

    })
    
    output$age_estimation5 <- DT::renderDataTable({
        
        # Make sure we update
        prepare()
        k = getK()
        
        age_table <- compute_age_estimation(method = 5)
        return(age_table)
        
    })
    
    output$age_estimation6 <- DT::renderDataTable({
        
        # Make sure we update
        prepare()
        k = getK()
        
        age_table <- compute_age_estimation(method = 6)
        return(age_table)
        
    })
    
    output$age_estimation6_cor <- DT::renderDataTable({
        
        # Make sure we update
        prepare()
        k = getK()
        
        age_table <- compute_age_estimation(method = "6_cor")
        return(age_table)
        
    })
    
    output$age_estimation8 <- DT::renderDataTable({
        
        # Make sure we update
        prepare()
        k = getK()
        
        age_table <- compute_age_estimation(method = 8)
        return(age_table)
        
    })

    output$age_estimation_old <- renderUI({
        # Make sure we update
        updatedField()
        #prepare()
        #fam = if (input$fam) "fam" else "individual"

        k = if (input$grouping == "Clustering") getK() else length(countries)

        if (input$grouping == "Clustering"){
            db_val = ""
            if (input$clustMethod == "DBSCAN"){
                db_val = paste0("-minPts_", input$minPts, "-eps_", input$eps)
            } else if (input$clustMethod == "HDBSCAN") {
                db_val = paste0("-minPts_", input$minPts)
            }

            filename = paste0(project_path, "/mutationAge_groups/", mut_name, "/", gene, "-", mut_name, "-MutationAge-",
                              "distMethod-", input$distMethod, "-clustMethod-", input$clustMethod, "-",
                              fam, "-k_", k, "-cutoff_", input$cutoff, db_val, ".txt")

            if (file.exists(filename)){
                age = read.table(filename, header = T, sep = "\t")
            } else {
                prepare()
                hc()
                dir.create(paste0(project_path, "/mutationAge_groups/", mut_name, "/"), showWarnings = F)
                age = computeMutationAge(gene, mut, fam=input$fam, filename = filename, cluster = T)
            }

        } else {
            quantile = 0.05
            countries <<- brca_data %>% group_by(Country) %>% summarise(n=n()/nrow(.)) %>% filter(n > quantile) %>% .$Country
            dir.create(paste0(project_path, "/mutationAge_country/", mut_name, "/"), showWarnings = F)
        }


        mutationAge_list <- lapply(1:nrow(age), function(i){
            if (input$grouping == "Clustering"){
                if (length(cluster_groups[cluster_groups==i]) < input$req_group_members){
                    list(tags$h3(paste("Age of Mutation: Group", i)), tags$html(
                        paste("Note: Group consists of less than", input$req_group_members, "members\n\n")),
                        HTML("<br><br><br>"))
                } else {
                    list(tagList(tags$h3(paste("Age of Mutation: Group", i))),
                         tags$h4(paste(" - Method 1:", round(age[i, "age_geneticDist"], digits = 2), "generations")),
                         tags$h4(paste(" - Method 2:", round(age[i, "age_physicalDist"], digits = 2), "generations")),
                         HTML("<br>"))
                }
            } else {
                list(tags$h3(""))
            }

        })

        do.call(tagList, mutationAge_list)

    })
}

shinyApp(ui, server)
