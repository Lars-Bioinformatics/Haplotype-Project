setwd("~/Documents/PhD/haplotype-project")
#setwd("Z:/Users/lars/Documents/PhD/haplotype-project-windows")
library(shiny)
library(shinycssloaders)
library(markdown)


source("Scripts/Haplotype-projekt-ver2.R")
source("Scripts/Clustering.R")

## Read input data into global variables
#readData()
#readFamData()

ui <- fluidPage(#theme = "bootstrap-darkly.css", # Does not seem to so anything
    titlePanel("Haplotypes"),
    
    sidebarLayout(
        sidebarPanel(
            h4("Choose mutation"),
            # selectInput("mut",
            #     label = "Select mutation:",
            #     #choices = unique(brca1_pheno_merged$Mut1HGVS),
            #     choices = list("BRCA1"=count(brca1_pheno_merged, Mut1HGVS) %>% filter(n>100) %>% arrange(n) %>% .$Mut1HGVS,
            #                    "BRCA2"=count(brca2_pheno_merged, Mut1HGVS) %>% filter(n>100) %>% arrange(n) %>% .$Mut1HGVS),
            #     selected = "c.5123C>A"),
            selectInput("gene",
                        label = "Select gene:",
                        choices = list("BRCA1",
                                       "BRCA2"),
                        #selected = "BRCA1"),
                        selected = "BRCA2"),
            conditionalPanel(
                condition = "input.gene=='BRCA1'",
                selectInput("mut1",
                            label = "Select mutation:",
                            #choices = unique(brca1_pheno_merged$Mut1HGVS),
                            choices = count(brca1_pheno_merged, Mut1HGVS) %>% filter(n>100) %>% arrange(n) %>% .$Mut1HGVS,
                            #choices = c("c.1016delA", "c.1016dupA"),
                            selected = "c.5123C>A")
                            #selected = "c.1016delA")
            ),
            conditionalPanel(
                condition = "input.gene=='BRCA2'",
                selectInput("mut2",
                            label = "Select mutation:",
                            choices = count(brca2_pheno_merged, Mut1HGVS) %>% filter(n>100) %>% arrange(n) %>% .$Mut1HGVS,
                            selected = "c.7617+1G>A")
                            #selected = "c.2808_2811delACAA")
                            #selected = "c.1310_1313delAAGA")
            ),
            h4("Grouping data"),
            radioButtons("grouping", 
                         label = "Group data by:", 
                         choices = list('Clustering', 'Country'), 
                         selected = 'Clustering'),
            selectInput("fam", 
                        label = "Combine family haplotypes or consider individually:",
                        choices = list('Family' = TRUE,
                                       'Individually' = FALSE),
                        selected = TRUE),
            h4("Clustering"),
            selectInput("distMethod", 
                        label = "Select distance metric:",
                        choices = list('firstBreakDist',
                                       "firstUpperBreakDist",
                                       "firstUnderBreakDist",
                                       "matchstates",
                                       "genpofad",
                                       "mrca",
                                       "nei",
                                       "euclidean weighted according to brca dist"="euclideanDist_BRCAdistWeighted"
                                       # 'euclidean',
                                       # "centered pearson"="correlation",
                                       # "pearson",
                                       # "canberra",
                                       # "spearman",
                                       # "manhattan",
                                       # "kendall"
                                       ),
                        selected = 'firstBreakDist'),
            selectInput("clustMethod", 
                        label = "Clustering method:",
                        choices = list('complete',
                                       'ward.D',
                                       'ward.D2',
                                       'single'),
                        selected = 'complete'),
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
                        label = "Define number of clusters (k):", 
                        choices = list("Automatic", "Automatic (extended)", 1,2,3,4,5,6,7,8,9,10,11,12),
                        selected = "Automatic"),
            # h4("Ancestral plots"),
            # numericInput("group_selection",
            #             label = "Select group to show plot:",
            #             min = 1, max = 10,
            #             value = 1)
            h4("Miscellaneous"),
            selectInput("cutoff",
                        label = "Select frequency cutoff for ancestral haplotype",
                        choices = seq(0.5, 1, 0.1),
                        selected = 0.5),
            numericInput("req_group_members",
                        label = "Minimum required members in group",
                        min = 1, max = 100,
                        value = 3)
        , width = 3 ),
        mainPanel(
            tabsetPanel(
                tabPanel("Hierachical Clustering",
                        h2("Hierachical clustering"),
                        textOutput("mutInfo"),
                        withSpinner(plotOutput("hc_plot", height = 600), type = 1, color.background = "#FFFFFF"),
                        h4("Clustering information:"),
                        htmlOutput("clust_selected"),
                        h5("Groups are read from right to left in the dendrogram above i.e. group 1 is the rightmost branch etc."),
                        verbatimTextOutput("clust_info")),
                tabPanel("Cluster Validity Checking",
                         h2("Estimate number of groups in clustering"),
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
                tabPanel("Breakpoints",
                        #h2("Ancestral plots"),
                        h2("Visualizing breakpoints"),
                        h5("Plots all homozygous breakpoints found on either side of the BRCA gene"),
                        br(),
                        withSpinner(uiOutput("breakpoints"))
                ),
                tabPanel("Nearest Breakpoints",
                         h2("Visualizing nearest breakpoints"),
                         h5("Plots only the nearest homozygous breakpoint on either side of the BRCA gene"),
                         #HTML("<pre><b>Plot description:</b>\n    Left:   Ordered by group\n    Middle: Ordered by positive genomic position\n    Right:  Ordered by negative genomic position</pre>"),
                         HTML("<pre><b>Plot descriptions:</b>\n    Left-top:     Ordered by group\n    Right-top:    Ordered by unbroken haplotype length\n    Left-bottom:  Ordered by positive genomic position\n    Right-bottom: Ordered by negative genomic position</pre>"),
                         br(),
                         withSpinner(uiOutput("nearest_breakpoint"))
                ),
                tabPanel("Age Estimation",
                    h2("Coming soon...!"),
                    tableOutput("nearestBreaksStatistics")
                ),
                tabPanel("How To Use",
                    includeMarkdown("instructions.md")
                )
            )
            
        , width = 9 )
    )
)

server <- function(input, output) {
    
    ###########################################################################
    ### Prepare data
    ###########################################################################
    prepare <- reactive({
        gene <<- input$gene
        mut <<- if (input$gene == "BRCA1") input$mut1 else input$mut2
        prepare_dataframe(mut, gene)
        get_matched_plink()
        if (input$fam==T){
            extractFamPlink()
        }
    })
    
    ###########################################################################
    ### Specify number of clusters - either manually or automatically
    ###########################################################################
    getK <- reactive({
        if (input$k == "Automatic"){
            dst = dst()
            fam = if (input$fam) "fam" else "individual"
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
        }
        #print(k)
        return(k)
    })
    
    ###########################################################################
    ### Clustering tab
    ###########################################################################
    
    output$mutInfo <- renderText({
        prepare()
        paste("Number of samples:", nrow(matched_single), "- Number of families (with more than 2 members):", nrow(matched_fam))
    })
    
    getDist <- function(distMethod){
        fam = if (input$fam) "fam" else "individual"
        filename = paste0("cache/distMatrices/", gene, "-", mut_name, "-distMethod-", 
                          distMethod, "-", fam, ".rds")
        
        if (distMethod == "firstBreakDist"){
            dst <- firstBreakDist(input$fam, filename)
        } else if (distMethod %in% c("firstUpperBreakDist", "firstUnderBreakDist", "matchstates", "genpofad", "mrca", "nei")){
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
        #gene = input$gene
        #mut <- if (input$gene == "BRCA1") input$mut1 else input$mut2
        distMethod = input$distMethod
        getDist(distMethod)
    })
    
    dst2 <- reactive({
        distMethod = input$distMethod2
        getDist(distMethod)
    })
    
    hc <- reactive({
        print(2)
        dst <- dst()
        hiearchical_clustering(dst, clustMethod = input$clustMethod, k = getK(), doHapPlot = F)
    })
    
    output$hc_plot <- renderPlot({
        #mut <<- if (input$gene == "BRCA1") input$mut1 else input$mut2
        print(1)
        #prepare()
        
        print(input$mut1)
        print(input$mut2)
        print(input$gene)
        print(input$fam)
        print(input$k)
        
        # if (input$fam==T){
        #     famData()
        # }
        
        #matched_plink <<- prepare()
        
        print(3)
        hc <- hc()
        k = getK()
        
        plot_dendrogram(hc, k)
      
    }, height = 600)
    
    output$clust_selected <- renderText({
        paste("Number of clusters: <b>", getK(), "</b>")
    })
    
    output$clust_info <- renderPrint({
        # Reactive dependencies
        mut <- if (input$gene == "BRCA1") input$mut1 else input$mut2
        distMethod = input$distMethod
        k = getK()
        fam <- input$fam
        clustMethod <- input$clustMethod
        # Print group information
        groups <- list()
        for (i in 1:max(cluster_groups)){
            #print(subset(brca_data.filtered, cluster_groups == i) %>% select(FamCode, Onc_ID, Country, cluster_groups))
            cat(paste0("Group ", i, ":\n"))
            #print(subset(brca_data.filtered, cluster_groups == i) %>% count(Country))
            groups[[i]] <- as.data.frame(subset(brca_data.filtered, cluster_groups == i) %>% count(Country))# %>% select(Country, number=n))
            print(groups[[i]], row.names = F)
            if (i < max(cluster_groups)) cat("\n")
        }
        #print(groups)
    })
    
    
    ###########################################################################
    ### Breakpoints tab
    ###########################################################################
    
    colors_country <- reactive({
        country_freq <- table(fct_lump(brca_data$Country, prop = 0.05))
        if (length(country_freq) < 10){
            cols <- brewer.pal(max(length(country_freq), 3), "Set1")[1:(length(country_freq))]
        } else {
            cols <- colorRampPalette(brewer.pal(9, "Set1"))(length(country_freq))
        }
        names(cols) <- names(country_freq)
        cols <<- cols
    })
    
    colors_clustering <- reactive({
        hc()
        Country <- brca_data.filtered$Country
        if (length(unique(Country)) < 10){
            cols <- brewer.pal(max(length(unique(Country)), 3), "Set1")[1:(length(unique(Country)))]
        } else {
            cols <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(Country)))
        }
        country_freq <- sort(table(as.character(Country)), decreasing = T)
        names(cols) <- names(country_freq)
        cols <<- cols
    })
    
    plot_country_breakpoints <- reactive({
        prepare()
        k = getK()
        colors_country()
        fam = if (input$fam) "fam" else "individual"
        quantile = 0.05
        countries <<- brca_data %>% group_by(Country) %>% summarise(n=n()/nrow(.)) %>% filter(n > quantile) %>% .$Country
        
        withProgress(message = "Computing breakpoints and creating plots", value = 0, {
            i = 0
            k = length(countries)
            for (country in countries){
                i = i + 1
                local({
                    filename2 = paste0("cache/nearestBreakpointPlotsCountry/", gene, "-", mut_name, 
                                       "-BreakpointPlot-Country_", country, "-cutoff_", input$cutoff, "-", fam, ".png")
                    filename3 = paste0("cache/nearestBreakpointPlotsCountry/", gene, "-", mut_name, 
                                       "-NearestBreakpointPlot-Country_", country, "-cutoff_", input$cutoff, "-", fam, ".png")
                    plotname = paste0("BreakpointsPlot-", country)
                    plotname2 = paste0("NearestBreakpointsPlot-",country)
                    
                    if (!file.exists(filename2)){
                        incProgress(amount = 0, detail = paste(i, "of", k))
                        
                        chr_pos = haplotype_mutation(mut, prepare = F, country = country, cutoff = input$cutoff)
                        incProgress(amount = 1/(k*3), detail = paste(i, "of", k))
                        
                        p = plot_haplotype(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, country = country)
                        ggsave(filename2, scale = 1, width = 17.5, height = 10)
                        incProgress(amount = 1/(k*3), detail = paste(i, "of", k))
                        
                        nearestBreakPlots = plot_nearest_brca_break_country(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, country = country)
                        p2 = ggarrange(nearestBreakPlots[[1]], nearestBreakPlots[[2]], nearestBreakPlots[[3]], ncol=3, nrow=1)
                        
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
                        list(src = filename3, width = "100%", height = 400)
                    }, deleteFile = F)
                })
            }
        })
    })
    
    plot_cluster_breakpoints <- reactive({
        # Preparing data - could be avoided if plots already cached, however no huge performance penalty
        prepare()
        #mut <- if (input$gene == "BRCA1") input$mut1 else input$mut2
        k = getK()
        fam = if (input$fam) "fam" else "individual"
        #distMethod = input$distMethod
        hc()
        brca_data <<- brca_data.filtered
        colors_clustering()
        #req_group_members <<- 3
        
        withProgress(message = "Computing breakpoints and creating plots", value = 0, {
            for (i in 1:k){
                if (length(cluster_groups[cluster_groups==i]) < input$req_group_members) {
                    incProgress(amount = 1/k, detail = paste(i, "of", k))
                    next
                }
                local({
                    i <- i
                    # filename = paste0("cache/breakpointPlots/", gene, "-", mut_name, "-BreakpointPlot-", 
                    #                   "distMethod-", input$distMethod, "-clustMethod-", input$clustMethod, "-",
                    #                   fam, "-k_", k, "-group_", i, ".rds")
                    filename2 = paste0("cache/breakpointPlots/", gene, "-", mut_name, "-BreakpointPlot-", 
                                      "distMethod-", input$distMethod, "-clustMethod-", input$clustMethod, "-",
                                      fam, "-k_", k, "-cutoff_", input$cutoff, "-group_", i, ".pdf")
                    filename3 = paste0("cache/nearestBreakpointPlots/", gene, "-", mut_name, "-NearestBreakpointPlot-", 
                                       "distMethod-", input$distMethod, "-clustMethod-", input$clustMethod, "-",
                                       fam, "-k_", k, "-cutoff_", input$cutoff, "-group_", i, ".pdf")
                    plotname = paste0("BreakpointsPlot-Group_",i)
                    plotname2 = paste0("NearestBreakpointsPlot-Group_",i)
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
                    if (!file.exists(filename2)){
                        incProgress(amount = 0, detail = paste(i, "of", k))
                        
                        chr_pos = haplotype_mutation(mut, prepare = F, group = i, cutoff = input$cutoff)
                        incProgress(amount = 1/(k*3), detail = paste(i, "of", k))
                        
                        p = plot_haplotype_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group = i)
                        ggsave(filename2, scale = 1, width = 17.5, height = 10)
                        incProgress(amount = 1/(k*3), detail = paste(i, "of", k))
                        
                        nearestBreakPlots = plot_nearest_brca_break_group(chr_pos, brca_chr, brca_name, 
                                                                          brca_start, brca_stop, mut, group = i)
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
                    # output[[plotname2]] <- renderPlot({
                    #     p2
                    # })
                })
            }
        })
        #plot_cluster_breakpoints
    })
    
    header <- function(i){
        if (input$grouping == "Clustering"){
            paste("Group", i, "consensus")
        } else {
            paste(countries[i], "consensus")
        }
    }
    
    output$breakpoints <- renderUI({
        
        plotname <- function(i){
            if (input$grouping == "Clustering"){
                plotname <- paste0("BreakpointsPlot-Group_",i)
            } else {
                plotname <- paste0("BreakpointsPlot-", countries[i])
            }
        }
        
        if (input$grouping == "Clustering"){
            plots <- plot_cluster_breakpoints()
        } else {
            plots <- plot_country_breakpoints()
        }
            
        
        plot_output_list <- lapply(1:getK(), function(i) {
            if (input$grouping == "Clustering" && length(cluster_groups[cluster_groups==i]) < input$req_group_members){
                list(tags$h4(header(i)), tags$html(
                    paste("Note: Group consists of less than", input$req_group_members, "members\n\n")), 
                    HTML("<br><br><br>"))
            } else if (input$grouping == "Country" && !(countries[i] %in% names(table(fct_lump(brca_data$Country, prop = 0.05))))){
                list(tags$h4(header(i)), tags$html(
                    paste("Note: Group consists of less than", input$req_group_members, "members\n\n")), 
                    HTML("<br><br><br>"))
            } else {
                plotname = plotname(i)
                list(tags$h4(header(i)), plotOutput(plotname, height = 720))
            }
        })
        
        # Convert the list to a tagList - this is necessary for the list of items
        # to display properly.
        do.call(tagList, plot_output_list)
    })
    
    output$nearest_breakpoint <- renderUI({
        
        plotname <- function(i){
            if (input$grouping == "Clustering"){
                plotname <- paste0("NearestBreakpointsPlot-Group_",i)
            } else {
                plotname <- paste0("NearestBreakpointsPlot-", countries[i])
            }
        }
        
        if (input$grouping == "Clustering"){
            plots <- plot_cluster_breakpoints()
        } else {
            plots <- plot_country_breakpoints()
        }
        
        plot_output_list <- lapply(1:getK(), function(i) {
            if (input$grouping == "Clustering" && length(cluster_groups[cluster_groups==i]) < input$req_group_members){
                list(tags$h4(header(i)), tags$html(
                    paste("Note: Group consists of less than", input$req_group_members, "members\n\n")), 
                    HTML("<br><br><br>"))
            } else if (input$grouping == "Country" && !(countries[i] %in% names(table(fct_lump(brca_data$Country, prop = 0.05))))){
                list(tags$h4(header(i)), tags$html(
                    paste("Note: Group consists of less than", input$req_group_members, "members\n\n")), 
                    HTML("<br><br><br>"))
            } else {
                plotname = plotname(i)
                list(tags$h4(header(i)), plotOutput(plotname, height = 820))
            }
        })
        
        # Convert the list to a tagList - this is necessary for the list of items
        # to display properly.
        do.call(tagList, plot_output_list)
    })
    
    # output$breakpoints <- renderPlot({
    #     
    #     plots <- plot_cluster_breakpoints()
    #     
    #     print(length(plots))
    #     print(input$group_selection)
    #     print(as.integer(input$group_selection))
    #     
    #     group_selection = if (input$group_selection > input$k) input$k else input$group_selection
    #     print(group_selection)
    #     
    #     #p = plots[[as.integer(group_selection)]] 
    #     #p = plots[[1]] 
    #     plots
    # }, height = 700)
    
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
            filename = paste0("cache/ClusterValidationPlots/", gene, "-", mut_name, "-ClusterValidationPlot-", index[i],
                              "_score-distMethod-", input$distMethod, "-clustMethod-", input$clustMethod, "-", fam, ".pdf")
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
        # for (i in 1:k){
        #     chr_pos = haplotype_mutation(mut, prepare = F, group = i, cutoff = input$cutoff)
        #     nearestBreaksStatistics(chr_pos, mut, country = F, group = i, input = input)
        # }
        i=2
        chr_pos = haplotype_mutation(mut, prepare = F, group = i, cutoff = input$cutoff)
        nearestBreaksStatistics(chr_pos, mut, country = F, group = i, input = input)
    })
    
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
}

shinyApp(ui, server) 
