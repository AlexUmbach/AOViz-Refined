## LIBRARIES:-
library(shiny)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ecodist)
library(tidyr)
library(forcats)
library(vegan)
library(ape)
library(ade4)
library(gtools)
library(SRS)
library(vroom)
library(reactable)
library(shinycssloaders)
library(shinydashboard)
library(stringi)
library(shinyWidgets)
library(fresh)
library(htmltools)
library(tidytext)
library(ggh4x)
library(shinyjs)
library(patchwork)
library(viridis)
library(GUniFrac)

#main page
Rversion = R.Version()$version.string
bloop_welcome = "
AOViz contains several R-based visualization scripts wrapped within an R Shiny UI, allowing users to quickly explore their short-read amplicon data.
Figures can be customized using a provided metadata file and exported as PDFs.
"

bloop_req = "AOViz requires an ASV/OTU table produced from sequence data processing pipelines (e.g., QIIME2) using the QIIME2/SILVA taxonomy formatting.
Any phylogenetic marker gene can be used (e.g., cpn60, amoA) so long as the taxonomy formatting is identical. All uploaded files must be in txt format. Your ASV/OTU table should be formatted as detailed in the example picture below:"
# asv_example
# meta_example


## This is some custom HTML settings for colour, background, textsize, and whatnot:
ui <- navbarPage(id = "navbarID",
                 setBackgroundColor(color = "#E7E7E7"),
                 # setBackgroundColor(color = "#E7E7E7"),
                 tags$style(HTML("
                                 
        .navbar-default .navbar-brand {
        color:white;
        font-size:28px
        }
        
        .navbar-default .navbar-brand:hover {
        color:white;
        }
        
        .navbar {
        background-color:#3C8DBC;
        }
        
        .navbar-default .navbar-nav > li > a {
        color:white;
        font-size:22px;
        black
        }
        
        .navbar-default .navbar-nav > .active > a,
        
        .navbar-default .navbar-nav > .active > a:focus,
        
        .navbar-default .navbar-nav > .active > a:hover {
        color:black;
        background-color:white;
        }
        .navbar-default .navbar-nav > li > a:hover {
        color:black;
        background-color:white;
        text-decoration
        }
        
        .well {
        background:white;
        }
        
        body {
        background-color:white;
        }
        
        #MetaTableOut {
          zoom: 0.80;
        }
        
        #MainTableOut {
          zoom: 0.80;
        }
        
        #ContamTableOut {
          zoom: 0.80;
        }
        
        #FinalProcessedTable {
          zoom: 0.80; nowrap;
        }
        
        #ReadTableOut {
          zoom: 0.80; nowrap;
        }
        
        #BarTableOut {
          zoom: 0.80;
        }
        
        #BubbleTableOut {
          zoom: 0.80;
        }
        
                                 ")),
                 
                 title = "AOViz v3.0",
                 
                 #title = img(src="ALEXIOME_logo.png", height = "100%", width = "100%"),
                 
                 #### Main Page ####
                 tabPanel("WELCOME",
                          sidebarLayout(
                            sidebarPanel(Rversion,
                                         p("Leave any comments or suggestions on", a("the git page", href = "https://github.com/AlexUmbach/AOViz")),
                                         img(src="MainLogo.svg", height = "100%", width = "100%"),
                                         style = "height: 500px; position:relative; border-color:#000000",
                                         width = 3
                            ),
                            mainPanel(width = 9,
                                      fluidRow(
                                        box(h1("Welcome to AOViz!"),
                                            p(bloop_welcome),
                                            p(h1("What kind of data are required?")),
                                            p(bloop_req),
                                            width = 12,
                                            style = "background-color:#FFFFFF; border-color:#ffffff; border-style: solid; border-width: 1.5px; margin-left:0px; margin-right:30px; padding: 10px",
                                        ),
                                      ),
                                      box(
                                        p("This is how your ASV/OTU table should be formatted. The column",strong("Feature.ID"),"is required and must contain unique identifiers for each row."),
                                        p("ASV/OTU tables can also be in a collapsed format, where the 'Feature.ID' column contains taxonomy. The",strong("Consensus.Lineage"),"and",strong("ReprSequence"),"headers must be as shown. The sample column names should not contain special characters (e.g., *, -, (), [], /). If you're encountering errors, try checking this first."),
                                        img(src="ASV_example.png",height = "60%", width = "60%"),
                                        width = 13,
                                        style = "background-color:#FFFFFF; border-color:#ffffff; border-style: solid; border-width: 1.5px; margin-left:0px; margin-right:30px; padding: 10px",
                                      ),
                                      br(),
                                      box(
                                        p("Your metadata file should appear as below. There is a single required feature of your metadata table: a column labelled",strong("SampleName."),
                                          "The 'SampleName' column should include a list of sample names that are",strong("identical"), "to the samples in your ASV/OTU table. Any additional columns containing sample information (e.g., temperature, location, group) can be included and will be incorporated as options during analyses."),
                                        img(src="meta_example.png",height = "50%", width = "50%"),
                                        width = 13,
                                        style = "background-color:#FFFFFF; border-color:#ffffff; border-style: solid; border-width: 1.5px; margin-left:0px; margin-right:30px; padding: 10px",
                                      ),
                                      br(),
                                      box(
                                        p("If interested in analyzing taxonomically collapsed table (such as a genus-collapsed), you need only rename the taxonomy column to 'Feature.ID', as below:"),
                                        img(src="collapsed_Example.png",height = "50%", width = "50%"),
                                        width = 13,
                                        style = "background-color:#FFFFFF; border-color:#ffffff; border-style: solid; border-width: 1.5px; margin-left:0px; margin-right:30px; padding: 10px",
                                      ),
                                      box(
                                        p(h1("Column name blacklist")),
                                        p("AOViz requires priority when assigning column names. Because of this, you",strong("must not use"),"the following names in your metadata:"),
                                        p("TaxaName"),
                                        p("input"),
                                        p("Feature ID (or variants"),
                                        p("rowID (or variants)"),
                                        p("Taxonomy (or variants)"),
                                        width = 13,
                                        style = "background-color:#FFFFFF; border-color:#ffffff; border-style: solid; border-width: 1.5px; margin-left:0px; margin-right:30px; padding: 10px",
                                      )
                                      
                            )
                          )
                 ),
                 
                 
                 
                 #### Upload Data and review ####
                 tabPanel("Data upload",
                          sidebarLayout(
                            sidebarPanel(h5("Please upload your data"),
                                         fileInput("MainFile","Primary ASV table"),
                                         checkboxInput("IsMainCollapsed","Is this a collapsed table?",value = FALSE),
                                         fileInput("MetaFile","Metadata table"),
                                         fileInput("ContamFile","ASV contaminant list"),
                                         width = 3,
                                         style = "overflow-y:scroll; max-height: 850px; position:relative; border-color:#000000"
                            ),
                            mainPanel(width = 9,
                                      fluidRow(
                                        box(h3("Upload your data"),
                                            p("This is where you can upload your ASV tables, metadata table, and",em("optional"),"contaminant list. Once you have uploaded your ASV table, indicate whether it is a true ASV table or a collapsed taxon table (e.g., genus-level). If you encounter errors, please check your tables to ensure proper column names and formatting."
                                            ),
                                            width = 12,
                                        ),
                                        style = "background-color:#FFFFFF; border-color:#ffffff; border-style: solid; border-width: 1.5px; margin-left:0px; margin-right:30px; padding: 10px",
                                        width = 6,        
                                      ),
                                      textOutput("maintext"),
                                      dataTableOutput("MainTableOut"),
                                      textOutput("metatext"),
                                      dataTableOutput("MetaTableOut"),
                                      textOutput("contamtext"),
                                      dataTableOutput("ContamTableOut"),
                                      style = "overflow-y:scroll; max-height: 850px; position:relative;"
                                      
                            )
                          )
                 ),
                 
                 
                 #### Processed Data #####
                 tabPanel("Processed data",
                          sidebarLayout(
                            sidebarPanel(h4("Select preferences"),
                                         # radioButtons("RemovePrefix","Do you want to remove prefixes?",c("Yes","No"),selected = "Yes", inline = TRUE),
                                         # radioButtons("TruncateTaxa","Do you want to truncate taxa to the most resolved taxon?",c("Yes","No"), selected = "Yes",inline = TRUE),
                                         radioButtons("ContamChoice","Do you want to remove or analyze contaminants?", c("Remove","Analyze","No"),selected = "No",inline = TRUE),
                                         checkboxInput("RemoveLowReads","Do you want to remove low abundance reads?", value = FALSE),
                                         
                                         conditionalPanel(
                                           condition = "input.RemoveLowReads == true",
                                           numericInput("ReadThreshold","Set a threshold",value = 100),
                                         ),
                                         
                                         width = 3,
                                         style = "overflow-y:scroll; max-height: 900px; position:relative;border-color:#000000"
                            ),
                            mainPanel(width = 9,
                                      fluidRow(
                                        box(h3("Processed data review"),
                                            p("This page is where you can review your processed data to ensure that everything is looking as it should. The table is interactable and searchable and should allow you to cross reference samples and metadata with your original tables if you have any concerns
                                        about whether the data has been processed properly."),
                                            width = 12)),
                                      br(),
                                      textOutput("proc_new_text"),
                                      dataTableOutput("FinalProcessedTable"),
                                      style = "overflow-y:scroll; max-height: 800px; position:relative;"
                                      
                                      
                            )
                          )
                 ),
                 
                 
                 
                 
                 #### Total read plot ####
                 tabPanel("Read plot", value = "readtab",
                          sidebarLayout(
                            sidebarPanel(h4("Select preferences"),
                                         radioButtons("BoxSelect","What kind of plot?",choices = c("Box","Bar"),selected = "Bar", inline = TRUE),
                                         hr(style = "border-width: 3px; border-color:#A9A9A9"),
                                         numericInput("ReadYAxisLimit","Change the y-axis limit",value = 100000),
                                         selectInput("ReadSortByAxis","How do you want to group your data?",choices = NULL),
                                         hr(style = "border-width: 3px; border-color:#A9A9A9"),
                                         hr(style = "border-width: 3px; border-color:#A9A9A9"),
                                         radioButtons("ReadPanel","Do you want panel borders?",c("Yes","No"),selected = "Yes", inline = TRUE),
                                         sliderInput("ReadPanelSpacing","Modify panel spacing",min = 0,max = 20,step = 1,value = 0),
                                         sliderInput("ReadWidth","Change bar width",min=0, max=1, step = 0.1, value = 0.9),
                                         width = 3,
                                         style = "overflow-y:scroll; max-height: 850px; position:relative;border-color:#000000"
                                         
                            ),
                            mainPanel("",
                                      width = 9,
                                      fluidRow(
                                        box(h3("Sequence read plot"),
                                            p("This is an ASV or sequence read plot. It summarizes your total read counts in your dataset, organized by sample name or chosen metadata category. The script uses base R functions to sum the columns for each sample in an ASV table and reports those totals as 'total reads'. Samples with read counts higher than the chosen y-axis limit will be removed from the dataset, so make sure you set an appropriate limit."),
                                            width = 12)),
                                      br(),
                                      dataTableOutput("ReadTableOut") %>% withSpinner(type = 8, color.background = "white"),
                                      fluidRow(
                                        box(
                                          sliderInput("ReadPlotOutW","Plot width",min = 0, max = 2000,step = 100,value = 600),
                                          sliderInput("ReadPlotOutH","Plot height",min = 0, max = 2000,step = 100,value = 400),
                                          actionButton("ReadStartButton",label = "Start!"),
                                          downloadButton("ReadDownload","Save figure"),
                                          downloadButton("ReadTableDownload", "Download table"),
                                        ),
                                      ),
                                      br(),
                                      br(),
                                      plotOutput("ReadPlotOut") %>% withSpinner(type = 1,color.background = "white",size = 3),
                                      style = "overflow-y:scroll; max-height: 850px; position:relative;"
                            )
                          )
                 ),
                 
                 
                 #### Taxa bar plot analysis ####
                 tabPanel("Taxonomy relative abundance",
                          sidebarLayout(
                            sidebarPanel(
                              numericInput("BarCutOff","Select your cutoff",value = 10),
                              actionButton("BarStartButton",label = "Start!",width = 334),
                              hr(style = "border-width: 3px; border-color:#A9A9A9"),
                              selectInput("BarSortByAxis", label = "How do you want to group your samples?",choice = "Updating"),
                              checkboxInput("BarRenameCheck", "Do you want to relabel your samples with a SampleShort column?", FALSE),
                              selectInput("BarTaxonLevel", label = "Colour by which level?", choice = c("Species","Genus","Family","Order","Class","Phylum","Sub1","Sub2","Sub3","Sub4")),
                              checkboxInput("BarFacet","Do you want a second facet?",value = FALSE),
                              
                              conditionalPanel(
                                condition = "input.BarFacet == true",
                                selectInput("BarSecondFacet","Choose the second ordering",choices = "Updating"),
                              ),
                              hr(style = "border-width: 3px; border-color:#A9A9A9"),
                              radioButtons("BarPanelBorder","Do you want panel borders?",choices = c("Yes","No"),selected = "Yes", inline = TRUE),
                              sliderInput("BarPanelSpacing","Modify panel spacing",min = 0, max = 30, value = 0, step = 1),
                              sliderInput("BarAlpha","Select your bar alpha",min= 0,max= 1,step= 0.1,value= 0.4),
                              width = 3,
                              style = "overflow-y:scroll; max-height: 850px; position:relative;border-color:#000000"
                              
                            ),
                            mainPanel(width = 9,
                                      fluidRow(
                                        box(h3("Taxonomy plot"),
                                            p("This is a taxonomy bar plot. It represents the relative abundances of each taxon within and among samples. The 'cut-off' value represents the relative abundance required for the taxon to be represented by name within the plot. For example, if set to '20', only taxa that have relative abundances equal to or greater than 20%, in any sample, will be presented in the plot. All other taxa are lumped into the 'other' category."),
                                            p("To create this barplot, the ASV table is first transformed into a proportion table and all taxa below the cut-off are removed. Each column is summed and substracted from 100, representing the total proportion of taxa not present at greater than the cut-off threshold (i.e., the 'other' category). These 'other' taxa are assigned the taxonomy identifer 'ZZOther' and then plotted."),
                                            p("You may also filter your table for specific taxa at any taxonomy level (e.g., Staphylococcus, Pseudomonadota)."),
                                            width = 12)),
                                      dataTableOutput("BarTableOut") %>% withSpinner(type = 1,color.background = "white"),
                                      fluidRow(
                                        box(
                                          sliderInput("TaxaPlotOutW","Plot width",min = 0, max = 2000,step = 100,value = 600),
                                          sliderInput("TaxaPlotOutH","Plot height",min = 0, max = 2000,step = 100,value = 400),
                                          downloadButton("BarDownload","Save figure"),
                                          downloadButton("BarTableDownload","Save table"),
                                          width = 3
                                        )
                                      ),
                                      br(),
                                      br(),
                                      plotOutput("TaxaPlotOut") %>% withSpinner(type = 1,color.background = "white",size = 3),
                                      style = "overflow-y:scroll; max-height: 850px; position:relative;"
                                      
                            )
                          )
                 ),
                 
                 
                 #### Bubble plot ####
                 tabPanel("Relative Abundance Bubble plot",
                          sidebarLayout(
                            sidebarPanel(
                              textInput("BubbleAbundThresh","Relative abundace threshold (%)", value = 5),
                              actionButton("BubbleStartButton",label = "Start!", width = 335),
                              hr(style = "border-width: 3px; border-color:#A9A9A9"),
                              radioButtons("BubbleInclPercent","Do you want to include percent abundance numbers?",c("Yes","No"), selected = "Yes", inline = TRUE),
                              numericInput("BubbleDec","Set the number of decimals (0 to 5)",value = 2,min = 0, max = 5),
                              hr(style = "border-width: 3px; border-color:#A9A9A9"),
                              selectInput("BubbleSortXAxis","How do you want to order your individual samples?",choices = sort("Updating")),
                              selectInput("BubbleFacet","How do you want to group samples (i.e., faceting)?", choices = "Updating"),
                              selectInput("BubbleColour","How do you want to colour your bubbles?", choices = "Updating"),
                              selectInput("BubbleTaxSort","By what taxonomic level will the y-axis be sorted?",choice = c("Species","Genus","Family","Order","Class","Phylum","Sub1","Sub2","Sub3","Sub4")),
                              checkboxInput("BubbleRename","Do you want to rename your samples on the plot using a SampleShort column?", FALSE),
                              hr(style = "border-width: 3px; border-color:#A9A9A9"),
                              checkboxInput("BubbleInclRead","Do you want to include sample counts in a read plot?", FALSE),
                              checkboxInput("BubbleInclTaxa", "Do you want to include taxon read proportions?", FALSE),
                              checkboxInput("BubbleFakeTaxon","Do you want to show all samples regardless of present taxa?", value = FALSE),
                              checkboxInput("BubbleFactorData","Do you want to factorize data?", value = FALSE),
                              hr(style = "border-width: 3px; border-color:#A9A9A9"),
                              checkboxInput("BubbleSecondFacet", "Do you want a second facet?", FALSE),
                              
                              conditionalPanel(
                                condition = "input.BubbleSecondFacet == true",
                                selectInput("BubbleSecondFacetMeta","Choose the second ordering",choices = "Updating"),
                                checkboxInput("BubbleThirdFacet", "Do you want a third facet?", FALSE),
                              ),
                              
                              conditionalPanel(
                                condition = "input.BubbleThirdFacet == true",
                                selectInput("BubbleThirdFacetMeta","Choose the third ordering",choices = "Updating"),
                              ),
                              hr(style = "border-width: 3px"),
                              textInput("BubbleTaxKeywords","Select specific or multiple taxa (example: staph,coryn). Must be separated by a comma; do not include spaces"),
                              textInput("BubbleTaxRemove", "Remove specific or multiple taxa (example: staph,coryn). Must be separated by a comma; do not include spaces"),
                              selectInput("BubbleMetaFilt","Select a metadata category for data filtering",choices = "Updating"),
                              textInput("BubbleMetaKeywords","Provide a specific keyword for the metadata filtering. For multiple, use a comma with no spaces (Example: canada,japan"),
                              hr(style = "border-width: 3px; border-color:#A9A9A9"),
                              radioButtons("BubblePanelBorder","Do you want panel borders?",choices = c("Yes","No"),selected = "No", inline = TRUE),
                              radioButtons("BubbleFacetSideX","Where do you want the sample bar?",c("Top","Bottom"),selected = "Bottom",inline = TRUE),
                              radioButtons("BubbleFacetSideY","Where do you want the taxon bar?",c("Left","Right"),selected = "Left",inline = TRUE),
                              sliderInput("BubblePanelSpacing","Modify the facet spacing",min = 0, max = 20, value = 0),  
                              width = 3,
                              style = "overflow-y:scroll; max-height: 850px; position:relative;border-color:#000000"
                            ),
                            mainPanel(width = 9,
                                      fluidRow(
                                        box(h3("Relative abundance bubble plot"),
                                            p("The relative abundance bubble plot shows the distribution of taxa among samples and their relative proportions in the dataset. Bubble size and values represent the proportion of reads associated with a specific taxon within a given sample."),
                                            p("To construct the bubble plot, an ASV table (or a collapsed table) is transformed into a proportion table by converting each read count into a proportion of total reads for a given sample. This table is then filtered using a defined relative abundance threshold, only keeping taxa if they are present at or above the threshold in any sample. Taxa are assigned uniqe identifiers, metadata is appended, and data is plotted as desired by the user."),
                                            p("The table can be filtered or organized as desired. Faceting can be done based on any supplied metadata, and specific taxa or samples can be emphasized (all other samples removed) using the sidebar."),
                                            p("A read plot can be appended to the plot showing the total read counts for each sample. A bar plot showing the proportion of taxa across all samples can also be added. In this case, read counts associated with each taxon are summed and then divided by the total read count of the entire dataset"),
                                            width = 12)),
                                      br(),
                                      dataTableOutput("BubbleTableOut"),
                                      fluidRow(box(
                                        sliderInput("BubblePlotOutW","Plot width",min = 0, max = 4000,step = 10,value = 600),
                                        sliderInput("BubblePlotOutH","Plot height",min = 0, max = 4000,step = 10,value = 600),
                                        downloadButton("BubblePlotDownload","Save figure"),
                                        downloadButton("BubbleTableDownload","Save table"),
                                        width = 3
                                      )
                                      ),
                                      br(),
                                      br(),
                                      plotOutput("BubblePlotOut") %>% withSpinner(type = 1,color.background = "white"),
                                      style = "overflow-y:scroll; max-height: 850px; position:relative;"
                                      
                                      
                            )
                          )
                 ),
                 
                 
                 #### Bray curtis PCOA ####
                 tabPanel("Bray-Curtis PCoA Triplot",
                          sidebarLayout(
                            sidebarPanel(
                              numericInput("SrsDepth","Select your SRS depth.",min = 0, value = 1000),
                              actionButton("PStartButton",label = "Start!",width = 335),
                              hr(style = "border-width: 3px"),
                              numericInput("PEnvPThresh", "Select your p-value threshold (0 - 1)", min = 0, max = 1, value = 0.5),
                              numericInput("PEnvRThresh", "Select your R-value threshold (0 - 1)", min = 0, max = 1, value = 0.5),
                              numericInput("PTaxaThresh", "Select your taxon abundance (0 - 100)", min = 0, max = 100, value = 5),
                              hr(style = "border-width: 3px"),
                              selectInput("PFillCol","Select your fill colour",choices = "Updating"),
                              sliderInput("PSizeSelect", "Change your point size", value = 5,min = 0, max = 15),
                              checkboxInput("PGradient", "Do you want a colour gradient instead?", value = FALSE),
                              selectInput("PPalletSelect","Select a specific colour pallet", choices = c("viridis","magma","plasma","inferno","cividis","mako","rocket","turbo")),
                              hr(style = "border-width: 3px"),
                              checkboxInput("PSampleLabel", "Do you want to label your samples?", value = FALSE),
                              checkboxInput("PShapeChoice","Do you want to add a shape variable?",value = FALSE),
                              
                              # Conditional for shape selection
                              conditionalPanel(
                                condition = "input.PShapeChoice == true",
                                selectInput("PShape","Select your shape",choices = "Updating"),
                              ),
                              checkboxInput("PElips","Do you want to include statistically generated elipses?", value = FALSE),
                              width = 3,
                              style = "overflow-y:scroll; max-height: 850px; position:relative;border-color:#000000"
                            ),
                            mainPanel(width = 9,
                                      fluidRow(
                                        box(h3("Bray-Curtis PCoA triplot"),
                                            p("A Bray-Curtis PCoA triplot displays sample dissimilarity along with associated influential taxa and environmental data. Samples that group together are more similar to eachother than samples further apart in ordination space. Thresholds for taxon relative abundance (proportion of reads) and p-values for significance of environmental variables can be selected in the left panel. Arrow length is proportional to the magnitude of the effect for a given environmental parameter."),
                                            p("To generate this triplot, an ASV table is first 'rarefied' to a desired number of reads using scaling with ranked subsampling (SRS), using the", em("SRS")," R package. This normalized table is then converted to a proportion table, then into a Bray-Curtis dissimilarity matrix and PCoA matrix using the",em("vegan"),"and",em("ape")," packages. Environmental variables are fit to the PCoA coordinates using the",em("envfit")," command at 10,000 permutations, and taxa are additionally mapped using normalized weighted scores based on taxon abundance."),
                                            p("All samples below the provided SRS depth will be removed prior to PCoA and envfit analyses and will not be plotted (so be careful)"),
                                            width = 12)),
                                      dataTableOutput("PTableOut"),
                                      fluidRow(
                                        box(
                                          sliderInput("PPlotOutW","Plot width",min = 0, max = 3000,step = 100,value = 1000),
                                          sliderInput("PPlotOutH","Plot height",min = 0, max = 3000,step = 100,value = 600),
                                          downloadButton("PEnvFitTableSave","Download triplot stats table"),
                                          downloadButton("PPlotSave","Save figure"),
                                          width = 6
                                        ),
                                      ),
                                      dataTableOutput("PStatsTableOut"),
                                      br(),
                                      br(),
                                      plotOutput("PPlotOut") %>% withSpinner(type = 1,color.background = "white"),
                                      style = "overflow-y:scroll; max-height: 850px; position:relative;"
                                      
                            )
                          )
                 ),
                 
                 #---- UniFrac Plot ----
                 tabPanel("UniFrac PCoA Triplot",
                          sidebarLayout(
                            sidebarPanel(
                              fileInput("UniTree","Please upload your rooted phylogenetic tree (in Newick format)"),
                              hr(style = "border-width: 3px"),
                              numericInput("UniDepth", label = "Select your sampling depth.", value = 5000),
                              selectInput("UniDissSelect",label = "Select your type of UniFrac", choices = c("unweighted","weighted")),
                              actionButton("UniStartButton","Start!", width = 335),
                              hr(style = "border-width: 3px"),
                              numericInput("UniEnvPThresh", label = "Select your p-value threshold (0 to 1)",value = 0.5, min = 0, max = 1),
                              numericInput("UniEnvRThresh", label = "Select your R-value threshold (0 to 1)", value = 0.5, min = 0, max = 1),
                              numericInput("UniTaxaThresh",label = "Select your taxon abundance (0 - 100)", value = 5),
                              hr(style = "border-width: 3px"),
                              selectInput("UniFilCol",label = "Select your fill colour", choices = "Updating"),
                              sliderInput("UniSize", "Change your point size", value = 5,min = 0, max = 15),
                              checkboxInput("UniGradient", "Do you want a colour gradient instead?", value = FALSE),
                              selectInput("UniPalletSelect","Select a specific colour pallet. To change, reselect the checkbox above.", choices = c("viridis","magma","plasma","inferno","cividis","mako","rocket","turbo")),
                              checkboxInput("UniSampleLabel", "Do you want to label your samples?", value = FALSE),
                              checkboxInput("UniShapeSelect",label = "Do you want to differentiate by shape as well?", value = FALSE),
                              conditionalPanel(
                                condition = "input.UniShapeSelect == true",
                                selectInput("UniShape", label = "Select your shape", choice = "Updating")
                              ),
                              checkboxInput("UniElips",label = "Do you want statistically generated elipses?", value = FALSE),
                              width = 3,
                              style = "overflow-y:scroll; max-height: 850px; position:relative; border-color:#000000"
                            ),
                            mainPanel(width = 9,
                                      fluidRow(
                                        box(h3("UniFrac PCoA triplot"),
                                            p("A UniFrac PCoA triplot displays sample dissimilarity along with associated influential taxa and environmental data. The UniFrac metric includes measures of taxon phylogenetic relatedness along with taxon presence/absence (unweighted) and presence/absence and abundance (weighted). This analysis requires an additional datafile: a rooted phylogenetic tree generated from the representative sequences in your dataset. This tree is generated in QIIME2 as part of the", em("core-metrics-phylogenetics")," workflow (export your rooted tree). Samples that group together are more similar to eachother than samples further apart in ordination space. Thresholds for taxon relative abundance (proportion of reads) and p-values for significance of environmental variables can be selected in the left panel. Arrow length is proportional to the magnitude of the effect for a given environmental parameter."),
                                            p("To generate this triplot, an ASV table is first 'rarefied' to a desired number of reads using scaling with ranked subsampling (SRS), using the", em("SRS")," R package. This normalized table is then converted to a proportion table, then into a UniFrac dissimilarity using the proportion table and provided phylogenetic tree using the", em("GUniFrac"),"package. The resulting distance matrix is converted into a PCoA matrix using the",em("vegan"),"and",em("ape")," packages. Environmental variables are fit to the PCoA coordinates using the",em("envfit")," command at 10,000 permutations, and taxa are additionally mapped using normalized weighted scores based on taxon abundance."),
                                            p("All samples below the provided SRS depth will be removed prior to PCoA and envfit analyses and will not be plotted (so be careful). Please be patient: it can take upwards of 30 seconds for the plot to appear."),
                                            width = 12)),
                                      fluidRow(
                                        box(
                                          
                                          sliderInput("UniPlotOutW","Plot width",min = 0, max = 3000,step = 100,value = 1000),
                                          sliderInput("UniPlotOutH","Plot height",min = 0, max = 3000,step = 100,value = 600),
                                          downloadButton("UniPlotDownload","Save figure"),
                                          downloadButton("UniStatsFullDownload","Download full stats table"),
                                          width = 6
                                        )
                                      ),
                                      dataTableOutput("UniTableOut"),
                                      br(),
                                      br(),
                                      plotOutput("UniPlotOut") %>% withSpinner(type = 1,color.background = "white"),
                                      style = "overflow-y:scroll; max-height: 850px; position:relative;"
                                      
                            )
                          )
                 ),
                 
                 #### Next Plot ####
                 
                 # #### FAQ PAGE ####
                 tabPanel(title = "FAQ",
                          sidebarLayout(NULL,
                                        mainPanel(h1("Frequently asked questions and common problems"),
                                                  uiOutput("helpme")
                                        )
                          )
                 )
                 
                 
)

