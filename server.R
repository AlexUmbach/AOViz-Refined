# ---- Initialize the server ----
  server <- function(input, output, session) {
  options(shiny.maxRequestSize = 1000 * 1024 ^ 2)
  
# ---- Loading animations ----    
  # Provide null data to prevent swirling animations within the plot space
  output$ReadTableOut <- renderDT(NULL)
  output$ReadPlot <- renderDT(NULL)
  output$B1TableOut <- renderDT(NULL)
  output$BubbleOut <- renderDT(NULL)
  output$TaxonomyBar <- renderDT(NULL)
  output$BarTableOut <- renderDT(NULL)
  output$PcoaPlotOut <- renderDT(NULL)
  output$PcoaTableOut <- renderDT(NULL)
  output$RankedTableOut <- renderDT(NULL)
  output$RankedPlotOut <- renderDT(NULL)
  output$UniPlotOut <- renderDT(NULL)
  
  
# ---- Load the data files ----
  # This is where the main ASV table, metadata table, and additional tables can be uploaded
  
  # Upload the main ASV table
  MainDataFileUpload <- reactive({
    read.table(
      file = input$MainFile$datapath,
      fill = TRUE,
      header = TRUE, 
      sep = "\t"
    )
  })
  
  # Output to the main panel of the page
  output$MainTableOut <- renderDataTable({
    req(input$MainFile)
    Table <- read.table(
      file = input$MainFile$datapath,
      fill = TRUE,
      header = TRUE,
      sep = "\t"
    )
    head(Table, 3)
  })
  
  # Upload the main metadata table
  MetaDataFileUpload <- reactive({
    read.table(
      file = input$MetaFile$datapath,
      fill = TRUE,
      header = TRUE,
      sep = "\t"
    )
  })
  # Output to the main panel of the page
  output$MetaTableOut <- renderDataTable({
    req(input$MetaFile)
    Table <- read.table(
      file = input$MetaFile$datapath,
      fill = TRUE,
      header = TRUE,
      sep = "\t"
    )
    head(Table, 6)
  })
  
  # Upload the optional ASV list file
  ContamDataFileUpload <- reactive({
    read.table(
      file = input$ContamFile$datapath,
      fill = TRUE,
      header = TRUE,
      sep = "\t"
    )
  })
  
  # Output to the main panel of the page
  output$ContamTableOut <- renderDataTable({
    req(input$ContamFile)
    Table <- read.table(
      file = input$ContamFile$datapath,
      fill = TRUE,
      header = TRUE,
      sep = "\t"
    )
    head(Table, 6)
  })
  
# ---- Data pre processing ----
  # Now that the files are uploaded and stored, they have to be filtered to eliminate missing data in both the primary ASV table and
  # the metadata file. This means filtering samples from both data tables. 
  
  # Start by removing samples in the ASV table that are not present in the metadata tables
  MainASVTable <- reactive({
    ASVTable <- MainDataFileUpload()
    MetaTable <- MetaDataFileUpload()
    
    # Fill any missing cells with zeros
    ASVTable[is.na(ASVTable)] <- 0
    
    MetaNames <-
      c(
        MetaTable$SampleName,
        "Consensus.Lineage",
        "rowID",
        "Feature.ID",
        "ReprSequence"
      )
    ASVTable <- ASVTable[, names(ASVTable) %in% MetaNames]
    rownames(ASVTable) <- ASVTable$Feature.ID
    ASVTable
  })
  
  # Now remove metadata columns associated with samples that are not present in the ASV table
  MainMetaTable <- reactive({
    ASVTable <- MainDataFileUpload()
    MetaTable <- MetaDataFileUpload()
    
    MetaNames <-
      c(
        MetaTable$SampleName,
        "Consensus.Lineage",
        "rowID",
        "Feature.ID",
        "ReprSequence"
      )
    ASVTable <- ASVTable[, names(ASVTable) %in% MetaNames]
    MetaTable <- MetaTable %>% filter(SampleName %in% colnames(ASVTable))
    MetaTable <- replace(MetaTable, is.na(MetaTable), "NA")
    MetaTable
  })
  
  # The main contaminant table
  MainContamTable <- reactive({
    Table <- ContamDataFileUpload()
    rownames(Table) <- Table$Feature.ID
    Table
  })
  
# ---- Data processing ----
  # We now process the main ASV table to format taxonomy labels and unique identifiers necessary for downstream visualizations and processing. 
  
  # The first requirement is to format taxonomically collapsed tables, if directed to do so
  TransDataCollapsed <- reactive({
    req(input$MainFile)
    req(input$MetaFile)
    ASVTable <- MainASVTable()
    
    # If collapse table checkbox is checked, then transform. Otherwise, leave unaltered.
    if (input$IsMainCollapsed == TRUE) {
      ASVTable$Consensus.Lineage <- ASVTable$Feature.ID # Set the feature ID to imitate the taxonomy in an uncollapsed table
      ASVTable$rowID <- 1:nrow(ASVTable) # Add a column of unique numbers to imitate unique feature IDs
      ASVTable$ReprSequence <- 1:nrow(ASVTable)
      ASVTable <- ASVTable %>% mutate(ReprSequence = "Not Applicable")
    } else {
      if (input$IsMainCollapsed == FALSE) {
          ASVTable$rowID <- 1:nrow(ASVTable) # Add a column called rowID; if present it will be overwritten, so user caution is advised
        }
    }
    ASVTable
  })

  # Now we make a key file that will keep Feature ids, row ids, taxonomy, and representative sequence information 
  # properly associated with samples. This is the only location in which taxonomic identifiers, labels, or otherwise 
  # should be modified so that everything is consistent.
  FeatureDataKey <- reactive({
    
    ASVTable <- TransDataCollapsed()
    
    FeatureKey <- data.frame(FeatureID = ASVTable$Feature.ID, 
                             rowID = ASVTable$rowID, 
                             OriginalTaxonomy = ASVTable$Consensus.Lineage,
                             ReprSequence = ASVTable$ReprSequence,
                             AppendedTaxonomy = paste(ASVTable$Consensus.Lineage, "_", ASVTable$rowID))
    
    # Now we transform the the taxonomy labels in the FeatureKey, so they're formatted for visuals
    FeatureKey$Taxon <- FeatureKey$OriginalTaxonomy

    # Remove numbers and special characters from lineages, and substitutes undefined ambiguous taxa
    labels <- FeatureKey$OriginalTaxonomy
    labels <- gsub(" ", "", labels)
    labels <- gsub("_[0-9]*$", "", labels)
    labels <- gsub(" ", "", labels)
    labels <- gsub("(;Ambiguous__taxa)", ";s__Ambiguous_taxa", labels)
    labels <- gsub("(;Ambiguous_taxa)", ";s__Ambiguous_taxa", labels)
    
    # Truncate taxonomic lineages for readability (default is "Yes")
    # if (input$TruncateTaxa == "Yes") {
    labels <- paste(";", sep = "", labels)
    labels <- gsub("(;\\s*Ambiguous_taxa)", "", labels)
    labels <- gsub("(uncultured.*)", "", labels)
    labels <- gsub("(__uncultured.*)", "", labels)
    labels <- gsub("(unidenti.*)", "", labels)
    labels <- gsub("(__unidenti.*)", "", labels)
    labels <- gsub("(;.__Ambiguous_taxa)", "", labels)
    labels <- gsub("(;._Ambiguous_taxa)", "", labels)
    labels <- gsub("(;s__$)", "", labels)
    labels <- gsub("(;g__$)", "", labels)
    # }
    
    # Remove the prefixes associated with SILVA classifiers (e.g., D_*__ or p_, f_)
    labels <- gsub("(D_.__)", "", labels)
    labels <- gsub(";$", "", labels)
    labels <- gsub("(D_.__$)", "", labels)
    
    # Option to include this exclude this step; leaving here incase someone requests it back
    # if (input$remove_prefix == "Yes") {
    #   labels <- gsub("(D_.__)", "", labels)
    #   labels <- gsub(";$", "", labels)
    # }
    # labels <- gsub("(D_.__$)", "", labels)
    
    labels <- gsub("(d__)", "", labels)
    labels <- gsub("(p__)", "", labels)
    labels <- gsub("(c__)", "", labels)
    labels <- gsub("(o__)", "", labels)
    labels <- gsub("(f__)", "", labels)
    labels <- gsub("(g__)", "", labels)
    labels <- gsub("(s__)", "", labels)
    labels <- gsub(" ", "", labels)
    labels <- gsub("(;metagenome$)", "", labels)
    labels <- gsub("(__.$)", "", labels)
    labels <- gsub("(;__)", "", labels)
    labels <- gsub("(;$)", "", labels)
    labels <- gsub("(;$)", "", labels) # This is not a duplicate. Leave here.
    labels <- gsub("^;","", labels)
    
    # Retrieves the terminal taxon string
    FeatureKey$Labels <- labels
    FeatureKey$Taxon <- paste(gsub(".*;", "", FeatureKey$Labels), FeatureKey$rowID, sep = "_")
    
    FeatureKey
  })
  
  # Generate a table containing only count data, using the FeatureIDs as rownames to track data
  # Also remove low-abundance sequences, if so desired
  FilteredTable <- reactive({
    ASVTable <- MainASVTable()
    
    ColsToFilter <- c("OTU.ID",
                      "Feature.ID",
                      "rowID",
                      "Consensus.Lineage",
                      "ReprSequence"
                      )
    ASVTable <- ASVTable[,!(names(ASVTable) %in% ColsToFilter)]
    
    # Set all reads below a threshold as 0; any taxa with zero reads will be removed downstream
    if (input$RemoveLowReads == TRUE) {
      ASVTable[ASVTable < input$ReadThreshold] <- 0
      ASVTable
    }
    
    # Set contaminant reads to zero
    if (input$ContamChoice == "Remove"){
      Contaminants <- MainContamTable()
      Matches <- intersect(rownames(ASVTable),rownames(Contaminants))
      ASVTable[Matches, !(colnames(ASVTable) %in% c("OTU.ID",
                                                    "Feature.ID",
                                                    "rowID",
                                                    "Consensus.Lineage",
                                                    "ReprSequence"))] <- 0 
      }
    
    ASVTable
    
  })
  
  # Show this table in the main panel
  output$FinalProcessedTable <- renderDataTable({
    Table <- FilteredTable()
    output$proc_maintext <- renderText("This is your processed data")
    Table
  })
  
  
# ---- Read Plot Visualization ----
  # This is our first plot: the read plot. It shows sequence depth for each sample or groups of samples
  
  # First we update 
  observeEvent(input$MetaFile,{
    req(input$MetaFile)
    MetaData <- MainMetaTable()
    MetaColNames <- colnames(MetaData)
    MetaColNames <- MetaColNames[MetaColNames != "SampleName"]
    updateSelectInput(session, "ReadSortByAxis", choices = sort(MetaColNames))
    updateSelectInput(session, "ReadMetaGroup", choices = sort(MetaColNames))
    updateSliderInput(session, "ReadPlotOutW")
    updateSliderInput(session, "ReadPlotOutH")
    # })
  })
  
  # Now we transform the raw ASV counts into totals
  ReadCountsTotal <- reactive({
    
    req(input$ReadStartButton)
    ASVTable <- FilteredTable()
    MetaData <- MainMetaTable()

    
    # If analyzing contaminants, match between ASV table and ASV list, and keep only those that match
    if (input$ContamChoice == "Analyze"){
      Contams <- MainContamTable()
      Matches <- intersect(rownames(ASVTable),rownames(Contams))
      ASVTable <- ASVTable[Matches, ]
      
      ASVTable
    }
    
    # Sum each column and add a total, then select only that row
    ASVTable <- rbind(ASVTable, Total = colSums(ASVTable))
    ReadTotal <- ASVTable["Total", ]
    
    # Now convert to long form
    ReadTotal <- reshape2::melt(ReadTotal, 
                                      variable.name = as.character("SampleName"), 
                                      value.name = "Total")
    
    # Add metadata
    ReadTotal <-left_join(ReadTotal,
                          MetaData,
                          by = "SampleName",
                          copy = FALSE)
    
    # Not sure this if this still required
    # ## Reorder x-axis to follow metadata category in data frame
    # data_read_total$SampleName <-
    #   as.character(data_read_total$SampleName)
    # data_read_total$SampleName <-
    #   factor(data_read_total$SampleName,
    #          levels = unique(data_read_total$SampleName))
    # 
    # ## Make a list of unique metadata
    # read_meta_list <- unique(data_read_total$input$read_meta_group)
    # 
    # read_meta_list
    ReadTotal
    
    })
  
  # Generate the read plot
  ReadPlotVisual <- reactive({
    
    req(input$ReadStartButton)
    ReadCounts <- ReadCountsTotal()
    MetaData <- MainMetaTable()
    
    # Determine whether it will be a bar or boxplot
    if (isolate(input$BoxSelect) == "Bar") {
      ReadPlot <- ggplot(data = ReadCounts, aes(x = SampleName,
                                                      y = Total,
                                                      width = isolate(input$ReadWidth)
                                                      ))
      
      ReadPlot <- ReadPlot + geom_bar(fill = "#D3D3D3",
                                      colour = "black",
                                      size = 0.4,
                                      alpha = 0.8,
                                      stat = "identity",
                                      position = "stack"
      )
      
      ReadPlot
      
    }
    
    if (isolate(input$BoxSelect) == "Box") {
      ReadPlot <- ggplot(data = ReadCounts, 
                         aes(x = as.factor(isolate(input$ReadSortByAxis)),
                             y = Total,
                             width = 2
                          ))
      
      ReadPlot <- ReadPlot + geom_boxplot(fill = "#D3D3D3",
                                          position = "dodge2",
                                          alpha = 0.8,
                                          # size = input$read_width,
                                          size = 0.4,
                                          outlier.shape = NA
      )
      
      ReadPlot <- ReadPlot + stat_summary(fun = mean,
                                          geom = "point",
                                          shape = 15,
                                          size = 1,
                                          colour = "red"
      )
      
      ReadPlot <- ReadPlot + geom_jitter()
      
      ReadPlot
    }
    
    # setting the graph so that it begins at the x-axis and there is no gap. Also sets the limits of the y-axis.
    ReadPlot <- ReadPlot + scale_y_continuous(
      expand = c(0, 0),
      limit = (c(0, isolate(
        input$ReadYAxisLimit)
        )
      )
    )
        
    
    # Add faceting for sorting
    ReadPlot <- ReadPlot +
      facet_grid( ~ eval(
        parse(
          text = isolate(
            input$ReadSortByAxis)
                          )
                        ),
        space = "free",
        scales = "free",
        switch = "both"
      )
    
    # Now set options for borders
    if (isolate(input$ReadPanel) == "Yes") {
      ReadPlot <- ReadPlot + theme_bw() +
        theme(
          panel.grid = element_blank(),
          text = element_text(colour = "black"),
          #axis.line = element_line(colour = "black"),
          axis.line = element_blank(),
          axis.text = element_text(colour = "black", size = 12),
          axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size = 12
          ),
          legend.text = element_text(face = "italic", size = 12),
          legend.title = element_text(size = 14),
          panel.spacing = isolate(unit(as.numeric(input$ReadPanelSpacing), "points")),
          #legend.position = "none",
          axis.title = element_text(size = 14, face = NULL),
          axis.text.y = element_text(size = 16),
          strip.text.x = element_text(size = 10, face = "bold"),
          strip.background = element_rect(fill = "white"),
        )
    }
    
    if (isolate(input$ReadPanel) == "No") {
      ReadPlot <- ReadPlot + theme_bw() +
        theme(
          panel.grid = element_blank(),
          text = element_text(colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 12),
          axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size = 12
          ),
          legend.text = element_text(face = "italic", size = 12),
          legend.title = element_text(size = 14),
          panel.border = element_blank(),
          panel.spacing = isolate(unit(as.numeric(input$ReadPanelSpacing), "points")),
          #legend.position = "none",
          axis.title = element_text(size = 14, face = NULL),
          axis.text.y = element_text(size = 16),
          strip.background = element_rect(fill = "white"),
          strip.text.x = element_text(size = 10, face = "bold")
        )
    }
    
    ReadPlot <- ReadPlot + labs(fill = "Sample category")
    ReadPlot <- ReadPlot + xlab("Samples")
    ReadPlot <- ReadPlot + ylab("Total reads following DADA2")
    
    ReadPlot
    
  }) # End of ReadPlotVisual
  
  # Output to the Read Plot panel
  output$ReadTableOut <- renderDataTable({
    req(input$ReadStartButton)
    ReadTable <- ReadCountsTotal()
    ReadTable
  })
  
  ## You must define the input for width/height within a reactive context, then call it in the output
  ReadPlotHeight = reactive(input$ReadPlotOutH)
  ReadPlotWidth = reactive(input$ReadPlotOutW)
  
  output$ReadPlotOut = renderPlot({
    ReadPlot = ReadPlotVisual()
    ReadPlot
    
  },
  width = ReadPlotWidth,
  height = ReadPlotHeight)
  
  # Download output for read plot
  output$ReadDownload = downloadHandler(
    filename = "read_plot.pdf",
    contentType = ".pdf",
    content = function(ReadFile) {
      ggsave(
        ReadFile,
        plot = ReadPlotVisual(),
        device = "pdf",
        height = as.numeric(input$ReadPlotOutH),
        width = as.numeric(input$ReadPlotOutW),
        units = "px",
        scale = 4
      )
    }
  )
  
  # Download the table for the read plot
  output$ReadTableDownload <- downloadHandler(
    filename = "read_table.csv",
    content = function(save) {
      write.csv(ReadCountsTotal(), save)
    }) # End of read plot
  
  
# ---- Taxonomic Bar Plot Visualization ----
  # This plot shows the distribution of taxa within a sample, using a cut-off wherein any taxa below said cut-off are relegated to an "other" category
  
  # Add events to monitor and change the dropdown menus and inputs
  observeEvent(input$MetaFile,{
    req(input$MetaFile)
    MetaData <- MainMetaTable()
    MetaColNames <- colnames(MetaData)
    MetaColNames <- MetaColNames[MetaColNames != "SampleName"]
    updateSelectInput(session, "BarSortByAxis", choices = sort(MetaColNames))
    updateSelectInput(session, "BarTaxonLevel")
    updateTextInput(session, "BarPlotHeight")
    updateTextInput(session, "BarPlotWidth")
    updateSelectInput(session, "BarSecondFacet", choices = sort(MetaColNames))
  })
  
  # Generate the dataframe necessary for plotting
  BarPlotDataTran <- reactive({
    
    req(input$BarStartButton)
    ASVTable <- FilteredTable()
    FeatureKey <- FeatureDataKey()
    
    # Produce a proportion table 
    BarProp <- prop.table(as.matrix(ASVTable), 2) * 100
    BarProp <- data.frame(BarProp)
    
    # Remove columns with zero total counts, then filter the table using a given user-defined cut-off. If any sample contains 
    # taxa with a proportion greater than a cut-off, that taxon is retained in all samples
    BarProp <- BarProp[, colSums(is.na(BarProp)) == 0]
    BarProp <- BarProp %>% filter_all(any_vars(. >= as.numeric(isolate(input$BarCutOff))))
    
    # Reassign taxonomy to the filtered table using the rownames and FeatureDataKey reactive file
    BarProp$FeatureID <- rownames(BarProp)
    BarProp <- left_join(BarProp, FeatureKey, by = "FeatureID")
    rownames(BarProp) <- BarProp$FeatureID
    
    # If analyzing contaminants, match between ASV table and ASV list, and keep only those that match
    if (input$ContamChoice == "Analyze"){
      Contams <- MainContamTable()
      Matches <- intersect(rownames(BarProp),rownames(Contams))
      BarProp <- BarProp[Matches, ]
      BarProp
  
    }
    
    BarProp
  })
  
  # Generate a long form dataframe; keep for output
  BarPlotDataLong <- reactive({

    req(input$BarStartButton)
    BarTable <- BarPlotDataTran()
    MetaData <- MainMetaTable()
    FeatureKey <- FeatureDataKey()
    
    BarTable <- reshape2::melt(BarTable,
                               id.vars = c("Taxon",
                                           "OriginalTaxonomy",
                                           "rowID",
                                           "FeatureID",
                                           "AppendedTaxonomy",
                                           "ReprSequence",
                                           "Labels"),
                               variable.name = as.character("SampleName"),
                               value.name = "Percentage"
    )
    
    # Remove zero-abundance taxa
    BarTable <- dplyr::filter(BarTable, Percentage > 0)
    
    # Separate the full lineage into individual taxonomic levels
    BarTable <- separate(BarTable,
                         Labels,
                         c("Domain",
                           "Phylum",
                           "Class",
                           "Order",
                           "Family",
                           "Genus",
                           "Species",
                           "Sub1",
                           "Sub2",
                           "Sub3",
                           "Sub4"),
                         sep = ";",
                         remove = TRUE,
                         convert = FALSE
                         )
    
    # Replace any NA with unresolved
    BarTable <- BarTable %>%
      mutate(across(c(Domain,
                      Phylum,
                      Class,
                      Order,
                      Family,
                      Genus,
                      Species,
                      Sub1,
                      Sub2,
                      Sub3,
                      Sub4),
                    ~replace_na(.x, "ZZZ_No_Taxon_Info")))
    
    
    
    # Add an identifier so the data can be mapped afterwards
    BarTable$Identifier <- 1:nrow(BarTable)
  
    BarTable
  })
  
  # Generate a summarized dataframe for plotting; includes setting filtered taxa as "other"
  BarPlotDataSum <- reactive({
    req(input$BarStartButton)
    BarTable <- BarPlotDataLong()
    MetaData <- MainMetaTable()
    
    # Now sort and combine proportions for samples with matching taxonomies
    BarTableSum <- BarTable %>% 
      group_by(SampleName, !!sym(isolate(input$BarTaxonLevel))) %>% 
      summarize(Percentage = sum(Percentage)) %>%
      rename(TaxonLevel = !!sym(input$BarTaxonLevel))
    
    # Filter samples to "other" category based on a cut-off
    IncludedSamples <- aggregate(BarTableSum$Percentage,
                                 by = list(BarTableSum$SampleName),
                                 FUN = sum
    )
    
    IncludedSamples <- cbind(IncludedSamples,
                             100 - IncludedSamples$x)
    colnames(IncludedSamples) <- c("SampleName", "FiltSum", "Percentage")
    IncludedSamples$TaxonLevel <- "ZZOther"
    
    IncludedSamples <- select(IncludedSamples,
                              -c("FiltSum"))
    BarTableSumFinal <- bind_rows(BarTableSum, IncludedSamples)
    
    # Add metadata. Because this collapses features, all additional feature data and full lineages are lost
    BarTableSumFinal <- left_join(BarTableSumFinal, MetaData, by = "SampleName")

    
    BarTableSumFinal
  })
  
  BarPlotVisual <- reactive({
    
    req(input$BarStartButton)
    BarTable <- BarPlotDataSum()
    MetaData <- MainMetaTable()
    
    BarPlot <- ggplot(data = BarTable,
                       aes(x = SampleName,
                           y = Percentage,
                           width = 1,
                           fill = TaxonLevel
                       ))
    
    BarPlot <- BarPlot + geom_bar(colour = "black",
                                  size = 0.5,
                                  alpha = isolate(input$BarAlpha),
                                  stat = "identity",
                                  position = "stack"
                                  )
    
    if (input$BarRenameCheck == TRUE){
      BarPlot <- BarPlot + scale_x_discrete(breaks = BarTable$SampleName,
                                            labels = BarTable$SampleShort)
    }
    
    ## Add faceting for sorting
    if (isolate(input$BarFacet) == FALSE){
      BarPlot <- BarPlot +
        facet_grid(
          ~ eval(parse(text = isolate(input$BarSortByAxis))),
          space = "free",
          scales = "free",
          switch = "both"
        )
    }
    
    else if (isolate(input$BarFacet) == TRUE){
      BarPlot <- BarPlot +
        facet_nested(~ eval(parse(text = isolate(input$BarSortByAxis))) + eval(parse(text = isolate(input$BarSecondFacet))),
                     space = "free",
                     scales = "free",
                     switch = "both"
        )
    }
    
    if (isolate(input$BarPanelBorder == "Yes")) {
      BarPlot <- BarPlot + theme_bw() +
        theme(panel.grid = element_blank(),
              text = element_text(colour = "black"),
              # axis.line = element_line(colour = "black"),
              axis.line = element_blank(),
              axis.text = element_text(colour = "black", size = 12),
              axis.text.x = element_text(angle = 90,
                                         hjust = 1,
                                         vjust = 0.5,
                                         size = 12),
              legend.text = element_text(face = "italic", size = 12),
              legend.title = element_text(size = 14),
              panel.spacing = unit(as.numeric(isolate(input$BarPanelSpacing)),
                                   "points"),
              #legend.position = "none",
              axis.title = element_text(size = 10, face = NULL),
              axis.text.y = element_text(size = 16),
              strip.background = element_rect(fill = "white"),
              strip.text.x = element_text(size = 10, face = "bold"),
        )
      }
    
    if (isolate(input$BarPanelBorder) == "No") {
      BarPlot <- BarPlot + theme_bw() +
        theme(panel.grid = element_blank(),
              text = element_text(colour = "black"),
              axis.line = element_line(colour = "black"),
              axis.text = element_text(colour = "black",
                                       size = 12),
              axis.text.x = element_text(angle = 90,
                                         hjust = 1,
                                         vjust = 0.5,
                                         size = 12
                                         ),
              legend.text = element_text(face = "italic", size = 12),
              legend.title = element_text(size = 14),
              panel.border = element_blank(),
              panel.spacing = unit(as.numeric(isolate(input$BarPanelSpacing)),
                                   "points"),
              #legend.position = "none",
              axis.title = element_text(size = 10, face = NULL),
              axis.text.y = element_text(size = 16),
              strip.background = element_rect(fill = "white"),
              strip.text.x = element_text(size = 10, face = "bold"),
        )
      }
    
    # setting the graph so that it begins at the x-axis and there is no gap. Also sets the limits of the y-axis.
    BarPlot <- BarPlot + scale_y_continuous(expand = c(0, 0))
    
    ColourNumber <- length(unique(BarTable[["TaxonLevel"]]))
    
    MyColours <- colorRampPalette(brewer.pal(8, "Paired"))(ColourNumber)
    ## colour fill using Brewers
    BarPlot <- BarPlot + scale_fill_manual(values = MyColours)
    #bar_plot <- bar_plot + scale_fill_brewer(palette = "Paired")
    
    BarPlot = BarPlot + ylab(
      paste0(
        "Proportion of affiliated taxa at >",
        isolate(input$BarCutOff),
        "% relative abundance"
      )
    )
      
    BarPlot <- BarPlot + xlab("Samples")
    BarPlot <- BarPlot + labs(fill = "Taxonomy")
    
    #bar_plot<-bar_plot + ylab(input$bar_yaxis)
    BarPlot
  })
  
  output$BarTableOut <- renderDataTable({
    BarTable <- BarPlotDataLong()
    BarTable
  })
  
  ## You must define the input for width/height within a reactive context, then call it in the output.
  TaxaPlotWidth <- reactive(input$TaxaPlotOutW)
  TaxaPlotHeight <- reactive(input$TaxaPlotOutH)
  
  output$TaxaPlotOut <- renderPlot({
    BarPlot <- BarPlotVisual()
    BarPlot
  },
  width = TaxaPlotWidth,
  height = TaxaPlotHeight)
  
  # Download options for the plot
  output$BarDownload <- downloadHandler(
    filename = "bar_plot.pdf",
    contentType = ".pdf",
    content = function(BarFile) {
      ggsave(
        BarFile,
        plot = BarPlotVisual(),
        device = "pdf",
        height = as.numeric((input$TaxaPlotOutH)),
        width = as.numeric((input$TaxaPlotOutW)),
        units = "px",
        scale = 4
      )
    }
  )
  
  # Download the table for the bar plot
  output$BarTableDownload <- downloadHandler(
    filename = "bar_table.csv",
    content = function(BarTable) {
      write.csv(BarPlotDataLong(), BarTable)
    }
  )
  
# ---- Bubble Plot Visualization ----
  # Formatting and visualizations for the bubble plot
  
  ## Update selection - Bubble plot
  observeEvent(input$MetaFile,{
    req(input$MetaFile)
    # Incorporate the meta file column names into every appropriate input
    MetaData <- MainMetaTable()
    MetaColNames <- colnames(MetaData)
    MetaColNames <- c(colnames(MetaData), "Taxon")
    MetaColNames <- MetaColNames[MetaColNames != "SampleName"]
    updateSelectInput(session, "BubbleSortXAxis", choices = sort(MetaColNames))
    updateSelectInput(session, "BubbleColour", choices = sort(MetaColNames))
    updateSelectInput(session, "BubbleMetaFilt", choices = sort(MetaColNames))
    updateSelectInput(session, "BubbleFacet", choices = sort(MetaColNames))
    updateTextInput(session, "BubbleWidth")
    updateTextInput(session, "BubbleHeight")
    updateTabItems(session, "BubbleMetaKeywords")
    updateSelectInput(session, "BubbleSecondFacetMeta", choices = sort(MetaColNames))
    updateSelectInput(session, "BubbleThirdFacetMeta", choices = sort(MetaColNames))
    updateSelectInput(session, "BubbleSortYAxis", choices = sort(MetaColNames))
    updateSelectInput(session, "BubbleRenameXAxis", choices = sort(MetaColNames))
  })
  
  # We must regenerate the proportion tables and associated data
  BubbleDataTran <- reactive({
    req(input$BubbleStartButton)
    BubbleTable <- FilteredTable()
    FeatureKey <- FeatureDataKey()
    
    # Create a proportion table
    BubbleProp <- prop.table(as.matrix(BubbleTable),2) * 100
    BubbleProp <- as.data.frame(BubbleProp)
    
    # Remove zeros and filter based on a user-defined cut-off
    BubbleProp <- BubbleProp[, colSums(is.na(BubbleProp)) == 0]
    BubblePropFilt <- BubbleProp %>% filter_all(any_vars(. >= as.numeric(isolate(input$BubbleAbundThresh))))
    
    # Reinsert FeatureKey information
    BubblePropFilt$FeatureID <- rownames(BubblePropFilt)
    BubblePropFilt <- left_join(BubblePropFilt, FeatureKey, by = "FeatureID")
    rownames(BubblePropFilt) <- BubblePropFilt$FeatureID
    
    # If analyzing contaminants, match between ASV table and ASV list, and keep only those that match
    if (input$ContamChoice == "Analyze"){
      Contams <- MainContamTable()
      Matches <- intersect(rownames(BubblePropFilt),rownames(Contams))
      BubblePropFilt <- BubblePropFilt[Matches, ]
      
      BubblePropFilt
      
    }

    BubblePropFilt
  })
  
  # Create a long form dataframe
  BubbleDataLong <- reactive({
    
    req(input$BubbleStartButton)
    BubbleTable <- BubbleDataTran()
    FeatureKey <- FeatureDataKey()
    MetaData <- MainMetaTable()
    
    # Create a long format
    BubbleTable <- reshape2::melt(BubbleTable,
                                  id.vars = c("Taxon",
                                              "OriginalTaxonomy",
                                              "rowID",
                                              "FeatureID",
                                              "AppendedTaxonomy",
                                              "ReprSequence",
                                              "Labels"),
                                  variable.name = as.character("SampleName"),
                                  value.name = "Percentage"
                                  )
    
    # Remove samples containing zero counts. Check this later
    BubbleTable <- filter(BubbleTable, Percentage > 0)
    
    # If desired, add a fake taxon (post abundance calculations) so all samples are plotted regardless of presence. 
    # Useful for image formatting and later editing (missing samples can be annoying)
    ## Add in a fake taxonomy to show all samples
    if (input$BubbleFakeTaxon == TRUE){
        FakeTaxonDf <- data.frame(SampleName = unique(BubbleTable$SampleName),
                                  Percentage = 1, 
                                  Taxon = "AAA_FAKE_TAXON_999999", 
                                  Labels ="FakeTaxon", 
                                  ReprSequence = "AAAAA",
                                  rowID = 999999,
                                  OriginalTaxonomy = "d__FakeTaxon",
                                  AppendedTaxonomy = "d__FakeTaxon_999999",
                                  FeatureID = "AAAAAAAAAAAAA"
        )
        
        # Now add the fake taxon to the table
        BubbleTable <- rbind(BubbleTable, FakeTaxonDf)
    }
    
    # Separate the full lineage into individual taxonomic levels
    BubbleTable <- separate(BubbleTable,
                         Labels,
                         c("Domain",
                           "Phylum",
                           "Class",
                           "Order",
                           "Family",
                           "Genus",
                           "Species",
                           "Sub1",
                           "Sub2",
                           "Sub3",
                           "Sub4"),
                         sep = ";",
                         remove = TRUE,
                         convert = FALSE
    )
    
    BubbleTable <- BubbleTable %>%
      mutate(across(c(Domain,
                      Phylum,
                      Class,
                      Order,
                      Family,
                      Genus,
                      Species,
                      Sub1,
                      Sub2,
                      Sub3,
                      Sub4),
                    ~replace_na(.x, "ZZZNoTaxonInfo")))

    
    # Now adjust percentages decimals
    BubbleTable$Percentage <- round(BubbleTable$Percentage,
                                    digits = as.numeric(isolate(input$BubbleDec))
                                    )
    
    # Add the metadata
    BubbleTable <- left_join(BubbleTable, MetaData, by = "SampleName")
    
    # Below are functions to remove or include specific taxa. These are performed after abundance calculations
    # so the subsequent plot abundances are not affected by their removal or isolation
    
    # Filtering to show specific taxa (removing all others). Splits the keywords into lists using "," separator
    if (!is.na(isolate(input$BubbleTaxKeywords))){
      # Split the input$b1_tax_keyword into a list of keywords using ","
      TaxKeywords <- strsplit(isolate(input$BubbleTaxKeywords), ",")[[1]]
      if (length(TaxKeywords) > 0){
        # Create a pattern by pasting the keywords with "|" for OR condition
        pattern <- paste(TaxKeywords, collapse = "|")
        
        # Perform pattern matching using grepl on the full lineage
        TaxonHits <- grepl(
          pattern = pattern,
          ignore.case = TRUE,
          x = BubbleTable$OriginalTaxonomy
        )
        
        # Subset the data_long_bubble dataset based on TaxonHits
        BubbleTable <- BubbleTable[TaxonHits, ]
        
        # Display a warning message indicating taxa filtering is selected
        warning("Taxa filtering selected.")
      }
    }
    
    # Removing specific taxa from the list
    if (!is.na(isolate(input$BubbleTaxRemove))){
      # Split the input$b1_tax_keyword into a list of keywords using ","
      TaxRemKeywords <- strsplit(isolate(input$BubbleTaxRemove), ",")[[1]]
      if (length(TaxRemKeywords) > 0){
        
        # Create a pattern by pasting the keywords with "|" for OR condition
        pattern <- paste(TaxRemKeywords, collapse = "|")
        
        # Perform pattern matching using grepl
        TaxonHits <- grepl(
          pattern = pattern,
          ignore.case = TRUE,
          x = BubbleTable$OriginalTaxonomy
        )
        
        # Subset the data_long_bubble dataset based on TaxonHits
        BubbleTable <- BubbleTable[!TaxonHits, ]
        
        # Display a warning message indicating taxa filtering is selected
        warning("Taxa filtering selected.")
      }
    }
    
    # Filter samples based on metadata categories
    if (is.na(isolate(input$BubbleMetaFilt))) {
      warning("No metadata category selected. Script will continue without filtering by groups.")
    } else {
      # Split the isolate(input$BubbleMetaKeywords) into a list of keywords using ","
      MetaKeywords <- strsplit(isolate(input$BubbleMetaKeywords), ",")[[1]]
      
      # Create a pattern by pasting the keywords with "|" for OR condition
      Pattern <- paste(MetaKeywords, collapse = "|")
      
      # Use any() to check if any keyword matches in the metadata group
      SampleHits <- sapply(BubbleTable[, isolate(input$BubbleMetaFilt)], function(GroupValue) {
        any(grepl(Pattern, GroupValue, ignore.case = TRUE))
      })
      
      # Subset the data_long_bubble dataset based on sample_hits
      BubbleTable <- BubbleTable[SampleHits, ]
      
      # Display a warning message indicating metadata filtering is selected
      warning("Metadata filtering selected.")
    }
    
    
    # Set taxonomy as factors so that the plot colours correctly, otherwise it isn't a gradient and repeats
    BubbleTable <-
      BubbleTable[with(BubbleTable,
                       order(eval(parse(
                              text = isolate(input$BubbleTaxSort)
                            )), Taxon, decreasing = TRUE)), ]
    BubbleTable$Taxon <- as.character(BubbleTable$Taxon)
    BubbleTable$Taxon <- factor(BubbleTable$Taxon,
                                levels = unique(BubbleTable$Taxon))
    BubbleTable
    
    BubbleTable
    
  })
  
  # Now create the bubbleplot visual
  BubblePlotVisual <- reactive({
    
    req(input$BubbleStartButton)
    BubbleTable <- BubbleDataLong()
    FeatureKey <- FeatureDataKey()
    MetaData <- MainMetaTable()
    FilteredTable <- FilteredTable()
    
    # The primary bubble plot
    BubblePlot <- ggplot(BubbleTable,
                         aes(
                           x = reorder(SampleName,
                                       isolate(!!sym(input$BubbleSortXAxis))),
                           y = Taxon,
                           fill = as.factor(isolate(!!sym(input$BubbleColour))),
                           colour = isolate(!!sym(input$BubbleColour)),
                           size = Percentage,
                           )
                         ) +
      guides(fill = "none", colour = "none") +
      labs(size = "Relative abundance")
    
    
    # Rename using SampleShort column
    if (input$BubbleRename == TRUE){
      BubblePlot <- BubblePlot + scale_x_discrete(breaks = BubbleTable$SampleName, labels = BubbleTable$SampleShort)
    }
    
    # Add faceting to the plot based on chosen metadata categories. This will be revisited later, but right now it works.
    
    # Adding second faceting
    if (isolate(input$BubbleSecondFacet) == TRUE) {
      BubblePlot <- BubblePlot +
        
        # First must check the facet options. Not sure how else to do this without running through all possible cases:
        if ((input$BubbleFacetSideX == "Top") &
            (input$BubbleFacetSideY == "Right")) {
          facet_nested(
            eval(parse(text = input$BubbleTaxSort)) ~ eval(parse(text = isolate(input$BubbleFacet))) +
              eval(parse(text = isolate(input$BubbleSecondFacetMeta))),
            space = "free",
            scales = "free"
          )}
      
      else if ((input$BubbleFacetSideX == "Top") &
               (input$BubbleFacetSideY == "Left")) {
        facet_nested(
          eval(parse(text = input$BubbleTaxSort)) ~ eval(parse(text = isolate(input$BubbleFacet))) +
            eval(parse(text = isolate(input$BubbleSecondFacetMeta))),
          space = "free",
          scales = "free",
          switch = "y"
        )}
      
      else if ((input$BubbleFacetSideX == "Bottom") &
               (input$BubbleFacetSideY == "Left")) {
        facet_nested(
          eval(parse(text = input$BubbleTaxSort)) ~ eval(parse(text = isolate(input$BubbleFacet))) +
            eval(parse(text = isolate(input$BubbleSecondFacetMeta))),
          space = "free",
          scales = "free",
          switch = "both"
        )}
      
      else if ((input$BubbleFacetSideX == "Bottom") &
               (input$BubbleFacetSideY == "Right")) {
        facet_nested(
          eval(parse(text = input$BubbleTaxSort)) ~ eval(parse(text = isolate(input$BubbleFacet))) +
            eval(parse(text = isolate(input$BubbleSecondFacetMeta))),
          space = "free",
          scales = "free",
          switch = "x"
        )}
    }
    
    # If no additional faceting is chosen:
    if (isolate(input$BubbleSecondFacet) == FALSE) {
      BubblePlot <- BubblePlot +
        
        # First must check the facet options. Not sure how else to do this without running through all possible cases:
        
        if ((input$BubbleFacetSideX == "Top") &
            (input$BubbleFacetSideY == "Right")) {
          facet_nested(
            eval(parse(text = input$BubbleTaxSort)) ~ eval(parse(text = isolate(input$BubbleFacet))),
            space = "free",
            scales = "free"
          )}
      
      else if ((input$BubbleFacetSideX == "Top") &
               (input$BubbleFacetSideY == "Left")) {
        facet_nested(
          eval(parse(text = input$BubbleTaxSort)) ~ eval(parse(text = isolate(input$BubbleFacet))),
          space = "free",
          scales = "free",
          switch = "y"
        )}
      
      else if ((input$BubbleFacetSideX == "Bottom") &
               (input$BubbleFacetSideY == "Left")) {
        facet_nested(
          eval(parse(text = input$BubbleTaxSort)) ~ eval(parse(text = isolate(input$BubbleFacet))),
          space = "free",
          scales = "free",
          switch = "both"
        )}
      
      else if ((input$BubbleFacetSideX == "Bottom") &
               (input$BubbleFacetSideY == "Right")) {
        facet_nested(
          eval(parse(text = input$BubbleTaxSort)) ~ eval(parse(text = isolate(input$BubbleFacet))),
          space = "free",
          scales = "free",
          switch = "x"
        )}
    }
    

    # Third level faceting. This code must come after all others. Not sure why. But if it comes first it doesn't remain dynamic
    if (isolate(input$BubbleThirdFacet) == TRUE) {
      BubblePlot <- BubblePlot +
        
        # First must check the facet options. Not sure how else to do this without running through all possible cases:
        if ((input$BubbleFacetSideX == "Top") &
            (input$BubbleFacetSideY == "Right")) {
          facet_nested(
            eval(parse(text = isolate(input$BubbleTaxSort))) ~ eval(parse(text = isolate(input$BubbleFacet))) *
              eval(parse(text = isolate(input$BubbleSecondFacetMeta))) * 
              eval(parse(text = isolate(input$BubbleThirdFacetMeta))),
            space = "free",
            scales = "free"
          )}
      
      else if ((input$BubbleFacetSideX == "Top") &
               (input$BubbleFacetSideY == "Left")) {
        facet_nested(
          eval(parse(text = isolate(input$BubbleTaxSort))) ~ eval(parse(text = isolate(input$BubbleFacet))) *
            eval(parse(text = isolate(input$BubbleSecondFacetMeta))) * 
            eval(parse(text = isolate(input$BubbleThirdFacetMeta))),
          space = "free",
          scales = "free",
          switch = "y"
        )}
      
      else if ((input$BubbleFacetSideX == "Bottom") &
               (input$BubbleFacetSideY == "Left")) {
        facet_nested(
          eval(parse(text = isolate(input$BubbleTaxSort))) ~ eval(parse(text = isolate(input$BubbleFacet))) *
            eval(parse(text = isolate(input$BubbleSecondFacetMeta))) * 
            eval(parse(text = isolate(input$BubbleThirdFacetMeta))),
          space = "free",
          scales = "free",
          switch = "both"
        )}
      else if ((input$BubbleFacetSideX == "Bottom") &
               (input$BubbleFacetSideY == "Right")) {
        facet_nested(
          eval(parse(text = isolate(input$BubbleTaxSort))) ~ eval(parse(text = isolate(input$BubbleFacet))) *
            eval(parse(text = isolate(input$BubbleSecondFacetMeta))) * 
            eval(parse(text = isolate(input$BubbleThirdFacetMeta))),
          space = "free",
          scales = "free",
          switch = "x"
        )}
    }
    
    
    # Add the bubbles and percentage labels to the plot:
    
    if (input$BubbleInclPercent == "Yes") {
      BubblePlot <- BubblePlot +
        geom_point(shape = 21,
                   alpha = 0.8) +
        geom_text(aes(label = Percentage),
                  colour = "black",
                  size = 3.0) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_size_area(max_size = 15) +
        ggtitle("") + 
        xlab("") +
        scale_size_area(max_size = 15) +
        ylab(paste0(
          "Proportion of affiliated taxa at >",
          isolate(input$BubbleAbundThresh),
          "% relative abundance"
          )
        )
      } else if (input$BubbleInclPercent == "No") {
        BubblePlot <- BubblePlot +
          geom_point(shape = 21,
                     alpha = 0.8) +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_size_area(max_size = 15) +
          ggtitle("") + 
          xlab("") +
          scale_size_area(max_size = 15) +
          ylab(paste0(
            "Proportion of affiliated taxa at >",
            isolate(input$BubbleAbundThresh),
            "% relative abundance"
          )
          )
      }
    
    ## Modify the general theme, including panel borders
    if (input$BubblePanelBorder == "Yes") {
      BubblePlot <- BubblePlot +
        theme_bw() + theme(
          axis.text = element_text(colour = "black", size = 10),
          axis.line = element_blank(),
          strip.background.y = element_rect(fill = "white"),
          strip.background.x = element_rect(fill = "white"),
          panel.spacing = unit(as.numeric(input$BubblePanelSpacing), "points"),
          legend.position = "bottom",
          #panel.grid = element_line(colour = "grey"),
          #axis.line.y = element_line(colour="black",size=0.5),
          #panel.border = element_blank(),
          #text = element_text(size=10),
          axis.text.x = element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1
          ),
          # plot.margin=unit(c(-0.30,0,0,0), "null")
          
          #legend.position = "none")
        )
    }
    
    ## modify the general theme, removing panel borders
    if (input$BubblePanelBorder == "No") {
      BubblePlot <- BubblePlot +
        theme_bw() + theme(
          axis.text = element_text(colour = "black", size = 10),
          strip.background.y = element_rect(fill = "white"),
          strip.background.x = element_rect(fill = "white"),
          panel.spacing = unit(as.numeric(input$BubblePanelSpacing), "points"),
          legend.position = "bottom",
          #panel.grid = element_line(colour = "grey"),
          #axis.line.y = element_line(colour="black",size=0.5),
          panel.border = element_blank(),
          #text = element_text(size=10),
          axis.text.x = element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1
          ),
          # plot.margin=unit(c(-0.30,0,0,0), "null")
          #legend.position = "none")
        )
    }
    BubblePlot
    
    if(input$BubbleInclRead == TRUE && input$BubbleInclTaxa){
      BubbleReadPlot <- BubbleReadVisual()
      BubbleTaxaPlot <- BubbleTaxonVisual()
      layout <- "AA##
                 BBCC"
      BubblePlot <- BubbleReadPlot / BubblePlot + BubbleTaxaPlot + 
        plot_layout(design = layout, heights = c(0.25,1,0.10), widths = c(0.25,1,0.01))
      BubblePlot
    }
    
    else if (input$BubbleInclRead == TRUE){
      BubbleReadPlot <- BubbleReadVisual()
      BubblePlot <- BubbleReadPlot + BubblePlot +
        plot_layout(ncol = 1, heights = c(0.1,1))
      BubblePlot
    }
    
    else if (input$BubbleInclTaxa == TRUE){
      BubbleTaxaPlot <- BubbleTaxonVisual()
      layout <- "AAC"
      BubblePlot <- BubblePlot + BubbleTaxaPlot +
        plot_layout(design = layout, nrow = 1, heights = c(0.10,1))
      BubblePlot
    }
    else {
      BubblePlot <- BubblePlot
    }
    
  })
  
  # Settings for rendering the bubble plot output
  BubblePlotHeight <- reactive(input$BubblePlotOutH)
  BubblePlotWidth <- reactive(input$BubblePlotOutW)
  
  output$BubblePlotOut <- renderPlot({
    BubbleReadPlot <- BubbleReadVisual()
    BubbleTaxaPlot <- BubbleTaxonVisual()
    BubblePlot <- BubblePlotVisual()
    BubblePlot
    },
  width = BubblePlotWidth,
  height = BubblePlotHeight)
  
  # Bubble plot download
  output$BubblePlotDownload <- downloadHandler(
    filename = "bubble_plot.pdf",
    contentType = ".pdf",
    content = function(Bubble) {
      ggsave(
        Bubble,
        plot = BubblePlotVisual(),
        device = "pdf",
        height = as.numeric(input$BubblePlotOutH),
        width = as.numeric(input$BubblePlotOutW),
        units = "px",
        scale = 4
      )
    }
  )
  
  # Bubble plot table download
  output$BubbleTableDownload <- downloadHandler(
    filename = "bubble_data_table.csv",
    content = function(BubbleTable) {
      write.csv(BubbleDataLong(), BubbleTable)
    }
  ) # End of bubble plot
  
# ---- Bubble Plot Read Plot ----
  # An optional addition to the bubble plot that adds read counts above each samples. Read counts remain constant regardless
  # of changes to added/removed taxa through filtering, and are only removed when their respecitve samples are filtered.
  
  BubbleReadVisual <- reactive({
    
    req(input$BubbleStartButton)
    ASVTable <- FilteredTable()
    MetaData <- MainMetaTable()
    BubbleTable <- BubbleDataLong()
    
    # Tally reads and add metadata
    FeatureCounts <- rbind(ASVTable, Total = colSums(ASVTable))
    FeatureCounts <- FeatureCounts["Total", ]
    FeatureCounts <- reshape2::melt(FeatureCounts,
                                    variable.name = as.character("SampleName"),
                                    value.name = "Total"
    )
    FeatureCounts <- left_join(FeatureCounts, MetaData, by = "SampleName")
    
    # Now we need to ratify any missing samples in the bubble plot by filtering the read plot based on bubble plot sample names
    UniqueSamples <- unique(BubbleTable$SampleName)
    FilteredFeatureCounts <- FeatureCounts %>%
      filter(SampleName %in% UniqueSamples)
    
    Plot <- ggplot(FilteredFeatureCounts,
                   aes(x = reorder(SampleName,
                                   isolate(!!sym(input$BubbleSortXAxis))),
                       y = Total,
                       width = 0.9
                   )
    )
    
    Plot <- Plot + geom_bar(aes(fill = "grey"),
                            colour = "black",
                            size = 0.5,
                            alpha = 0.8,
                            stat = "identity",
                            position = "stack"
                            )
    
    Plot <- Plot + scale_fill_manual(values = "grey")
    Plot <- Plot + ylab("Reads")
    
    Plot <- Plot + theme_bw() +
      theme(
        panel.grid = element_blank(),
        text = element_text(colour = "black"),
        #axis.line = element_line(colour = "black"),
        axis.line = element_blank(),
        axis.text = element_text(colour = "black", size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(face = "italic", size = 12),
        legend.title = element_text(size = 10),
        panel.spacing = unit(as.numeric(input$BubblePanelSpacing), "points"),
        legend.position = "none",
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        axis.title = element_text(size = 10, face = NULL),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_blank(),
        # plot.margin=unit(c(-0.30,0,0,0), "null")
      )
    
    if (isolate(input$BubbleSecondFacet) == TRUE){
      Plot <- Plot + facet_nested(~ eval(parse(text = isolate(input$BubbleFacet))) + 
                                    eval(parse(text = isolate(input$BubbleSecondFacetMeta))),
                                  space = "free",
                                  scales = "free",
                                  switch = "both"
      )
      }
    
    # Default faceting - single level
    if (isolate(input$BubbleSecondFacet) == FALSE){
      Plot <- Plot + facet_nested(~ eval(parse(text = isolate(input$BubbleFacet))),
                     space = "free",
                     scales = "free",
                     switch = "both"
      )
      }
    
    if (isolate(input$BubbleThirdFacet) == TRUE){
      Plot <- Plot + facet_nested( ~ eval(parse(text = isolate(input$BubbleFacet))) * eval(parse(text = isolate(input$BubbleSecondFacetMeta))) * eval(parse(text = isolate(input$BubbleThirdFacetMeta))),
                    space = "free",
                    scales = "free",
                    switch = "both"
        )
    }
    
    Plot
    
  })
  
# ---- Bubble Taxon Proportion Plot ----
  # An option plot that calculates the proportion of taxa among all samples in the dataset. Pretty niche use case
  BubbleTaxonVisual <- reactive({
    
    req(input$BubbleStartButton)
    
    FilteredTable <- FilteredTable()
    MetaData <- MainMetaTable()
    FeatureKey <- FeatureDataKey()
    BubbleTable <- BubbleDataLong()
    
    # Count total number of reads
    TotalReads <- sum(FilteredTable)
    
    # Generate total counts for each taxon among the entire dataset
    FeatureSums <- data.frame(Sums = rowSums(FilteredTable))
    FeatureSums$FeatureID <- rownames(FeatureSums)
    
    # Join the taxon labels
    FeatureSums <- left_join(FeatureSums, FeatureKey, by = "FeatureID")
    
    
    # Calculate a proportion
    FeatureSums$Proportion <- FeatureSums$Sums / TotalReads * 100
    
    # Now remove any taxa not present in the main bubble plot
    UniqeTaxa <- unique(BubbleTable$Taxon)
    FilteredFeatureSums <- FeatureSums %>% filter(Taxon %in% UniqeTaxa)
    
    # Separate the full lineage into individual taxonomic levels
    FilteredFeatureSums <- separate(FilteredFeatureSums,
                                    Labels,
                                    c("Domain",
                                      "Phylum",
                                      "Class",
                                      "Order",
                                      "Family",
                                      "Genus",
                                      "Species",
                                      "Sub1",
                                      "Sub2",
                                      "Sub3",
                                      "Sub4"),
                                    sep = ";",
                                    remove = TRUE,
                                    convert = FALSE
    )
    
    # Now plot the visualization
    Plot <- ggplot(data = FilteredFeatureSums,
                   aes(x = reorder(Taxon, desc(Taxon)),
                       y = Proportion))
    
    Plot <- Plot + geom_bar(aes(),
                            position = position_dodge2(),
                            stat = "identity",
                            colour = "black",
                            size = 0.6,
                            alpha = 0.7,
                            width = 0.9
    )
    
    Plot <- Plot + scale_fill_manual(values = "grey")
    
    Plot <- Plot + facet_nested(get(input$BubbleTaxSort)~., scales = "free", space = "free")
    
    Plot <- Plot + coord_flip()
    

    
    Plot <- Plot + theme_bw() + theme(
      panel.grid = element_blank(),
      text = element_text(colour = "black"),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black",size=0),
      axis.line.x.bottom = element_line(size=-0),
      axis.text = element_text(colour = "black",size=12),
      axis.text.x = element_text(angle = 60, hjust =1.4, vjust=1.2,size=12,face = "plain"),
      legend.text = element_text(face = "plain",size = 16),
      legend.title = element_text(size=16),
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      strip.text.x = element_text(size=10,face="bold"),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.spacing = unit(as.numeric(input$BubblePanelSpacing), "points"),
      panel.border = element_rect(colour="black",size=0,fill=NA),
      axis.ticks = element_line(colour = "black"),
      axis.ticks.y = element_blank(),
      #panel.grid.major.y = element_line(colour = "black"),
      axis.line.y = element_line(colour="black"),
      strip.text.y = element_blank(),
      plot.margin=unit(c(0,0,0,0), "null")
    )
    
    Plot
    
  })
  
  
# ---- Bray-Curtis Triplot Visualization ----
  
  # Populate the drop down menus and update inputs
  observeEvent(input$MetaFile, {
    req(input$MetaFile)
    MetaData <- MainMetaTable()
    MetaColNames <- c(colnames(MetaData), "Taxon")
    MetaColNames <- MetaColNames[MetaColNames != "SampleName"]
    updateSelectInput(session, "PFillCol", choices = sort(MetaColNames))
    updateSelectInput(session, "PElipCol", choices = sort(MetaColNames))
    updateSelectInput(session, "PShape", choices = sort(MetaColNames))
    updateSelectInput(session, "PPalletSelect")
    updateSelectInput(session, "PGradient")
  })
  
  # Create an SRS table
  PSrsTable <- reactive({
    
    req(input$PStartButton)
    
    SrsTable <- FilteredTable()
    FeatureKey <- FeatureDataKey()
    SrsDepth <- isolate(input$SrsDepth)

    SrsTable <- SRS(SrsTable, Cmin = SrsDepth, seed = 123)
    rownames(SrsTable) <- FeatureKey$FeatureID
    
    SrsTable
  })
  
  # Generate a proportion table
  PPropTable <- reactive({
    
    req(input$PStartButton)
    
    SrsTable <- PSrsTable()
    MetaData <- MainMetaTable()
    
    # Create a proportion table
    SrsPropTable <- prop.table(as.matrix(SrsTable),2) * 100
    SrsPropTable <- as.data.frame(SrsPropTable)
    
    # Remove zeros
    # SrsTable <- SrsTable[, colSums(is.na(SrsTable)) == 0]
    
    # If analyzing contaminants, match between ASV table and ASV list, and keep only those that match
    if (input$ContamChoice == "Analyze"){
      Contams <- MainContamTable()
      Matches <- intersect(rownames(SrsPropTable),rownames(Contams))
      SrsPropTable <- SrsPropTable[Matches, ]
      
      SrsPropTable
      
    }
    SrsPropTable
  })
  
  #Generate the PCoA table
  PPcoa <- reactive({
    
    req(input$PStartButton)
    PPropTable <- PPropTable()
    
    TPropTable <- t(PPropTable)
    BrayCurtis <- as.matrix(vegdist(TPropTable, method = "bray"))
    
    # Perform the PCoA
    PPcoa <- pcoa(as.dist(BrayCurtis))
    PPcoa
  })
  
  # Isolate the axis1 and axis2 coordinates
  PCoords <- reactive({
    
    req(input$PStartButton)
    PPcoa <- PPcoa()
    
    PCoords <- data.frame(PCoA1 = PPcoa$vectors[, 1],
                              PCoA2 = PPcoa$vectors[, 2],
                              row.names = rownames(PPcoa$vectors)
                              )
    PCoords
  })
  
  #Isolate the relative corrected eigen values
  PEigenValues <- reactive({
    
    req(input$PStartButton)
    PPcoa <- PPcoa()
    
    # Extract relative eigenvalues
    EigenValues <- PPcoa$values[3]
    
    # Extract the first two columns (axis 1 and axis 2)
    PEigenValues <- data.frame(Axis1 = EigenValues[1,] * 100, Axis2 = EigenValues[2,] * 100)
    PEigenValues
  })
  
  # Filter the metadata for missing samples after SRS rarefaction
  PMetaFilt <- reactive({
    
    req(input$PStartButton)
    Pcoa <- PPcoa()
    MetaData <- MainMetaTable()
    SrsTable <- PSrsTable()
    
    # Collect the column names after SRS rarefaction
    PColNames <- colnames(SrsTable)
    
    # Filter the metadata table
    MetaData <-
      MetaData %>% filter(SampleName %in% PColNames)
    MetaData
    MetaData
  })
  
  # Fit the environment variables to the PCoA coordinates
  PEnvFit <- reactive({
    
    req(input$PStartButton)
    Pcoa <- PPcoa()
    MetaData <- PMetaFilt()
    SrsTable <- PSrsTable()
    
    # convert to a dataframe
    PVectors <- as.data.frame(Pcoa$vectors)
    PEnvFit <- envfit(PVectors, MetaData, perm = 10000)
    
    ## Scales the arrow vectors so they aren't huge
    PEnvFitFinal <- as.data.frame(PEnvFit$vectors$arrows * sqrt(PEnvFit$vectors$r))
    PEnvFitFinal <- cbind(PEnvFitFinal, PEnvFit$vectors$r)
    PEnvFitFinal <- cbind(PEnvFitFinal, PEnvFit$vectors$pvals)
    colnames(PEnvFitFinal) <- c("Axis1", "Axis2","R2", "pvalue")
    PEnvFitFinal$R2_rounded <- round(PEnvFitFinal$R2, 5)
    PEnvFitFinal$pvalue <- round(PEnvFitFinal$pvalue, 5)
    PEnvFitFinal
  })
  
  # Filter the data below a given R or p-value threshold
  PEnvFitFilt <- reactive({
    
    req(input$PStartButton)
    PEnvFit <- PEnvFit()
    
    # Filter
    PEnvFitFilt <- filter(PEnvFit,
                          PEnvFit$pvalue < (input$PEnvPThresh) &
                            PEnvFit$R2 > (input$PEnvRThresh))
    PEnvFitFilt
  })
  
  #Fit taxonomy abundances to PCoA data
  PTaxonScores <- reactive({
    
    req(input$PStartButton)
    MetaData <- PMetaFilt()
    FeatureKey <- FeatureDataKey()
    SrsProp <- PPropTable()
    Pcoa <- PPcoa()
  
    # Transpose
    TSrsProp <- t(SrsProp)
    
    # Use the wascores function to calculate weighted average scores 
    TaxonWeightedScores <- wascores(Pcoa$vectors[, 1:3], TSrsProp)

    
    # Remove NA values
    TaxonWeightedScores[is.na(TaxonWeightedScores)] <- 0
    
    # Calculate normalized scores, based on total abundance of each taxon
    NormWeightedScores <- data.frame(Abundance = (apply(TSrsProp, 2, sum)) / sum(TSrsProp))
    
    # Combine weighted scores with normalized abundance
    TaxonWeightedScores <- cbind(TaxonWeightedScores, NormWeightedScores$Abundance)
    colnames(TaxonWeightedScores) <- c("Axis1","Axis2","Axis3","Abundance")
    TaxonWeightedScores <- as.data.frame(TaxonWeightedScores)
    
    # Filter based on user-defined taxonomy relative abundance thresholds
    TaxonWeightedScores <- filter(TaxonWeightedScores,
                                  TaxonWeightedScores$Abundance > input$PTaxaThresh / 100
    )
    
    # Set a FeatureID column, then append taxonomic information using the FeatureKey
    TaxonWeightedScores$FeatureID <- rownames(TaxonWeightedScores)
    TaxonWeightedScores <- left_join(TaxonWeightedScores, FeatureKey, by = "FeatureID")
    
    TaxonWeightedScores
    
  })
  
  # Generate the PCoA plot
  PPlotVisual <- reactive({
    req(input$PStartButton)
    Pcoa <- PCoords()
    MetaData <- PMetaFilt()
    PEnvFitFilt <- PEnvFitFilt()
    Eigen <- PEigenValues()
    TaxonScores <- PTaxonScores()
    
    # Insert a SampleName column then merge the metadata
    Pcoa$SampleName <- rownames(Pcoa)
    Pcoa <- left_join(Pcoa, MetaData, by = "SampleName")
    
    # Change to categorical data 
    Pcoa <- Pcoa %>% mutate_if(!names(.) %in% c("PCoA1", "PCoA2"), factor)
    
    # merged_df <-
    #     merged_df %>% mutate_if(!names(.) %in% c("PCoA1", "PCoA2"), factor)
      
      # Define the available shapes and colors
      # I need to add the option to include a gradient
      AvailShapes <- c(21, 22, 23, 24, 25, 14, 13:1)
      AvailColours <- 2:27
      AvailFill <- 2:27
      
      
      
      #If the user selects shape options
      if ((input$PShapeChoice) == TRUE){
        PPlot <- ggplot(Pcoa, aes(x = PCoA1, y = PCoA2)) +
          geom_point(size = input$PSizeSelect, aes(
            shape = get(isolate(input$PShape)),
            colour = get(input$PFillCol),
            fill = get(input$PFillCol))) + 
          
          labs(x = paste0(
            c("Axis1","(",round(Eigen$Axis1,1),"%)")),
            y = paste0(
            c("Axis1","(",round(Eigen$Axis2,1),"%)")),
            fill = input$PFillCol,
            colour = input$PFillCol,
            shape = input$PShape
          ) +
          scale_shape_manual(values = AvailShapes)
        
        # Uses a colour gradient if selected

        if (input$PGradient == TRUE){
          PPlot <- PPlot + scale_colour_viridis(option = input$PPalletSelect,discrete = TRUE, direction = -1)
          PPlot <- PPlot + scale_fill_viridis(option = input$PPalletSelect,discrete = TRUE, direction = -1)
        } else {
          PPlot <- PPlot + scale_fill_manual(values = AvailFill) +
            scale_colour_manual(values = AvailFill)
        }
      }
      
      #If the user does not select shapes
      if ((input$PShapeChoice) == FALSE){
        
        PPlot <- ggplot(data = Pcoa, aes(x = PCoA1, y = PCoA2)) +
          geom_point(size = input$PSizeSelect, aes(
            colour = !!sym(input$PFillCol),
            fill = !!sym(input$PFillCol)
          )) + 
        labs(x = paste0(
          c("Axis1","(",round(Eigen$Axis1,1),"%)")),
          y = paste0(
            c("Axis1","(",round(Eigen$Axis2,1),"%)")),
          fill = input$PFillCol,
          colour = input$PFillCol,
          shape = input$PShape
        )
        
        ## If colou r gradient is selected
        if (input$PGradient == TRUE){
          PPlot <- PPlot + scale_colour_viridis(option = input$PPalletSelect,discrete = TRUE, direction = -1)
          PPlot <- PPlot + scale_fill_viridis(option = input$PPalletSelect,discrete = TRUE, direction = -1)
        } else {
          PPlot <- PPlot + scale_fill_manual(values = AvailFill) +
            scale_colour_manual(values = AvailFill)
        }
      }
      
      #Add coordinates and line segments for the environmental data, but only if it is present 
      if (dim(PEnvFitFilt)[1] != 0) {
        PPlot <- PPlot + geom_segment(
          data = PEnvFitFilt,
          aes(
            x = 0,
            y = 0,
            xend = PEnvFitFilt$Axis1,
            yend = PEnvFitFilt$Axis2,
          ),
          show.legend = FALSE,
          arrow = arrow(ends = "last")
        ) +
          geom_label(
            data = PEnvFitFilt,
            aes(
              label = rownames(PEnvFitFilt),
              x = PEnvFitFilt$Axis1 / 2,
              y = PEnvFitFilt$Axis2 / 2
            ),
            size = 4
          )}
      
      #If sample labels are selected:
      if (input$PSampleLabel == TRUE){
        PPlot <- PPlot +
          geom_text(aes(label = SampleName))
      }
      
      
      #Add elipses to the data
      if ((input$PElips) == TRUE) {
        
        PPlot <-PPlot + stat_ellipse(aes(color = get(input$PFillCol)),
                                   show.legend = FALSE)
      }
      
      # Add taxon abundance data, but only if it is present
      if (dim(TaxonScores)[1] != 0) {
        PPlot = PPlot + geom_point(data = TaxonScores,
                                   aes(Axis1, Axis2, size = round(Abundance * 100, digits = 0)),
                                   inherit.aes = FALSE,
                                   shape = 21,
                                   fill = NA,
                                   colour = "black",
                                   show.legend = TRUE
        ) +
          labs(size = "Relative abundance")+
          scale_size_area(max_size = 15)
      }
      
      if (dim(TaxonScores)[1] != 0) {
        # Add taxon annotation
        PPlot <- PPlot + geom_text(
          data = TaxonScores,
          aes(Axis1, Axis2, label = Taxon),
          inherit.aes = FALSE,
          size = 4
        )
      }
      
      # Customize plot aesthetics
      PPlot <- PPlot +
        theme(
          panel.grid = element_blank(),
          text = element_text(colour = "black"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 12),
          axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size = 14,
            face = "plain"
          ),
          axis.text.y.left = element_text(size = 14, face = "plain"),
          legend.text = element_text(face = "italic", size = 16),
          legend.title = element_text(size = 16),
          axis.title = element_text(size = 16, face = NULL),
          axis.text.y = element_text(size = 14),
          strip.text.x = element_text(size = 10, face = "bold"),
          panel.spacing = unit(0, "lines"),
          panel.border = element_rect(
            colour = "black",
            size = 1,
            fill = NA
          ),
          axis.ticks = element_line(colour = "black"),
          axis.ticks.y = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black")
        ) +
        guides(fill = guide_legend(override.aes = list(shape = 22)))
      
      # Plot PCoA
      PPlot
    
  })
  
  
  # Downloads the stats table
  output$PEnvFitTableSave <- downloadHandler(
    filename = "BrayCurtisStats.csv",
    content = function(Table) {
      write.csv(PEnvFit(), Table)
    })
  
  output$PPlotSave <- downloadHandler(
    filename = "BrayCurtisTriplot.pdf",
    contentType = ".pdf",
    content = function(Pcoa) {
      ggsave(
        Pcoa,
        plot = PPlotVisual(),
        device = "pdf",
        height = as.numeric(input$PPlotOutH),
        width = as.numeric(input$PPlotOutW),
        units = "px",
        scale = 4
      )
    })
  
  # Render the stats table
  output$PTableOut <- renderTable({
    Table <- PEnvFit()
    Table
  })
  
  # Render the PCoA plot
  PPlotHeight <- reactive(input$PPlotOutH)
  PPlotWidth <- reactive(input$PPlotOutW)
  
  output$PPlotOut <- renderPlot({
    PPlot <- PPlotVisual()
    PPlot
  },
  width = PPlotWidth,
  height = PPlotHeight)
  
  
# ---- UniFrac Plot ----
  # Herein lies the code of the fabled UniFrac triplot. Some of the code may seem redundant, but it isn't. I've written it so that
  # it collects all the necessary data with one pass through the UniFrac function, because it is time consuming so I don't want
  # to repeat it for minor changes
  
  ## Update selection - UniFrac PCoA plot
  observe({
    req(input$MetaFile)
    
    MetaData <- MainMetaTable()
    MetaColNames <- c(colnames(MetaData), "Taxon")
    MetaColNames <- MetaColNames[MetaColNames != "SampleName"]
    
    updateSelectInput(session, "UniFilCol", choices = sort(MetaColNames))
    updateSelectInput(session, "UniShape", choices = sort(MetaColNames))
    updateSelectInput(session, "UniPalletSelect")
    updateSelectInput(session, "UniGradient")
  })
  
  # Import the tree
  MainUniTree <- reactive({
    
    req(input$UniStartButton)
    req(input$UniTree)
    # req(input$main_file)
    # req(input$meta_file)
    read.tree(
      file = input$UniTree$datapath
    )
  })
  
  # SRS the table
  UniSrsTable <- reactive({
    
    req(input$UniStartButton)
    req(input$UniTree)
    
    ASVTable <- FilteredTable()
    MetaData <- MainMetaTable()
    FeatureKey <- FeatureDataKey()

    SrsTable <- SRS(ASVTable, Cmin = isolate(input$UniDepth), seed = 123)
    rownames(SrsTable) <- FeatureKey$FeatureID
    SrsTable
  })
  
  # Create a proportion table
  UniPropTable <- reactive({
    req(input$UniStartButton)
    
    SrsTable <- UniSrsTable()
    
    SrsPropTable <- prop.table(as.matrix(SrsTable),2) * 100
    SrsPropTable <- as.data.frame(SrsPropTable)
    
    # If analyzing contaminants, match between ASV table and ASV list, and keep only those that match
    if (input$ContamChoice == "Analyze"){
      Contams <- MainContamTable()
      Matches <- intersect(rownames(SrsPropTable),rownames(Contams))
      SrsPropTable <- SrsPropTable[Matches, ]
      
      SrsPropTable
      
    }
    
    SrsPropTable
    
  })
  
  #Filter the metadata table
  UniMetaFilt <- reactive({
    
    req(input$UniStartButton)
    
    SrsTable <- UniSrsTable()
    MetaData <- MainMetaTable()
    
    # Collect the column names after SRS rarefaction
    FilColNames <- colnames(SrsTable)
    
    # # Filter the metadata table to remove samples no longer present after SRS
    MetaData <- MetaData %>% filter(SampleName %in% FilColNames)
    MetaData
  })
  
  #Generate the UniFrac distances; must use a transposed proportion table
  UniDiss <- reactive({
    
    req(input$UniStartButton)
    PropTable <- UniPropTable()
    UniTree <- MainUniTree()
    
    PropTable <- t(PropTable)
    
    UniDiss <- GUniFrac(PropTable, UniTree, alpha = c(0,0.5,1))$unifracs
    UniDiss
  })
  
  #Separate the weighted UniFrac; d_1 is weighted
  UniDissWeighted <- reactive({
    req(input$UniStartButton)
    UniDiss <- UniDiss()
    UniDiss <- UniDiss[, , "d_1"]
    UniDiss
  })
  
  #Separate the unweighted UniFrac; d_UW unweighted
  UniDissUnweighted <- reactive({
    req(input$UniStartButton)
    UniDiss <- UniDiss()
    UniDiss <- UniDiss[, , "d_UW"]
    UniDiss
  })
  
  #Generate the PcoA for the weighted UniFrac
  UniWeightedPcoa <- reactive({
    req(input$UniStartButton)
    UniDiss <- UniDissWeighted()
    
    Pcoa <- ape::pcoa(UniDiss, correction = "cailliez")
    Pcoa
  })
  
  #Generate the PcoA for the unweighted UniFrac
  UniUnweightedPcoa <- reactive({
    req(input$UniStartButton)
    UniDiss <- UniDissUnweighted()
    
    Pcoa <- ape::pcoa(UniDiss, correction = "cailliez")
    Pcoa
  })
  
  # Capture the first and second axis coordinates for the weighted UniFrac
  UniWeightedCoords <- reactive({
    
    req(input$UniStartButton)
    Pcoa <- UniWeightedPcoa()
    
    Coords <- data.frame(Axis1 = Pcoa$vectors[, 1],
                         Axis2 = Pcoa$vectors[, 2],
                         row.names = row.names(Pcoa$vectors))
    Coords
  })
  
  # Capture the first and second axis coordinates for the weighted UniFrac
  UniUnweightedCoords <- reactive({
    
    req(input$UniStartButton)
    Pcoa <- UniUnweightedPcoa()
    
    Coords <- data.frame(Axis1 = Pcoa$vectors[, 1],
                         Axis2 = Pcoa$vectors[, 2],
                         row.names = row.names(Pcoa$vectors))
    Coords
  })
  
  # Capture the Eigen values
  UniWeightedEigen <- reactive({
    
    req(input$UniStartButton)
    Pcoa <- UniWeightedPcoa()
    
    EigenValues <- Pcoa$values[3]
    
    EigenValues <- data.frame(Axis1 = EigenValues[1,] * 100,
                              Axis2 = EigenValues[2,] * 100)
    EigenValues
  })
  
  # Capture the Eigen values
  UniUnweightedEigen <- reactive({
    
    req(input$UniStartButton)
    Pcoa <- UniUnweightedPcoa()
    
    EigenValues <- Pcoa$values[3]
    
    EigenValues <- data.frame(Axis1 = EigenValues[1,] * 100,
                              Axis2 = EigenValues[2,] * 100)
    EigenValues
  })
  
  #Fit the environmental variables to the PCoA for the unweighted UniFrac. 
  UniEnvFit <- reactive({
    
    req(input$UniStartButton)
    MetaData <- UniMetaFilt()
    
    # meta_data_table <- meta_datafile()
    if (input$UniDissSelect == "unweighted"){
      Pcoa <- UniUnweightedPcoa()
      PcoaVectors <- as.data.frame(Pcoa$vectors)
      EnvFit <- envfit(PcoaVectors, MetaData, perm = 10000)
      ## Scales the arrow vectors so they aren't huge
      EnvFitDf <- as.data.frame(EnvFit$vectors$arrows * sqrt(EnvFit$vectors$r))
      EnvFitDf <- cbind(EnvFitDf, EnvFit$vectors$r)
      EnvFitDf <- cbind(EnvFitDf, EnvFit$vectors$pvals)
      colnames(EnvFitDf) <- c("Axis1", "Axis2", "R2", "pvalue")
      EnvFitDf$R2Rounded <- round(EnvFitDf$R2, 2)
      EnvFitDf$pvalue <- round(EnvFitDf$pvalue, 4)
      # pcoa_envfit_df$pvaluecorr <- pcoa_envfit_df$pvalue*100
      # colnames(pcoa_envfit_df) <- c("axis1", "axis2", "R", "pvalue","pcorrected")
      EnvFitDf
    }
    #Fit the environmental variables to the PCoA for the weighted UniFrac
    else if (input$UniDissSelect == "weighted"){
      Pcoa <- UniWeightedPcoa()
      PcoaVectors <- as.data.frame(Pcoa$vectors)
      EnvFit <- envfit(PcoaVectors, MetaData, perm = 10000)
      ## Scales the arrow vectors so they aren't huge
      EnvFitDf <- as.data.frame(EnvFit$vectors$arrows * sqrt(EnvFit$vectors$r))
      EnvFitDf <- cbind(EnvFitDf, EnvFit$vectors$r)
      EnvFitDf <- cbind(EnvFitDf, EnvFit$vectors$pvals)
      colnames(EnvFitDf) <- c("Axis1", "Axis2", "R2", "pvalue")
      EnvFitDf$R2Rounded <- round(EnvFitDf$R2, 2)
      EnvFitDf$pvalue <- round(EnvFitDf$pvalue, 4)
      # pcoa_envfit_df$pvaluecorr <- pcoa_envfit_df$pvalue*100
      # colnames(pcoa_envfit_df) <- c("axis1", "axis2", "R", "pvalue","pcorrected")
      EnvFitDf
    }
  })
  
  #Filter the environmental statistics based on R or p-values.
  UniEnvFitFilt <- reactive({
    
    req(input$UniStartButton)
    EnvFit <- UniEnvFit()
    
    EnvFitFilt <- filter(EnvFit,
                         EnvFit$pvalue < input$UniEnvPThresh &
                           EnvFit$R2 > input$UniEnvRThresh)
    EnvFitFilt
  })
  
  
  
  #Fit taxonomy abundances to PCoA data
  UniTaxonScores <- reactive({
    
    req(input$UniStartButton)
    MetaData <- UniMetaFilt()
    FeatureKey <- FeatureDataKey()
    PropTable <- UniPropTable()
    
    if (input$UniDissSelect == "unweighted"){
      Pcoa <- UniUnweightedPcoa()
      Coords <- UniUnweightedCoords()
      Pcoa
    } 
    
    else if (input$UniDissSelect == "weighted"){
      Pcoa <- UniWeightedPcoa()
      Coords <- UniWeightedCoords()
      Pcoa
    }
    
    TPropTable <- t(PropTable)
    
    # Use the wascores function to calculate weighted average scores 
    TaxonWeightedScores <- wascores(Pcoa$vectors[, 1:3], TPropTable)
    
    # Remove NA values
    TaxonWeightedScores[is.na(TaxonWeightedScores)] <- 0
    
    # Calculate normalized scores, based on total abundance of each taxon
    NormWeightedScores <- data.frame(Abundance = (apply(TPropTable, 2, sum)) / sum(TPropTable))
    
    # Combine weighted scores with normalized abundance
    TaxonWeightedScores <- cbind(TaxonWeightedScores, NormWeightedScores$Abundance)
    colnames(TaxonWeightedScores) <- c("Axis1","Axis2","Axis3","Abundance")
    TaxonWeightedScores <- as.data.frame(TaxonWeightedScores)
    
    # Filter based on user-defined taxonomy relative abundance thresholds
    TaxonWeightedScores <- filter(TaxonWeightedScores,
                                  TaxonWeightedScores$Abundance > input$UniTaxaThresh / 100
    )
    
    # Set a FeatureID column, then append taxonomic information using the FeatureKey
    TaxonWeightedScores$FeatureID <- rownames(TaxonWeightedScores)
    TaxonWeightedScores <- left_join(TaxonWeightedScores, FeatureKey, by = "FeatureID")
    
    TaxonWeightedScores
    
  })
  
  # Generate the UniFrac triplot visual
  UniPlotVisual <- reactive({
    
    req(input$UniStartButton)
    MetaData <- UniMetaFilt()
    EnvFit <- UniEnvFitFilt()
    TaxonScores <- UniTaxonScores()

    # Call the specific UniFrac data
    if(input$UniDissSelect == "unweighted"){
      Pcoa <- UniUnweightedCoords()
      Eigen <- UniUnweightedEigen()
    }
    
    else if (input$UniDissSelect == "weighted"){
      Pcoa <- UniWeightedCoords()
      Eigen <- UniWeightedEigen()
    }
    
    # Insert a SampleName column then merge the metadata
    Pcoa$SampleName <- rownames(Pcoa)
    Pcoa <- left_join(Pcoa, MetaData, by = "SampleName")
    
    # Change to categorical data 
    Pcoa <- Pcoa %>% mutate_if(!names(.) %in% c("Axis1", "Axis2"), factor)
    
    # merged_df <-
    #     merged_df %>% mutate_if(!names(.) %in% c("PCoA1", "PCoA2"), factor)
    
    # Define the available shapes and colors
    # I need to add the option to include a gradient
    AvailShapes <- c(21, 22, 23, 24, 25, 14, 13:1)
    AvailColours <- 2:27
    AvailFill <- 2:27
    
    
    
    #If the user selects shape options
    if ((input$UniShapeSelect) == TRUE){
      UniPlot <- ggplot(Pcoa, aes(x = Axis1, y = Axis2)) +
        geom_point(size = input$UniSize, aes(
          shape = get(isolate(input$UniShape)),
          colour = get(input$UniFilCol),
          fill = get(input$UniFilCol))) + 
        
        labs(x = paste0(
          c("Axis1","(",round(Eigen$Axis1,1),"%)")),
          y = paste0(
            c("Axis1","(",round(Eigen$Axis2,1),"%)")),
          fill = input$UniFilCol,
          colour = input$UniFilCol,
          shape = input$UniShape
        ) +
        scale_shape_manual(values = AvailShapes)
      
      # Uses a colour gradient if selected
      
      if (input$UniGradient == TRUE){
        UniPlot <- UniPlot + scale_colour_viridis(option = input$UniPalletSelect,discrete = TRUE, direction = -1)
        UniPlot <- UniPlot + scale_fill_viridis(option = input$UniPalletSelect,discrete = TRUE, direction = -1)
      } else {
        UniPlot <- UniPlot + scale_fill_manual(values = AvailFill) +
          scale_colour_manual(values = AvailFill)
      }
    }
    
    #If the user does not select shapes
    if ((input$UniShapeSelect) == FALSE){
      
      UniPlot <- ggplot(data = Pcoa, aes(x = Axis1, y = Axis2)) +
        geom_point(size = input$UniSize, aes(
          colour = !!sym(input$UniFilCol),
          fill = !!sym(input$UniFilCol)
        )) + 
        labs(x = paste0(
          c("Axis1","(",round(Eigen$Axis1,1),"%)")),
          y = paste0(
            c("Axis1","(",round(Eigen$Axis2,1),"%)")),
          fill = input$UniFilCol,
          colour = input$UniFilCol,
          shape = input$UniShape
        )
      
      ## If colou r gradient is selected
      if (input$UniGradient == TRUE){
        UniPlot <- UniPlot + scale_colour_viridis(option = input$UniPalletSelect,discrete = TRUE, direction = -1)
        UniPlot <- UniPlot + scale_fill_viridis(option = input$UniPalletSelect,discrete = TRUE, direction = -1)
      } else {
        UniPlot <- UniPlot + scale_fill_manual(values = AvailFill) +
          scale_colour_manual(values = AvailFill)
      }
    }
    
    #Add coordinates and line segments for the environmental data, but only if it is present 
    if (dim(EnvFit)[1] != 0) {
      UniPlot <- UniPlot + geom_segment(
        data = EnvFit,
        aes(
          x = 0,
          y = 0,
          xend = EnvFit$Axis1,
          yend = EnvFit$Axis2,
        ),
        show.legend = FALSE,
        arrow = arrow(ends = "last")
      ) +
        geom_label(
          data = EnvFit,
          aes(
            label = rownames(EnvFit),
            x = EnvFit$Axis1 / 2,
            y = EnvFit$Axis2 / 2
          ),
          size = 4
        )}
    
    #If sample labels are selected:
    if (input$UniSampleLabel == TRUE){
      UniPlot <- UniPlot +
        geom_text(aes(label = SampleName))
    }
    
    
    #Add elipses to the data
    if ((input$UniElips) == TRUE) {
      
      UniPlot <- UniPlot + stat_ellipse(aes(color = get(input$UniFilCol)),
                                   show.legend = FALSE)
    }
    
    # Add taxon abundance data, but only if it is present
    if (dim(TaxonScores)[1] != 0) {
      UniPlot <- UniPlot + geom_point(data = TaxonScores,
                                 aes(Axis1, Axis2, size = round(Abundance * 100, digits = 0)),
                                 inherit.aes = FALSE,
                                 shape = 21,
                                 fill = NA,
                                 colour = "black",
                                 show.legend = TRUE
      ) +
        labs(size = "Relative abundance")+
        scale_size_area(max_size = 15)
    }
    
    if (dim(TaxonScores)[1] != 0) {
      # Add taxon annotation
      UniPlot <- UniPlot + geom_text(
        data = TaxonScores,
        aes(Axis1, Axis2, label = Taxon),
        inherit.aes = FALSE,
        size = 4
      )
    }
    
    # Customize plot aesthetics
    UniPlot <- UniPlot +
      theme(
        panel.grid = element_blank(),
        text = element_text(colour = "black"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 12),
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = 14,
          face = "plain"
        ),
        axis.text.y.left = element_text(size = 14, face = "plain"),
        legend.text = element_text(face = "italic", size = 16),
        legend.title = element_text(size = 16),
        axis.title = element_text(size = 16, face = NULL),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 10, face = "bold"),
        panel.spacing = unit(0, "lines"),
        panel.border = element_rect(
          colour = "black",
          size = 1,
          fill = NA
        ),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.y = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")
      ) +
      guides(fill = guide_legend(override.aes = list(shape = 22)))
    
    # Plot PCoA
    UniPlot
  })
  
  
  # Downloads for two tables and the plot
  output$UniStatsFullDownload <- downloadHandler(
    filename = "UniFractStatsTable.csv",
    content = function(Table) {
      write.csv(UniEnvFit(), Table)
    })

  output$UniPlotDownload <- downloadHandler(
    filename = "UniFracPlot.pdf",
    contentType = ".pdf",
    content = function(Plot) {
      ggsave(
        Plot,
        plot = UniPlotVisual(),
        device = "pdf",
        height = as.numeric(input$UniPlotOutH),
        width = as.numeric(input$UniPlotOutW),
        units = "px",
        scale = 4
      )
    })
  
  
  # Output the stats table
  output$UniTableOut <- renderDataTable({
    Table <- UniEnvFit()
    Table
  })
  
  # Render the PCoA plot
  UniPlotHeight <- reactive(input$UniPlotOutH)
  UniPlotWidth <- reactive(input$UniPlotOutW)

  output$UniPlotOut <- renderPlot({
    UniPlot <- UniPlotVisual()
    UniPlot
  },
  width = UniPlotWidth,
  height = UniPlotHeight)

    
  
  
  
 } # End of server