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
  output$UnifracPlotOut <- renderDT(NULL)
  
  
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
    read.table(
      file = input$MainFile$datapath,
      fill = TRUE,
      header = TRUE,
      sep = "\t"
    )
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
  output$MetaTable <- renderDataTable({
    req(input$MetaFile)
    read.table(
      file = input$MetaFile$datapath,
      fill = TRUE,
      header = TRUE,
      sep = "\t"
    )
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
  output$ContamTable <- renderDataTable({
    req(input$ContamFile)
    read.table(
      file = input$ContamFile$datapath,
      fill = TRUE,
      header = TRUE,
      sep = "\t"
    )
  })
  
# ---- Data pre processing ----
  # Now that the files are uploaded and stored, they have to be filtered to eliminate missing data in both the primary ASV table and the metadata file.     # This means filtering samples from both data tables. 
  
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
  
# ---- Data processing ----
  # We now process the main ASV table to format taxonomy labels and unique identifiers necessary for downstream visualizations and processing. 
  
  # The first requirement is to format taxonomically collapsed tables, if directed to do so
  TransDataCollapsed <- reactive({
    ASVTable <- MainASVTable()
    
    # If collapse table checkbox is checked, then transform. Otherwise, leave unaltered.
    if (input$IsMainCollapsed == TRUE) {
      ASVTable$Consensus.Lineage <- ASVTable$Feature.ID # Set the feature ID to imitate the taxonomy in an uncollapsed table
      ASVTable$rowID <- 1:nrow(ASVTable) # Add a column of unique numbers to imitate unique feature IDs
    } else {
      if (input$IsMainCollapsed == FALSE) {
          ASVTable$rowID <- 1:nrow(ASVTable) # Add a column called rowID; if present it will be overwritten, so user caution is advised
        }
    }
    ASVTable
  })
  
  # Show this table in the main panel
  output$ProcessedMainTable <- renderDataTable({
    Table <- TransDataCollapsed()
    output$proc_maintext <- renderText("This is your processed data")
    Table
  })
  
  # Now we make a key file that will keep Feature ids, row ids, taxonomy, and representative sequence information properly associated with samples
  # This is the only location in which taxonomic identifiers, labels, or otherwise should be modified so that everything is consistent.
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
    if (input$TruncateTaxa == "Yes") {
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
    }
    
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
    
    # FullLineage <- data.frame(paste(FeatureKey$Consensus.Lineage, data_tran$rowID, sep = "_"))
    # lineage_OTU <- as.data.frame(paste(gsub(".*;", "", labels), data_tran$rowID, sep = "_"))
    # colnames(lineage_OTU) <- "TaxaName"
    # colnames(full_lineage) <- "TaxaName"
    # rownames(data_tran) <- full_lineage$TaxaName
    # data_tran$appended_taxonomy <- full_lineage$TaxaName
    # 
    # ## Keep this unfiltered
    # unfiltered_table <- data_tran
    # unfiltered_table
    FeatureKey
  })
  
  # Show this the FeatureKey in the main panel
  output$ProcessedFeatureTable <- renderDataTable({
    Table <- FeatureDataKey()
    output$proc_maintext <- renderText("This is your processed data")
    Table
  })
  
  # Generate a table containing only count data, using the FeatureIDs as rownames to track data
  FilteredTable <- reactive({
    ASVTable <- MainASVTable()
    
    ColsToFilter <- c("OTU.ID",
                      "Feature.ID",
                      "rowID",
                      "Consensus.Lineage",
                      "ReprSequence"
                      )
    ASVTable <- ASVTable[,!(names(ASVTable) %in% ColsToFilter)]
  })
  
  
# ---- Read Plot Visualization ----
  # This is our first plot: the read plot. It shows sequence depth for each sample or groups of samples
  
  # First we update 
  observeEvent(input$MetaFile,{
    req(input$MetaFile)
    ## Incorporate the meta file column names into every appropriate input
    MetaData <- MainMetaTable()
    MetaColNames <- colnames(MetaData)
    # meta_colnames <- c(colnames(meta_datafile), "TaxaName")
    MetaColNames <- MetaColNames[MetaColNames != "SampleName"]
    ## Update the selections - read plot
    # isolate({
    updateSelectInput(session, "ReadSortByAxis", choices = sort(MetaColNames))
    updateSelectInput(session, "ReadMetaGroup", choices = sort(MetaColNames))
    # updateSelectInput(session,"read_meta_key",choices = read_meta_list)
    updateSliderInput(session, "ReadPlotOutW")
    # updateSelectInput(session,"read_colour",choices = colnames(meta_datafile()))
    # updateSelectInput(session, "read_sortby_axis", choices = sort(MetaColNames))
    # updateSelectInput(session, "read_meta_group", choices = sort(MetaColNames))
    updateSliderInput(session, "ReadPlotOutH")
    # })
  })
  
  # Now we transform the raw ASV counts into totals
  ReadCountsTotal <- reactive({
    
    req(input$ReadStartButton)
    ASVTable <- FilteredTable()
    MetaData <- MainMetaTable()
    
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
                           "Species"),
                         sep = ";",
                         remove = TRUE,
                         convert = FALSE
                         )
    
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
    
    # Add metadata. Because this collapses features, all feature data is lost
    BarTableSumFinal <- left_join(BarTableSumFinal, MetaData, by = "SampleName")
    
    # Now remove all columns excep SampleName? Why?
    # BarFilt <- colnames(MetaData)
    # BarFilt <- BarFilt[BarFilt != "SampleName"]
    # BarTableSumFinal <- select(BarTableSumFinal, -BarFilt)
    
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
                           "Species"),
                         sep = ";",
                         remove = TRUE,
                         convert = FALSE
    )

    
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
    
    # # Set the taxonomy as factors for later ordering. This doesn't do much, but if you remove it the colouring settings change
    # data_long_bubble <-
    #   data_long_bubble[with(data_long_bubble,
    #                         order(eval(parse(
    #                           text = isolate(input$BubbleTaxSort)
    #                         )), TaxaName, decreasing = TRUE)), ]
    # data_long_bubble$TaxaName <-
    #   as.character(data_long_bubble$TaxaName)
    # data_long_bubble$TaxaName <-
    #   factor(data_long_bubble$TaxaName,
    #          levels = unique(data_long_bubble$TaxaName))
    # data_long_bubble
    
    BubbleTable
  })
  
  # Now create the bubbleplot visual
  BubblePlotVisual <- reactive({
    
    req(input$BubbleStartButton)
    BubbleTable <- BubbleDataLong()
    FeatureKey <- FeatureDataKey()
    MetaData <- MainMetaTable()
    
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
    
    if (isolate(input$BubbleInclPercent) == "Yes") {
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
      } else {
        BubblePlot <- BubblePlot +
        geom_point(shape = 21,
                   alpha = 0.8) +
        ggtitle("")
        + xlab("") +
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
          plot.margin=unit(c(-0.30,0,0,0), "null")
          
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
          plot.margin=unit(c(-0.30,0,0,0), "null")
          #legend.position = "none")
        )
    }
  
    BubblePlot
    
  })
  
  
  
  
  
  BubblePlotHeight <- reactive(input$BubblePlotOutH)
  BubblePlotWidth <- reactive(input$BubblePlotOutW)
  
  output$BubblePlotOut <- renderPlot({
    BubblePlot <- BubbleDataTran()
    Test <- BubbleDataLong()
    BubblePlot <- BubblePlotVisual()
    BubblePlot
  },
  width = BubblePlotWidth,
  height = BubblePlotHeight)
  



 } # End of server