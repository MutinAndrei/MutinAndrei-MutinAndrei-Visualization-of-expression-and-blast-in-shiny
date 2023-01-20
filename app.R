# needed packages
library(DT)
library(thematic)
library(bslib)
library(XML)
library(plyr)
library(dplyr)
library(rBLAST)
library(ggplot2)
library(tidyr)
library(xlsx)
library(shinyalert)
library(Biostrings)
library(shinycssloaders)
library(ggtext)
library(plotly)
library(glue)
library(rclipboard)
library(BiocManager)



custom_db <- c("Eulimnogammarus verrucosus")
custom_db1 <- c("Eulimnogammarus cyaneus")
custom_db2 <- c("Gammarus lacustris")

options(repos = BiocManager::repositories())


#The Blast code was taken from "https://github.com/ScientistJake/Shiny_BLAST"

ui <- fluidPage(theme = bs_theme(version = 5, bootswatch = "united"),
                tagList(
                  tags$head(
                    tags$link(rel="stylesheet", type="text/css",href="style.css"),
                    tags$script(type="text/javascript", src = "busy.js")
                  )
                ),
                
                #This block gives us all the inputs:
                
                  titlePanel(h1("GenexploT")),
                  tabsetPanel(
                    tabPanel(h2('Blast'),
                             br(),
                             useShinyalert(force = TRUE),  # Set up shinyalert
                             actionButton("preview", "Instructions for using this app"),
                             br(),
                             br(),
                             br(),
                             helpText("Insert your nucleotide / amino acid sequence"),
                             textAreaInput('query', 'Input sequence:', value = "", placeholder = "", width = "600px", height="200px"),
                             helpText("Choose a database and parametrs Blast"),
                             br(),
                             selectInput("db" , "Data base:", choices=c(custom_db, "Eulimnogammarus verrucosus", custom_db1, "Eulimnogammarus cyaneus", 
                                                                       custom_db2, "Gammarus lacustris"), width="350px"),
                             div(style="display:inline-block",
                                 selectInput("program", "Program:", choices=c("blastn","tblastn"), width="200px")),
                             div(style="display:inline-block",
                                 selectInput("eval", "e-value:", choices=c(0.001,1e-4,1e-5,1e-10, 1), width="150px")),
                             br(),
                             helpText("Click on the button 'BLAST!' and wait for the result"),
                             br(),
                             br(),
                             actionButton("blast", "BLAST!"),
                             
                             
                             mainPanel(
                               h4("Results"),
                               rclipboardSetup(), 
                               DT::dataTableOutput("blastResults"),
                               p("Alignment:", tableOutput("clicked") ),
                               verbatimTextOutput("alignment")
                             ),
                             
              
                             ),
                    
                    tabPanel(h2("Fasta"),
                             br(),
                             actionButton("fasta", "Instructions for using the 'Fasta' tab"),
                             br(),
                             br(),
                             helpText("Select the data base to continue the work"),
                             selectInput("fa" , "Data base:", choices=c("Eulimnogammarus verrucosus", "Eulimnogammarus cyaneus", "Gammarus lacustris"), 
                                         width="350px"),
                             shinycssloaders::withSpinner(
                               DT::dataTableOutput("mytable2"), image = "https://i.gifer.com/ALdI.gif", 
                               image.width = 200, image.height = 200)),
                    
                    
                    tabPanel(h2("Table of Genes"),
                             br(),
                             actionButton("tog", "Instructions for using the 'Table of Genes' tab"),
                             br(),
                             br(),
                             shinycssloaders::withSpinner(
                               DT::dataTableOutput("mytable"), 
                               image = "https://i.gifer.com/ALdI.gif", image.width = 200, image.height = 200)),
                    
                    
                    tabPanel(h2("Expression plot"),
                             br(),
                             actionButton("expr", "Instructions for using the 'Expression plot' tab"),
                             br(),
                             br(),
                             helpText("Select the data base to continue the work"),
                             selectInput("exp" , "Data base:", choices=c("Eulimnogammarus verrucosus", "Eulimnogammarus cyaneus", "Gammarus lacustris"), 
                                         width="350px"),
                             shinycssloaders::withSpinner(DT::dataTableOutput("mytable1"), 
                                                          image = "https://i.gifer.com/ALdI.gif", image.width = 200, image.height = 200),
                             br(),
                             br(),
                             
                             shinycssloaders::withSpinner(plotlyOutput('plot1', height = 800, width = "100%"), 
                                                          image = "https://i.gifer.com/ALdI.gif", image.width = 200, image.height = 200))
                    
                    ),
)



thematic::thematic_shiny()

#The Blast code was taken from "https://github.com/ScientistJake/Shiny_BLAST"
server <- function(input, output, session){
  
  
  #page 1
  
  
  custom_db <- c("Eulimnogammarus verrucosus")
  custom_db1 <- c("Eulimnogammarus cyaneus")
  custom_db2 <- c("Gammarus lacustris")
  custom_db_path <- c("E.verr_TSA.blastdb")
  custom_db_path1 <- c("allcyaneusTSA1.blastdb")
  custom_db_path2 <- c("alllacustrisTSA1.blastdb")
  
  blastresults <- eventReactive(input$blast, {
    
    #gather input and set up temp file
    query <- input$query
    tmp <- tempfile(fileext = ".fa")
    
    #if else chooses the right database
    if (input$db == custom_db){
      db <- custom_db_path
      remote <- c("")
    } else if (input$db == custom_db1) {
      db <- custom_db_path1
      remote <- c("")
    }
      else if (input$db == custom_db2){
       db <- custom_db_path2
       remote <- c("")
    }
      else if (input$db != custom_db){
       db <- c("E.verr_TSA")
       #add remote option for nr since we don't have a local copy
       remote <- c("-remote")
    }
      else if (input$db != custom_db1) {
        db <- c("allcyaneusTSA1")
        #add remote option for nr since we don't have a local copy
        remote <- c("-remote")
    }
      else if (input$db != custom_db2) {
        db <- c("alllacustrisTSA1")
    #add remote option for nr since we don't have a local copy
    remote <- c("-remote")
      }
    
    #this makes sure the fasta is formatted properly
    if (startsWith(query, ">")){
      writeLines(query, tmp)
    } else {
      writeLines(paste0(">Query\n",query), tmp)
    }
    
    #calls the blast
    data <- system(paste0(input$program," -query ",tmp," -db ",db," -evalue ",input$eval," -outfmt 5 -max_hsps 1 -max_target_seqs 10 ", remote), intern = T)
    xmlParse(data)
  }, ignoreNULL= T)
  
  #Now to parse the results...
  parsedresults <- reactive({
    if (is.null(blastresults())){}
    else {
      xmltop = xmlRoot(blastresults())
      
      #the first chunk is for multi-fastas
      results <- xpathApply(blastresults(), '//Iteration',function(row){
        query_ID <- getNodeSet(row, 'Iteration_query-def') %>% sapply(., xmlValue)
        hit_IDs <- getNodeSet(row, 'Iteration_hits//Hit//Hit_id') %>% sapply(., xmlValue)
        hit_length <- getNodeSet(row, 'Iteration_hits//Hit//Hit_len') %>% sapply(., xmlValue)
        bitscore <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_bit-score') %>% sapply(., xmlValue)
        eval <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_evalue') %>% sapply(., xmlValue)
        cbind(query_ID,hit_IDs,hit_length,bitscore,eval)
      })
      #this ensures that NAs get added for no hits
      results <-  rbind.fill(lapply(results,function(y){as.data.frame((y),stringsAsFactors=FALSE)}))
    }
  })
  
  #makes the datatable
  output$blastResults <- renderDataTable({
    if (is.null(blastresults())){
    } else {
      parsedresults()
    }
  }, selection="single")
  
  #this chunk gets the alignemnt information from a clicked row
  output$clicked <- renderTable({
    if(is.null(input$blastResults_rows_selected)){}
    else{
      xmltop = xmlRoot(blastresults())
      clicked = input$blastResults_rows_selected
      tableout<- data.frame(parsedresults()[clicked,])
      
      tableout <- t(tableout)
      names(tableout) <- c("")
      rownames(tableout) <- c("Query ID","Hit ID", "Length", "Bit Score", "e-value")
      colnames(tableout) <- NULL
      data.frame(tableout)
    }
  },rownames =T,colnames =F)
  
  #this chunk makes the alignments for clicked rows
  output$alignment <- renderText({
    if(is.null(input$blastResults_rows_selected)){}
    else{
      xmltop = xmlRoot(blastresults())
      
      clicked = input$blastResults_rows_selected
      
      #loop over the xml to get the alignments
      align <- xpathApply(blastresults(), '//Iteration',function(row){
        top <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_qseq') %>% sapply(., xmlValue)
        mid <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_midline') %>% sapply(., xmlValue)
        bottom <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_hseq') %>% sapply(., xmlValue)
        rbind(top,mid,bottom)
      })
      
      #split the alignments every 40 carachters to get a "wrapped look"
      alignx <- do.call("cbind", align)
      splits <- strsplit(gsub("(.{40})", "\\1,", alignx[1:3,clicked]),",")
      
      #paste them together with returns '\n' on the breaks
      split_out <- lapply(1:length(splits[[1]]),function(i){
        rbind(paste0("Q-",splits[[1]][i],"\n"),paste0("M-",splits[[2]][i],"\n"),paste0("H-",splits[[3]][i],"\n"))
      })
      unlist(split_out)
    }
  })
  47
  #Data for plot and table
  

  
  #Loads the expression data and converts the data frame
  
  
  

  #page 2
  sfas <- read.csv(file = "Blastsfas.csv")
  fastaEC <- read.csv(file = "fastaEC.csv")
  fastaGL <- read.csv(file = "fastaGL.csv")
  

  
  observe({
    if (input$fa == "Eulimnogammarus verrucosus") {output$mytable2 = renderDataTable(
    sfas, options = list(pageLength = 1,
                         server = FALSE,
                         dom = 'Bfrtip',
                         buttons = c('copy', 'csv')
    ),
    extensions = 'Buttons')
  } else if (input$fa == "Eulimnogammarus cyaneus") {output$mytable2 = renderDataTable(
    fastaEC, options = list(pageLength = 1,
                         server = FALSE,
                         dom = 'Bfrtip',
                         buttons = c('copy', 'csv')
    ),
    extensions = 'Buttons')
  } else if (input$fa == "Gammarus lacustris") {output$mytable2 = renderDataTable(
    fastaGL, options = list(pageLength = 1,
                            server = FALSE,
                            dom = 'Bfrtip',
                            buttons = c('copy', 'csv')
    ),
    extensions = 'Buttons')
  } else if (input$fa == "") {return(NULL)}
  })
  
  
  
  
  #make table 
  tableoggenes <- read.csv("ToGa.csv")
  output$mytable <- renderDataTable(
    tableoggenes, options = list(paging = TRUE,    ## paginate the output
                                 pageLength = 6,  ## number of rows to output for each page
                                 server = FALSE,   ## use client-side processing
                                 dom = 'Bfrtip',
                                 buttons = c('csv', 'excel'),
                                 columnDefs = list(list(targets = '_all', className = 'dt-center'),
                                                   list(targets = c(0, 8, 9), visible = FALSE))
    ),
    extensions = 'Buttons')
  

  
  
  
  #page 3 
  expressionEC <- read.csv("expressionplotEC.csv")
  expressionEC$Temperature <- as.factor(expressionEC$Temperature)
  expressionEC[,1] <- factor(expressionEC[,1], levels=unique(expressionEC[,1]))
  
  expression <- read.csv("expressionplot.csv")
  expression$Temperature <- as.factor(expression$Temperature)
  expression[,1] <- factor(expression[,1], levels=unique(expression[,1]))
  
  expressionGL <- read.csv("expressionplotGL.csv")
  expressionGL$Temperature <- as.factor(expressionGL$Temperature)
  expressionGL[,1] <- factor(expressionGL[,1], levels=unique(expressionGL[,1]))
  
  #make table2
observe({
if (input$exp == "Eulimnogammarus verrucosus"){
    output$mytable1 = DT::renderDataTable({
    datatable(expression, options = list(pageLength = 3,
                                        searchDelay	= 3000,
                                        server = FALSE,
                                        dom = 'Bfrtip',
                                        buttons = c('csv', 'excel')
    ), extensions = 'Buttons', filter = 'top')
  })

  #filtred table
    filtered_table <- reactive({
    req(input$mytable1_rows_all)
    expression[input$mytable1_rows_all, ]
  })
  #make plot and quantitative error
output$plot1 <- renderPlotly({
  if (nrow(filtered_table()) == 2662011) {
    return(NULL)}
  else if (nrow(filtered_table()) == 1943298) {
    return(NULL)
  }
  else if (nrow(filtered_table()) == 2120534) {
    return(NULL)
  }
  else if (nrow(filtered_table()) >= 1000) {
    return(shinyalert("Very much data, specify the parameters of the post :')", type = "error"))}
  else{
      g <- ggplot(filtered_table(), aes(x = filtered_table()$Temperature, y = filtered_table()$Expression, 
                                        col = Characteristic)) + geom_boxplot() + 
        ggtitle("Visualise the expression level of your chosen transcript/characteristic in the table above.") +
        xlab("Exposure temperature, °C") + 
        ylab("Expression level") +
        theme(text=element_text(size=14, color = "Black", family = "Helvetica"),
              legend.title = element_text(face = "bold"), 
              axis.title.x = element_text(face = "bold"),
              axis.title.y = element_text(face = "bold"),
              title = element_text(face = "bold", size = 16))}
    
    g <- ggplotly(g)
    g <- g %>% layout(boxmode = "group")
    g <- g %>% layout(legend = list(x = 0, y = -0.2))
    g <- g %>% layout(legend = list(orientation = 'h'))})}
  
else if (input$exp == "Eulimnogammarus cyaneus") {
    output$mytable1 = DT::renderDataTable({
      datatable(expressionEC, options = list(pageLength = 3,
                                             searchDelay	= 3000,
                                             server = FALSE,
                                             dom = 'Bfrtip',
                                             buttons = c('csv', 'excel')
        ), extensions = 'Buttons', filter = 'top')})

      #filtred table
      filtered_table <- reactive({
        req(input$mytable1_rows_all)
        expressionEC[input$mytable1_rows_all, ]
      })
      #make plot and quantitative error
output$plot1 <- renderPlotly({
  if (nrow(filtered_table()) == 2662011) {
    return(NULL)}
  else if (nrow(filtered_table()) == 1943298) {
    return(NULL)
  }
  else if (nrow(filtered_table()) == 2120534) {
    return(NULL)
  }
  else if (nrow(filtered_table()) >= 1000) {
    return(shinyalert("Very much data, specify the parameters of the post :')", type = "error"))}
  else{
    g <- ggplot(filtered_table(), aes(x = filtered_table()$Temperature, y = filtered_table()$Expression, 
                                            col = Characteristic)) + geom_boxplot() + 
      ggtitle("Visualise the expression level of your chosen transcript/characteristic in the table above.") +
      xlab("Exposure temperature, °C") + 
      ylab("Expression level") +
      theme(text=element_text(size=14, color = "Black", family = "Helvetica"),
                  legend.title = element_text(face = "bold"), 
                  axis.title.x = element_text(face = "bold"),
                  axis.title.y = element_text(face = "bold"),
                  title = element_text(face = "bold", size = 16))}
        
        g <- ggplotly(g)
        g <- g %>% layout(boxmode = "group")
        g <- g %>% layout(legend = list(x = 0, y = -0.2))
        g <- g %>% layout(legend = list(orientation = 'h'))})} 
else if (input$exp == "Gammarus lacustris") {
      output$mytable1 = DT::renderDataTable({
        datatable(expressionGL, options = list(pageLength = 3,
                                             searchDelay	= 3000,
                                             server = FALSE,
                                             dom = 'Bfrtip',
                                             buttons = c('csv', 'excel')
        ), extensions = 'Buttons', filter = 'top')})
    
      #filtred table
      filtered_table <- reactive({
        req(input$mytable1_rows_all)
        expressionGL[input$mytable1_rows_all, ]
      })
      #make plot and quantitative error
output$plot1 <- renderPlotly({
  if (nrow(filtered_table()) == 2662011) {
    return(NULL)}
  else if (nrow(filtered_table()) == 1943298) {
    return(NULL)
  }
  else if (nrow(filtered_table()) == 2120534) {
    return(NULL)
  }
  else if (nrow(filtered_table()) >= 1000) {
    return(shinyalert("Very much data, specify the parameters of the post :')", type = "error"))}
  else{
    g <- ggplot(filtered_table(), aes(x = filtered_table()$Temperature, y = filtered_table()$Expression, 
                                            col = Characteristic)) + geom_boxplot() + 
      ggtitle("Visualise the expression level of your chosen transcript/characteristic in the table above.") +
      xlab("Exposure temperature, °C") + 
      ylab("Expression level") +
      theme(text=element_text(size=14, color = "Black", family = "Helvetica"),
                  legend.title = element_text(face = "bold"), 
                  axis.title.x = element_text(face = "bold"),
                  axis.title.y = element_text(face = "bold"),
                  title = element_text(face = "bold", size = 16))}
        
        g <- ggplotly(g)
        g <- g %>% layout(boxmode = "group")
        g <- g %>% layout(legend = list(x = 0, y = -0.2))
        g <- g %>% layout(legend = list(orientation = 'h'))})} 
else {return(NULL)}})
  
  

  #make user guide 
  
  observeEvent(input$preview, {
    # Show a modal when the button is pressed
    shinyalert("Instructions for using the Blast tab!", "1. This tab allows you to search for transcript homologues and proteins. 
2. To use the blast in GenexploT you need to insert your nucleotide or amino acid sequence in the 'Input sequence:' window in .fasta format
3. Once your sequence has been placed in the 'Input sequence:' you need to select your search parameters.
3.1 The 'Program:' parameter allows you to choose the type of sequence to search for homologues, amino acid or nucleotide.
3.2 The 'e-value:' parameter allows you to choose the similarity criterion of the sequence being searched for.
4. After all the parameters have been configured, you need to press the 'BLAST!
5. After the search is done you will get the results of the blasta, namely the alignment and the sequence table.

Mail for requests and suggestions: andreimutin97@gmail.com

This application was developed by researchers of the Institute of Biology at Irkutsk State University with the financial support of Grant № 091-22-318", type = "info", size = "l", showCancelButton = TRUE, showConfirmButton = FALSE)
  })
  
  observeEvent(input$fasta, {
    # Show a modal when the button is pressed
    shinyalert("Instructions for using the Fasta tab", "1. This tab allows you to search for transcripts and helps you create a .fasta file. 
2. First, you need to select a database from the 'Data base:' pop-up list
3. After that you need to wait for a while until the table with sequences is loaded.
4. After the table is loaded you can search the transcript by its ID.
5. The final step is to click on the 'Copy' button which will save the row or rows of the table.

Mail for requests and suggestions: andreimutin97@gmail.com

This application was developed by researchers of the Institute of Biology at Irkutsk State University with the financial support of Grant № 091-22-318", type = "info", size = "l", showCancelButton = TRUE, showConfirmButton = FALSE)
  })
  
  observeEvent(input$tog, {
    # Show a modal when the button is pressed
    shinyalert("Instructions for using the Table of gens tab!", "1. The page where you are now is BLAST (if you need to find similar sequences to your reference sequence, go here).
               1.1. For BLAST you need to insert the sequence you are interested in in the blast field Input sequence 
               1.2. Once the sequence has been placed in the field, you need to select the database you want to search for.
1.3. After that you need to define search parameters Program: (define search type, by amino acid (tblastn) or nucleotide (blastn) sequence and e-value: (sequence similarity) P.S. the lower value of e-value: the higher similarity.
               1.4. Press the BLAST button and obtain the result 
2. On the second page called Table of Genes you can find the sequence of interest or what this sequence encodes, and also here you can perform a generalized search of gene expression depending on the incubation temperature, for example how the expression of the enzyme family dehydrogenase depends on the temperature.
2.1. To do this you need to enter the characteristics of interest (transcript ID or product name) in the search field
2.2. Also, on this page you can save the table you have created as .csv or .xlsx 
3. On page three, Expression plot, you can plot the expression plot of the sequences you are interested in, either by ID or by the name of the substance you are investigating.
3.1. To do this you need to enter either the ID or the name of the substance in the search box.

Mail for requests and suggestions: andreimutin97@gmail.com

This application was developed by researchers of the Institute of Biology at Irkutsk State University with the financial support of Grant № 091-22-318", type = "info", size = "l", showCancelButton = TRUE, showConfirmButton = FALSE)
  })
  
  observeEvent(input$expr, {
    # Show a modal when the button is pressed
    shinyalert("Instructions for using the Expression plot tab!", "1. This tab allows you to visualize the expression of a gene or genes by performing either a general search by gene ID or protein characteristic, or a pointwise search by either gene ID or protein characteristic.
2.  To visualize the expression you need to select the 'Data base:' database
3. After selecting the database you have to wait a little while for the table with the expression levels to load
4. Once the table is loaded you can search for the nucleotide sequence or protein you want by its name from the database https://www.uniprot.org/uniprotkb?query=*.
5. Once you have your graph you can save it using the tools in the upper right corner of the graph

Mail for requests and suggestions: andreimutin97@gmail.com

This application was developed by researchers of the Institute of Biology at Irkutsk State University with the financial support of Grant № 091-22-318", type = "info", size = "l", showCancelButton = TRUE, showConfirmButton = FALSE)
  })


  
}
      
      
      
shinyApp(ui = ui, server = server)