#-------------------------------------------------------------------------
#  Roseric Azondekon,
#  December 17, 2018
#  Last Update: February 11, 2019
#  ZSPH, Milwaukee, WI, USA
#-------------------------------------------------------------------------

options(shiny.maxRequestSize=50*1024^2) # Increasing file upload size to 50MB

source('./server/functions.R')

# #############################################
# ### FUNCTIONS
# ############################################
## Capture selected folder
setDir <- function(session, input){
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$directory
    },
    handlerExpr = {
      if (input$directory > 0) {
        # condition prevents handler execution on initial app launch

        # launch the directory selection dialog with initial path read from the widget
        default <- readDirectoryInput(session, 'directory')
        path <- choose.dir(default)

        # update the widget value
        updateDirectoryInput(session, 'directory', value = path)
      }
    }
  )
}

## Split path function
split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))

#
chooseFile <- function(){
  if(!is.null(input$Btn_GetFile)){
    file_selected<-parseFilePaths(volumes, input$Btn_GetFile)
    output$txt_file <- renderText(as.character(file_selected$datapath))
  }
  if(length(file_selected$datapath) != 0){
    output$create_index_btn <- renderUI({
      actionButton("create_index", "Create Index!")
    })
  }
}

#Server side programming...
server = function(input, output, session) {
  # output$session<-renderPrint(session_info())
  # values <- reactiveValues()
  # setDir(session, input)
  searchURL <- "http://useast.ensembl.org/Multi/Search/Results?q="#"https://www.ncbi.nlm.nih.gov/gene/?term="
  setDir(session, input)
  volumes <- getVolumes()
  smpDir <- "" #Global sample Directory variable
  wkd <- "" # Global working directory variable
  smp <- c() # Global samples vector variable
  rpt <- "" # Global QC report path
  QUANT_LOGS <- list()
  QUANT_SF <- list()
  conditions<-factor()
  results_dge <- list()
  edgeR_fit_mod <- list()

  local(shinyFileChoose(input, "Btn_GetFile", roots = volumes, session = session, filetypes=c('fa', 'fna'),defaultPath=getwd()))

  # Observe user interactions with Shiny Dashboard
  observe({

    if(!is.null(input$Btn_GetFile)){
      file_selected<-parseFilePaths(volumes, input$Btn_GetFile)
      output$txt_file <- renderText(as.character(file_selected$datapath))
    }
    if(length(file_selected$datapath) != 0){
      output$create_index_btn <- renderUI({
        actionButton("create_index", "Create Index!")
      })
    }

    # Get working directory
    workDir <- session$input[[sprintf("%s__chosen_dir", 'directory')]]
    home <- system("echo $HOME",intern=T)
    if(workDir != home){
      workDir <<- workDir
      setwd(workDir)
      files <- file.path(getwd(), list.files(pattern = "\\.fastq.gz$", recursive = TRUE))
      samples <- c()

      for(f in files){
        split <- split_path(f)
        samples<-c(samples,split[2])
      }
      smp <<- samples
      sampleDir<-""
      for(i in (length(split)-1):3){
        sampleDir <- paste0(sampleDir,split[i],'/')
      }
      sampleDir<-paste0("/",sampleDir)
      smpDir <- sampleDir

      if(length(files)>0){
        output$fastqc <- renderUI({
          actionButton("perform_qc","Run Quality Control!")
        })
      }

      output$summary <- renderUI({
        multiInput(
          inputId = "selSamples",
          label = "Found Samples:",
          choices = NULL,
          selected = unique(samples),
          choiceNames = unique(samples),
          choiceValues = unique(samples),
          options = list(
            enable_search = FALSE,
            non_selected_header = "Choose between:",
            selected_header = "You have selected:"
          )
        )
      })
    }

    # Display working Directory
    output$workDir <- renderPrint({
      workDir
    })

    report <- file.path(workDir,"fastqc_results/QC","multiqc_report.html")
    rpt <- report

    #################################################
    ### FASTQ QC check
    #################################################

    ## perform_qc_process()
    perform_qc_process <- function(){
      fq_reports <- list.files(file.path(workDir,"fastqc_results"),pattern = ".html$", recursive = TRUE)
      # if(file.exists(file.path(workDir,"fastqc_results/QC/multiqc_report.html"))){
      if(file.exists(file.path(workDir,"fastqc_results/QC/multiqc_report.html")) && length(fq_reports)>0){
        sendSweetAlert(
          session = session,
          title = "Information",
          text = "QC report already exists!",
          type = "info"
        )
        output$multiQc_report <- renderUI({
          actionButton("mqc", "Click to View MultiQC Report!")
        })
        output$delQc_report <- renderUI({
          actionButton("delQC", "Click to delete QC Report!")
        })
      }
      # else if(!file.exists(report) && length(input$selSamples)<=0){
      #   sendSweetAlert(
      #     session = session,
      #     title = "Error...",
      #     text = "Please choose at least one sample!",
      #     type = "error"
      #   )
      # }
      else if(!file.exists(file.path(workDir,"fastqc_results/QC/multiqc_report.html"))){
        # report <- rpt
        # workDir <<- getwd()
        samp <<- input$selSamples
        sampleDir<<-smpDir
        system(paste("cd",workDir,"&& mkdir -p fastqc_results"),wait = F)
        for(i in 1:length(samp)){
          cmd <- paste0("cd ",workDir," && mkdir -p fastqc_results/",samp[i]
                        ," && fastqc ",sampleDir,samp[i],"/*fastq.gz -o fastqc_results/",samp[i])
          system(cmd, wait=F)
          # updateProgressBar(session = session, id = "pb1", value = i*100/length(samp))
        }
        # system("echo -e '${LCYAN}QC check completed!${NC}'", wait=F)

        output$multiQc_report <- renderUI({
          actionButton("mqc", "Click to View MultiQC Report!")
        })
        output$delQc_report <- renderUI({
          actionButton("delQC", "Click to delete QC Report!")
        })
        local(sendSweetAlert(
          session = session,
          title = "Processing...",
          text = "QC report in progress... \nPlease wait and check progress in your terminal!",
          type = "info"
        ))
      }
    }

    # QC Check...
    observeEvent(input$perform_qc, {
      local(perform_qc_process())
    })

    #################################################################
    ## When "Click to View MultiQC Report!" button is clicked
    #################################################################
    ##mqc_process
    mqc_process <- function(){
      # report <- rpt
      if(file.exists(file.path(workDir,"fastqc_results/QC/multiqc_report.html"))){
        system(paste("xdg-open",file.path(workDir,"fastqc_results/QC/multiqc_report.html")),wait = F)
      }
      else{
        system(paste("cd",file.path(workDir,"fastqc_results")," && multiqc . --dirs -o QC/ ",
                     "&& xdg-open",file.path(workDir,"fastqc_results/QC/multiqc_report.html")),wait = F)
      }
    }

    observeEvent(input$mqc, {
      local(mqc_process())
    })

    #################################################################
    ## When "Click to delete QC Report!" button is clicked
    #################################################################
    ## delQC_process()
    delQC_process <- function(){
      # report <- rpt
      if(file.exists(file.path(workDir,"fastqc_results/QC/multiqc_report.html"))){
        confirmSweetAlert(
          session = session,
          inputId = "delQCconfirm",
          type = "warning",
          title = "Want to confirm ?",
          danger_mode = TRUE
        )
      }else if(!file.exists(file.path(workDir,"fastqc_results/QC/multiqc_report.html"))){
        sendSweetAlert(
          session = session,
          title = "Error...",
          text = "QC report does not exist!",
          type = "error"
        )
      }
    }

    observeEvent(input$delQC, {
      local(delQC_process())
    })


    ##############################
    observeEvent(input$delQCconfirm,{
      if(input$delQCconfirm){
        system(paste("cd .. && rm -r",file.path(getwd(),"fastqc_results")),wait = F)
        local(sendSweetAlert(
          session = session,
          title = "Information",
          text = "QC report successfully deleted!",
          type = "info"
        ))
      }
    })

    observeEvent(input$get_files_btn,{
      dt <- filedata(input$datafile)
      dat <- dt[,-1]
      rownames(dat) <- dt[,1]
      txi <<- list("counts"=dat)
      txi.header <<- colnames(dat)
      output$files_data_report <- renderPrint({
        paste(dim(dat)[1],"genes and",dim(dat)[2],"samples have been loaded!")
      })
      DT_output <- dat
      row.names(DT_output) <- paste0('<a target="_blank" href="',searchURL,rownames(DT_output),'">',rownames(DT_output),'</a>')
      output$files_data <- renderDataTable({DT_output} ,options = list(scrollX = TRUE), escape = F)

      output$sel_conds <- renderUI({
        multiInput(
          inputId = "sel_samp",
          label = "Select samples for each condition:",
          choices = NULL,
          selected = txi.header[1:round(length(txi.header)/2,0)],
          choiceNames = txi.header,
          choiceValues = txi.header,
          options = list(
            enable_search = FALSE,
            non_selected_header = "Condition A:",
            selected_header = "Condition B:"
          )
        )
      })

      output$set_cond <- renderUI({
        actionButton("set_cond_btn2", "Set conditions!")
      })
    })

    #################################################
    ### Transcriptome Indexing
    #################################################
    ## When 'Create Index!' button is clicked!
    indexing <- eventReactive(input$create_index, {
      if(workDir %in% c(home)){
        out<-"Please, set working directory!"
        out
      }
      else if(length(file_selected$datapath)==0){
        out<-"Please, select a reference transcriptome file!"
        out
      }
      else{
        ref <- as.character(file_selected$datapath)
        cmd <- paste("cd",workDir,"&&","if [ ! -d reference_index ]; then salmon index -t",ref,"-i reference_index; fi")
        cmd <- paste(cmd,"&& cat",file.path(workDir,"reference_index/indexing.log"))
        # cmd<-"cat /media/New_Volume/Roseric/salmon_tutorial/reference_index1/indexing.log"
        # capture.output(system(cmd,wait = F,intern=F))
        # cmd
        output$quantSamples <- renderUI({
          multiInput(
            inputId = "qSamples",
            label = "Select sample for read quantification:",
            choices = NULL,
            selected = unique(samples),
            choiceNames = unique(samples),
            choiceValues = unique(samples),
            options = list(
              enable_search = FALSE,
              non_selected_header = "Choose between:",
              selected_header = "You have selected:"
            )
          )
        })
        output$quant_read <- renderUI({
          actionButton("quantRead", "Click to Quantify Read!")
        })
        output$nCores <- renderUI({
          sliderInput("cores" ,"Number of threads:", min = 1, step = 1, round = 1, max = detectCores(), value = detectCores())
        })
        system(cmd,intern = T)
      }
    })

    output$index_output <- renderPrint({
      indexing()
    })

    #################################################################
    ## When "Click to Quantify Read!" button is clicked
    #################################################################
    ## quantRead_process
    quandRead_process <- function(){
      quants_files <- list.files(file.path(workDir,"quants"),pattern = ".sf$", recursive = TRUE)
      if(dir.exists(file.path(workDir,"quants")) && length(quants_files)){
        confirmSweetAlert(
          session = session,
          inputId = "myconfirmation",
          type = "warning",
          title = "Reads Quantification already exists! Overwrite?",
          danger_mode = TRUE
        )
        # rv$started <- TRUE
      }
      else{
        s <- input$qSamples
        for(sample in s){
          bashcmd <- "echo ''"
          fq_files <- list.files(file.path(smpDir,sample),pattern = ".fastq.gz$", recursive = TRUE)
          if(length(fq_files)==1){
            bashcmd <- paste0("cd ",workDir," && salmon quant -i reference_index -l A -r "
                              ,smpDir,sample,"/*.fastq.gz -p "
                              ,input$cores," -o quants/", sample,"_quant")
          }else if(length(fq_files)==2){
            bashcmd <- paste0("cd ",workDir," && salmon quant -i reference_index -l A -1 "
                              ,smpDir,sample,"/*_1.fastq.gz -2 "
                              ,smpDir,sample,"/*_2.fastq.gz -p "
                              ,input$cores," -o quants/", sample,"_quant")
          }

          system(bashcmd,wait = F)
          # rv$started <- TRUE
        }
        system("Read Quantification is completed!",wait = F)
      }
    }

    observeEvent(input$quantRead, {
      local(quandRead_process())
    })

    observeEvent(input$myconfirmation, {
      if(input$myconfirmation){
        s <- input$qSamples
        for(sample in s){
          bashcmd <- "echo ''"
          fq_files <- list.files(file.path(smpDir,sample),pattern = ".fastq.gz$", recursive = TRUE)
          if(length(fq_files)==1){
            bashcmd <- paste0("cd ",workDir," && salmon quant -i reference_index -l A -r "
                              ,smpDir,sample,"/*.fastq.gz -p "
                              ,input$cores," -o quants/", sample,"_quant")
          }else if(length(fq_files)==2){
            bashcmd <- paste0("cd ",workDir," && salmon quant -i reference_index -l A -1 "
                              ,smpDir,sample,"/*_1.fastq.gz -2 "
                              ,smpDir,sample,"/*_2.fastq.gz -p "
                              ,input$cores," -o quants/", sample,"_quant")
          }

          system(bashcmd,wait = F)
          # rv$started <- TRUE
        }
        # rv$started <- TRUE
      }
      else{
        # rv$started <- TRUE
      }
    })

    #################################################################
    ## When "Get TCGA RNASeq data!" button is clicked
    #################################################################
    observeEvent(input$getTCGA_btn,{
      withProgress(message=paste(input$dataType,"data will be imported! This may take some time!"),{
      tcga_DATA <- getTCGA(disease=input$disease_select, data.type=input$dataType)
      setProgress(value = 1/2,detail = "Preparing TCGA data...")
      if(!is.null(tcga_DATA)){
      txi <<- list("counts"=tcga_DATA$dat)
      # txi <<- list("counts"=tcga_data$dat, "lengths"=tcga_data$dat+1, "countsFromAbundance"="no")
      txi.header <<- colnames(tcga_DATA$dat)
      output$tcga_report <- renderPrint({
        n_genes <- dim(tcga_DATA$dat)[1]
        n_samples <- dim(tcga_DATA$dat)[2]
        paste(n_genes,"genes and",n_samples,"samples have been imported")
      })

      DT_output <- txi$counts
      row.names(DT_output) <- paste0('<a target="_blank" href="',searchURL,rownames(DT_output),'">',rownames(DT_output),'</a>')
      output$tcga_data <- renderDataTable({DT_output},options = list(scrollX = TRUE), escape = F)
      output$tcga_data_dld <- renderUI({
        downloadButton("tcga_data_dld_btn","Download TCGA data!")
      })

      output$sel_conds <- renderUI({
        multiInput(
          inputId = "sel_samp",
          label = "Select samples for each condition:",
          choices = NULL,
          selected = txi.header[1:round(length(txi.header)/2,0)],
          choiceNames = txi.header,
          choiceValues = txi.header,
          options = list(
            enable_search = FALSE,
            non_selected_header = "Condition A:",
            selected_header = "Condition B:"
          )
        )
      })
      setProgress(value = 1,detail = "Done!")
      output$set_cond <- renderUI({
        actionButton("set_cond_btn2", "Set conditions!")
      })}else{
        output$tcga_report <- renderPrint(
          "Error: No data available for download. Please ensure the data is available from TCGA."
        )
      }
      })

    })

    output$tcga_data_dld_btn <- downloadHandler(
      filename = paste0(input$disease_select,"_TCGA_data.csv"),
      content = function(file){
        write.csv(tcga_data, file)
      }
    )

    #############################################################
    ## When show_detail button is clicked
    #############################################################

    showDetails <- function(){
      if(dir.exists(file.path(workDir,"quants")) && length(file_selected$datapath)>0){
        quant_sf <-list()
        quant_logs <- list()

        sf_files <- file.path(workDir,"quants",list.files(path=file.path(workDir,"quants"), pattern = ".sf$", recursive = T))
        log_files <- file.path(workDir,"quants",list.files(path=file.path(workDir,"quants"), pattern = ".log$", recursive = T))
        for(sf in sf_files){
          quant_sf[paste0(split_path(sf)[2],".sf")]<-sf
        }
        for(log in log_files){
          quant_logs[paste0(split_path(log)[3],".log")] <- log
        }
        # log_files <- file.path(workDir,list.files(path=file.path(workDir,"quants"), pattern = ".sf$", recursive = T))
        # rv$quant_logs <- paste(readLines(file.path(workDir,"shinylogs/quantify_reads.txt")), collapse = "<br/>")
        # rv$quant_sf <- file.path(workDir,list.files(path=file.path(workDir,"quants"), pattern = ".log$", recursive = T))
      }
      if(exists("quant_logs") && length(quant_logs)>0){
        QUANT_LOGS <<- quant_logs # updting global variable QUANT_LOGS
        output$quant_log_files <- renderUI({
          pickerInput(
            inputId = "log",
            label = "Select a log",
            choices = names(QUANT_LOGS),
            # options = list(
            #   `actions-box` = TRUE,
            #   size = 10,
            #   `selected-text-format` = "count > 3"
            # ),
            multiple = F
          )
        })
        output$view_log <- renderUI({
          actionButton("view_log_btn", "Click to View Log!")
        })
      }

      if(exists("quant_sf") && length(quant_sf)>0){
        QUANT_SF <<- quant_sf # updting global variable QUANT_SF
        output$quant_sf_files <- renderUI({
          pickerInput(
            inputId = "quant",
            label = "Select a quantification file",
            choices = names(QUANT_SF),
            # options = list(
            #   `actions-box` = TRUE,
            #   size = 10,
            #   `selected-text-format` = "count > 3"
            # ),
            multiple = F
          )
        })
        output$view_quant <- renderUI({
          actionButton("view_quant_btn", "Display Read quant. table!")
        })

        output$annot <- renderUI({
          tagList(fluidRow(column(4,
                                  radioButtons("tx2gene",
                                               label = h3("Select a gene annotation for:"),
                                               choices = list("Homo sapiens" = "org.Hs.eg.db", "Mouse" = "org.Mm.eg.db",
                                                              "Zebrafish" = "org.Dr.eg.db", "Other"="other"),
                                               inline=T, selected = "org.Hs.eg.db")
          )
          ,column(4,
                  radioButtons("refType",
                               label = h3("Select the reference source:"),
                               choices = list("RefSeq NCBI Transcripts" = "REFSEQ", "Ensembl transcripts" = "ENSEMBLTRANS"),
                               inline=T, selected = "REFSEQ")
          )
          ,column(4,
                  radioButtons("transOut",
                               label = h3("Level estimate:"),
                               choices = list("Gene level" = "FALSE", "Transcript level" = "TRUE"),
                               inline=T, selected = "FALSE")
          )
          ))
        })

        output$cond_annot <- renderUI({
          conditionalPanel(
            condition = "input.tx2gene == 'other'",
            textInput("other_annot", "Specify a genome wide annotation database", value = "",
                     placeholder = "e.g: Enter 'org.Bt.eg.db' for Bovine Genome wide annotation
                     (see https://bioconductor.org/packages/release/BiocViews.html#___AnnotationData)")
          )
        })

        output$tr_output <- renderUI({
          actionButton("import","Import Estimates!")
        })

        output$import_output <- renderUI({
          withSpinner(verbatimTextOutput("tr_counts"),size=1,proxy.height=200,type=sample(c(1,4,5:8), 1))
        })
      }
    }
    ###end of showDetails()

    observeEvent(input$show_logs, {
      d_output <- showDetails()
    })

    observeEvent(input$view_quant_btn, {
      # if(length(QUANT_SF)>0){
      f_quant <- QUANT_SF[[input$quant]]
      system(paste("gedit",f_quant),wait=F)
      # }
    })
    # output$quant_output<-renderPrint({
    #   out_quant()
    # })

    observeEvent(input$view_log_btn, {
      # if(length(QUANT_LOGS)>0){
      f_log <- QUANT_LOGS[[input$log]]
      system(paste("gedit",f_log),wait=FALSE)
      # }
    })

    ###############################################################
    # When 'import transcript level estimates' button is clicked
    ###############################################################
    import_transcript <- eventReactive(input$import, {
      # Get qunt.sf files
      files <- file.path(workDir,"quants",list.files(file.path(workDir,"quants"),pattern = "\\.sf$", recursive = TRUE))
      header<-c()
      for(f in files){header<-c(header,gsub("_quant","",split_path(f)[2]))}
      txi.header <<- header
      # Annotation file
      annot_lib <- input$tx2gene
      if(input$tx2gene == 'other'){
        annot_lib <- input$other_annot
        }
      if(!(annot_lib %in% installed.packages()[,"Package"])){
        install.packages(annot_lib, repos = repositories)
        lapply(annot_lib, require, character.only = TRUE)
      }else{
        lapply(annot_lib, require, character.only = TRUE)
      }
      db<-get(annot_lib)

      # Get text2gene dataset
      txt2gene <- select(db,
                         keys = keys(db),
                         columns=c(input$refType,"SYMBOL")#columns=c("REFSEQ","SYMBOL")#,"GENENAME")
                         # ,keytype="ENTREZID"
                         )

      # Cleaning dataset
      txt2gene<-txt2gene[,-1]
      colnames(txt2gene)<-c("TXNAME","GENEID")#,"GENENAME")
      # txt2gene<-txt2gene[,-3]
      txt2gene<-txt2gene[complete.cases(txt2gene),]

      # transcript_out <- FALSE
      #
      # if(input$refType == "ENSEMBLTRANS"){transcript_out <- TRUE}

      # Import quantification reads into R
      txi <<- tximport(files, type = "salmon", tx2gene = txt2gene,ignoreTxVersion=T,txOut = as.logical(input$transOut))

      output$readCounts_dld <- renderUI({
        downloadButton(outputId = "counts_download",label = "Download transcripts data!")
      })

      output$tbChoice <- renderUI({
        selectInput("tb_select", "Select a table"
                    ,list("abundance"="abundance","counts"="counts","length"="length"),
                    selected = "counts", multiple = FALSE,
                    selectize = TRUE, width = "200px", size = NULL)
      })

      output$tb <- renderUI({
        dataTableOutput('table')
      })

      output$sel_conds <- renderUI({
        multiInput(
          inputId = "sel_samp",
          label = "Select samples for each condition:",
          choices = NULL,
          selected = txi.header[1:round(length(txi.header)/2,0)],
          choiceNames = txi.header,
          choiceValues = txi.header,
          options = list(
            enable_search = FALSE,
            non_selected_header = "Condition A:",
            selected_header = "Condition B:"
          )
        )
      })

      output$set_cond <- renderUI({
        actionButton("set_cond_btn", "Set conditions!")
      })

      # output$perform_DESeq2 <- renderUI({
      #   actionButton("perform_DESeq2_btn", "Run DESeq2!")
      # })

      list("countsFromAbundance"=txi$countsFromAbundance,'samples'=header)
    })
    ################ end of import_transcript() ###################

    output$counts_download <- downloadHandler(
      filename = paste0("read_",input$tb_select,".csv"),
      content = function(file){
        down_table <- txi[[input$tb_select]]
        colnames(down_table) <- txi.header
        write.csv(down_table, file)
      }
    )

    output$tr_counts <- renderPrint({
      import_transcript()
    })

    output$table <- renderDataTable({
      t <- txi[[input$tb_select]]
      colnames(t) <- txi.header
      row.names(t) <- paste0('<a target="_blank" href="',searchURL,rownames(t),'">',rownames(t),'</a>')
      t
    },options = list(scrollX = TRUE), escape = F)

    ######print_conds
    get_conds <- eventReactive(input$set_cond_btn,{
      vs <- list()
      for(s in txi.header){vs[[s]]<-"A"} #set condition "A" by default
      for(s in input$sel_samp){vs[[s]]<-"B"}
      conditions <<- factor(as.vector(factor(unlist(vs))))
      output$get_DESeq2_data <- renderUI({
        actionButton("get_DESeq2_data_btn", "Get dataset!")
      })
      conditions
    })

    get_conds2 <- eventReactive(input$set_cond_btn2,{
      vs <- list()
      for(s in txi.header){vs[[s]]<-"A"} #set condition "A" by default
      for(s in input$sel_samp){vs[[s]]<-"B"}
      conditions <<- factor(as.vector(factor(unlist(vs))))
      output$get_DESeq2_data <- renderUI({
        actionButton("get_DESeq2_data_btn2", "Get dataset!")
      })
      conditions
    })
    #######################################

    output$print_conds <- renderPrint({
      get_conds()
    })

    output$print_conds2 <- renderPrint({
      get_conds2()
    })

    ##################################################
    ### When "Get DESeq2 dataset!" button is clicked
    #################################################
    get_DESeq_data <- function(){
      library(magrittr)
      sampleTable <- data.frame(condition = conditions)
      rownames(sampleTable) <- txi.header
      dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
      # dds <- DESeqDataSetFromMatrix(txi, sampleTable, ~condition)
      dds <- dds[ rowSums(counts(dds)) > 0, ] # remove  genes  without  any  counts
      dds <<- dds

      #Output DESeq2 dataset summary
      output$conditions_info <- renderPrint({head(colData(dds))})
      DT_output <- assay(dds)
      row.names(DT_output) <- paste0('<a target="_blank" href="',searchURL,rownames(DT_output),'">',rownames(DT_output),'</a>')
      output$assay_info <- renderDataTable({DT_output},options = list(scrollX = TRUE), escape = F)
      output$assay_info_dld <- renderUI({
        downloadButton(outputId = "assay_info_dld_btn",label = "Download raw counts data!")
      })
      output$library_size <- renderPrint(colSums(counts(dds)))
      output$library_counts <- renderPrint(colSums(txi$counts))

      output$normalization <- renderUI({
        actionButton("normalization_btn", "Normalize dataset!")
      })

      output$bar_counts <- renderPlotly({
        p <- plotly::plot_ly(y=matrixStats::colSums2(counts(dds)), x=names(colSums(counts(dds))), type="bar")
        p
      })
    }
    ### end of get_DESeq_data()

    get_DESeq_data2 <- function(){
      library(magrittr)
      sampleTable <- data.frame(condition = conditions)
      rownames(sampleTable) <- txi.header
      # dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
      dds <- DESeqDataSetFromMatrix(txi$counts, sampleTable, ~condition)
      dds <- dds[ rowSums(counts(dds)) > 0, ] # remove  genes  without  any  counts
      dds <<- dds

      DT_output <- assay(dds)
      row.names(DT_output) <- paste0('<a target="_blank" href="',searchURL,rownames(DT_output),'">',rownames(DT_output),'</a>')
      #Output DESeq2 dataset summary
      output$conditions_info <- renderPrint({head(colData(dds))})
      output$assay_info <- renderDataTable({DT_output},options = list(scrollX = TRUE), escape = F)
      output$assay_info_dld <- renderUI({
        downloadButton(outputId = "assay_info_dld_btn",label = "Download raw counts data!")
      })
      output$library_size <- renderPrint(colSums(counts(dds)))
      output$library_counts <- renderPrint(colSums(txi$counts))

      output$normalization <- renderUI({
        actionButton("normalization_btn", "Normalize dataset!")
      })

      output$bar_counts <- renderPlotly({
        p <- plotly::plot_ly(y=matrixStats::colSums2(counts(dds)), x=names(colSums(counts(dds))), type="bar")
        p
      })
    }
    ### end of get_DESeq_data2()

    output$assay_info_dld_btn <- downloadHandler(
      filename = "RawCounts.csv",
      content = function(file){write.csv(assay(dds),file)}
    )

    #################################################

    observeEvent(input$get_DESeq2_data_btn,{
      withProgress(message="Preprocessing DESeq2 dataset...",{
        get_DESeq_data()
      })
      # get_DESeq_data()
    })

    observeEvent(input$get_DESeq2_data_btn2,{
      withProgress(message="Preprocessing DESeq2 dataset...",{
        get_DESeq_data2()
      })
      # get_DESeq_data2()
    })

    ##################################################
    ### When "Get DESeq2 dataset!" button is clicked
    #################################################
    norm_DESeq_data <- function(){
      setProgress(value = 1/5, detail = "Estimating size factors...")
      dds_norm <- estimateSizeFactors(dds)
      setProgress(value = 2/5, detail = "Count normalization...")
      counts.sf_normalized  <<- counts(dds_norm, normalized = TRUE)
      setProgress(value = 3/5, detail = "Count log-normalization...")
      log.norm.counts  <<- log2(counts.sf_normalized + 1)
      # DESeq.rlog  <<- rlog(dds_norm , blind = TRUE)
      if(dim(dds_norm)<1000){
        DESeq.rlog  <<- varianceStabilizingTransformation(dds_norm, blind = F) # added line for TCGA data
      }else{
        # DESeq.rlog  <<- varianceStabilizingTransformation(dds_norm, blind = F)
        DESeq.rlog  <<- vst(dds_norm) # added line for TCGA data
      }
      setProgress(value = 4/5, detail = "Count regularized log-normalization...")
      rlog.norm.counts  <<- assay(DESeq.rlog)
      dds <<- dds_norm

      DT_out <- counts.sf_normalized
      row.names(DT_out) <- paste0('<a target="_blank" href="',searchURL,rownames(DT_out),'">',rownames(DT_out),'</a>')
      output$norm_assay_info <- renderDataTable({round(DT_out,3)},options = list(scrollX = TRUE), escape = F)
      output$norm_assay_info_dld <- renderUI({
        downloadButton(outputId = "norm_assay_info_dld_btn",label = "Download Normalized RC data!")
      })

      DT_out2 <- log.norm.counts
      row.names(DT_out2) <- paste0('<a target="_blank" href="',searchURL,rownames(DT_out2),'">',rownames(DT_out2),'</a>')
      output$lnorm_assay_info <- renderDataTable({round(DT_out2,3)},options = list(scrollX = TRUE), escape = F)
      output$lnorm_assay_info_dld <- renderUI({
        downloadButton(outputId = "lnorm_assay_info_dld_btn",label = "Download log2-transformed RC!")
      })

      DT_out3 <- rlog.norm.counts
      row.names(DT_out3) <- paste0('<a target="_blank" href="',searchURL,rownames(DT_out3),'">',rownames(DT_out3),'</a>')
      output$rlnorm_assay_info <- renderDataTable({round(DT_out3,3)},options = list(scrollX = TRUE), escape = F)

      output$rlnorm_assay_info_dld <- renderUI({
        downloadButton(outputId = "rlnorm_assay_info_dld_btn",label = "Download rlog2-transformed RC!")
      })
      # output$norm_counts <- renderDataTable(rlog.norm.counts)
      # boxplots  of non -transformed  read  counts
      output$plot_un_read_counts <- renderPlot(
        boxplot(counts.sf_normalized, notch = TRUE, main = "untransformed  read  counts", ylab = "read  counts")
        )
      # box plots of log2 -transformed read counts
      output$plot_tr_read_counts <- renderPlot(
        boxplot(log.norm.counts, notch = TRUE, main = "log2 -transformed  read  counts", ylab = "log2(read  counts)")
      )
      # Normalized log2(read counts)
      output$plot_lnorm_counts <- renderPlot(
        plot(log.norm.counts [,1:2], cex=.1, main = "Normalized log2(read counts) before variance shrinkage")
      )
      # Normalized log2(read counts) after variance shrinkage
      output$plot_rlnorm_counts <- renderPlot(
        plot(rlog.norm.counts  [,1:2], cex=.1, main = "Normalized log2(read counts) after variance shrinkage")
      )
      # sequencing  depth  normalized  log2(read  counts) before variance shrinkage
      output$vsn_lnorm_counts <- renderPlot({
        msd_plot  <- meanSdPlot(log.norm.counts, ranks=FALSE, plot = FALSE)
        msd_plot$gg +
          ggtitle("log2(read  counts) before variance shrinkage") +
          ylab("standard  deviation")
      })
      # sequencing  depth  normalized  log2(read  counts) after variance shrinkage
      output$vsn_rlnorm_counts <- renderPlot({
        msd_plot  <- meanSdPlot(rlog.norm.counts ,ranks=FALSE, plot = FALSE)
        msd_plot$gg +
          ggtitle("rlog -transformed  read  counts after variance shrinkage") +
          ylab("standard  deviation")
      })

      output$global_rc <- renderUI({
        actionButton("global_rc_btn", "Explore global read counts!")
      })

      output$DESeq2_runDGE <- renderUI({
        actionButton("DESeq2_runDGE_btn", "Run DESeq2 DGE Analysis!")
      })

      output$edgeR_runDGE <- renderUI({
        actionButton("edgeR_runDGE_btn", "Run edgeR DGE Analysis!")
      })
      output$edgeR_type <- renderUI({
        selectInput("edgeR_model", "Select a model type:"
                    ,c("negative binomial GLM"="glmFit", "Quasi-likelihood GLM"="glmQLFit")
                    ,selected = "glmFit", multiple = FALSE,
                    selectize = TRUE, width = NULL, size = NULL)
      })

      output$limma_runDGE <- renderUI({
        actionButton("limma_runDGE_btn", "Run limma DGE Analysis!")
      })
      setProgress(value = 1, detail = "Done!")
    }
    ###################### end of norm_DESeq_data ######################

    output$norm_assay_info_dld_btn <- downloadHandler(
      filename = "NormalizedRawCounts.csv",
      content = function(file){
        write.csv(counts.sf_normalized,file)
      }
    )

    output$lnorm_assay_info_dld_btn <- downloadHandler(
      filename = "logTransfNormRawCounts.csv",
      content = function(file){
        write.csv(log.norm.counts,file)
      }
    )

    output$rlnorm_assay_info_dld_btn <- downloadHandler(
      filename = "reglogTransfNormRawCounts.csv",
      content = function(file){
        write.csv(rlog.norm.counts,file)
      }
    )
    #################################################

    observeEvent(input$normalization_btn,{
      withProgress(message = "Data Normalization in progress... please wait!",{
        norm_DESeq_data()
      })
      # norm_DESeq_data()
    })


    ##########################################################
    ### When "Explore global read counts!" button is clicked
    #########################################################
    global_rc_process <- function(){
      # cor()  calculates  the  correlation  between  columns  of a matrix
      setProgress(value = 1/4, detail = "Calculating Pearson correlation  between  columns  of a matrix...")
      distance.m_rlog  <- as.dist(1 - cor(rlog.norm.counts , method = "pearson" ))
      setProgress(value = 2/4, detail = "Principal Components Analysis...")
      pc <- prcomp(t(rlog.norm.counts))

      # output Pearson Correlation table
      output$DESeq2_pcorr <- renderDataTable({
        round(cor(counts.sf_normalized, method = "pearson"),3)
      },options = list(scrollX = TRUE))

      output$DESeq2_pcorr_dld <- renderUI({
        downloadButton("DESeq2_pcorr_dld_btn", "Download Correlation data!")
      })

      output$pcorr_mat <- renderPlotly({
        heatmaply::heatmaply(as.matrix(cor(counts.sf_normalized, method = "pearson")))
      })

      # plot Hierarchical clustering tree
      setProgress(value = 3/4, detail = "Hierarchical clustering tree...")
      output$DESeq2_hr_tree <- renderPlot(
        plot(hclust(distance.m_rlog), labels = colnames(rlog.norm.counts),
             main = "rlog  transformed  read  counts\ndistance: Pearson  correlation")
      )

      # plot PCA
      output$DESeq2_pca <- renderPlotly({
        P <- plotPCA2(pc,colorby=conditions, title="Rlog  transformed  counts")
        print(P)
      })
      setProgress(value = 1, detail = "Done!")
    }
    ############# end of global_rc_process() ###############

    output$DESeq2_pcorr_dld_btn <- downloadHandler(
      filename = "pairwiseCorr_results.csv",
      content = function(file){
        write.csv(cor(counts.sf_normalized, method = "pearson"),file)
      }
    )
    #########################################################

    observeEvent(input$global_rc_btn,{
      withProgress(message = "Processing in progress.. Please wait!",{
        global_rc_process()
      })
    })

    ##########################################################
    ### When "Run DESeq2 DGE Analysis!" button is clicked
    #########################################################
    observeEvent(input$DESeq2_runDGE_btn,{
      withProgress(message = "Performing DGE Analysis with DESeq2... Please wait!", value = 0,{
        n <- 16
        colData(dds)$condition <- relevel(colData(dds)$condition , "A")
        setProgress(value = 4/n, #message = "Fitting model, please wait...",
                    detail = "Sequencing  depth  normalization between  the  samples...")
        # dds2 <- DESeq(dds)
        # sequencing  depth  normalization between  the  samples
        dds2 <- estimateSizeFactors(dds)

        setProgress(value = 6/n, #message = "Fitting model, please wait...",
                    detail = "Gene-wise dispersion  estimates across  all  samples...")
        # gene -wise  dispersion  estimates across  all  samples
        dds2 <- estimateDispersions(dds2)

        setProgress(value = 9/n, #message = "Fitting model, please wait...",
                    detail = "Fitting a negative  binomial  GLM...")
        # this  fits a negative  binomial  GLM and applies  Wald  statistics  to each  gene
        dds2 <- nbinomWaldTest(dds2)

        setProgress(value = 10/n, detail = "Results filtering...")
        DGE.results  <- results(dds2, independentFiltering = TRUE , alpha = 0.05)

        setProgress(value = 11/n, detail = "Print Rendering...")
        output$DESeq2_pvalues <- renderPlot(
          hist(DGE.results$pvalue, col = "grey", border = "white", xlab = "", ylab = "", main = "frequencies  of p-values")
        )

        output$DESeq2_sumout <- renderPrint({
          summary(dds2)
        })

        setProgress(value = 12/n, detail = "MA plot processing...")
        output$DESeq2_plotMA <-renderPlot(
          DESeq2::plotMA(DGE.results , alpha = 0.05,  main = "A vs. B conditions", ylim = c(-4,4))
        )

        setProgress(value = 13/n, detail = "Volcano plot processing...")
        DGE_data <<- as.data.frame(DGE.results)

        results_dge[["DESeq2"]] <<- DGE_data

        d_choices <- colnames(DGE.results)

        output$DESeq2_volcano <- renderPlot({
          apval<-input$DESeq2_vp_pval
          lfc<-input$DESeq2_vp_lfc
          lim<-input$DESeq2_vp_limit
          with(DGE_data, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=lim#c(-10,10)
          ))
          with(subset(DGE_data, padj<apval ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
          with(subset(DGE_data, padj<apval & abs(log2FoldChange)>lfc), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
        })

        setProgress(value = 14/n, detail = "Table rendering...")
        DT_output<-DGE_data
        feature <- paste0('<a target="_blank" href="',searchURL,rownames(DT_output),'">',rownames(DT_output),'</a>')
        row.names(DT_output)<-feature

        output$DESeq2_dge_res <- renderDataTable({DT_output},options = list(scrollX = TRUE), escape = F)

        output$DESeq2_dge_res_dld <- renderUI({
          downloadButton(outputId = "DESeq2_dge_download",label = "Download DGE results Data!")
        })

        output$DESeq2_dge_res2 <- renderDataTable({DT_output},options = list(scrollX = TRUE), escape = F)

        setProgress(value = 15/n, detail = "UI rendering...")
        output$filterGenes <- renderUI({
          tagList(fluidRow(column(5
                                  ,selectInput("gprof_select", "Filter features by:",colnames(DGE_data),
                                               selected = "padj", multiple = FALSE,
                                               selectize = TRUE, width = "200px", size = NULL))
                           ,column(5
                                   ,conditionalPanel(condition='input.gprof_select=="log2FoldChange"',
                                                     sliderInput("gp_log2FoldChange","Define a range:",min = 0#floor(min(DGE_data$log2FoldChange, na.rm = T))
                                                                 ,max = abs(ceiling(max(DGE_data$log2FoldChange, na.rm = T)))
                                                                 ,step = 0.1#round(max(DGE_data$log2FoldChange, na.rm = T)/100,2)
                                                                 ,value = 2#c(floor(min(DGE_data$log2FoldChange, na.rm = T)),floor(min(DGE_data$log2FoldChange, na.rm = T))
                                                                            #+max(DGE_data$log2FoldChange, na.rm = T)/100)
                                                                 ,round = -2))
                                   ,conditionalPanel(condition='input.gprof_select=="lfcSE"',
                                                     sliderInput("gp_lfcSE","Define a range:"
                                                                 ,min = floor(min(DGE_data$lfcSE, na.rm = T))
                                                                 ,max = ceiling(max(DGE_data$lfcSE, na.rm = T))
                                                                 ,step = round(max(DGE_data$lfcSE, na.rm = T)/100,2)
                                                                 ,value = c(floor(min(DGE_data$lfcSE, na.rm = T)),floor(min(DGE_data$lfcSE, na.rm = T))
                                                                            +max(DGE_data$lfcSE, na.rm = T)/100)
                                                                 ,round = -2))
                                   ,conditionalPanel(condition='input.gprof_select=="baseMean"',
                                                     sliderInput("gp_baseMean","Define a range:"
                                                                 ,min = floor(min(DGE_data$baseMean, na.rm = T))
                                                                 ,max = ceiling(max(DGE_data$baseMean, na.rm = T))
                                                                 ,step = round(max(DGE_data$baseMean, na.rm = T)/100,2)
                                                                 ,value = c(floor(min(DGE_data$baseMean, na.rm = T)),floor(min(DGE_data$baseMean, na.rm = T))
                                                                            +max(DGE_data$baseMean, na.rm = T)/100)
                                                                 ,round = -2))
                                   ,conditionalPanel(condition='input.gprof_select=="pvalue"',
                                                     sliderInput("gp_pvalue","Define a range:",min = floor(min(DGE_data$pvalue, na.rm = T))
                                                                 ,max = ceiling(max(DGE_data$pvalue, na.rm = T))
                                                                 ,step = round(max(DGE_data$pvalue, na.rm = T)/100,2)
                                                                 ,value = c(0,0.05)
                                                                 ,round = -2))
                                   ,conditionalPanel(condition='input.gprof_select=="padj"',
                                                     sliderInput("gp_padj","Define a range:",min = floor(min(DGE_data$padj, na.rm = T))
                                                                 ,max = ceiling(max(DGE_data$padj, na.rm = T))
                                                                 ,step = round(max(DGE_data$padj, na.rm = T)/100,2)
                                                                 ,value = c(0,0.05)
                                                                 ,round = -2))
                                   ,conditionalPanel(condition='input.gprof_select=="stat"',
                                                     sliderInput("gp_stat","Define a range:"
                                                                 ,min = floor(min(DGE_data$stat, na.rm = T))
                                                                 ,max = ceiling(max(DGE_data$stat, na.rm = T))
                                                                 ,step = round(max(DGE_data$stat, na.rm = T)/100,2)
                                                                 ,value = c(floor(min(DGE_data$stat, na.rm = T)),floor(min(DGE_data$stat, na.rm = T))
                                                                            +max(DGE_data$stat, na.rm = T)/100)
                                                                 ,round = -2))
                           )
          )
          ,actionButton("filt_genes_btn","Filter features!")
          ,verbatimTextOutput("filt_genes")
          )
        })

        output$DESeq2_gprofile_par <- renderUI({
          tagList(
            fluidRow(hr(),
                     column(5,selectInput("gp_organism", "Choose organism"
                                          ,list("Homo sapiens"="hsapiens", "Mus musculus"="mmusculus", "Gorilla gorilla"="ggorilla"
                                                ,"Bonobo"="ppaniscus", "Ovis aries"="oaries","Other"="other"),
                                          selected = "hsapiens", multiple = FALSE,
                                          selectize = TRUE, width = "200px", size = NULL)
                            ,conditionalPanel(
                              condition = "input.gp_organism == 'other'",
                              textInput("other_gpo", "Specify a Profiler organism", value = "aaegypti",
                                        placeholder = "e.g: Enter 'aaegypti' for Aedes aegypti (https://biit.cs.ut.ee/gprofiler/help.cgi?help_id=64)"))
                     )
                     ,column(5,selectInput("correction_method", "Correction method"
                                           ,list("Analytical"="analytical", "gSCS"="gSCS", "fdr"="fdr"
                                                 ,"bonferroni"="bonferroni"),
                                           selected = "fdr", multiple = FALSE,
                                           selectize = TRUE, width = "200px", size = NULL))
                     ,column(5,selectInput("hier_filtering", "Hierarchical filtering strength"
                                           ,list("None"="none", "Moderate"="moderate", "Strong"="strong"),
                                           selected = "none", multiple = FALSE,
                                           selectize = TRUE, width = "200px", size = NULL))
                     ,column(5,selectInput("domain_size", "Statistical domain size"
                                           ,list("Annotated"="annotated", "Known"="known"),
                                           selected = "annotated", multiple = FALSE,
                                           selectize = TRUE, width = "200px", size = NULL))
                     ,column(5, materialSwitch("sort_by_structure","Sort by structure",value=T
                                               ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
                     ,column(5, materialSwitch("ordered_query","Ordered Query",value=F
                                               ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
                     ,column(5, materialSwitch("significant","Significant Only",value=T
                                               ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
                     ,column(5, materialSwitch("exclude_iea","Exclude electronic annotations (IEA)",value=F
                                               ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
            )
            ,uiOutput("perf_gprof_btn"),hr()
            ,uiOutput("gprofileLink"), hr()
            ,uiOutput("webgestaltLink"), hr()
            ,dataTableOutput("gprofile_table"), hr()
            ,uiOutput("gprofile_res_dld")
          )
        })

        output$DESeq2_settings <- renderUI({
          tagList(
            selectInput("DESeq2_sort_by","Sort by:", d_choices, selected="padj",width = "200px")
            ,conditionalPanel(condition='input.DESeq2_sort_by=="log2FoldChange"',
                              fluidRow(column(6
                                              ,sliderInput("d_log2FoldChange","Define a range:",min = 0#floor(min(DGE_data$log2FoldChange))
                                                           ,max = abs(ceiling(max(DGE_data$log2FoldChange)))
                                                           ,step = 0.1#round(max(DGE_data$log2FoldChange)/100,2)
                                                           ,value = 2#c(floor(min(DGE_data$log2FoldChange)),floor(min(DGE_data$log2FoldChange))
                                                                      #+max(DGE_data$log2FoldChange)/100)
                                                           ,round = -2))
                              )
            )
            ,conditionalPanel(condition='input.DESeq2_sort_by=="lfcSE"',
                              fluidRow(column(6
                                              ,sliderInput("d_lfcSE","Define a range:",min = floor(min(DGE_data$lfcSE))
                                                           ,max = ceiling(max(DGE_data$lfcSE))
                                                           ,step = round(max(DGE_data$lfcSE)/100,2)
                                                           ,value = c(floor(min(DGE_data$lfcSE)),floor(min(DGE_data$lfcSE))
                                                                      +max(DGE_data$lfcSE)/100)
                                                           ,round = -2))
                              )
            )
            ,conditionalPanel(condition='input.DESeq2_sort_by=="baseMean"',
                              fluidRow(column(6
                                              ,sliderInput("d_DGE_data","Define a range:",min = floor(min(DGE_data$baseMean))
                                                           ,max = ceiling(max(DGE_data$baseMean))
                                                           ,step = round(max(DGE_data$baseMean)/100,2)
                                                           ,value = c(floor(min(DGE_data$baseMean)),floor(min(DGE_data$baseMean))
                                                                      +max(DGE_data$baseMean)/100)
                                                           ,round = -2))
                              )
            )
            ,conditionalPanel(condition='input.DESeq2_sort_by=="pvalue"',
                              fluidRow(column(6
                                              ,sliderInput("d_pvalue","Define a range:",min = 0#floor(min(DGE.results_limma$P.Value))
                                                           ,max = 1#ceiling(max(DGE.results_limma$P.Value))
                                                           ,step = 0.01#round(max(DGE.results_limma$P.Value)/100,2)
                                                           ,value = c(0,0.05)#c(floor(min(DGE.results_limma$P.Value)),floor(min(DGE.results_limma$P.Value))
                                                           #+max(DGE.results_limma$P.Value)/100)
                                                           ,round = -2))
                              )
            )
            ,conditionalPanel(condition='input.DESeq2_sort_by=="padj"',
                              fluidRow(column(6
                                              ,sliderInput("d_padj","Define a range:",min = 0#floor(min(DGE.results_limma$adj.P.Val))
                                                           ,max = 1#ceiling(max(DGE.results_limma$adj.P.Val))
                                                           ,step = 0.01#round(max(DGE.results_limma$adj.P.Val)/100,2)
                                                           ,value = c(0,0.05)#c(floor(min(DGE.results_limma$adj.P.Val)),floor(min(DGE.results_limma$adj.P.Val))
                                                           #+max(DGE.results_limma$adj.P.Val)/100)
                                                           ,round = -2))
                              )
            )
            ,conditionalPanel(condition='input.DESeq2_sort_by=="stat"',
                              fluidRow(column(6
                                              ,sliderInput("d_stat","Define a range:",min = floor(min(DGE_data$stat))
                                                           ,max = ceiling(max(DGE_data$stat))
                                                           ,step = round(max(DGE_data$stat)/100,2)
                                                           ,value = c(floor(min(DGE_data$stat)),floor(min(DGE_data$stat))
                                                                      +max(DGE_data$stat)/100)
                                                           ,round = -2))
                              )
            )
            ,actionButton("d_filter","Filter features!"),br(), hr()
            ,verbatimTextOutput("DESeq2_num_dge"), hr()
            ,uiOutput("d_heatmap_settings"), br()
            ,hr()
          )
      })
        setProgress(value = 1, detail = "Done!")
      })
    })
    ##################################################################################################################

    output$DESeq2_dge_download <- downloadHandler(
        filename = "DESeq2_DGE_results.csv",
        content = function(file){
          write.csv(DGE_data,
                    file)
        }
      )

    observeEvent(input$d_filter,{
      if(input$DESeq2_sort_by=="log2FoldChange"){
        lfc <- eval(parse(text = paste0("input$d_",input$DESeq2_sort_by)))
        DGEgenes_DESeq2  <<- rownames(subset(DGE_data, abs(eval(parse(text = input$DESeq2_sort_by))) > lfc))
      }else{
        # DGE.results_lima.sorted  <- DGE.results_limma[order(DGE.results_limma$adj.P.val), ]
        minn <- eval(parse(text = paste0("input$d_",input$DESeq2_sort_by)))[1]
        maxx <- eval(parse(text = paste0("input$d_",input$DESeq2_sort_by)))[2]
        # list(minn,maxx)
        DGEgenes_DESeq2  <<- rownames(subset(DGE_data, eval(parse(text = input$DESeq2_sort_by))>=minn &
                                               eval(parse(text = input$DESeq2_sort_by))<=maxx))
      }
      if(length(DGEgenes_DESeq2)>0){
        output$d_heatmap_settings <- renderUI({
          tagList(fluidRow(column(2,materialSwitch("d_rdl", label="Row_dend_left:", value=F
                                                   ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
                           ,column(2,materialSwitch("d_Rowv", label="Rowv", value=T
                                                    ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
                           ,column(2,materialSwitch("d_Colv", label="Colv", value=T
                                                    ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
                           , column(2,materialSwitch("d_colorbar", label="Colorbar:", value=F
                                                     ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
          )
          ,fluidRow(column(2, pickerInput(inputId = "d_colors", label = "Colors:"
                                          ,choices=c("BrBG","Spectral","RdYlBu","RdYlGn","RdGy","RdBu","PuOr","PRGn","PiYG")
                                          , selected = "PiYG"))
                    ,column(2, pickerInput(inputId = "d_dist_method", label = "Distance method:"
                                           ,choices = c("euclidean", "maximum","manhattan","canberra","binary","minkowski")
                                           ,selected="euclidean"))
                    ,column(2, pickerInput(inputId = "d_hclust_method", label = "Hclust method:"
                                           ,choices = c("complete", "ward.D","ward.D2","single","average","mcquitty","median","centroid")
                                           ,selected="average"))
                    ,column(2, pickerInput(inputId = "d_dendogram", label = "Dendogram:", choices = c("both", "row", "column", "none")
                                           ,selected="both"))
          )
          ,fluidRow(column(2, pickerInput(inputId = "d_plot_method", label = "Plot method:", choices = c("plotly", "ggplot")
                                           ,selected="ggplot"))
                    ,column(2,numericInput("d_k_col", "k_col:", 1, min = 1, max = 100, step = 1))
                    ,column(2,numericInput("d_k_row", "k_row:", 1, min = 1, max = 100, step = 1))
          )
          ,hr()
          ,actionButton("d_heatmap","Generate Heatmap!")
          )
        })
      }
      output$DESeq2_num_dge <- renderPrint({
        list("DGE"=paste(length(DGEgenes_DESeq2),"differentially expressed features"),"Differentially Expressed features"=DGEgenes_DESeq2)
      })
    })

    observeEvent(input$d_heatmap,{
      # extract  the  normalized  read  counts  for DE genes  into a matrix
      hm.mat_DGEgenes  <- log.norm.counts[DGEgenes_DESeq2, ]
      output$DESeq2_heatmap <- renderPlotly({
        withProgress(message = "Generating heatmap...Please wait!", value = 0.8,{
          isolate(heatmaply::heatmaply(as.matrix(hm.mat_DGEgenes), xlab = "Samples", ylab = "Features",
                    scale = "row", row_dend_left = input$d_rdl, plot_method = input$d_plot_method, dendrogram=input$d_dendogram,
                    Rowv = input$d_Rowv, Colv = input$d_Colv, #distfun = input$d_distfun, #hclustfun = "average",
                    colorbar=input$d_colorbar, dist_method = input$d_dist_method, hclust_method = input$d_hclust_method,
                    colors = eval(parse(text=input$d_colors)),
                    k_col = input$d_k_col, k_row = input$d_k_row,
                    margins = c(60,100,40,20)))
        })
      })
    })

    ######################################

    observeEvent(input$filt_genes_btn,{
      output$filt_genes <- renderPrint({
        if(input$gprof_select=="log2FoldChange"){
          gprofile_genes <<- rownames(subset(DGE_data
                                             ,abs(eval(parse(text = input$gprof_select))) >
                                               eval(parse(text = paste0("input$gp_",input$gprof_select)))))
        }else{
          gprofile_genes <<- rownames(subset(DGE_data
                                             ,eval(parse(text = input$gprof_select)) <= eval(parse(text = paste0("input$gp_",input$gprof_select)))[2] &
                                               eval(parse(text = input$gprof_select)) >= eval(parse(text = paste0("input$gp_",input$gprof_select)))[1]))
        }
        if(any(startsWith(as.character(gprofile_genes),"ENS"))==T)
          gprofile_genes <<- unique(apply(as.data.frame(gprofile_genes),MARGIN = 1,
                   FUN = function(x){if(startsWith(as.character(x),"ENS")){return(strsplit(as.character(x),fixed = T,split = '.')[[1]][1])}
                   }))
        list("Number of filtered features:"=paste(length(gprofile_genes),"filtered features!"),"Features:"=gprofile_genes)
      })
      output$perf_gprof_btn <- renderUI(actionButton("gprofile_btn","Perform gProfile pathway Analysis!"))
    })

    #########################################################
    observeEvent(input$gprofile_btn,{
      withProgress(message = "Gene Pathway analysis in progress...",{
        gpr_organism <- input$gp_organism
        if(gpr_organism=='other'){gpr_organism<-input$other_gpo}
        setProgress(value = 2/6, detail = "gProfiler pathway Analysis...")
        gprof_table <<- as.data.frame(gprofiler(gprofile_genes
                                                ,organism=gpr_organism, sort_by_structure = input$sort_by_structure, ordered_query = input$ordered_query
                                                ,significant = input$significant, exclude_iea = input$exclude_iea, correction_method = input$correction_method
                                                ,hier_filtering = input$hier_filtering, domain_size = input$domain_size
        ))
        query=paste(gprofile_genes,collapse = "%0A")
        setProgress(value = 3/6, detail = "gProfiler link generation...")
        gprof_url <- a("Click me to be directed to gprofile Website!"
                       ,href=paste0("https://biit.cs.ut.ee/gprofiler/gost"#"https://biit.cs.ut.ee/gprofiler/index.cgi"
                                    ,'?query=',query,'&organism=',gpr_organism
                                    ,"&significant=",as.numeric(input$significant)
                                    ,"&ordered=",tolower(as.character(input$ordered_query))#"&ordered_query=",as.numeric(input$ordered_query)
                                    ,"&no_iea=",as.numeric(input$exclude_iea)
                                    ,"&domain_scope=",input$domain_size#"&domain_size_type=",input$domain_size
                                    # ,"&sort_by_structure=",as.numeric(input$sort_by_structure)
                                    ,"&significance_threshold_method=",input$correction_method#"&threshold_algo=",input$correction_method
                       )
                       ,target="_blank")

        setProgress(value = 4/6, detail = "WebGestalt link generation...")
        wgst_url <- a("Click me to be directed to WebGestalt Website!"
                      ,href=paste0("http://www.webgestalt.org/option.php",'?gene_list=',query,'&organism=',gpr_organism
                                   ,"&enrich_method=ORA&fdr_method=BH&enriched_database_category=geneontology"
                                   ,"&enriched_database_name=Biological_Process_noRedundant&sig_method=top&sig_value=5"
                                   ,"&max_num=200&id_type=genesymbol"
                                   ,"&ref_set=genome")
                      ,target="_blank")
        setProgress(value = 5/6, detail = "gProfiler Table rendering...")
        output$gprofile_table <- renderDataTable({gprof_table},options = list(scrollX = TRUE))

        output$gprofileLink <- renderUI({
          tagList("gprofile Link:", gprof_url)
        })

        output$webgestaltLink <- renderUI({
          tagList("WebGestalt Link:", wgst_url)
        })

        output$gprofile_res_dld<-renderUI({
          downloadButton(outputId = "gprofiler_download",label = "Download gProfile results!")
        })
        setProgress(value = 1, detail = "Done!")
      })
    })
    #########################################################

    output$gprofiler_download <- downloadHandler(
      filename = "gProfileR_results.csv",
      content = function(file){
        write.csv(gprof_table,
                  file)
      }
    )

    ##########################################################
    ### When "Run limma DGE Analysis!" button is clicked
    #########################################################
    observeEvent(input$limma_runDGE_btn,{
      withProgress(message = "Performing DGE Analysis with limma... Please wait!",{
        output$limma_sumout <-renderUI({
          withSpinner(verbatimTextOutput("limma_DGE_results_sum"),size=1,proxy.height=200,type=sample(c(1,4,5:8),1))
        })

        sample_info.edger <- relevel(conditions , ref = "A")
        # limma  also  needs a design  matrix , just  like  edgeR
        design  <- model.matrix(~sample_info.edger)
        # transform  the  count  data to log2 -counts -per -million  and  estimate
        # the mean -variance  relationship , which is used to  compute  weights
        # for  each  count  -- this is  supposed  to make  the  read  counts
        # amenable  to be used  with  linear  models
        design  <- model.matrix(~sample_info.edger)
        edgeR.DGElist <- DGEList(counts = txi$counts , group = sample_info.edger)
        # keep <- rowSums( cpm(edgeR.DGElist) >= 1) >= 5
        # edgeR.DGElist <- edgeR.DGElist[keep ,]
        edgeR.DGElist <- calcNormFactors(edgeR.DGElist , method = "TMM")
        rownames(design) <- colnames(edgeR.DGElist)

        setProgress(value = 1/8, detail = "Transform count data to log2-counts per million...")
        voomTransformed  <- voom(edgeR.DGElist , design , plot=F)

        # fit a linear  model  for  each  gene
        setProgress(value = 3/8, detail = "Fit a linear  model  for  each  gene...")
        voomed.fitted  <- lmFit(voomTransformed , design = design)

        output$limma_DGE_results_sum <- renderPrint({
          summary(voomed.fitted)
        })

        # compute  moderated t-statistics , moderated F-statistics ,
        # and log -odds of  differential  expression
        setProgress(value = 4/11, detail = "Computing statistics for differential expression...")
        voomed.fitted  <- eBayes(voomed.fitted)

        # extract  gene  list  with  logFC  and  statistical  measures
        colnames(design) # check  how  the  coefficient  is named

        setProgress(value = 5/8, detail = "Extracting gene list...")
        DGE.results_limma  <<- topTable(voomed.fitted , coef = "sample_info.edgerB",
                                        number = Inf , adjust.method = "BH",
                                        sort.by = "logFC")

        results_dge[["limma"]] <<- DGE.results_limma

        setProgress(value = 6/8, detail = "Rendering results table...")
        DT_output <- DGE.results_limma
        feature <- paste0('<a target="_blank" href="',searchURL,rownames(DT_output),'">',rownames(DT_output),'</a>')
        row.names(DT_output)<-feature

        output$limma_dge_res <- renderDataTable({DT_output},options = list(scrollX = TRUE), escape = F)
        output$limma_dge_res2 <- renderDataTable({DT_output},options = list(scrollX = TRUE), escape = F)

        setProgress(value = 7/8, detail = "UI rendering...")
        output$limma_pvalues <- renderPlot(
          hist(DGE.results_limma$P.Value, col = "grey", border = "white", xlab = "", ylab = "", main = "frequencies  of p-values")
        )

        output$limma_plotMA <-renderPlot(
          limma::plotMA(DGE.results_limma , alpha = 0.05,  main = "A vs. B conditions")
        )

        output$limma_volcano <- renderPlot({
          apval<-input$limma_vp_pval
          lfc<-input$limma_vp_lfc
          lim<-input$limma_vp_limit
          with(DGE.results_limma, plot(logFC, -log10(adj.P.Val), pch=20, main="Volcano plot", xlim=lim#c(-10,10)
          ))
          with(subset(DGE.results_limma, adj.P.Val<apval ), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
          with(subset(DGE.results_limma, adj.P.Val<apval & abs(logFC)>lfc), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
        })

        output$limma_dge_res_dld <- renderUI({
          downloadButton(outputId = "limma_dge_download",label = "Download limma DGE results!")
        })

        l_choices <- colnames(DGE.results_limma)

        output$limma_settings <- renderUI({
          tagList(
            selectInput("limma_sort_by","Sort by:", l_choices, selected="adj.P.Val",width = "200px")
            ,conditionalPanel(condition='input.limma_sort_by=="logFC"',
                              fluidRow(column(6
                                              ,sliderInput("s_logFC","Define a range:",min = 0#floor(min(DGE.results_limma$logFC))
                                                           ,max = abs(ceiling(max(DGE.results_limma$logFC)))
                                                           ,step = 0.1#round(max(DGE.results_limma$logFC)/100,2)
                                                           ,value = 2#c(floor(min(DGE.results_limma$logFC)),floor(min(DGE.results_limma$logFC))
                                                                      #+max(DGE.results_limma$logFC)/100)
                                                           ,round = -2))
                              )
            )
            ,conditionalPanel(condition='input.limma_sort_by=="AveExpr"',
                              fluidRow(column(6
                                              ,sliderInput("s_AveExpr","Define a range:",min = floor(min(DGE.results_limma$AveExpr))
                                                           ,max = ceiling(max(DGE.results_limma$AveExpr))
                                                           ,step = round(max(DGE.results_limma$AveExpr)/100,2)
                                                           ,value = c(floor(min(DGE.results_limma$AveExpr)),floor(min(DGE.results_limma$AveExpr))
                                                                      +max(DGE.results_limma$AveExpr)/100)
                                                           ,round = -2))
                              )
            )
            ,conditionalPanel(condition='input.limma_sort_by=="t"',
                              fluidRow(column(6
                                              ,sliderInput("s_t","Define a range:",min = floor(min(DGE.results_limma$t))
                                                           ,max = ceiling(max(DGE.results_limma$t))
                                                           ,step = round(max(DGE.results_limma$t)/100,2)
                                                           ,value = c(floor(min(DGE.results_limma$t)),floor(min(DGE.results_limma$t))
                                                                      +max(DGE.results_limma$t)/100)
                                                           ,round = -2))
                              )
            )
            ,conditionalPanel(condition='input.limma_sort_by=="P.Value"',
                              fluidRow(column(6
                                              ,sliderInput("s_P.Value","Define a range:",min = 0#floor(min(DGE.results_limma$P.Value))
                                                           ,max = 1#ceiling(max(DGE.results_limma$P.Value))
                                                           ,step = 0.01#round(max(DGE.results_limma$P.Value)/100,2)
                                                           ,value = c(0,0.05)#c(floor(min(DGE.results_limma$P.Value)),floor(min(DGE.results_limma$P.Value))
                                                           #+max(DGE.results_limma$P.Value)/100)
                                                           ,round = -2))
                              )
            )
            ,conditionalPanel(condition='input.limma_sort_by=="adj.P.Val"',
                              fluidRow(column(6
                                              ,sliderInput("s_adj.P.Val","Define a range:",min = 0#floor(min(DGE.results_limma$adj.P.Val))
                                                           ,max = 1#ceiling(max(DGE.results_limma$adj.P.Val))
                                                           ,step = 0.01#round(max(DGE.results_limma$adj.P.Val)/100,2)
                                                           ,value = c(0,0.05)#c(floor(min(DGE.results_limma$adj.P.Val)),floor(min(DGE.results_limma$adj.P.Val))
                                                           #+max(DGE.results_limma$adj.P.Val)/100)
                                                           ,round = -2))
                              )
            )
            ,conditionalPanel(condition='input.limma_sort_by=="B"',
                              fluidRow(column(6
                                              ,sliderInput("s_B","Define a range:",min = floor(min(DGE.results_limma$B))
                                                           ,max = ceiling(max(DGE.results_limma$B))
                                                           ,step = round(max(DGE.results_limma$B)/100,2)
                                                           ,value = c(floor(min(DGE.results_limma$B)),floor(min(DGE.results_limma$B))
                                                                      +max(DGE.results_limma$B)/100)
                                                           ,round = -2))
                              )
            )
            ,actionButton("l_filter","Filter features!"),br(), hr()
            ,verbatimTextOutput("limma_num_dge"), hr()
            ,uiOutput("l_heatmap_settings"), br()
            ,hr()
          )
        })

        ##############################limma_filter_genes(on limma_gprofile)
        output$limma_filterGenes <- renderUI({
          tagList(fluidRow(column(5
                                  ,selectInput("limma_gprof_select", "Filter features by:",colnames(DGE.results_limma),
                                               # ,list("P-value"="pvalue", "Adjusted P-value"="padj", "Fold Change"="log2FoldChange"
                                               #       ,"Base Mean"="baseMean", "lfcSE"="lfcSE", "stat"="stat"),
                                               selected = "logFC", multiple = FALSE,
                                               selectize = TRUE, width = "200px", size = NULL))
                           ,column(5
                                   ,conditionalPanel(condition='input.limma_gprof_select=="logFC"',
                                                     sliderInput("limma_gp_logFC","Define a range:"
                                                                 ,min = 0#floor(min(DGE.results_limma$logFC, na.rm = T))
                                                                 ,max = abs(ceiling(max(DGE.results_limma$logFC, na.rm = T)))
                                                                 ,step = 0.1#round(max(DGE.results_limma$logFC, na.rm = T)/100,2)
                                                                 ,value = 2#c(floor(min(DGE.results_limma$logFC, na.rm = T)),floor(min(DGE.results_limma$logFC, na.rm = T))
                                                                            #+max(DGE.results_limma$logFC, na.rm = T)/100)
                                                                 ,round = -2))
                                   ,conditionalPanel(condition='input.limma_gprof_select=="AveExpr"',
                                                     sliderInput("limma_gp_AveExpr","Define a range:"
                                                                 ,min = floor(min(DGE.results_limma$AveExpr, na.rm = T))
                                                                 ,max = ceiling(max(DGE.results_limma$AveExpr, na.rm = T))
                                                                 ,step = round(max(DGE.results_limma$AveExpr, na.rm = T)/100,2)
                                                                 ,value = c(floor(min(DGE.results_limma$AveExpr, na.rm = T)),floor(min(DGE.results_limma$AveExpr, na.rm = T))
                                                                            +max(DGE.results_limma$AveExpr, na.rm = T)/100)
                                                                 ,round = -2))
                                   ,conditionalPanel(condition='input.limma_gprof_select=="t"',
                                                     sliderInput("limma_gp_t","Define a range:",min = floor(min(DGE.results_limma$t, na.rm = T))
                                                                 ,max = ceiling(max(DGE.results_limma$t, na.rm = T))
                                                                 ,step = round(max(DGE.results_limma$t, na.rm = T)/100,2)
                                                                 ,value = c(floor(min(DGE.results_limma$t, na.rm = T)),floor(min(DGE.results_limma$t, na.rm = T))
                                                                            +max(DGE.results_limma$t)/100)
                                                                 ,round = -2))
                                   ,conditionalPanel(condition='input.limma_gprof_select=="P.Value"',
                                                     sliderInput("limma_gp_P.Value","Define a range:"
                                                                 ,min = floor(min(DGE.results_limma$P.Value, na.rm = T))
                                                                 ,max = ceiling(max(DGE.results_limma$P.Value, na.rm = T))
                                                                 ,step = round(max(DGE.results_limma$P.Value, na.rm = T)/100,2)
                                                                 ,value = c(0,0.05)
                                                                 ,round = -2))
                                   ,conditionalPanel(condition='input.limma_gprof_select=="adj.P.Val"',
                                                     sliderInput("limma_gp_adj.P.Val","Define a range:"
                                                                 ,min = floor(min(DGE.results_limma$adj.P.Val, na.rm = T))
                                                                 ,max = ceiling(max(DGE.results_limma$adj.P.Val, na.rm = T))
                                                                 ,step = round(max(DGE.results_limma$adj.P.Val, na.rm = T)/100,2)
                                                                 ,value = c(0,0.05)
                                                                 ,round = -2))
                                   ,conditionalPanel(condition='input.limma_gprof_select=="B"',
                                                     sliderInput("limma_gp_B","Define a range:",min = floor(min(DGE.results_limma$B))
                                                                 ,max = ceiling(max(DGE.results_limma$B))
                                                                 ,step = round(max(DGE.results_limma$B)/100,2)
                                                                 ,value = c(floor(min(DGE.results_limma$B)),floor(min(DGE.results_limma$B))
                                                                            +max(DGE.results_limma$B)/100)
                                                                 ,round = -2))
                           )
          )
          ,actionButton("limma_filt_genes_btn","Filter features!")
          ,verbatimTextOutput("limma_filt_genes")
          )
        })
        #############################


        output$limma_gprofile_par <- renderUI({
          tagList(
            fluidRow(hr(),
                     column(5,selectInput("lgp_organism", "Choose organism"
                                          ,list("Homo sapiens"="hsapiens", "Mus musculus"="mmusculus", "Gorilla gorilla"="ggorilla"
                                                ,"Bonobo"="ppaniscus", "Ovis aries"="oaries","Other"="other"),
                                          selected = "hsapiens", multiple = FALSE,
                                          selectize = TRUE, width = "200px", size = NULL)
                            ,conditionalPanel(
                              condition = "input.lgp_organism == 'other'",
                              textInput("l_other_gpo", "Specify a Profiler organism", value = "aaegypti",
                                        placeholder = "e.g: Enter 'aaegypti' for Aedes aegypti (https://biit.cs.ut.ee/gprofiler/help.cgi?help_id=64)"))
                     )
                     ,column(5,selectInput("l_correction_method", "Correction method"
                                           ,list("Analytical"="analytical", "gSCS"="gSCS", "fdr"="fdr"
                                                 ,"bonferroni"="bonferroni"),
                                           selected = "fdr", multiple = FALSE,
                                           selectize = TRUE, width = "200px", size = NULL))
                     ,column(5,selectInput("l_hier_filtering", "Hierarchical filtering strength"
                                           ,list("None"="none", "Moderate"="moderate", "Strong"="strong"),
                                           selected = "none", multiple = FALSE,
                                           selectize = TRUE, width = "200px", size = NULL))
                     ,column(5,selectInput("l_domain_size", "Statistical domain size"
                                           ,list("Annotated"="annotated", "Known"="known"),
                                           selected = "annotated", multiple = FALSE,
                                           selectize = TRUE, width = "200px", size = NULL))
                     ,column(5, materialSwitch("l_sort_by_structure","Sort by structure",value=T
                                               ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
                     ,column(5, materialSwitch("l_ordered_query","Ordered Query",value=F
                                               ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
                     ,column(5, materialSwitch("l_significant","Significant Only",value=T
                                               ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
                     ,column(5, materialSwitch("l_exclude_iea","Exclude electronic annotations (IEA)",value=F
                                               ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
            )
            ,uiOutput("l_perf_gprof_btn"),hr()
            ,uiOutput("l_gprofileLink"), hr()
            ,uiOutput("l_webgestaltLink"), hr()
            ,dataTableOutput("l_gprofile_table"), hr()
            ,uiOutput("l_gprofile_res_dld")
          )
        })

        setProgress(value = 1, detail = "Done!")
      })
    })
    ##################################################################################################################

    output$limma_dge_download <- downloadHandler(
      filename = "limma_DGE_results.csv",
      content = function(file){
        write.csv(DGE.results_limma, file)
      }
    )

    observeEvent(input$l_filter,{
      if(input$limma_sort_by=="logFC"){
        lfc <- eval(parse(text = paste0("input$s_",input$limma_sort_by)))
        DGEgenes_lima  <<- rownames(subset(DGE.results_limma, abs(eval(parse(text = input$limma_sort_by))) > lfc))
      }else{
        # DGE.results_lima.sorted  <- DGE.results_limma[order(DGE.results_limma$adj.P.val), ]
        minn <- eval(parse(text = paste0("input$s_",input$limma_sort_by)))[1]
        maxx <- eval(parse(text = paste0("input$s_",input$limma_sort_by)))[2]
        # list(minn,maxx)
        DGEgenes_lima  <<- rownames(subset(DGE.results_limma, eval(parse(text = input$limma_sort_by))>=minn &
                                             eval(parse(text = input$limma_sort_by))<=maxx))
      }
      if(length(DGEgenes_lima)>0){
        output$l_heatmap_settings <- renderUI({
          tagList(fluidRow(column(2, materialSwitch("l_rdl", label="Row_dend_left:", value=F
                                                   ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
              ,column(2, materialSwitch("l_Rowv", label="Rowv", value=T
                                       ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
              ,column(2, materialSwitch("l_Colv", label="Colv", value=T
                                       ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
              , column(2, materialSwitch("l_colorbar", label="Colorbar:", value=F
                                        ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
            )
            ,fluidRow(column(2, pickerInput(inputId = "l_colors", label = "Colors:"
                                           ,choices=c("BrBG","Spectral","RdYlBu","RdYlGn","RdGy","RdBu","PuOr","PRGn","PiYG")
                                           , selected = "PiYG"))
                      # ,column(2, pickerInput(inputId = "l_distfun", label = "Distance function:", choices = c("pearson", "spearman","kendall")
                      #                        ,selected="pearson"))
                      ,column(2, pickerInput(inputId = "l_dist_method", label = "Distance method:"
                                             ,choices = c("euclidean", "maximum","manhattan","canberra","binary","minkowski")
                                             ,selected="euclidean"))
                      ,column(2, pickerInput(inputId = "l_hclust_method", label = "Hclust method:"
                                             ,choices = c("complete", "ward.D","ward.D2","single","average","mcquitty","median","centroid")
                                             ,selected="average"))
                      ,column(2, pickerInput(inputId = "l_dendogram", label = "Dendogram:", choices = c("both", "row", "column", "none")
                                             ,selected="both"))
            )
            ,fluidRow(column(2, pickerInput(inputId = "l_plot_method", label = "Plot method:", choices = c("plotly", "ggplot")
                                             ,selected="ggplot"))
                      ,column(2,numericInput("l_k_col", "k_col:", 1, min = 1, max = 100, step = 1))
                      ,column(2,numericInput("l_k_row", "k_row:", 1, min = 1, max = 100, step = 1))
            )
            ,hr()
            ,actionButton("l_heatmap","Generate Heatmap!")
          )
        })
      }
      output$limma_num_dge <- renderPrint({
        list("DGE"=paste(length(DGEgenes_lima),"differentially expressed features"),"Differentially Expressed features"=DGEgenes_lima)
      })
    })

    observeEvent(input$l_heatmap,{
      # extract  the  normalized  read  counts  for DE genes  into a matrix
      hm.mat_DGEgenes.lima  <- log.norm.counts[DGEgenes_lima, ]
      output$limma_heatmap <- renderPlotly({
        withProgress(message = "Generating heatmap...Please wait!", value = 0.8,{
          isolate(heatmaply::heatmaply(as.matrix(hm.mat_DGEgenes.lima), xlab = "Samples", ylab = "Features",
                    scale = "row", row_dend_left = input$l_rdl, plot_method = input$l_plot_method, dendrogram=input$l_dendogram,
                    Rowv = input$l_Rowv, Colv = input$l_Colv, #distfun = input$l_distfun, #hclustfun = "average",
                    colorbar=input$l_colorbar, dist_method = input$l_dist_method, hclust_method = input$l_hclust_method,
                    colors = eval(parse(text=input$l_colors)),
                    k_col = input$l_k_col, k_row = input$l_k_row,
                    margins = c(60,100,40,20)))
        })
      })
    })
    ##################################################################################################################

    output$gprofiler_download <- downloadHandler(
      filename = "gProfileR_results.csv",
      content = function(file){
        write.csv(gprof_table, file)
      }
    )

    ##################################################################################################################
    observeEvent(input$limma_filt_genes_btn,{
      output$limma_filt_genes <- renderPrint({
        if(input$limma_gprof_select=="logFC"){
          limma_gprofile_genes <<- rownames(subset(DGE.results_limma
                                                   ,abs(eval(parse(text = input$limma_gprof_select))) >
                                                     eval(parse(text = paste0("input$limma_gp_",input$limma_gprof_select)))))
        }else{
          limma_gprofile_genes <<- rownames(subset(DGE.results_limma
                                                   ,eval(parse(text = input$limma_gprof_select)) <= eval(parse(text = paste0("input$limma_gp_",input$limma_gprof_select)))[2] &
                                                     eval(parse(text = input$limma_gprof_select)) >= eval(parse(text = paste0("input$limma_gp_",input$limma_gprof_select)))[1]))
        }
        if(any(startsWith(as.character(limma_gprofile_genes),"ENS"))==T)
          limma_gprofile_genes<<- unique(apply(as.data.frame(limma_gprofile_genes),MARGIN = 1,
                                     FUN = function(x){if(startsWith(as.character(x),"ENS")){return(strsplit(as.character(x),fixed = T,split = '.')[[1]][1])}
                                     }))
        list("Number of filtered features:"=paste(length(limma_gprofile_genes),"filtered features!"),"Features:"=limma_gprofile_genes)
      })
      output$l_perf_gprof_btn <- renderUI(actionButton("limma_gprofile_btn","Perform gProfile pathway Analysis!"))
    })

    #########################################################
    observeEvent(input$limma_gprofile_btn,{
      withProgress(message = "Gene Pathway analysis in progress...",{
        lgpr_organism <- input$lgp_organism
        if(lgpr_organism=='other'){lgpr_organism<-input$l_other_gpo}
        setProgress(value = 2/6, detail = "gProfiler pathway Analysis...")
        limma_gprof_table <<- as.data.frame(gprofiler(limma_gprofile_genes
                                                      ,organism=lgpr_organism, sort_by_structure = input$l_sort_by_structure, ordered_query = input$l_ordered_query
                                                      ,significant = input$l_significant, exclude_iea = input$l_exclude_iea, correction_method = input$l_correction_method
                                                      ,hier_filtering = input$l_hier_filtering, domain_size = input$l_domain_size
        ))
        query=paste(limma_gprofile_genes,collapse = "%0A")
        setProgress(value = 3/6, detail = "gProfiler link generation...")
        gprof_url <- a("Click me to be directed to gprofile Website!"
                       ,href=paste0("https://biit.cs.ut.ee/gprofiler/gost"#"https://biit.cs.ut.ee/gprofiler/index.cgi"
                                    ,'?query=',query,'&organism=',lgpr_organism
                                    ,"&significant=",as.numeric(input$l_significant)
                                    ,"&ordered=",tolower(as.character(input$l_ordered_query))#"&ordered_query=",as.numeric(input$l_ordered_query)
                                    ,"&no_iea=",as.numeric(input$l_exclude_iea)
                                    ,"&domain_scope=",input$l_domain_size#"&domain_size_type=",input$l_domain_size
                                    # ,"&sort_by_structure=",as.numeric(input$l_sort_by_structure)
                                    ,"&significance_threshold_method=",input$l_correction_method#"&threshold_algo=",input$l_correction_method
                       )
                       ,target="_blank")

        setProgress(value = 4/6, detail = "WebGestalt link generation...")
        wgst_url <- a("Click me to be directed to WebGestalt Website!"
                      ,href=paste0("http://www.webgestalt.org/option.php",'?gene_list=',query,'&organism=',lgpr_organism
                                   ,"&enrich_method=ORA&fdr_method=BH&enriched_database_category=geneontology"
                                   ,"&enriched_database_name=Biological_Process_noRedundant&sig_method=top&sig_value=5"
                                   ,"&max_num=200&id_type=genesymbol"
                                   ,"&ref_set=genome")
                      ,target="_blank")
        setProgress(value = 5/6, detail = "gProfiler Table rendering...")
        output$l_gprofile_table <- renderDataTable({limma_gprof_table},options = list(scrollX = TRUE))

        output$l_gprofileLink <- renderUI({
          tagList("gprofile Link:", gprof_url)
        })

        output$l_webgestaltLink <- renderUI({
          tagList("WebGestalt Link:", wgst_url)
        })

        output$l_gprofile_res_dld<-renderUI({
          downloadButton(outputId = "limma_gprofiler_download",label = "Download gProfile results!")
        })
        setProgress(value = 1, detail = "Done!")
      })
    })
    #########################################################

    output$limma_gprofiler_download <- downloadHandler(
      filename = "gProfileR_results.csv",
      content = function(file){
        write.csv(limma_gprof_table,
                  file)
      }
    )


    ##########################################################
    ### When "Run edgeR DGE Analysis!" button is clicked
    #########################################################
    observeEvent(input$edgeR_runDGE_btn,{
      withProgress(message = "Performing DGE Analysis with edgeR... Please wait!",{
        n <- 12
        setProgress(value = 1/n, detail = "DGE List constructor...")
        sample_info.edger <- relevel(conditions , ref = "A")
        design  <- model.matrix(~sample_info.edger)
        design  <- model.matrix(~sample_info.edger)
        edgeR.DGElist <- DGEList(counts = txi$counts , group = sample_info.edger)
        # keep <- rowSums( cpm(edgeR.DGElist) >= 1) >= 5
        # edgeR.DGElist <- edgeR.DGElist[keep ,]
        setProgress(value = 2/n, detail = "Calculating normalization factors...")
        edgeR.DGElist <- calcNormFactors(edgeR.DGElist , method = "TMM")
        rownames(design) <- colnames(edgeR.DGElist)
        # estimate  the  dispersion  for all  read  counts  across  all  samples
        setProgress(value = 3/n, detail = "Estimating dispersion...")
        edgeR.DGElist  <- estimateDisp(edgeR.DGElist , design)
        setProgress(value = 4/n, detail = "Estimating dispersion...")
        edgeR.DGElist  <- estimateCommonDisp(edgeR.DGElist)
        setProgress(value = 5/n, detail = "Estimating dispersion...")
        edgeR.DGElist  <- estimateTrendedDisp(edgeR.DGElist)
        setProgress(value = 6/n, detail = "Estimating dispersion...")
        edgeR.DGElist  <- estimateTagwiseDisp(edgeR.DGElist)
        if(input$edgeR_model=="glmFit"){
          setProgress(value = 7/n, detail = "Fitting a negative binomial generalized log-linear model...")
          edger_fit  <- glmFit(edgeR.DGElist , design)
          setProgress(value = 8/n, detail = "Conducting likelihood ratio tests...")
          edger_lrt  <- glmLRT(edger_fit)
          # extract  results  from  edger_lrt$table  plus  adjusted p-values
          setProgress(value = 9/n, detail = "Extracting DGE results...")
          DGE.results_edgeR  <<- as.data.frame(topTags(edger_lrt , n = Inf , # to  retrieve  all  genes
                                                       sort.by = "PValue", adjust.method = "BH"))
          edgeR_fit_mod[["glmFit"]] <<- DGE.results_edgeR
          results_dge[["edgeR"]] <<- DGE.results_edgeR
        }else if(input$edgeR_model=="glmQLFit"){
          setProgress(value = 7/n, detail = "Fitting a quasi-likelihood negative binomial log-linear model...")
          edger_fit <- glmQLFit(edgeR.DGElist, design)
          setProgress(value = 8/n, detail = "Conducting empirical Bayes quasi-likelihood F-tests...")
          edger_qlf <- glmQLFTest(edger_fit,coef=2)
          setProgress(value = 9/n, detail = "Extracting DGE results...")
          DGE.results_edgeR  <<- as.data.frame(topTags(edger_qlf , n = Inf , # to  retrieve  all  genes
                                                       sort.by = "PValue", adjust.method = "BH"))
          edgeR_fit_mod[["glmQLFit"]] <<- DGE.results_edgeR
          results_dge[["edgeR_QL"]] <<- DGE.results_edgeR
        }

        setProgress(value = 10/n, detail = "Rendering DGE results table...")
        output$edgeR_DGE_results_sum <- renderPrint({
          summary(edger_fit)
        })

        output$edgeR_type2 <- renderUI({
          selectInput("edgeR_type_sel", "Select an edgeR results data", names(edgeR_fit_mod)
                      ,selected = names(DGE.results_edgeR)[1]
                      ,multiple = FALSE
                      ,selectize = TRUE, width = NULL, size = NULL)
        })

        output$edgeR_getData <- renderUI({
          actionButton("edgeR_getData_btn","Get edgeR data results!")
        })

        DT_output <- DGE.results_edgeR
        feature <- paste0('<a target="_blank" href="',searchURL,rownames(DT_output),'">',rownames(DT_output),'</a>')
        row.names(DT_output)<-feature
        # DT_output$Feature <- feature

        output$edgeR_dge_res <- renderDataTable({DT_output},options = list(scrollX = TRUE), escape = F)
        #output$edgeR_dge_res2 <- renderDataTable({DGE.results_edgeR},options = list(scrollX = TRUE))

        output$edgeR_pvalues <- renderPlot(
          hist(DGE.results_edgeR$PValue, col = "grey", border = "white", xlab = "", ylab = "", main = "frequencies  of p-values")
        )

        setProgress(value = 11/n, detail = "UI rendering...")
        output$edgeR_sumout <-renderUI({
          # plotMDS()
          withSpinner(verbatimTextOutput("edgeR_DGE_results_sum"),size=1,proxy.height=200,type=sample(c(1,4,5:8),1))
        })

        output$edgeR_volcano <- renderPlot({
          apval<-input$edgeR_vp_pval
          lfc<-input$edgeR_vp_lfc
          lim<-input$edgeR_vp_limit
          with(DGE.results_edgeR, plot(logFC, -log10(FDR), pch=20, main="Volcano plot", xlim=lim#c(-10,10)
          ))
          with(subset(DGE.results_edgeR, FDR<apval ), points(logFC, -log10(FDR), pch=20, col="blue"))
          with(subset(DGE.results_edgeR, FDR<apval & abs(logFC)>lfc), points(logFC, -log10(FDR), pch=20, col="red"))
        })

        output$edgeR_plotMA <-renderPlot({
          plotBCV(edgeR.DGElist, main="Biological Coefficient of Variation")
        })

        output$edgeR_dge_res_dld <- renderUI({
          downloadButton(outputId = "edgeR_dge_download",label = "Download limma DGE results!")
        })

        e_choices <- colnames(DGE.results_edgeR)

        output$edgeR_settings <- renderUI({
          tagList(
            selectInput("edgeR_sort_by","Sort by:", e_choices, selected="FDR",width = "200px")
            ,conditionalPanel(condition='input.edgeR_sort_by=="logFC"',
                              fluidRow(column(6
                                              ,sliderInput("e_logFC","Define a range:",min = 0#floor(min(DGE.results_edgeR$logFC, na.rm = T))
                                                           ,max = abs(ceiling(max(DGE.results_edgeR$logFC, na.rm = T)))
                                                           ,step = 0.1#round(max(DGE.results_edgeR$logFC, na.rm = T)/100,2)
                                                           ,value = 2#c(floor(min(DGE.results_edgeR$logFC, na.rm = T))
                                                                      #,floor(min(DGE.results_edgeR$logFC, na.rm = T))
                                                                      #+max(DGE.results_edgeR$logFC, na.rm = T)/100)
                                                           ,round = -2))
                              )
            )
            ,conditionalPanel(condition='input.edgeR_sort_by=="logCPM"',
                              fluidRow(column(6
                                              ,sliderInput("e_logCPM","Define a range:",min = floor(min(DGE.results_edgeR$logCPM, na.rm = T))
                                                           ,max = ceiling(max(DGE.results_edgeR$logCPM, na.rm = T))
                                                           ,step = round(max(DGE.results_edgeR$logCPM, na.rm = T)/100,2)
                                                           ,value = c(floor(min(DGE.results_edgeR$logCPM, na.rm = T))
                                                                      ,floor(min(DGE.results_edgeR$logCPM, na.rm = T))
                                                                      +max(DGE.results_edgeR$logCPM, na.rm = T)/100)
                                                           ,round = -2))
                              )
            )
            ,conditionalPanel(condition='input.edgeR_sort_by=="F"',
                              fluidRow(column(6
                                              ,sliderInput("e_F","Define a range:",min = floor(min(DGE.results_edgeR$F, na.rm = T))
                                                           ,max = ceiling(max(DGE.results_edgeR$F, na.rm = T))
                                                           ,step = round(max(DGE.results_edgeR$F, na.rm = T)/100,2)
                                                           ,value = c(floor(min(DGE.results_edgeR$F, na.rm = T))
                                                                      ,floor(min(DGE.results_edgeR$F, na.rm = T))
                                                                      +max(DGE.results_edgeR$F, na.rm = T)/100)
                                                           ,round = -2))
                              )
            )
            ,conditionalPanel(condition='input.edgeR_sort_by=="LR"',
                              fluidRow(column(6
                                              ,sliderInput("e_LR","Define a range:",min = floor(min(DGE.results_edgeR$LR, na.rm = T))
                                                           ,max = ceiling(max(DGE.results_edgeR$LR, na.rm = T))
                                                           ,step = round(max(DGE.results_edgeR$LR, na.rm = T)/100,2)
                                                           ,value = c(floor(min(DGE.results_edgeR$LR, na.rm = T))
                                                                      ,floor(min(DGE.results_edgeR$LR, na.rm = T))
                                                                      +max(DGE.results_edgeR$LR, na.rm = T)/100)
                                                           ,round = -2))
                              )
            )
            ,conditionalPanel(condition='input.edgeR_sort_by=="PValue"',
                              fluidRow(column(6
                                              ,sliderInput("e_PValue","Define a range:",min = 0#floor(min(DGE.results_edgeR$P.Value))
                                                           ,max = 1#ceiling(max(DGE.results_edgeR$PValue))
                                                           ,step = 0.01#round(max(DGE.results_edgeR$PValue)/100,2)
                                                           ,value = c(0,0.05)#c(floor(min(DGE.results_edgeR$PValue)),floor(min(DGE.results_edgeR$PValue))
                                                           #+max(DGE.results_edgeR$PValue)/100)
                                                           ,round = -2))
                              )
            )
            ,conditionalPanel(condition='input.edgeR_sort_by=="FDR"',
                              fluidRow(column(6
                                              ,sliderInput("e_FDR","Define a range:",min = 0#floor(min(DGE.results_edgeR$FDR))
                                                           ,max = 1#ceiling(max(DGE.results_edgeR$FDR))
                                                           ,step = 0.01#round(max(DGE.results_edgeR$FDR)/100,2)
                                                           ,value = c(0,0.05)#c(floor(min(DGE.results_edgeR$FDR)),floor(min(DGE.results_edgeR$FDR))
                                                           #+max(DGE.results_edgeR$FDR)/100)
                                                           ,round = -2))
                              )
            )
            ,actionButton("e_filter","Filter features!"),br(), hr()
            ,verbatimTextOutput("edgeR_num_dge"), hr()
            ,uiOutput("e_heatmap_settings"), br()
            ,hr()
          )
        })
        setProgress(value = 1, detail = "Done!")
      })
    })
    ##################################################################################################################

    observeEvent(input$edgeR_getData_btn,{
      DT_output<-edgeR_fit_mod[[input$edgeR_type_sel]]
      feature <- paste0('<a target="_blank" href="',searchURL,rownames(DT_output),'">',rownames(DT_output),'</a>')
      row.names(DT_output)<-feature
      output$edgeR_dge_res2 <- renderDataTable({DT_output},options = list(scrollX = TRUE), escape = F)

      output$edgeR_filterGenes <- renderUI({
        edgeR_data <- as.data.frame(edgeR_fit_mod[[input$edgeR_type_sel]])
        tagList(fluidRow(column(5
                                ,selectInput("edgeR_gprof_select", "Filter features by:",colnames(DGE.results_edgeR),
                                             selected = "logFC", multiple = FALSE,
                                             selectize = TRUE, width = "200px", size = NULL))
                         ,column(5
                                 ,conditionalPanel(condition='input.edgeR_gprof_select=="logFC"',
                                                   sliderInput("edgeR_gp_logFC","Define a range:"
                                                               ,min = 0#floor(min(edgeR_data$logFC, na.rm = T))
                                                               ,max = abs(ceiling(max(edgeR_data$logFC, na.rm = T)))
                                                               ,step = 0.1#round(max(edgeR_data$logFC, na.rm = T)/100,2)
                                                               ,value = 2#c(floor(min(edgeR_data$logFC, na.rm = T))
                                                                          #,floor(min(edgeR_data$logFC, na.rm = T))
                                                                          #+max(edgeR_data$logFC, na.rm = T)/100)
                                                               ,round = -2))
                                 ,conditionalPanel(condition='input.edgeR_gprof_select=="logCPM"',
                                                   sliderInput("edgeR_gp_logCPM","Define a range:"
                                                               ,min = floor(min(edgeR_data$logCPM, na.rm = T))
                                                               ,max = ceiling(max(edgeR_data$logCPM, na.rm = T))
                                                               ,step = round(max(edgeR_data$logCPM, na.rm = T)/100,2)
                                                               ,value = c(floor(min(edgeR_data$logCPM, na.rm = T)),floor(min(edgeR_data$logCPM, na.rm = T))
                                                                          +max(edgeR_data$logCPM, na.rm = T)/100)
                                                               ,round = -2))
                                 ,conditionalPanel(condition='input.edgeR_gprof_select=="F"',
                                                   sliderInput("edgeR_gp_F","Define a range:",min = floor(min(edgeR_data$F, na.rm = T))
                                                               ,max = ceiling(max(edgeR_data$F, na.rm = T))
                                                               ,step = round(max(edgeR_data$F, na.rm = T)/100,2)
                                                               ,value = c(floor(min(edgeR_data$F, na.rm = T)),floor(min(edgeR_data$F, na.rm = T))
                                                                          +max(edgeR_data$F)/100)
                                                               ,round = -2))
                                 ,conditionalPanel(condition='input.edgeR_gprof_select=="LR"',
                                                   sliderInput("edgeR_gp_LR","Define a range:",min = floor(min(edgeR_data$LR, na.rm = T))
                                                               ,max = ceiling(max(edgeR_data$LR, na.rm = T))
                                                               ,step = round(max(edgeR_data$LR, na.rm = T)/100,2)
                                                               ,value = c(floor(min(edgeR_data$LR, na.rm = T)),floor(min(edgeR_data$LR, na.rm = T))
                                                                          +max(edgeR_data$LR)/100)
                                                               ,round = -2))
                                 ,conditionalPanel(condition='input.edgeR_gprof_select=="PValue"',
                                                   sliderInput("edgeR_gp_PValue","Define a range:"
                                                               ,min = floor(min(edgeR_data$PValue, na.rm = T))
                                                               ,max = ceiling(max(edgeR_data$PValue, na.rm = T))
                                                               ,step = round(max(edgeR_data$PValue, na.rm = T)/100,2)
                                                               ,value = c(0,0.05)
                                                               ,round = -2))
                                 ,conditionalPanel(condition='input.edgeR_gprof_select=="FDR"',
                                                   sliderInput("edgeR_gp_FDR","Define a range:"
                                                               ,min = floor(min(edgeR_data$FDR, na.rm = T))
                                                               ,max = ceiling(max(edgeR_data$FDR, na.rm = T))
                                                               ,step = round(max(edgeR_data$FDR, na.rm = T)/100,2)
                                                               ,value = c(0,0.05)
                                                               ,round = -2))
                         )
        )
        ,actionButton("edgeR_filt_genes_btn","Filter features!")
        ,verbatimTextOutput("edgeR_filt_genes")
        )
      })

      output$edgeR_gprofile_par <- renderUI({
        tagList(
          fluidRow(hr(),
                   column(5,selectInput("egp_organism", "Choose organism"
                                        ,list("Homo sapiens"="hsapiens", "Mus musculus"="mmusculus", "Gorilla gorilla"="ggorilla"
                                              ,"Bonobo"="ppaniscus", "Ovis aries"="oaries","Other"="other"),
                                        selected = "hsapiens", multiple = FALSE,
                                        selectize = TRUE, width = "200px", size = NULL)
                          ,conditionalPanel(
                            condition = "input.egp_organism == 'other'",
                            textInput("e_other_gpo", "Specify a Profiler organism", value = "aaegypti",
                                      placeholder = "e.g: Enter 'aaegypti' for Aedes aegypti (https://biit.cs.ut.ee/gprofiler/help.cgi?help_id=64)"))
                   )
                   ,column(5,selectInput("e_correction_method", "Correction method"
                                         ,list("Analytical"="analytical", "gSCS"="gSCS", "fdr"="fdr"
                                               ,"bonferroni"="bonferroni"),
                                         selected = "fdr", multiple = FALSE,
                                         selectize = TRUE, width = "200px", size = NULL))
                   ,column(5,selectInput("e_hier_filtering", "Hierarchical filtering strength"
                                         ,list("None"="none", "Moderate"="moderate", "Strong"="strong"),
                                         selected = "none", multiple = FALSE,
                                         selectize = TRUE, width = "200px", size = NULL))
                   ,column(5,selectInput("e_domain_size", "Statistical domain size"
                                         ,list("Annotated"="annotated", "Known"="known"),
                                         selected = "annotated", multiple = FALSE,
                                         selectize = TRUE, width = "200px", size = NULL))
                   ,column(5, materialSwitch("e_sort_by_structure","Sort by structure",value=T
                                             ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
                   ,column(5, materialSwitch("e_ordered_query","Ordered Query",value=F
                                             ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
                   ,column(5, materialSwitch("e_significant","Significant Only",value=T
                                             ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
                   ,column(5, materialSwitch("e_exclude_iea","Exclude electronic annotations (IEA)",value=F
                                             ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
          )
          ,uiOutput("e_perf_gprof_btn"),hr()
          ,uiOutput("e_gprofileLink"), hr()
          ,uiOutput("e_webgestaltLink"), hr()
          ,dataTableOutput("e_gprofile_table"), hr()
          ,uiOutput("e_gprofile_res_dld")
        )
      })
    })

    output$edgeR_dge_download <- downloadHandler(
      filename = "edgeR_DGE_results.csv",
      content = function(file){
        write.csv(DGE.results_edgeR, file)
      }
    )

    observeEvent(input$e_filter,{
      if(input$edgeR_sort_by=="logFC"){
        lfc <- eval(parse(text = paste0("input$e_",input$edgeR_sort_by)))
        DGEgenes_edgeR  <<- rownames(subset(abs(DGE.results_edgeR), eval(parse(text = input$edgeR_sort_by))>lfc))
      }else{
        minn <- eval(parse(text = paste0("input$e_",input$edgeR_sort_by)))[1]
        maxx <- eval(parse(text = paste0("input$e_",input$edgeR_sort_by)))[2]
        DGEgenes_edgeR  <<- rownames(subset(DGE.results_edgeR, eval(parse(text = input$edgeR_sort_by))>=minn &
                                              eval(parse(text = input$edgeR_sort_by))<=maxx))
      }
      if(length(DGEgenes_edgeR)>0){
        output$e_heatmap_settings <- renderUI({
          tagList(fluidRow(column(2,materialSwitch("e_rdl", label="Row_dend_left:", value=F
                                                   ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
                           ,column(2,materialSwitch("e_Rowv", label="Rowv", value=T
                                                    ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
                           ,column(2,materialSwitch("e_Colv", label="Colv", value=T
                                                    ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
                           , column(2,materialSwitch("e_colorbar", label="Colorbar:", value=F
                                                     ,status = sample(c('primary', 'info', 'success', 'warning', 'danger'), 1)))
          )
          ,fluidRow(column(2, pickerInput(inputId = "e_colors", label = "Colors:"
                                          ,choices=c("BrBG","Spectral","RdYlBu","RdYlGn","RdGy","RdBu","PuOr","PRGn","PiYG")
                                          , selected = "PiYG"))
                    ,column(2, pickerInput(inputId = "e_dist_method", label = "Distance method:"
                                           ,choices = c("euclidean", "maximum","manhattan","canberra","binary","minkowski")
                                           ,selected="euclidean"))
                    ,column(2, pickerInput(inputId = "e_hclust_method", label = "Hclust method:"
                                           ,choices = c("complete", "ward.D","ward.D2","single","average","mcquitty","median","centroid")
                                           ,selected="average"))
                    ,column(2, pickerInput(inputId = "e_dendogram", label = "Dendogram:", choices = c("both", "row", "column", "none")
                                           ,selected="both"))
          )
          ,fluidRow(column(2, pickerInput(inputId = "e_plot_method", label = "Plot method:", choices = c("plotly", "ggplot")
                                           ,selected="ggplot"))
                    ,column(2,numericInput("e_k_col", "k_col:", 1, min = 1, max = 100, step = 1))
                    ,column(2,numericInput("e_k_row", "k_row:", 1, min = 1, max = 100, step = 1))
          )
          ,hr()
          ,actionButton("e_heatmap","Generate Heatmap!")
          )
        })
      }
      output$edgeR_num_dge <- renderPrint({
        list("DGE"=paste(length(DGEgenes_edgeR),"differentially expressed features"),"Differentially Expressed features"=DGEgenes_edgeR)
      })
    })

    observeEvent(input$e_heatmap,{

      # extract  the  normalized  read  counts  for DE genes  into a matrix
      hm.mat_DGEgenes.edgeR  <- log.norm.counts[DGEgenes_edgeR, ]
      output$edgeR_heatmap <- renderPlotly({
        withProgress(message = "Generating heatmap...Please wait!", value = 0.8,{
          isolate(heatmaply::heatmaply(as.matrix(hm.mat_DGEgenes.edgeR), xlab = "Samples", ylab = "Features",
                    scale = "row", row_dend_left = input$e_rdl, plot_method = input$e_plot_method, dendrogram=input$e_dendogram,
                    Rowv = input$e_Rowv, Colv = input$e_Colv, #distfun = input$e_distfun, #hclustfun = "average",
                    colorbar=input$e_colorbar, dist_method = input$e_dist_method, hclust_method = input$e_hclust_method,
                    colors = eval(parse(text=input$e_colors)),#PiYG,
                    k_col = input$e_k_col, k_row = input$e_k_row,
                    margins = c(60,100,40,20)))
        })
      })
    })

    ##################################################################################################################
    observeEvent(input$edgeR_filt_genes_btn,{
      output$edgeR_filt_genes <- renderPrint({
        if(input$edgeR_gprof_select=="logFC"){
          edgeR_gprofile_genes <<- rownames(subset(as.data.frame(edgeR_fit_mod[[input$edgeR_type_sel]])
                                                   ,abs(eval(parse(text = input$edgeR_gprof_select))) >
                                                     eval(parse(text = paste0("input$edgeR_gp_",input$edgeR_gprof_select)))))
        }else{
          edgeR_gprofile_genes <<- rownames(subset(as.data.frame(edgeR_fit_mod[[input$edgeR_type_sel]])
                                                   ,eval(parse(text = input$edgeR_gprof_select)) <= eval(parse(text = paste0("input$edgeR_gp_",input$edgeR_gprof_select)))[2] &
                                                     eval(parse(text = input$edgeR_gprof_select)) >= eval(parse(text = paste0("input$edgeR_gp_",input$edgeR_gprof_select)))[1]))
        }
        if(any(startsWith(as.character(edgeR_gprofile_genes),"ENS"))==T)
          edgeR_gprofile_genes <<- unique(apply(as.data.frame(edgeR_gprofile_genes),MARGIN = 1,
                   FUN = function(x){if(startsWith(as.character(x),"ENS")){return(strsplit(as.character(x),fixed = T,split = '.')[[1]][1])}
                   }))
        list("Number of filtered features:"=paste(length(edgeR_gprofile_genes),"filtered features!"),"Features:"=edgeR_gprofile_genes)
      })
      output$e_perf_gprof_btn <- renderUI(actionButton("edgeR_gprofile_btn","Perform gProfile pathway Analysis!"))
    })

    #########################################################
    observeEvent(input$edgeR_gprofile_btn,{
      withProgress(message = "Gene Pathway analysis in progress...",{
        egpr_organism <- input$egp_organism
        if(egpr_organism=='other'){egpr_organism<-input$e_other_gpo}
        setProgress(value = 2/6, detail = "gProfiler pathway Analysis...")
        edgeR_gprof_table <<- as.data.frame(gprofiler(edgeR_gprofile_genes
                                                      ,organism=egpr_organism, sort_by_structure = input$e_sort_by_structure, ordered_query = input$e_ordered_query
                                                      ,significant = input$e_significant, exclude_iea = input$e_exclude_iea, correction_method = input$e_correction_method
                                                      ,hier_filtering = input$e_hier_filtering, domain_size = input$e_domain_size
        ))
        query=paste(edgeR_gprofile_genes,collapse = "%0A")
        setProgress(value = 3/6, detail = "gProfiler link generation...")
        gprof_url <- a("Click me to be directed to gprofile Website!"
                       ,href=paste0("https://biit.cs.ut.ee/gprofiler/gost"#"https://biit.cs.ut.ee/gprofiler/index.cgi"
                                    ,'?query=',query,'&organism=',egpr_organism
                                    ,"&significant=",as.numeric(input$e_significant)
                                    ,"&ordered=",tolower(as.character(input$e_ordered_query))#"&ordered_query=",as.numeric(input$e_ordered_query)
                                    ,"&no_iea=",as.numeric(input$e_exclude_iea)
                                    ,"&domain_scope=",input$e_domain_size#"&domain_size_type=",input$e_domain_size
                                    # ,"&sort_by_structure=",as.numeric(input$e_sort_by_structure)
                                    ,"&significance_threshold_method=",input$e_correction_method#"&threshold_algo=",input$e_correction_method
                       )
                       ,target="_blank")

        setProgress(value = 4/6, detail = "WebGestalt link generation...")
        wgst_url <- a("Click me to be directed to WebGestalt Website!"
                      ,href=paste0("http://www.webgestalt.org/option.php",'?gene_list=',query,'&organism=',egpr_organism
                                   ,"&enrich_method=ORA&fdr_method=BH&enriched_database_category=geneontology"
                                   ,"&enriched_database_name=Biological_Process_noRedundant&sig_method=top&sig_value=5"
                                   ,"&max_num=200&id_type=genesymbol"
                                   ,"&ref_set=genome")
                      ,target="_blank")
        setProgress(value = 5/6, detail = "gProfiler Table rendering...")
        output$e_gprofile_table <- renderDataTable({edgeR_gprof_table},options = list(scrollX = TRUE))

        output$e_gprofileLink <- renderUI({
          tagList("gprofile Link:", gprof_url)
        })

        output$e_webgestaltLink <- renderUI({
          tagList("WebGestalt Link:", wgst_url)
        })

        output$e_gprofile_res_dld<-renderUI({
          downloadButton(outputId = "edgeR_gprofiler_download",label = "Download gProfile results!")
        })
        setProgress(value = 1, detail = "Done!")
      })
    })
    #########################################################

    output$edgeR_gprofiler_download <- downloadHandler(
      filename = "gProfileR_results.csv",
      content = function(file){
        write.csv(edgeR_gprof_table, file)
      }
    )

    ########################################################
    # VENN DIAGRAM
    ########################################################
    observeEvent(input$venn_data_btn,{
      if(length(names(results_dge))<=0){
        output$venn_message <- renderPrint({"No Venn Data found"
        })
      }else if(length(names(results_dge))==1){
        output$venn_message <- renderPrint({paste("Only one DGE model found:",names(results_dge))})
      }else{
        output$venn_message <- renderPrint({paste(length(names(results_dge))
                                                  ," Venn data table are found:"
                                                  ,paste(names(results_dge),collapse = ", "))})

        i=1
        if("DESeq2" %in% names(results_dge)){
          output[[paste0("sel",i)]] <- renderUI({
            selectInput(paste0("vsel_","DESeq2")
                        ,paste0("Select a column for DESeq2")
                        ,colnames(results_dge[["DESeq2"]])
                        ,selected = colnames(results_dge[["DESeq2"]])[1]
                        ,multiple = F
                        ,selectize = T, width = "200px", size = NULL)
          })
          i = i + 1
        }

        if("edgeR" %in% names(results_dge)){
          output[[paste0("sel",i)]] <- renderUI({
            selectInput(paste0("vsel_","edgeR")
                        ,paste0("Select a column for edgeR")
                        ,colnames(results_dge[["edgeR"]])
                        ,selected = colnames(results_dge[["edgeR"]])[1]
                        ,multiple = F
                        ,selectize = T, width = "200px", size = NULL)
          })
          i = i + 1
        }

        if("edgeR_QL" %in% names(results_dge)){
          output[[paste0("sel",i)]] <- renderUI({
            selectInput(paste0("vsel_","edgeR_QL")
                        ,paste0("Select a column for edgeR_QL")
                        ,colnames(results_dge[["edgeR_QL"]])
                        ,selected = colnames(results_dge[["edgeR_QL"]])[1]
                        ,multiple = F
                        ,selectize = T, width = "200px", size = NULL)
          })
          i = i + 1
        }

        if("limma" %in% names(results_dge)){
          output[[paste0("sel",i)]] <- renderUI({
            selectInput(paste0("vsel_","limma")
                        ,paste0("Select a column for limma")
                        ,colnames(results_dge[["limma"]])
                        ,selected = colnames(results_dge[["limma"]])[1]
                        ,multiple = F
                        ,selectize = T, width = "200px", size = NULL)
          })
          i = i + 1
        }
        output$spar <- renderUI({actionButton("setpar_btn", "Set Parameters!")})
      }
    })

    observeEvent(input$setpar_btn,{
      i = 1

      if("DESeq2" %in% names(results_dge)){
        output[[paste0("m_slider",i)]] <- renderUI({
          sliderInput("v_DESeq2_slider", paste("Define threshold for",input$vsel_DESeq2,"(DESeq2):")
                      ,min = 0#abs(round(min(eval(parse(text = paste0("results_dge$DESeq2$",input$vsel_DESeq2))), na.rm = T)))
                      ,max = abs(round(max(eval(parse(text = paste0("results_dge$DESeq2$",input$vsel_DESeq2))), na.rm = T)))
                      ,value = 0.05 ,step = 0.01)
        })
        i = i + 1
      }

      if("edgeR" %in% names(results_dge)){
        output[[paste0("m_slider",i)]] <- renderUI({
          sliderInput("v_edgeR_slider", paste("Define threshold for",input$vsel_edgeR,"(edgeR):")
                      ,min = 0#abs(round(min(eval(parse(text = paste0("results_dge$edgeR$",input$vsel_edgeR))), na.rm = T)))
                      ,max = abs(round(max(eval(parse(text = paste0("results_dge$edgeR$",input$vsel_edgeR))), na.rm = T)))
                      ,value = 0.05 ,step = 0.01)
        })
        i = i + 1
      }

      if("edgeR_QL" %in% names(results_dge)){
        output[[paste0("m_slider",i)]] <- renderUI({
          sliderInput("v_edgeR_QL_slider", paste("Define threshold for",input$vsel_edgeR_QL,"(edgeR_QL):")
                      ,min = 0#abs(round(min(eval(parse(text = paste0("results_dge$edgeR_QL$",input$vsel_edgeR_QL))), na.rm = T)))
                      ,max = abs(round(max(eval(parse(text = paste0("results_dge$edgeR_QL$",input$vsel_edgeR_QL))), na.rm = T)))
                      ,value = 0.05 ,step = 0.01)
        })
        i = i + 1
      }

      if("limma" %in% names(results_dge)){
        output[[paste0("m_slider",i)]] <- renderUI({
          sliderInput("v_limma_slider", paste("Define threshold for",input$vsel_limma,"(limma):")
                      ,min = 0#abs(round(min(eval(parse(text = paste0("results_dge$limma$",input$vsel_limma))), na.rm = T)))
                      ,max = abs(round(max(eval(parse(text = paste0("results_dge$limma$",input$vsel_limma))), na.rm = T)))
                      ,value = 0.05 ,step = 0.01)
        })
        i = i + 1
      }

      output$gen_venn <- renderUI({actionButton("venndiagram_btn", "Generate Venn Diagrams!")})
    })

    observeEvent(input$venndiagram_btn,{
      withProgress(message = "Generating Venn Diagram...",{
        DE_list <- list()

        if("DESeq2" %in% names(results_dge)){
          if(input$vsel_DESeq2=="log2FoldChange"){
            dsq <- rownames(subset(results_dge$DESeq2 , abs(eval(parse(text = input$vsel_DESeq2))) > input$v_DESeq2_slider))
          }else{
            dsq <- rownames(subset(results_dge$DESeq2 , abs(eval(parse(text = input$vsel_DESeq2))) <= input$v_DESeq2_slider))
          }

          if(any(startsWith(as.character(dsq),"ENS"))==T)
            dsq <- unique(apply(as.data.frame(dsq),MARGIN = 1,
                                FUN = function(x){if(startsWith(as.character(x),"ENS")){return(strsplit(as.character(x),fixed = T,split = '.')[[1]][1])}
                                }))
          DE_list[["DESeq2"]] <- dsq
        }
        if("edgeR" %in% names(results_dge)){
          if(input$vsel_edgeR=="logFC"){
            edgr <- rownames(subset(results_dge$edgeR , abs(eval(parse(text = input$vsel_edgeR))) > input$v_edgeR_slider))
          }else{
            edgr <- rownames(subset(results_dge$edgeR , abs(eval(parse(text = input$vsel_edgeR))) <= input$v_edgeR_slider))
          }

          if(any(startsWith(as.character(edgr),"ENS"))==T)
            edgr <- unique(apply(as.data.frame(edgr),MARGIN = 1,
                                 FUN = function(x){if(startsWith(as.character(x),"ENS")){return(strsplit(as.character(x),fixed = T,split = '.')[[1]][1])}
                                 }))
          DE_list[["edgeR"]] <- edgr
        }
        if("edgeR_QL" %in% names(results_dge)){
          if(input$vsel_edgeR_QL=="logFC"){
            edgrq <- rownames(subset(results_dge$edgeR_QL , abs(eval(parse(text = input$vsel_edgeR_QL))) > input$v_edgeR_QL_slider))
          }else{
            edgrq <- rownames(subset(results_dge$edgeR_QL , abs(eval(parse(text = input$vsel_edgeR_QL))) <= input$v_edgeR_QL_slider))
          }

          if(any(startsWith(as.character(edgrq),"ENS"))==T)
            edgrq <- unique(apply(as.data.frame(edgrq),MARGIN = 1,
                                  FUN = function(x){if(startsWith(as.character(x),"ENS")){return(strsplit(as.character(x),fixed = T,split = '.')[[1]][1])}
                                  }))
          DE_list[["edgeR_QL"]] <- edgrq
        }
        if("limma" %in% names(results_dge)){
          if(input$vsel_limma=="logFC"){
            lma <- rownames(subset(results_dge$limma , abs(eval(parse(text = input$vsel_limma))) > input$v_limma_slider))
          }else{
            lma <- rownames(subset(results_dge$limma , abs(eval(parse(text = input$vsel_limma))) <= input$v_limma_slider))
          }

          if(any(startsWith(as.character(lma),"ENS"))==T)
            lma <- unique(apply(as.data.frame(lma),MARGIN = 1,
                                FUN = function(x){if(startsWith(as.character(x),"ENS")){return(strsplit(as.character(x),fixed = T,split = '.')[[1]][1])}
                                }))
          DE_list[["limma"]] <- lma
        }

        DE_gns  <- UpSetR :: fromList(DE_list)
        output$venn_output <- renderPrint({
          out<-gplots::venn(DE_list, show.plot = F)
          out <<- gsub('"','',capture.output(out))
          out
        })
        output$dld_venn_output <- renderUI({
          downloadButton(outputId = "venn_output_btn",label = "Download Venn Output!")
        })
        output$venn_diagram <- renderPlot({
          gplots::venn(DE_list)
        })
        output$upset_plot <- renderPlot({
          UpSetR ::upset(DE_gns , order.by = "freq", text.scale = 2)
        })
      })
    })

    output$venn_output_btn <- downloadHandler(
      filename = "VennDiagram_output.txt",
      content = function(file){write(out,file)}
    )

  })

  }
