##readQuant: read Quantification
readQuant <- tabItem(tabName = "readQuant", br(), br(),
                     h2("Read Quantification"), br()
                     # ,radioButtons("quantType", h3(""),list("From FASTQ files"="fqc","From TCGA"="tcga","From files"="files"), selected = "fqc", inline=T)
                     ,radioGroupButtons(inputId = "quantType"
                                        ,status = sample(c("default","primary","danger","warning","success"),1)
                                        ,width = "600px"
                                        ,label = "", list("From FASTQ files"="fqc","From TCGA"="tcga","From files"="files")
                                        ,justified = TRUE, selected = "fqc", individual = TRUE
                                        ,checkIcon = list(yes = icon("ok",lib = "glyphicon")))
                     ,hr()
                     ,conditionalPanel(condition = "input.quantType == 'fqc'"
                                       ,fluidRow(
                                         column(width=12
                                                ,boxPlus(collapsible=T, closable=F, width = 8, title = "Transcriptome Reference File"
                                                         ,status = "primary", solidHeader = TRUE,
                                                         shinyFilesButton("Btn_GetFile", "Choose a reference file" ,
                                                                          title = "Please select a file:", multiple = FALSE,
                                                                          buttonType = "default", class = NULL),
                                                         verbatimTextOutput("txt_file"),
                                                         br(), br()
                                                         ,fluidRow(column(width=4, offset = 0, uiOutput("create_index_btn")))#align ="center",
                                                         ,withSpinner(verbatimTextOutput("index_output")
                                                                      ,size=1,proxy.height=200,type=sample(c(1,4,5:8), 1))
                                                )
                                         )
                                         ,column(width=12
                                                 ,boxPlus(collapsible=T, closable=F, width = 8, title = "Read Quantification"
                                                          ,status = "primary", solidHeader = TRUE
                                                          ,uiOutput('quantSamples')
                                                          ,fluidRow(column(width=4, offset = 0, uiOutput("quant_read"))
                                                                    ,column(width=3, offset = 2, uiOutput("nCores")))
                                                 )
                                                 ,boxPlus(collapsible=T, closable = F, width = 8, title = "Read Quantification Logs"
                                                          ,status = "primary", solidHeader = T, actionButton("show_logs","Show Details!")
                                                          ,fluidRow(column(width=4, uiOutput("quant_log_files"), uiOutput("view_log"))
                                                                    ,column(width=4, uiOutput("quant_sf_files"), uiOutput("view_quant"))
                                                          )
                                                 )
                                                 ,boxPlus(collapsible=T, closable = F, width = 8, title = "Import transcript-level estimates"
                                                          ,status = "primary", solidHeader = TRUE
                                                          ,uiOutput("annot"), uiOutput("cond_annot"), uiOutput("tr_output"),uiOutput("import_output")
                                                          ,uiOutput("tbChoice"), uiOutput("tb"), uiOutput("readCounts_dld")
                                                 )
                                         )
                                       )
                     )
                     ,conditionalPanel(condition = "input.quantType == 'tcga'"
                                       ,fluidRow(
                                         column(width=12
                                                ,boxPlus(collapsible=T, closable=F, width = 10, title = "Export TCGA data"
                                                         ,status = "primary", solidHeader = TRUE
                                                         ,fluidRow(column(6,selectInput("disease_select", "Select a table"
                                                                                        ,list("adrenocortical carcinoma (ACC)"="ACC", "bladder urothelial carcinom (BLCA)"="BLCA",
                                                                                              "breast invasive carcinoma (BRCA)"="BRCA", "cervical and endocervical cancers (CESC)"="CESC",
                                                                                              "cholangiocarcinoma (CHOL)"="CHOL", "colon adenocarcinoma (COAD)"="COAD",
                                                                                              "colorectal adenocarcinoma (COADREAD)"="COADREAD", "large B-cell lymphoma (DLBC)"="DLBC",
                                                                                              "esophageal carcinoma (ESCA)"="ESCA", "FFPE Pilot Phase II (FPPP)"="FPPP",
                                                                                              "glioblastoma multiforme (GBM)"="GBM", "glioma (GBMLGG)"="GBMLGG",
                                                                                              "head and neck squamous cell carcinoma (HNSC)"="HNSC", "kidney chromophobe (KICH)"="KICH",
                                                                                              "pan-kidney cohort (KIPAN = KICH+KIRC+KIRP)"="KIPAN", "kidney renal clear cell carcinoma (KIRC)"="KIRC",
                                                                                              "kidney renal papillary cell carcinoma (KIRP)"="KIRP", "acute myeloid leukemia (LAML)"="LAML",
                                                                                              "brain lower grade glioma (LGG)"="LGG", "liver hepatocellular carcinoma (LIHC)"="LIHC",
                                                                                              "lung adenocarcinoma (LUAD)"="LUAD", "lung squamous cell carcinoma (LUSC)"="LUSC",
                                                                                              "mesothelioma (MESO)"="MESO", "ovarian serous cystadenocarcinoma (OV)"="OV",
                                                                                              "pancreatic adenocarcinoma (PAAD)"="PAAD", "pheochromocytoma and paraganglioma (PCPG)"="PCPG",
                                                                                              "prostate adenocarcinoma (PRAD)"="PRAD", "rectum adenocarcinoma (READ)"="READ",
                                                                                              "sarcoma (SARC)"="SARC", "skin cutaneous melanoma (SKCM)"="SKCM", "stomach adenocarcinoma (STAD)"="STAD",
                                                                                              "testicular germ cell tumors (TGCT)"="TGCT", "thyroid carcinoma (THCA)"="THCA",
                                                                                              "thymoma (THYM)"="THYM", "uterine corpus endometrial carcinoma (UCEC)"="UCEC",
                                                                                              "uterine carcinosarcoma (UCS)"="UCS", "uveal melanoma (UVM)"="UVM"),
                                                                                        selected = "OV", multiple = FALSE,
                                                                                        selectize = TRUE, width = "400px", size = NULL))
                                                                   ,column(4,selectInput("dataType", "Select a data type"
                                                                                         ,c("RNASeq", "miRNASeq"#, "RNASeq2", "Mutation", "Methylation", "CNA_CGH", "mRNA_Array"
                                                                                         ), selected = "RNASeq", multiple = F
                                                                                         ,selectize = T, width = "200px", size = NULL)))
                                                         ,actionButton("getTCGA_btn","Get TCGA RNASeq data!")
                                                         ,hr()
                                                         ,verbatimTextOutput("tcga_report"), hr()#withSpinner(verbatimTextOutput("tcga_report"),size=1,proxy.height=200,type=sample(c(1,4,5:8), 1)), hr()
                                                         ,dataTableOutput("tcga_data", width="100%")
                                                         ,uiOutput("tcga_data_dld")
                                                )
                                         )
                                       )
                     )
                     ,conditionalPanel(condition = "input.quantType == 'files'"
                                       ,boxPlus(collapsible=T, closable=F, width = 10, title = "Load read counts dataset from file"
                                                ,status = "primary", solidHeader = T
                                                ,fileInput('datafile', 'Load a dataset (.xlsx, .xls, .txt or .csv)',
                                                           accept=c('text/csv','text/comma-separated-values,text/plain',".xls",".xlsx",".csv")
                                                           ,width = "800px")
                                                ,hr()
                                                ,actionButton("get_files_btn","Get data!"), br()
                                                ,verbatimTextOutput("files_data_report"), hr()
                                                # ,uiOutput("get_files_data"), hr()
                                                ,dataTableOutput("files_data")
                                       )
                     )
)
