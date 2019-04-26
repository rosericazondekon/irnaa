qcCheck <- tabItem(tabName = "qcCheck", br(), br(), 
                   h2("Quality Control of raw FASTQ data"), br(), 
                   boxPlus(collapsible=T, closable=F, width = 8, title = "Detected FASTQ samples"
                           ,status = "primary", solidHeader = TRUE, verbatimTextOutput("workDir"),
                           uiOutput("summary"), verbatimTextOutput("sampleMsg")
                           ,fluidRow(column(width=4, offset = 0, uiOutput("fastqc"))
                                     ,column(width=4, offset = 0, uiOutput("multiQc_report"))
                                     ,column(width=4, offset = 0, uiOutput("delQc_report")))
                   )
)