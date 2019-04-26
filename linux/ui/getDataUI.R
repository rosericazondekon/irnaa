## DESeq2
DESeq2 <- tabItem(tabName = "DESeq2", br(), br(),
                  h2("Set conditions and get datasets"), br(),
                  boxPlus(collapsible=T, closable=F, width = 10, title = "Define conditions:"
                          ,status = "primary", solidHeader = TRUE
                          ,uiOutput("sel_conds"), uiOutput("set_cond"), verbatimTextOutput("print_conds"), verbatimTextOutput("print_conds2")
                          ,uiOutput("get_DESeq2_data")
                  )
                  ,boxPlus(collapsible=T, closable=F, width = 10, title = "Dataset Overview"
                           ,status = "primary", solidHeader = TRUE
                           ,fluidRow(column(1), column(9, plotlyOutput("bar_counts", width = "90%", height = "600px")))
                           ,fluidRow(column(width=4, offset = 0, verbatimTextOutput("conditions_info"))
                           )
                           ,fluidRow(column(width=6, offset = 0, verbatimTextOutput("library_size"))
                                     ,column(width=6, offset = 0, verbatimTextOutput("library_counts"))
                           )
                  )
)