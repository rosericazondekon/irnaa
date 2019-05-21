################edgeR
edgeR_dge <- tabItem(tabName = "edgeR_dge", br(), br(),
                     h2("Differential Gene Analysis with edgeR"), br()
                     # ,fluidRow(column(3,uiOutput("edgeR_runDGE")),column(3,uiOutput("edgeR_type")))
                     ,uiOutput("edgeR_type")
                     ,uiOutput("edgeR_runDGE")
                     ,hr()
                     ,boxPlus(collapsible=T, closable=F, width = 10, title = "DGE analysis Summary with edgeR"
                              ,status = "primary", solidHeader = TRUE
                              ,fluidRow(column(width=4, plotOutput("edgeR_pvalues"))
                                        ,column(width=4, plotOutput("edgeR_plotMA"))
                                        ,column(width=4,uiOutput("edgeR_sumout"))
                              )
                     )
                     ,boxPlus(collapsible=T, closable=F, width = 10, title = "edgeR DGE results"
                              ,status = "primary", solidHeader = TRUE
                              ,dataTableOutput("edgeR_dge_res")
                              ,uiOutput("edgeR_dge_res_dld")
                     )
                     ,boxPlus(collapsible=T, closable=F, width = 10, title = "Volcano Plot"
                              ,status = "primary", solidHeader = TRUE
                              ,plotOutput("edgeR_volcano", height = "800px")#,plotlyOutput("edgeR_volcano", height = "800px")
                              ,footer = tagList(fluidRow(column(2)
                                                         ,column(3, sliderInput("edgeR_vp_pval", "Define threshold for FDR:"
                                                                                ,min = 0 ,max = 1 ,value = 0.05 ,step = 0.01))
                                                         ,column(3, sliderInput("edgeR_vp_lfc", "Define threshold for logFC:"
                                                                                ,min = 0 ,max = 10 ,value = 2 ,step = 0.1))
                                                         ,column(3, sliderInput("edgeR_vp_limit", "Define plot x-limit:"
                                                                                ,min = -15 ,max = 15 ,value = c(-3, 3) ,step = 1))
                              ))
                     )
                     ,boxPlus(collapsible=T, closable=F, width = 10, title = "Heatmap"
                              ,status = "primary", solidHeader = TRUE
                              ,uiOutput("edgeR_settings")
                              # ,verbatimTextOutput("edgeR_num_dge"), hr()
                              ,plotlyOutput("edgeR_heatmap", height = '800px')#,plotOutput("edgeR_heatmap", height = '800px')
                              # ,d3heatmapOutput("edgeR_heatmap")
                     )
)