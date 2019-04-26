#################### DESeq2_dge: DESeq2 Differential Gene Analysis
DESeq2_dge <- tabItem(tabName = "DESeq2_dge", br(), br(),
                      h2("Differential Gene Analysis with DESeq2"), br()
                      ,uiOutput("DESeq2_runDGE")
                      ,hr()
                      ,boxPlus(collapsible=T, closable=F, width = 10, title = "Exploratory plots following DGE analysis with DESeq2"
                               # ,loadingState()
                               ,status = "primary", solidHeader = TRUE
                               # ,verbatimTextOutput("DESeq2_DGE_results_sum")
                               
                               ,fluidRow(
                                 column(width=4, plotOutput("DESeq2_pvalues"))
                                 ,column(width=4, plotOutput("DESeq2_plotMA"))
                                 ,column(width=4, verbatimTextOutput("DESeq2_sumout"))
                               )
                      )
                      ,boxPlus(collapsible=T, closable=F, width = 10, title = "DGE results"
                               # ,loadingState()
                               ,status = "primary", solidHeader = TRUE
                               ,dataTableOutput("DESeq2_dge_res")
                               ,uiOutput("DESeq2_dge_res_dld")
                               # ,d3heatmapOutput("DESeq2_heatmap")
                      )
                      ,boxPlus(collapsible=T, closable=F, width = 10, title = "Volcano Plot"
                               # ,loadingState()
                               ,status = "primary", solidHeader = TRUE
                               ,plotOutput("DESeq2_volcano", height = "800px")
                               ,footer = tagList(fluidRow(column(2)
                                                          ,column(3, sliderInput("DESeq2_vp_pval", "Define threshold for adjusted p-value:"
                                                                                 ,min = 0 ,max = 1 ,value = 0.05 ,step = 0.01))
                                                          ,column(3, sliderInput("DESeq2_vp_lfc", "Define threshold for log2FC:"
                                                                                 ,min = 0 ,max = 10 ,value = 2 ,step = 0.1))
                                                          ,column(3, sliderInput("DESeq2_vp_limit", "Define plot x-limit:"
                                                                                 ,min = -15 ,max = 15 ,value = c(-3, 3) ,step = 1))
                               ))
                      )
                      ,boxPlus(collapsible=T, closable=F, width = 10, title = "Heatmap"
                               ,status = "primary", solidHeader = TRUE
                               ,uiOutput("DESeq2_settings")
                               # ,verbatimTextOutput("edgeR_num_dge"), hr()
                               ,plotlyOutput("DESeq2_heatmap", height = '800px')#,plotOutput("edgeR_heatmap", height = '800px')
                               # ,d3heatmapOutput("edgeR_heatmap")
                      )
)