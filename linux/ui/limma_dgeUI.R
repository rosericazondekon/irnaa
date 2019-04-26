################limma-voom
limma_dge <- tabItem(tabName = "limma_dge", br(), br(),
                     h2("Differential Gene Analysis with limma-voom"), br()
                     ,uiOutput("limma_runDGE")
                     ,hr()
                     ,boxPlus(collapsible=T, closable=F, width = 10, title = "DGE analysis Summary with limma-voom"
                              ,status = "primary", solidHeader = TRUE
                              ,fluidRow(column(width=4, plotOutput("limma_pvalues"))
                                        ,column(width=4, plotOutput("limma_plotMA"))
                                        ,column(width=4, uiOutput("limma_sumout"))
                              )
                     )
                     ,boxPlus(collapsible=T, closable=F, width = 10, title = "limma-voom DGE results"
                              ,status = "primary", solidHeader = TRUE
                              ,dataTableOutput("limma_dge_res")
                              ,uiOutput("limma_dge_res_dld")
                     )
                     ,boxPlus(collapsible=T, closable=F, width = 10, title = "Volcano Plot"
                              ,status = "primary", solidHeader = TRUE
                              ,plotOutput("limma_volcano", height = "800px")
                              ,footer = tagList(fluidRow(column(2)
                                                         ,column(3, sliderInput("limma_vp_pval", "Define threshold for adjusted p-value:"
                                                                                ,min = 0 ,max = 1 ,value = 0.05 ,step = 0.01))
                                                         ,column(3, sliderInput("limma_vp_lfc", "Define threshold for log2FC:"
                                                                                ,min = 0 ,max = 10 ,value = 2 ,step = 0.1))
                                                         ,column(3, sliderInput("limma_vp_limit", "Define plot x-limit:"
                                                                                ,min = -15 ,max = 15 ,value = c(-3, 3) ,step = 1))
                              ))
                     )
                     ,boxPlus(collapsible=T, closable=F, width = 10, title = "Heatmap"
                              ,status = "primary", solidHeader = TRUE
                              ,uiOutput("limma_settings")
                              # ,verbatimTextOutput("limma_num_dge"), hr()
                              ,plotlyOutput("limma_heatmap", height = '800px')#,plotOutput("limma_heatmap", height = '800px')
                              # ,d3heatmapOutput("DESeq2_heatmap")
                     )
)