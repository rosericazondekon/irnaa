################### DESeq2_grp : Exploring global Read Counts pattern
DESeq2_grp <- tabItem(tabName = "DESeq2_grp", br(), br(),
                      h2("Global Read Counts Pattern"), br()
                      ,uiOutput("global_rc")
                      ,hr()
                      ,boxPlus(collapsible=T, closable=F, width = 10, title = "Hierarchical Clustering"
                               # ,loadingState()
                               ,status = "primary", solidHeader = TRUE
                               ,uiOutput("clust_method")
                               ,plotOutput("DESeq2_hr_tree")
                      )
                      ,boxPlus(collapsible=T, closable=F, width = 10, title = "Principal Components Analysis"
                               # ,loadingState()
                               ,status = "primary", solidHeader = TRUE
                               ,plotlyOutput("DESeq2_pca", width = "75%", height = "600px")#,plotOutput("DESeq2_pca")
                      )
                      ,boxPlus(collapsible=T, closable=F, width = 10, title = "Pairwise Correlation"
                               # ,loadingState()
                               ,status = "primary", solidHeader = TRUE
                               ,dataTableOutput("DESeq2_pcorr", height="500px"),br()
                               ,uiOutput("DESeq2_pcorr_dld")
                      )
                      ,boxPlus(collapsible=T, closable=F, width = 10, title = "Correlation Heatmap"
                               # ,loadingState()
                               ,status = "primary", solidHeader = TRUE
                               ,plotlyOutput("pcorr_mat", height = '500px')
                      )
                      # ,fluidRow(column(10,flipBox(id=2,main_img = "corr.svg"
                      #                             ,header_img = "https://image.flaticon.com/icons/svg/119/119596.svg"
                      #                             ,front_title = "Pairwise Correlation Table"
                      #                             ,back_title = "Correlation matrix"
                      #                             ,front_btn_text = "View Correlation Matrix"
                      #                             ,back_btn_text = "Back to Correlation Table", width = 10
                      #                             ,tagList(dataTableOutput("DESeq2_pcorr", height="500px",width="99%"),br()
                      #                                      ,uiOutput("DESeq2_pcorr_dld")
                      #                             )
                      #                             ,back_content = tagList(
                      #                               plotlyOutput("pcorr_mat", height = '500px')
                      #                             )
                      # ))
                      # )
)
