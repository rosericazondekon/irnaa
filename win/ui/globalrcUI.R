################### DESeq2_grp : Exploring global Read Counts pattern
DESeq2_grp <- tabItem(tabName = "DESeq2_grp", br(), br(),
                      h2("Global Read Counts Pattern"), br()
                      ,uiOutput("global_rc")
                      ,hr()
                      ,boxPlus(collapsible=T, closable=F, width = 10, title = "Hierarchical Clustering"
                               # ,loadingState()
                               ,status = "primary", solidHeader = TRUE
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
)