# Normalization
DESeq2_norm <- tabPanel("Normalization of dataset", icon = icon("line-chart")
                        ,boxPlus(collapsible=T, closable=F, width = 10, title = "Normalization for sequencing depth differences"
                                 ,status = "primary", solidHeader = TRUE
                                 ,uiOutput('normalization')
                                 ,fluidRow(column(width=5, offset = 0, plotOutput("plot_un_read_counts"))
                                           ,column(width=5, offset = 0, plotOutput("plot_tr_read_counts"))
                                 )
                        )
                        ,boxPlus(collapsible=T, closable=F, width = 10, title = "Checking for heteroskedasticity"
                                 ,status = "primary", solidHeader = TRUE
                                 ,fluidRow(column(width=5, offset = 0, plotOutput("plot_lnorm_counts"))
                                           ,column(width=5, offset = 0, plotOutput("plot_rlnorm_counts"))
                                 )
                                 ,fluidRow(column(width=5, offset = 0, plotOutput("vsn_lnorm_counts"))
                                           ,column(width=5, offset = 0, plotOutput("vsn_rlnorm_counts"))
                                 )
                        )
)

################ Normalization datasets
DESeq2_norm_data <- tabPanel("Normalized datasets", icon = icon("database")
                             ,boxPlus(collapsible=T, closable=F, width = 10, title = "Raw read counts"
                                      # ,loadingState()
                                      ,status = "primary", solidHeader = TRUE
                                      ,dataTableOutput("assay_info", width="100%")
                                      ,uiOutput("assay_info_dld")
                             )
                             ,boxPlus(collapsible=T, closable=F, width = 10, title = "Normalized read counts"
                                      # ,loadingState()
                                      ,status = "primary", solidHeader = TRUE
                                      ,dataTableOutput("norm_assay_info", width="100%")
                                      ,uiOutput("norm_assay_info_dld")
                             )
                             ,boxPlus(collapsible=T, closable=F, width = 10, title = "log2-Transformed Normalized read counts"
                                      # ,loadingState()
                                      ,status = "primary", solidHeader = TRUE
                                      ,dataTableOutput("lnorm_assay_info", width="100%")
                                      ,uiOutput("lnorm_assay_info_dld")
                             )
                             ,boxPlus(collapsible=T, closable=F, width = 10, title = "Regularized log2-Transformed Normalized read counts"
                                      # ,loadingState()
                                      ,status = "primary", solidHeader = TRUE
                                      ,dataTableOutput("rlnorm_assay_info", width="100%")
                                      ,uiOutput("rlnorm_assay_info_dld")
                             )
)

dataNorm <- tabItem(tabName = "dataNorm", br(), br(), br()
                    #,h2("Data Normalization")
                    ,fluidRow(tabBox(width=12,height = "2000px",DESeq2_norm,DESeq2_norm_data),selected=DESeq2_norm))