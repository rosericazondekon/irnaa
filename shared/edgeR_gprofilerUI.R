edgeR_gprofiler <- tabItem(tabName = "edgeR_gprofiler", br(), br(),
                           h2("Pathway and Gene Ontology"), br()
                           ,fluidRow(column(3,uiOutput("edgeR_type2")),column(3,uiOutput("edgeR_getData")))
                           ,boxPlus(collapsible=T, closable=F, width = 10, title = "DGE results table"
                                    ,status = "primary", solidHeader = TRUE
                                    ,dataTableOutput("edgeR_dge_res2")
                           )
                           ,boxPlus(collapsible=T, closable=F, width = 10, title = "Pathway Analysis and Gene Ontology"
                                    ,status = "primary", solidHeader = TRUE
                                    ,fluidRow(column(10, uiOutput("edgeR_filterGenes"))
                                              ,column(10, uiOutput("edgeR_gprofile_par"))
                                    )
                           )
)