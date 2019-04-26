DESeq2_gprofiler <- tabItem(tabName = "DESeq2_gprofiler", br(), br(),
                            h2("Pathway and Gene Ontology"), br()
                            ,boxPlus(collapsible=T, closable=F, width = 10, title = "DGE results table"
                                     ,status = "primary", solidHeader = TRUE
                                     ,dataTableOutput("DESeq2_dge_res2")
                            )
                            ,boxPlus(collapsible=T, closable=F, width = 10, title = "Pathway Analysis and Gene Ontology"
                                     ,status = "primary", solidHeader = TRUE
                                     ,fluidRow(column(10, uiOutput("filterGenes"))
                                               ,column(10, uiOutput("DESeq2_gprofile_par"))
                                     )
                            )
)