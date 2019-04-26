## Venn Diagram
venn <- tabItem(tabName = "venn", br(), br(),
                h2("Venn Diagrams"), br()
                ,actionButton("venn_data_btn","Get Data!"), hr()
                ,fluidRow(column(12,verbatimTextOutput("venn_message")))
                # ,uiOutput("venn_settings")
                ,fluidRow(column(3,uiOutput("sel1"))
                          ,column(3,uiOutput("sel2"))
                          ,column(3,uiOutput("sel3"))
                          ,column(3,uiOutput("sel4"))
                )
                ,uiOutput("spar") ,hr()
                ,fluidRow(column(3,uiOutput("m_slider1"))
                          ,column(3,uiOutput("m_slider2"))
                          ,column(3,uiOutput("m_slider3"))
                          ,column(3,uiOutput("m_slider4"))
                )
                ,uiOutput("gen_venn"), hr()
                ,fluidRow(column(10,flipBox(id=1,main_img = "venn-diagram.svg"
                                            ,header_img = "https://image.flaticon.com/icons/svg/119/119598.svg"
                                            ,front_title = "Venn Diagram"
                                            ,back_title = "UpSetPlot"
                                            ,front_btn_text = "View UpSetPlot"
                                            ,back_btn_text = "Back to main Venn Diagram", width = 10
                                            ,plotOutput("venn_diagram", height = "600px")
                                            ,back_content = tagList(
                                              plotOutput("upset_plot", height = "600px")
                                            )
                ))
                ), hr(), br()
                ,verbatimTextOutput("venn_output")
                # ,boxPlus(collapsible=T, closable=F, width = 10, title = "Venn Diagram"
                #          ,status = "primary", solidHeader = TRUE
                #           ,plotOutput("venn_diagram1", height = "800px")
                #           ,footer = tagList(fluidRow(column(2)
                #                                      ,column(3, selectInput("v_", "Select a data type"
                #                                                             ,c("RNASeq", "RNASeq2"), selected = "RNASeq", multiple = F
                #                                                             ,selectize = T, width = "200px", size = NULL))
                #                                      ,column(3, sliderInput("limma_vp_lfc", "Define threshold for log2FC:"
                #                                                             ,min = 0 ,max = 10 ,value = 2 ,step = 0.1))
                #                                      ,column(3, sliderInput("limma_vp_limit", "Define plot x-limit:"
                #                                                             ,min = -15 ,max = 15 ,value = c(-3, 3) ,step = 1))
                #           ))
                # )
)