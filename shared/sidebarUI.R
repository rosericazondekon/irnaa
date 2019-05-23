################################
#Define sidebar
################################
sidebar = dashboardSidebar(
  hr(),
  sidebarMenu(id="tabs" #,style = "position: fixed; overflow: visible; width: 225px; white-space: nowrap;"
              # ,menuItem("QC Check", tabName="qcCheck", icon=icon("stethoscope"), selected=TRUE)
              ,menuItem("Read Quantification", tabName = "readQuant", icon=icon("flask"))
              ,menuItem("Datasets", tabName = "DESeq2", icon = icon("hdd-o"))
              ,menuItem("Normalization", tabName = "dataNorm", icon = icon("line-chart"))
              # ,menuItem("Normalization", tabName = "DESeq2_norm", icon = icon("line-chart"))
              # ,menuItem("Normalized datasets", tabName = "DESeq2_norm_data", icon = icon("database"))
              ,menuItem("Global RC patterns", tabName = "DESeq2_grp", icon = icon("globe"))
              ,menuItem("DGE with DESeq2", tabName = "DESeq2_all", icon=icon("microchip")
                        ,menuSubItem("DGE Analysis", tabName = "DESeq2_dge", icon = icon("flask"))
                        ,menuSubItem("gProfile/WebGestaalt", tabName = "DESeq2_gprofiler", icon=icon("gg"))
              )
              ,menuItem("DGE with edgeR", tabName = "edgeR_all", icon=icon("ship")
                        ,menuSubItem("DGE Analysis", tabName = "edgeR_dge", icon = icon("flask"))
                        ,menuSubItem("gProfile/WebGestaalt", tabName = "edgeR_gprofiler", icon=icon("gg"))
              )
              ,menuItem("DGE with limma-voom", tabName = "limma_all", icon=icon("magic")
                        ,menuSubItem("DGE Analysis", tabName = "limma_dge", icon = icon("flask"))
                        ,menuSubItem("gProfile/WebGestaalt", tabName = "limma_gprofiler", icon=icon("gg"))
              )
              ,menuItem("Venn Diagrams", tabName = "venn", icon=icon("toggle-off"))
              ,menuItem("ReadMe", tabName = "readme", icon=icon("book"))#icon("mortar-board")
              # ,menuItem("About", tabName = "about", icon = icon("question"))
              ,menuItem("R Session Info", tabName = "sessInfo", icon = icon("info-circle"))
  )
)

################################
#Define rightsidebar
################################
rightsidebar <- rightSidebar(
  background = "dark", style = "position: fixed; overflow: visible;",
  rightSidebarTabContent(
    id = 1,
    title = "Working Directory",
    icon = "desktop",
    active = F#TRUE,
    # directoryInput('directory', label = 'select working directory', value = '~')
    # sliderInput(
    #   "obs",
    #   "Number of observations:",
    #   min = 0, max = 1000, value = 500
    # )
  )
  # ,rightSidebarTabContent(
  #   id = 2,
  #   title = "Reference transcriptome"#,
  #   # textInput("caption", "Caption", "Data Summary"),
  #   # shinyFilesButton("Btn_GetFile", "Choose a reference file" ,
  #   #                  title = "Please select a file:", multiple = FALSE,
  #   #                  buttonType = "default", class = NULL),
  #   # textOutput("txt_file2",inline = T)
  # )
  # ,rightSidebarTabContent(
  #   id = 3,
  #   icon = "paint-brush",
  #   title = "Tab 3",
  #   numericInput("obs", "Observations:", 10, min = 1, max = 100)
  # )
)
