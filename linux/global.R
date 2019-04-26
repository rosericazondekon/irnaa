#-------------------------------------------------------------------------
#  Roseric Azondekon,
#  December 17, 2018
#  Last Update: February 11, 2019
#  ZSPH, Milwaukee, WI, USA
#-------------------------------------------------------------------------

system("source ~/.bashrc",wait = F)
system("LCYAN='\033[1;36m' && NC='\033[0m'",wait = F)

repositories <<- c("https://cloud.r-project.org", "https://bioconductor.org/packages/3.7/bioc",
           "https://bioconductor.org/packages/3.7/data/annotation", "https://bioconductor.org/packages/3.7/data/experiment",
           "https://www.stats.ox.ac.uk/pub/RWin", "http://www.omegahat.net/R", "https://R-Forge.R-project.org",
           "https://www.rforge.net", "https://cloud.r-project.org", "http://www.bioconductor.org",
           "http://www.stats.ox.ac.uk/pub/RWin")
packages <- c("shinydashboardPlus", "shiny", "shinydashboard", "devtools", "vsn", "hexbin", "UpSetR", "gplots", "ROCR", "NMF", "DT",
              "org.Hs.eg.db", "org.Mm.eg.db", "org.Dr.eg.db", "pheatmap", "tximport", "readr", "TxDb.Hsapiens.UCSC.hg19.knownGene", 
              "rmarkdown", "edgeR", "shinyWidgets", "shinycssloaders", "CNTools", "openxlsx", "VennDiagram", "plyr", #"shinyjs",
              "plotly", "heatmaply", "cowplot", "dplyr",
              "DESeq2","ggpubr", "ggplot2", "limma", "biomaRt", "htmlwidgets", "AnnotationDbi", "Biobase", "ensembldb", "shinyjqui",
              "styler", "shinyAce","shinyFiles", "d3heatmap", "rhandsontable","magrittr","shinyjs","gProfileR", "TCGA2STAT") # CRAN packages

#"shinyTree"
# Install and load missing packages
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  install.packages(new.packages, repos = repositories)
}

lapply(packages, require, character.only = TRUE)

if(!'shinyDirectoryInput' %in% installed.packages()){
  devtools::install_github('wleepang/shiny-directory-input')
  library('shinyDirectoryInput')
}else{
  library('shinyDirectoryInput')
}

# if(!'shinyIncubator' %in% installed.packages()){
#   devtools::install_github("rstudio/shiny-incubator")
#   library('shinyIncubator')
# }else{
#   library('shinyIncubator')
# }


# =============================================

# Set up a button to have an animated loading indicator and a checkmark
# for better user experience
# Need to use with the corresponding `withBusyIndicator` server function
# withBusyIndicatorUI <- function(button) {
#   id <- button[['attribs']][['id']]
#   div(
#     `data-for-btn` = id,
#     button,
#     span(
#       class = "btn-loading-container",
#       hidden(
#         img(src = "ajax-loader-bar.gif", class = "btn-loading-indicator"),
#         icon("check", class = "btn-done-indicator")
#       )
#     ),
#     hidden(
#       div(class = "btn-err",
#           div(icon("exclamation-circle"),
#               tags$b("Error: "),
#               span(class = "btn-err-msg")
#           )
#       )
#     )
#   )
# }
# 
# # Call this function from the server with the button id that is clicked and the
# # expression to run when the button is clicked
# withBusyIndicatorServer <- function(buttonId, expr) {
#   # UX stuff: show the "busy" message, hide the other messages, disable the button
#   loadingEl <- sprintf("[data-for-btn=%s] .btn-loading-indicator", buttonId)
#   doneEl <- sprintf("[data-for-btn=%s] .btn-done-indicator", buttonId)
#   errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
#   shinyjs::disable(buttonId)
#   shinyjs::show(selector = loadingEl)
#   shinyjs::hide(selector = doneEl)
#   shinyjs::hide(selector = errEl)
#   on.exit({
#     shinyjs::enable(buttonId)
#     shinyjs::hide(selector = loadingEl)
#   })
#   
#   # Try to run the code when the button is clicked and show an error message if
#   # an error occurs or a success message if it completes
#   tryCatch({
#     value <- expr
#     shinyjs::show(selector = doneEl)
#     shinyjs::delay(2000, shinyjs::hide(selector = doneEl, anim = TRUE, animType = "fade",
#                                        time = 0.5))
#     value
#   }, error = function(err) { errorFunc(err, buttonId) })
# }
# 
# # When an error happens after a button click, show the error
# errorFunc <- function(err, buttonId) {
#   errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
#   errElMsg <- sprintf("[data-for-btn=%s] .btn-err-msg", buttonId)
#   errMessage <- gsub("^ddpcr: (.*)", "\\1", err$message)
#   shinyjs::html(html = errMessage, selector = errElMsg)
#   shinyjs::show(selector = errEl, anim = TRUE, animType = "fade")
# }
# 
# appCSS <- "
# .btn-loading-container {
#   margin-left: 10px;
#   font-size: 1.2em;
# }
# .btn-done-indicator {
#   color: green;
# }
# .btn-err {
#   margin-top: 10px;
#   color: red;
# }
# "
