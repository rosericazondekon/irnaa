#-------------------------------------------------------------------------
#  Roseric Azondekon,
#  December 17, 2018
#  Last Update: May 21, 2019
#  ZSPH, Milwaukee, WI, USA
#-------------------------------------------------------------------------

# system("source ~/.bashrc",wait = F)
# system("LCYAN='\033[1;36m' && NC='\033[0m'",wait = F)

repositories <<- c("https://cloud.r-project.org", "https://bioconductor.org/packages/3.7/bioc",
           "https://bioconductor.org/packages/3.7/data/annotation", "https://bioconductor.org/packages/3.7/data/experiment",
           "https://www.stats.ox.ac.uk/pub/RWin", "http://www.omegahat.net/R", "https://R-Forge.R-project.org",
           "https://www.rforge.net", "https://cloud.r-project.org", "http://www.bioconductor.org",
           "http://www.stats.ox.ac.uk/pub/RWin")
packages <- c("shinydashboardPlus", "shiny", "shinydashboard", "devtools", "vsn", "hexbin", "UpSetR", "gplots", "ROCR", "NMF", "DT",
              "org.Hs.eg.db", "org.Mm.eg.db", "org.Dr.eg.db", "pheatmap", "tximport", "readr", "TxDb.Hsapiens.UCSC.hg19.knownGene", 
              "rmarkdown", "edgeR", "shinyWidgets", "shinycssloaders", "CNTools", "openxlsx", "VennDiagram", "plyr", #"shinyjs",
              "plotly", "heatmaply", "cowplot", "dplyr", "R6", "rlang", "httpud", "jsonlite", "later", "codetools",
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

