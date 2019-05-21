#-------------------------------------------------------------------------
#  Roseric Azondekon,
#  December 17, 2018
#  Last Update: May 21, 2019
#  ZSPH, Milwaukee, WI, USA
#-------------------------------------------------------------------------
tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css')
library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(shinycssloaders)
library('shinyDirectoryInput')

# Load UI sources...
source('../shared/sidebarUI.R')
source('./ui/readQuantUI.R')
source('../shared/getDataUI.R')
source('../shared/normUI.R')
source('../shared/globalrcUI.R')
source('../shared/DESeq2_dgeUI.R')
source('../shared/DESeq2_gprofilerUI.R')
source('../shared/edgeR_dgeUI.R')
source('../shared/edgeR_gprofilerUI.R')
source('../shared/limma_dgeUI.R')
source('../shared/limma_gprofilerUI.R')
source('../shared/vennUI.R')
source('../shared/readmeUI.R')
source('../shared/sessionInfoUI.R')
source('../shared/headerUI.R')

# use Sweet Alert...
useSweetAlert()



#Define menu...
menu <- tabItems(readQuant, DESeq2, dataNorm, DESeq2_grp, DESeq2_dge, DESeq2_gprofiler, 
                 edgeR_dge, edgeR_gprofiler, limma_dge, limma_gprofiler, venn, readme, #about, 
                 sessInfo
)

#Define body
body = dashboardBody(menu)
# body = dashboardBody(
#   setShadow(class = "dropdown-menu")
# )

#Define ui
ui = dashboardPagePlus(header, sidebar, body, rightsidebar, title = "DashboardPage")
