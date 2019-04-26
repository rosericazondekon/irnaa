#-------------------------------------------------------------------------
#  Roseric Azondekon,
#  December 17, 2018
#  Last Update: March 8, 2019
#  ZSPH, Milwaukee, WI, USA
#-------------------------------------------------------------------------
tags$link(rel = 'stylesheet', type = 'text/css', href = 'www/styles.css')

# Load packages
library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(shinycssloaders)
library('shinyDirectoryInput')

# useShinyjs()

# Load UI sources...
source('./ui/sidebarUI.R')
source('./ui/qcCheckUI.R')
source('./ui/readQuantUI.R')
source('./ui/getDataUI.R')
source('./ui/normUI.R')
source('./ui/globalrcUI.R')
source('./ui/DESeq2_dgeUI.R')
source('./ui/DESeq2_gprofilerUI.R')
source('./ui/edgeR_dgeUI.R')
source('./ui/edgeR_gprofilerUI.R')
source('./ui/limma_dgeUI.R')
source('./ui/limma_gprofilerUI.R')
source('./ui/vennUI.R')
source('./ui/readmeUI.R')
source('./ui/sessionInfoUI.R')
source('./ui/headerUI.R')

# use Sweet Alert...
useSweetAlert()



#Define menu...
menu <- tabItems(qcCheck, readQuant, DESeq2, dataNorm, DESeq2_grp, DESeq2_dge, DESeq2_gprofiler, 
                 edgeR_dge, edgeR_gprofiler, limma_dge, limma_gprofiler, venn, readme, #about, 
                 sessInfo
                 )

#Define body...
body = dashboardBody(menu)
# body = dashboardBody(
#   setShadow(class = "dropdown-menu")
# )


#Define ui...
ui = dashboardPagePlus(header, sidebar, body, rightsidebar, title = "DashboardPage")
