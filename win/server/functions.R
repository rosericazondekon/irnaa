#-------------------------------------------------------------------------
#  Roseric Azondekon,
#  December 17, 2018
#  Last Update: March 11, 2019
#  ZSPH, Milwaukee, WI, USA
#-------------------------------------------------------------------------

#############################################
### HELPER FUNCTIONS
############################################
# Loading file data...
filedata <- function(infile){
  # infile <- input$datafile
  if (is.null(infile)) {
    # User has not uploaded a file yet
    return(NULL)
  }
  if(endsWith(infile$name, ".csv")){
    read.csv(infile$datapath)
  } else if(endsWith(infile$name, ".txt")){
    read.table(infile$datapath,header = T)
  }else if(endsWith(infile$name, ".xls")){
    read_excel(infile$datapath)
  } else if(endsWith(infile$name, ".xlsx")){
    read_xlsx(infile$datapath,sheet=1)
  }
}


# plotPCA2 function
plotPCA2 <- function(
  my_pca,
  colorby = NULL,
  xaxs = "PC1", yaxs = "PC2",
  title = NULL,
  palette = function(x) rainbow(x, s = 0.6),
  continuous_color = FALSE,
  ...
) {
  
  # obtaining % of variance explained by each axis
  varxaxs <- (summary(my_pca)$importance["Proportion of Variance", xaxs] * 100) %>% round(digits = 1)
  varyaxs <- (summary(my_pca)$importance["Proportion of Variance", yaxs] * 100) %>% round(digits = 1)
  
  myplot <- data.frame(
    myxaxs = my_pca$x[, xaxs],
    myyaxs = my_pca$x[, yaxs],
    texts = rownames(my_pca$x),
    colors = if (is.null(colorby)) {NA} else {colorby}
  ) %>%
    ggplot(aes(x = myxaxs, y = myyaxs, color = colors, text = texts)) +
    geom_point() +
    labs(x = paste0(xaxs, " (", varxaxs, "%)"), y = paste0(yaxs, " (", varyaxs, "%)"), title = title, color = NULL) +
    background_grid()
  
  if (continuous_color) {
    myplot <- myplot + scale_color_gradientn(colors = palette(64))
    # } else {
    #   myplot <- myplot + scale_color_manual(values = palette(n_distinct(colorby)))
  }
  
  ggplotly(myplot, tooltip = "text", ...)
  
}

