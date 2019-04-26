if (!require("shiny")) install.packages("shiny")
shiny::runApp(system('echo $IRNAAPATH',intern=T),launch.browser=T)
