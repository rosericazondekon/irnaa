@echo off
title IRNAA shinyApp installation
rem IF exist %path%;%~dp0win (echo %path%;%~dp0win already exists) ELSE (setx path "%path%;%~dp0win")
setx path "%path%;%~dp0win"
setx IRNAAPATH "%~dp0win"
mklink "%userprofile%\Desktop\irnaa" "%~dp0win\irnaa.bat"
echo IRNAA shinyApp successfully installed
pause
