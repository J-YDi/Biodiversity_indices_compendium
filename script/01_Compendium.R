#_______________________________________________________________________________
# Title              : 01_Compendium.r
# Date               : 16/01/2025
# Object             : Script to determine a list of reliable R packages to 
#                      calculate diversity indices
# Authors            : Jean-Yves Dias
# R version          : 4.5.0
# Github link        : 
#_______________________________________________________________________________

#_______________________________ Packages_______________________________________####

# Function to install packages if not present and/or load them
loadpackages <- function(packages){
  for (pkg in packages){
    if (!requireNamespace(pkg,quietly = T)){
      install.packages(pkg)
      library(pkg,character.only = T)
    }else{
      library(pkg,character.only = T)
    }
  }
}

packages_needed <- c("readr","dplyr","tidyr","stringr")

loadpackages(packages_needed)

#_______________________________________________________________________________####

library(vegan)
data("BCI")
data("dune")
data("mite")
data("pyrifos")
data("sipoo")
data("varespec")
