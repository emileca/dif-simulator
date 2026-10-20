if (!requireNamespace("difR", quietly = TRUE) ||
    packageVersion("difR") < "6.1.0")
  install.packages("difR")

# VGAM : requis par difPolyLogistic (vglm)
if (!requireNamespace("VGAM", quietly = TRUE))
  install.packages("VGAM")

# ltm : requis par Lord et Raju
if (!requireNamespace("ltm", quietly = TRUE))
  install.packages("ltm")

library(shiny)
library(difR)
library(markdown)
library(shinythemes)