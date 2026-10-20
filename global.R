options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))

if (!requireNamespace("difR", quietly = TRUE) ||
    packageVersion("difR") < "6.1.0")
  install.packages("difR")

# VGAM : requis par difPolyLogistic (vglm) — installation binaire rapide
if (!requireNamespace("VGAM", quietly = TRUE))
  install.packages("VGAM")

if (!requireNamespace("ltm", quietly = TRUE))
  install.packages("ltm")

library(shiny)
library(difR)
library(markdown)
library(shinythemes)