## Load packages #################################################################################################
library(deSolve)
library(reshape2)
library(ggplot2)

## Compile and load C code #######################################################################################
dyn.unload("TB_model_v7.dll") # Unload - only really need to do this before recompiling
system("R CMD SHLIB TB_model_v7.c") # Compile
dyn.load("TB_model_v7.dll") # Load

# load logcurve function #########################################################################################
source("logcurve.R")

