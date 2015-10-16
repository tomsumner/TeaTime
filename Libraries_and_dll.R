## Load packages #################################################################################################
library(deSolve)
library(reshape2)
library(ggplot2)

## Compile and load C code #######################################################################################
if (is.loaded("derivsc")){
  dyn.unload("TB_model_v8.dll") # Unload - do this if currently loaded (only really necessary if recompiling)
}
system("R CMD SHLIB TB_model_v8.c") # Compile
dyn.load("TB_model_v8.dll") # Load

# load logcurve function #########################################################################################
source("logcurve.R")

