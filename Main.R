## Script to run the TB model as compiled C code
## Runs the model once for a single parameter set

## Set working directory - input folders (Demog, HIV, TB) need to be subfolders of this directory
setwd("C:/Users/TOM SUMMER/Documents/TIME_research/TeaTime")

## Load packages, compile model and load DLL
source("Libraries_and_dll.R")

# Define number of disease states in model, number of age groups and number of HIV/ART states
# These values are used in other functions (e.g. plotting and data load)
n_dis <- 15
n_age <- 81
n_HIV <- 7
n_ART <- 3

## Define Country (currently 1=South_Africa, 2=Vietnam, 3=Bangladesh) ##################################################
cn <- 1
##
c_list <- c("South_Africa","Vietnam","Bangladesh")
cn <- c_list[cn]

# Load external data sources and create additional forcing functions where necessary #############################
source("Data_load.R")

# Set up the forcing functions and parameters
source(paste("Para_",cn,".R",sep=""))

# Run the model (and time it) 
strt<-Sys.time()

source("Run_model.R",local=TRUE) 

print(Sys.time()-strt)







