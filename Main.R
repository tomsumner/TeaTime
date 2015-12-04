## Script to run TIME research model - runs the model for a single parameter set for a defined country
## Model can be run using 1 or 5 year age bins (5 year is faster but less accurate demographically)

## First set the working directory - input folders (Demog, HIV, TB) need to be subfolders of this directory
setwd("C:/Users/TOM SUMMER/Documents/TIME_research/TeaTime")

## Define Country (1=South_Africa, 2=Vietnam, 3=Bangladesh)
cn <- 1

## Define number of age groups in the model (17=5 year age bins; 81=1 year age bins)
n_age <- 17

## PLOTTING

##################################################################################################################################

## Load packages, compile model and load DLL
## Load external data sources and create additional forcing functions where necessary
if (n_age==17) {
  source("Libraries_and_dll_5yr.R")
  source("Data_load_5yr.R")
}
if (n_age==81){
  source("Libraries_and_dll.R")
  source("Data_load.R")
}

# Set up the forcing functions and parameters
source(paste("Para_",cn,".R",sep=""))

# Run the model (and time it) 
strt<-Sys.time()

if (n_age==17) source("Run_model_5yr.R") 
if (n_age==81) source("Run_model.R")

print(Sys.time()-strt)

if (n_age==17) source("Plots_5yr.R") 
if (n_age==81) source("Plots.R")



