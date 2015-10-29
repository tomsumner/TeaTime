## Script to run the TB model as compiled C code

## Set working directory - input folders (Demog, HIV, TB) need to be subfolders of this directory
setwd("C:/Users/TOM SUMMER/Documents/TIME_research/TeaTime")

## This script loads the required libraries and compiles the dll of the model (which is coded in C) ##############
source("Libraries_and_dll.R")

ptm <- proc.time()

## Define the country ############################################################################################
# Pick from "South_Africa", "Vietnam","Bangladesh","India" .... (Maybe set this up for 5 countries that have been fitted - China, India)
cn <- "South_Africa"

# Step size (in years) to report outputs at - solver currently uses an adaptive step size solver
ss <- 1

# Load external data sources and create additional forcing functions where necessary #############################
source("Data_load.R")

## Indicate whether to generate and save plots or not 1=yes, 0=no ################################################
plot_gen <- 1
plot_save <- 0

# Generate forcing functions and parameter values for current country
source(paste("Para_",cn,".R",sep=""))

# Run the model ##################################################################################################
# Script runs model with the inputs defined above. Runs an equilibrium phase then from 1970 to 2050

model_time <- system.time(source("Run_model.R"))

#### CHECK FOR NEGATIVE VALUES #################
neg_cols<- apply(out,2,function(x) any(x<0))
which(neg_cols)

neg_rows <- apply(out,1,function(x) any(x<0))
out[which(neg_rows),"time"]  

# Generate plots for comparison against TIME and other data ######################################################
if (plot_gen == 1){
plot_time <- system.time(source("Plots.R"))
}

proc.time() - ptm


