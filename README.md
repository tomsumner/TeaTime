# Time_R, Version 1.0
# Tom Sumner, 11/12/2015

# Research Version of TIME Impact model

# This version replicates the core dynamic transmission model of the TIME Impact module of the Spectrum suite developed by Avenir Health and the London School of Hygiene and Tropical Medicine
 
# The model is implemented in R, but uses compiled C to code the equations (for computational efficiency)
# The model can be run in two age structure configurations: single year age bins; 5 year age bins
 
# The model consists of the following files:
 
# Main.R - main script for running the model, the user must define the file path, the country and the age structure to use
# Libraries_and_dll.R - script to load libraries and compile the C code 
# TB_model.c - the equations implemented in C
# logcurve.R - function for defining generalised logisitic function
# Data_load.R - loads and processes required input data
# Para_cn.R (where cn is the country) - file to define the parameters for the current model run
# Run_model.R - calls R package desolve to solve the equations and return outputs
# Plots.R - plots 

# Libraries_and_dll.R, TB_model.c, Data_load.R, Run_model.R and Plots.R have seperate versions for the 5 year age bin model (postscript _5ry). The correct version is called based on the age structure defined in Main.R

# Folders (Demog, HIV and TB) containing input files must be in the directory defined in Main.R 