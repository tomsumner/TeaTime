## Script to run the TB model as compiled C code

## Set working directory
setwd("C:/Users/TOM SUMMER/Documents/TIME_research/TeaTime")

## Load packages #################################################################################################
library(deSolve)
library(reshape2)
library(ggplot2)

## Compile and load C code #######################################################################################
dyn.unload("TB_model_v6.dll") # Unload - only really need to do this before recompiling
system("R CMD SHLIB TB_model_v6.c") # Compile
dyn.load("TB_model_v6.dll") # Load

# load logcurve function #########################################################################################
source("logcurve.R")

## Define the country ############################################################################################
# Pick from "SA", "Vietnam", .... (Maybe set this up for 5 countries that have been fitted - China, India, Bangladesh?)
cn <- "SA"

# Create forcing functions - these are time varying inputs used in the C code ####################################

# They need to consist of x,y pairs where x is the year and y is the value in that year
# They can be entered as numeric values or using the "logcurve" function to simulate a logistic curve:

# logcurve(base, target, base_year, target_year, growth, shape)
# "base" and "target" are the initial and final values
# "shape" determines the symmetry of the curve
# with "shape=1" the curve is symmetric around the midpoint of "base_year" and "target_year"
# "growth" determines the steepness of the curve

# BCG coverage - currently assume 100% at all times
BCG_cov <- cbind(c(1972,1973,2050),c(1,1,1))

# Case detection rate by HIV (neg/pos)
kneg <- logcurve(50,90,1990,2020,0.3,1)
kpos <- logcurve(50,90,1990,2020,0.3,1) 
  
# Relative detection smear neg
rel_d <- cbind(seq(1970,2050),0.8)

# DST coverage among new and previously treated cases 
dstneg_n <- logcurve(0,95,1975,2010,1,2)
dstneg_p <- logcurve(0,95,1975,2010,1,2)
dstpos_n <- logcurve(0,95,1975,2010,1,2)
dstpos_p <- logcurve(0,95,1975,2010,1,2)

# Linkage to care
l_s <- cbind(seq(1970,2050),0.83)
l_m <- cbind(seq(1970,2050),0.7)

# Treatment success by HIV (neg, pos no ART, pos on ART) and susceptibility
tneg_s <- cbind(seq(1970,2050),0.76)
tpos_s <- cbind(seq(1970,2050),0.76)
tART_s <- cbind(seq(1970,2050),0.76)
tneg_m <- cbind(seq(1970,2050),0.50)
tpos_m <- cbind(seq(1970,2050),0.50)
tART_m <- cbind(seq(1970,2050),0.50)
           
# Set up TB parameters ###########################################################################################

# Fitness of MDR (fit_cost), used to calculate parameter for superinfections (g)
# Both of these are passed into "parms" together with the MDR acquisition rate (e)
fit_cost=0.72
g = fit_cost/(1+fit_cost) # superinfections 
e = 0.01

# beta = contact rate; a = proportion developing primary TB; p = protection due to previous infection; v = reactivation rate; 
# sig = proportion of cases developing smear positive TB; rel_inf = relative infectiousness of smear negative cases;
# theta = rate of conversion from smear negative to smear positive; r = self cure rate; 
# mu_N = mortality for smear negative cases; mu_I = mortality for smear positive cases; 
# eff_n/p = relative efficacy of first line treatment in new/previously treated MDR cases; 

# _H indicates values for HIV+ (mu_N, mu_I, sig, r and rel_inf)
# other natural history parameters for HIV+ (a,v,p) are adjusted using the rate ratio parameters RR1a etc
# ART modifies a,v,p, muN_H and muI_H by ART_TB1 etc 

# proportion primary (a), proportion smear pos (sig) and mortality rates (muN and muI) take different values for 
# adults (>15) (_a), 0-4 (_0), 5-9 (_5) and 10-14 (_10)

parms <- c(beta = 21, 
           a_a = 0.115, a0 = 0.2171, a5 = 0.1155, a10 = 0.046,  
           p = 0.56, v = 0.0015, 
           sig_a = 0.45, sig0 = 0.0684, sig5 = 0.0414, sig10 = 0.0846, rel_inf = 0.25, theta = 0.02, r = 0.25, 
           mu_N = 0.18, mu_N0 = 0.3067, mu_I = 0.25, mu_I0 = 0.4260, fit_cost = fit_cost, e = e, g=g,
           eff_n = 0.61, eff_p = 0.45, 
           muN_H = 0.4, muI_H = 0.5, RR1a = 3.2, RR2a = 1.42, RR1v = 3.2, RR2v = 1.42, RR1p = 0.6, RR2p = 1.1,
           ART_TB1 = 0.16, ART_TB2 = 0.45, ART_TB3 = 0.55, ART_mort1 = 0.25, ART_mort2 = 0.6, ART_mort3 = 0.7,
           BCG_eff = 0.7,
           sig_H = 0.25,r_H=0.15,rel_inf_H=0.15)

# Load external data sources and create additional forcing functions where necessary #############################
source("Data_load.R")

# Run the model ##################################################################################################
# Script runs model with the inputs defined above. Runs an equilibrium phase then from 1970 to 2050
source("Run_model.R")

# Generate plots for comparison against TIME and other data ######################################################
source("Plots.R")





