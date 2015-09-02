## Script to run the TB model as compiled C code

## Set working directory
setwd("C:/Users/TOM SUMMER/Documents/TIME_research/TeaTime")

## Load packages #################################################################################################
library(deSolve)
library(reshape2)
library(ggplot2)

## Function to set up logistic curves for care and control parameters ############################################
## "base" and "target" are the initial and final values
## "shape" determines the symmetry of the curve
## with "shape=1" the curve is symmetric around the midpoint of "base_year" and "target_year"
## "growth" determines the steepness of the curve
logcurve <- function(base,target,base_year,target_year,growth,shape){
  x <- seq(1970,2050)
  year <- mean(c(base_year,target_year))
  y <- (base + ((target-base)/((1+exp(-growth*(x-year)))^(1/shape))))/100
  z <- cbind(x,y)
  return(z)
}

## Compile and load C code #######################################################################################
dyn.unload("TB_model_v6.dll") # Unload - need to do this before recompiling
system("R CMD SHLIB TB_model_v6.c") # Compile
dyn.load("TB_model_v6.dll") # Load

## Define the country ############################################################################################
# Pick from "SA", "Vietnam", .... (Maybe set this up for 5 countries that have been fitted - China, India, Bangladesh?)
cn <- "SA"

# Load data and create forcing functions where necessary #########################################################
source("Data_load.R")

# Set up age structure ###########################################################################################
ages <- c(seq(1,80),100) # upper end of age classes
num_ages <- length(ages) # calculates the number of age classes

# Create some additional forcing functions - time varying inputs used in the C code ##############################
# These describe the care and control parameters in the model
# They need to consist of x,y pairs where x is the year and y is the value in that year
# They can be entered as numeric values or using the logcurve function to simulate a logistic curve

# BCG coverage - currently assume 90% at all times
BCG_cov <- cbind(c(1972,1973,2050),c(0.9,0.9,0.9))

# Case detection rate by HIV (neg/pos) - uses logistic curves (base,target,base_year,target_year,growth,shape)
kneg <- logcurve(50,90,1990,2020,0.3,1)
kpos <- logcurve(50,90,1990,2020,0.3,1) 
  
# Relative detection smear neg
rel_d <- cbind(seq(1970,2050),0.8)

# DST coverage among new and previously treated cases - uses logistic curves (base,target,base_year,target_year,growth,shape)
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
           
## Now put together all the forcing fucntions in a list to be passed to the C code ###############################
# s are mortality at age i       ]
# h are HIV incidence at age i   ] see "Data_load.R" for source of these
# A are ART inputs               ]

force <- list(birth_rate,
              s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,
              s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,s32,s33,s34,s35,s36,s37,s38,s39,s40,
              s41,s42,s43,s44,s45,s46,s47,s48,s49,s50,s51,s52,s53,s54,s55,s56,s57,s58,s59,s60,
              s61,s62,s63,s64,s65,s66,s67,s68,s69,s70,s71,s72,s73,s74,s75,s76,s77,s78,s79,s80,s81,
              h0,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h14,h15,h16,h17,h18,h19,h20,
              h21,h22,h23,h24,h25,h26,h27,h28,h29,h30,h31,h32,h33,h34,h35,h36,h37,h38,h39,h40,
              h41,h42,h43,h44,h45,h46,h47,h48,h49,h50,h51,h52,h53,h54,h55,h56,h57,h58,h59,h60,
              h61,h62,h63,h64,h65,h66,h67,h68,h69,h70,h71,h72,h73,h74,h75,h76,h77,h78,h79,h80,
              Ahigh,A500,A349,A249,A199,A99,A50,Athresh,
              BCG_cov,pop_ad,
              kneg,kpos,rel_d,dstneg_n,dstneg_p,dstpos_n,dstpos_p,l_s,l_m,tneg_s,tpos_s,tART_s,tneg_m,tpos_m,tART_m)

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


ptm <- proc.time()
# Run the model ##################################################################################################
# Script runs model with the inputs defined above
# Runs an equilibrium phase then from 1970 to 2050
for (irt in 1:100) source("Run_model.R")
proc.time() - ptm

# Generate plots for comparison against TIME and other data ######################################################
source("Plots.R")












