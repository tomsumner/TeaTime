## Script to run the TB model as compiled C code

## Set working directory
setwd("C:/Users/TOM SUMMER/Documents/TIME_research/TeaTime")

## Load packages
library(deSolve)
library(reshape2)
library(ggplot2)

## Compile and load the C code
system("R CMD SHLIB TB_model_v2.c") # Compile
dyn.load("TB_model.dll")

# Set up the model parameters

fit_cost=0.7
g = fit_cost/(1+fit_cost) # superinfections 

parms <- c(age1 = 1/5, age2 = 1/21, beta = 10, a = 0.115, p = 0.65, v = 0.0001, sig = 0.62, rel_inf = 0.22, theta = 0.015,   
           r = 0.2, mu_N = 0.2, mu_I = 0.3, fit_cost = 0.7, e = 0.014, g=g, k = 0, l_s = 0, l_m = 0, d = 0, tau_s = 0, tau_m = 0,
           eff_n = 0, eff_p = 0, dst_n = 0, dst_p = 0) 

# Define times to solve for
times <- seq(1950,2050 , by=1)

# Set up the initial conditions

UN_pop_age <- as.data.frame(read.table("SA_pop_age.txt",header=TRUE)) # Load UN Population data
temp <- c()
for (i in 1:num_ages){temp[i]<-UN_pop_age[1,i+1]}
xstart <- c(S=c(0.9999*temp),
           Lsn=rep(0,num_ages),Lsp=rep(0,num_ages),Lmn=rep(0,num_ages),Lmp=rep(0,num_ages),
           Nsn=rep(0,num_ages),Nsp=rep(0,num_ages),Nmn=rep(0,num_ages),Nmp=rep(0,num_ages),
           Isn=c(0.0001*temp),Isp=rep(0,num_ages),Imn=rep(0,num_ages),Imp=rep(0,num_ages))

# Set up the forcing functions

## Fertility
# Uses crude birth rate (per 1000 population) from UN pop for South Africa which has values for 5 year periods
# values post 2010 based on medium fertiliy
birth_rate <- cbind(c(0,1950,seq(1952.5,2047.5,5)),
                    c(30.52,30.52,44,42,41,38,38,36,34,31,27,25,24,22,21,20,18,17,17,16,15,14))

## Survival
Survive_age <- as.data.frame(read.table("SA_survival_age.txt",header=TRUE)) # Load survival proportions calculated from life tables
# Proportion surviving from age 0 to 1 - used to determine entry to first age group
s_birth <- cbind(Survive_age$Year,Survive_age$X1)
# and to other ages
s5 <- cbind(Survive_age$Year,Survive_age$X5)
s10 <- cbind(Survive_age$Year,Survive_age$X10)
s15 <- cbind(Survive_age$Year,Survive_age$X15)
s20 <- cbind(Survive_age$Year,Survive_age$X20)
s25 <- cbind(Survive_age$Year,Survive_age$X25)
s30 <- cbind(Survive_age$Year,Survive_age$X30)
s35 <- cbind(Survive_age$Year,Survive_age$X35)
s40 <- cbind(Survive_age$Year,Survive_age$X40)
s45 <- cbind(Survive_age$Year,Survive_age$X45)
s50 <- cbind(Survive_age$Year,Survive_age$X50)
s55 <- cbind(Survive_age$Year,Survive_age$X55)
s60 <- cbind(Survive_age$Year,Survive_age$X60)
s65 <- cbind(Survive_age$Year,Survive_age$X65)
s70 <- cbind(Survive_age$Year,Survive_age$X70)
s75 <- cbind(Survive_age$Year,Survive_age$X75)
s80 <- cbind(Survive_age$Year,Survive_age$X80)

# Combine forcing functions into a list
force <- list(birth_rate,s_birth,s5,s10,s15,s20,s25,s30,s35,s40,s45,s50,s55,s60,s65,s70,s75,s80)

# Run the model
time_C<-system.time(for(i in 1:10) out <- ode(y=xstart, times, func = "derivsc",
                               parms = parms, dllname = "TB_model",initforc = "forcc",
                               forcings=force, initfunc = "parmsc", nout = 13,
                               outnames = c("Total","Total_S","Total_L","Total_Ns","Total_Nm",
                                            "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM"), method = rkMethod("rk34f")))

# Unload the C code
dyn.unload("TB_model.dll")


