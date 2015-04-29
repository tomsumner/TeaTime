## Script to run the TB model as compiled C code

## Set working directory
setwd("C:/Users/TOM SUMMER/Documents/TIME_research/TeaTime")

## Load packages
library(deSolve)
library(reshape2)
library(ggplot2)

## Compile and load the C code
system("R CMD SHLIB TB_model_v4.c") # Compile
dyn.load("TB_model_v4.dll") # Load
dyn.unload("TB_model_v4.dll") # Unload - need to do this before recompiling

##############################################################################################################################

## Load UN population data
UN_pop_age <- as.data.frame(read.table("SA_pop_age.txt",header=TRUE)) # Load UN Population data
# add total to data
UN_pop_age_t <- cbind(UN_pop_age,rowSums(UN_pop_age[,2:18]))
colnames(UN_pop_age_t) <- c(colnames(UN_pop_age),"Total")

# Set up age structure
ages <- c(4,9,14,19,24,29,34,39,44,49,54,59,64,69,74,79,100) # upper end of age classes
num_ages <- length(ages) # calculates the number of age classes

# Set up the forcing functions for birth and death - all from 1970 onwards

## Fertility
# Uses crude birth rate (per 1000 population) from UN pop for South Africa which has values for 5 year periods
# values post 2010 based on medium fertiliy
birth_rate <- cbind(seq(1972.5,2047.5,5),
                    c(38,36,34,31,27,25,24,22,21,20,18,17,17,16,15,14))

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

# HIV Incidence by age and year - based on AIM output, but ignoring childhood infections (this is waht Carel does in TIME)
HIV_Inc_age <- as.data.frame(read.table("HIV_Inc_age.txt",header=TRUE)) # Load HIV incidence data taken from AIM                                       # Data from AIM is rate per 1000 
HIV_Inc_age[,2:18]=HIV_Inc_age[,2:18]*0

h0 <- 0*cbind(HIV_Inc_age$Year,HIV_Inc_age$X0/1000)
h5 <- 0*cbind(HIV_Inc_age$Year,HIV_Inc_age$X5/1000)
h10 <- 0*cbind(HIV_Inc_age$Year,HIV_Inc_age$X10/1000)
h15 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X15/1000)
h20 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X20/1000)
h25 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X25/1000)
h30 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X30/1000)
h35 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X35/1000)
h40 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X40/1000)
h45 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X45/1000)
h50 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X50/1000)
h55 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X55/1000)
h60 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X60/1000)
h65 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X65/1000)
h70 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X70/1000)
h75 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X75/1000)
h80 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X80/1000)

# Combine forcing functions into a list
force <- list(birth_rate,s_birth,s5,s10,s15,s20,s25,s30,s35,s40,s45,s50,s55,s60,s65,s70,s75,s80,
              h0,h5,h10,h15,h20,h25,h30,h35,h40,h45,h50,h55,h60,h65,h70,h75,h80)

# Set up TB parameters

# Fitness of MDR, used to calculate parameter for superinfections
# Both of these are passed into "parms"
fit_cost=0.7
g = fit_cost/(1+fit_cost) # superinfections 

# proportion primary (a), proportion smear pos (sig) and mortality rates (muN and mu_I) take different values for 
# adults (>15) (_a), 0-4 (_0), 5-9 (_5) and 10-14 (_10)

parms <- c(age1 = 1/5, age2 = 1/21, beta = 18, 
           a_a = 0.14, a0 = 0.26432, a5 = 0.14056, a10 = 0.056,  
           p = 0.65, v = 0.001, 
           sig_a = 0.45, sig0 = 0.0684, sig5 = 0.0414, sig10 = 0.0846, rel_inf = 0.25, theta = 0.02, r = 0.25, 
           mu_N = 0.25, mu_N0 = 0.426, mu_I = 0.35, mu_I0 = 0.59, fit_cost = fit_cost, e = 0, g=g, k = 0.3, l_s = 0.83, l_m = 0.0, d = 0.8, tau_s = 0.76, tau_m = 0.0,
           eff_n = 0.0, eff_p = 0.0, dst_n = 0.0, dst_p = 0.0) 

##############################################################################################################################
# Model initialisation
# run the model from 1970 pop with 1970 birth/death rates and care and control parameters for 400 years to get stable age structure
# Then rerun with 100 TB cases at time 0, no MDR or HIV (note no HIV pre 1975). Run for 400 years to get stable disease state

# Times to run model for
times <- seq(0,200 , by=1)

# Initial conditions - all susceptible
temp <- c()
for (i in 1:num_ages){temp[i]<-UN_pop_age[21,i+1]}
xstart <- c(S=c(temp),
            Lsn=rep(0,num_ages),Lsp=rep(0,num_ages),Lmn=rep(0,num_ages),Lmp=rep(0,num_ages),
            Nsn=rep(0,num_ages),Nsp=rep(0,num_ages),Nmn=rep(0,num_ages),Nmp=rep(0,num_ages),
            Isn=rep(0,num_ages),Isp=rep(0,num_ages),Imn=rep(0,num_ages),Imp=rep(0,num_ages),
            S_H=rep(0,num_ages*7),S_A=rep(0,num_ages*7*3))

# Run the model
time_1 <- system.time(out_pop <- ode(y=xstart, times, func = "derivsc",
            parms = parms, dllname = "TB_model_v4",initforc = "forcc",
            forcings=force, initfunc = "parmsc", nout = 31,
            outnames = c("Total","Total_S","Total_SH","Total_SA","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                         "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM",
                         "CD4500","CD4350_500","CD4250_349","CD4200_249","CD4100_199","CD450_99","CD450",
                         "ART500","ART350_500","ART250_349","ART200_249","ART100_199","ART50_99","ART50"), method = rkMethod("rk34f")))

# Update initial conditions based on end of last run and add 100 TB cases
temp <- c()
for (i in 1:num_ages){temp[i]<-out_pop[dim(out_pop)[1],i+1]}
xstart <- c(S=c(temp),
            Lsn=rep(0,num_ages),Lsp=rep(0,num_ages),Lmn=rep(0,num_ages),Lmp=rep(0,num_ages),
            Nsn=rep(0,num_ages),Nsp=rep(0,num_ages),Nmn=rep(0,num_ages),Nmp=rep(0,num_ages),
            Isn=c(rep(0,5),100,rep(0,11)),Isp=rep(0,num_ages),Imn=rep(0,num_ages),Imp=rep(0,num_ages),
            S_H=rep(0,num_ages*7),S_A=rep(0,num_ages*7*3))

# Run the model
time_2 <- system.time(out_TB <- ode(y=xstart, times, func = "derivsc",
                       parms = parms, dllname = "TB_model_v4",initforc = "forcc",
                       forcings=force, initfunc = "parmsc", nout = 31,
                       outnames = c("Total","Total_S","Total_SH","Total_SA","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                                    "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM",
                                    "CD4500","CD4350_500","CD4250_349","CD4200_249","CD4100_199","CD450_99","CD450",
                                    "ART500","ART350_500","ART250_349","ART200_249","ART100_199","ART50_99","ART50"), method = rkMethod("rk34f")))

# Adjust pop down to 1970 values and reassign initial conditions - model can now be run from 1970 with TB and HIV
temp <- out_TB[dim(out_TB)[1],2:698]
temp <- temp/(sum(temp)/22502) # 22502 is total pop from UN estimates in 1970)
xstart <- temp

##############################################################################################################################

# Now run the model for TB and HIV

# Set times to run for
times <- seq(1970,2070 , by=1)
# Run the model
time_3 <-system.time(out <- ode(y=xstart, times, func = "derivsc",
           parms = parms, dllname = "TB_model_v4",initforc = "forcc",
           forcings=force, initfunc = "parmsc", nout = 31,
           outnames = c("Total","Total_S","Total_SH","Total_SA","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                        "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM",
                        "CD4500","CD4350_500","CD4250_349","CD4200_249","CD4100_199","CD450_99","CD450",
                        "ART500","ART350_500","ART250_349","ART200_249","ART100_199","ART50_99","ART50"), method = rkMethod("rk34f")))


######## Some plots for testing things

# Plot of CD4 distribution #####################
# Get the CD4 outputs

temp <- as.data.frame(cbind(seq(1970,2070),out[,716:722]))
colnames(temp) <- c("Year",colnames(temp[,2:8]))
temp_CD4_4 <- melt(temp,id="Year")

plot_CD4 <- ggplot(temp_CD4_4,aes(x=Year,y=value,fill=variable))+
  geom_area(colour="black", size=.2, alpha=.4) +
  xlim(c(1970,2050))

# Plot pop against UN data ###################

# convert UN data to long format
temp_data <- melt(UN_pop_age_t,id="Year")

# sum up model outputs over age groups and turn into long format
tot<-mat.or.vec(101,17)
for(i in 1:17){
  tot[,i] <- apply(out,1,function(x) sum(x[seq(i+1,698,17)]))
}

temp_model1 <- as.data.frame(cbind(seq(1970,2070),tot,out[,"Total"]))
colnames(temp_model1) <- colnames(UN_pop_age_t)
temp_model1 <- melt(temp_model1,id="Year")

# and plot
plot_pop <- ggplot(temp_model,aes(x=Year,y=value))+
  geom_line(colour="red")+
  geom_point(data=temp_data,aes(x=Year,y=value))+
  geom_line(data=temp_model1,aes(x=Year,y=value),col="green",linetype=2)+
  facet_wrap(~variable,scales="free")+
  xlim(c(1970,2100))





