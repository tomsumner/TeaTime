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
UN_pop_age_low <- as.data.frame(read.table("SA_pop_age_low.txt",header=TRUE)) # Load UN Population data
UN_pop_age_high <- as.data.frame(read.table("SA_pop_age_high.txt",header=TRUE)) # Load UN Population data
# add total to data
UN_pop_age_t <- cbind(UN_pop_age,rowSums(UN_pop_age[,2:18]))
colnames(UN_pop_age_t) <- c(colnames(UN_pop_age),"Total")
UN_pop_age_low_t <- cbind(UN_pop_age_low,rowSums(UN_pop_age_low[,2:18]))
colnames(UN_pop_age_low_t) <- c(colnames(UN_pop_age_low),"Total")
UN_pop_age_high_t <- cbind(UN_pop_age_high,rowSums(UN_pop_age_high[,2:18]))
colnames(UN_pop_age_high_t) <- c(colnames(UN_pop_age_high),"Total")

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

# HIV Incidence by age and year - based on AIM output, but ignoring childhood infections (this is what Carel does in TIME)
HIV_Inc_age <- as.data.frame(read.table("HIV_Inc_age.txt",header=TRUE)) # Load HIV incidence data taken from AIM                                       # Data from AIM is rate per 1000 
#HIV_Inc_age[,2:18]=HIV_Inc_age[,2:18]*0

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

# ART coverage - based on AIM, use CD4 eligibility threshold and % of those in need on ART
ART_data <- as.data.frame(read.table("ART_data.txt",header=TRUE)) # Load data
# Create forcing function of threshold category
Athresh <- cbind(ART_data[,"Year"],ART_data[,"CD4_cat"])


# Create forcing functions which account for threshold and coverage
A50 <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>50,1,0)/100)
A99 <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>99,1,0)/100)
A199 <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>199,1,0)/100)
A249 <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>249,1,0)/100)
A349 <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>349,1,0)/100)
A500 <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>499,1,0)/100)
Ahigh <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>500,1,0)/100)

# Pop adjust - to turn off population adjust for TB/HIV deaths from 2015 onwards
pop_ad <- cbind(c(2014,2015,2016),c(1,0,0))

# BCG coverage - currently assume 90% at all times
BCG_cov <- cbind(c(1972,1973,2050),c(0.9,0.9,0.9))

# Case detection rate - generalised logistic function, don't think this quite matches TIME 
k <- cbind(seq(1970,2050),(30 + ((120-30)/((1+exp(-0.5*(seq(1970,2050)-2004)))^(1/2))))/100)

#k <- cbind(c(1970,1990,2050),c(0.30,0.30,1.20))

# DST coverage among new and previously treated cases
dst_n <- cbind(seq(1970,2050),(0 + ((95-0)/((1+exp(-1*(seq(1970,2050)-1993)))^(1/2))))/100)
dst_p <- cbind(seq(1970,2050),(0 + ((95-0)/((1+exp(-1*(seq(1970,2050)-1993)))^(1/2))))/100)

# Combine forcing functions into a list
force <- list(birth_rate,s_birth,s5,s10,s15,s20,s25,s30,s35,s40,s45,s50,s55,s60,s65,s70,s75,s80,
              h0,h5,h10,h15,h20,h25,h30,h35,h40,h45,h50,h55,h60,h65,h70,h75,h80,
              Ahigh,A500,A349,A249,A199,A99,A50,Athresh,
              BCG_cov,pop_ad,k,dst_n,dst_p)

# Set up TB parameters

# Fitness of MDR, used to calculate parameter for superinfections
# Both of these are passed into "parms" together with e, the MDR acquisition rate (we set this to zero to exclude MDR in equilibirum phase)
fit_cost=0.7
g = fit_cost/(1+fit_cost) # superinfections 
e = 0.01

# proportion primary (a), proportion smear pos (sig) and mortality rates (muN and mu_I) take different values for 
# adults (>15) (_a), 0-4 (_0), 5-9 (_5) and 10-14 (_10)

parms <- c(age1 = 1/5, age2 = 1/21, beta = 18, 
           a_a = 0.14, a0 = 0.26432, a5 = 0.14056, a10 = 0.056,  
           p = 0.65, v = 0.001, 
           sig_a = 0.45, sig0 = 0.0684, sig5 = 0.0414, sig10 = 0.0846, rel_inf = 0.25, theta = 0.02, r = 0.25, 
           mu_N = 0.25, mu_N0 = 0.426, mu_I = 0.35, mu_I0 = 0.59, fit_cost = fit_cost, e = e, g=g, l_s = 0.83, l_m = 0.7, d = 0.8, tau_s = 0.76, tau_m = 0.5,
           eff_n = 0.61, eff_p = 0.45, 
           muN_H = 0.45, muI_H = 0.6, RR1a = 2, RR2a = 1.288, RR1v = 3, RR2v = 3, RR1p = 0.5, RR2p = 1.1,
           ART_TB1 = 0.7, ART_TB2 = 0.5, ART_TB3 = 0.35, ART_mort1 = 0.5, ART_mort2 = 0.4, ART_mort3 = 0.3,
           #ART_TB1 = 1, ART_TB2 = 1, ART_TB3 = 1, ART_mort1 = 1, ART_mort2 = 1, ART_mort3 = 1,
           BCG_eff = 0.39,
           sig_H = 0.35,r_H=0.15)

##############################################################################################################################
# Run the model 

# First run the model from 1970 pop with 1970 birth/death rates and care and control parameters with 100 TB cases (no MDR or HIV) 
# for 100 years to get stable age structure and disease state

# Times to run model for
times <- seq(0,100, by=1)

# Initial conditions - all susceptible
temp <- c()
for (i in 1:num_ages){temp[i]<-UN_pop_age[21,i+1]}
xstart <- c(S=c(temp),
            Lsn=rep(0,num_ages),Lsp=rep(0,num_ages),Lmn=rep(0,num_ages),Lmp=rep(0,num_ages),
            Nsn=rep(0,num_ages),Nsp=rep(0,num_ages),Nmn=rep(0,num_ages),Nmp=rep(0,num_ages),
            Isn=c(rep(0,5),0,rep(0,11)),Isp=rep(0,num_ages),Imn=rep(0,num_ages),Imp=rep(0,num_ages),
            S_H=rep(0,num_ages*7),
            Lsn_H=rep(0,num_ages*7),Lsp_H=rep(0,num_ages*7),Lmn_H=rep(0,num_ages*7),Lmp_H=rep(0,num_ages*7),
            Nsn_H=rep(0,num_ages*7),Nsp_H=rep(0,num_ages*7),Nmn_H=rep(0,num_ages*7),Nmp_H=rep(0,num_ages*7),
            Isn_H=rep(0,num_ages*7),Isp_H=rep(0,num_ages*7),Imn_H=rep(0,num_ages*7),Imp_H=rep(0,num_ages*7),
            S_A=rep(0,num_ages*7*3),
            Lsn_A=rep(0,num_ages*7*3),Lsp_A=rep(0,num_ages*7*3),Lmn_A=rep(0,num_ages*7*3),Lmp_A=rep(0,num_ages*7*3),
            Nsn_A=rep(0,num_ages*7*3),Nsp_A=rep(0,num_ages*7*3),Nmn_A=rep(0,num_ages*7*3),Nmp_A=rep(0,num_ages*7*3),
            Isn_A=rep(0,num_ages*7*3),Isp_A=rep(0,num_ages*7*3),Imn_A=rep(0,num_ages*7*3),Imp_A=rep(0,num_ages*7*3))

# For initialisation run turn off MDR by setting e = 0
parms["e"]=0

# Run the model
time_eq <- system.time(out_eq <- ode(y=xstart, times, func = "derivsc",
            parms = parms, dllname = "TB_model_v4",initforc = "forcc",
            forcings=force, initfunc = "parmsc", nout = 44,
            outnames = c("Total","Total_S","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                         "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM",
                         "CD4500","CD4350_500","CD4250_349","CD4200_249","CD4100_199","CD450_99","CD450",
                         "ART500","ART350_500","ART250_349","ART200_249","ART100_199","ART50_99","ART50",
                         "v1","v2","v3","v4","v5","v6","v7","ART_tot","ART_need","ART_new","ART_on","TB_deaths",
                         "Cases_neg","Cases_pos","Cases_ART"), method = rkMethod("rk34f")))

# Adjust pop down to 1970 values and reassign initial conditions - model can now be run from 1970 with TB and HIV
temp <- out_eq[dim(out_eq)[1],2:6410]
temp <- temp/(sum(temp)/22502) # 22502 is total pop from UN estimates in 1970)
xstart <- temp

# Reset e to allow MDR
parms["e"]=e

# Set times to run for
times <- seq(1970,2050 , by=1)
# Run the model
time_run <-system.time(out <- ode(y=xstart, times, func = "derivsc",
           parms = parms, dllname = "TB_model_v4",initforc = "forcc",
           forcings=force, initfunc = "parmsc", nout = 44,
           outnames = c("Total","Total_S","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                        "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM",
                        "CD4500","CD4350_500","CD4250_349","CD4200_249","CD4100_199","CD450_99","CD450",
                        "ART500","ART350_500","ART250_349","ART200_249","ART100_199","ART50_99","ART50",
                        "v1","v2","v3","v4","v5","v6","v7","ART_tot","ART_need","ART_new","ART_on","TB_deaths",
                        "Cases_neg","Cases_pos","Cases_ART"), method = rkMethod("rk34f")))

######## Some plots for testing things against spectrum

# Plot of CD4 distribution #####################

temp1 <- as.data.frame(cbind(seq(1970,2050),out[,6426:6432],"No ART"))
colnames(temp1) <- c("Year",colnames(temp1[,2:8]),"Type")
temp2 <- as.data.frame(cbind(seq(1970,2050),out[,6433:6439],"ART"))
colnames(temp2) <- c("Year",colnames(temp1[,2:8]),"Type")
temp <- rbind(temp1,temp2)
temp_CD4 <- melt(temp,id=c("Year","Type"))

plot_CD4 <- ggplot(temp_CD4,aes(x=as.numeric(as.character(Year)),y=as.numeric(as.character(value)),fill=variable))+
  geom_area(colour="black", size=.2, alpha=.4) +
  facet_wrap(~Type)+
  xlim(c(1970,2050))

# Plot pop against UN data ###################

# convert UN data to long format
temp_data <- melt(UN_pop_age_t,id="Year")
temp_data_l <- melt(UN_pop_age_low_t,id="Year")
temp_data_h <- melt(UN_pop_age_high_t,id="Year")

# sum up model outputs over age groups and turn into long format
tot<-mat.or.vec(81,17)
for(i in 1:17){
  tot[,i] <- apply(out,1,function(x) sum(x[seq(i+1,6410,17)]))
}

temp_model <- as.data.frame(cbind(seq(1970,2050),tot,out[,"Total"]))
colnames(temp_model) <- colnames(UN_pop_age_t)
temp_model <- melt(temp_model,id="Year")

# and plot
plot_pop <- ggplot(temp_model,aes(x=Year,y=value))+
  geom_line(colour="red")+
  geom_line(data=temp_data,aes(x=Year,y=value),colour="black")+
  geom_line(data=temp_data_l,aes(x=Year,y=value),colour="black",linetype="dashed")+
  geom_line(data=temp_data_h,aes(x=Year,y=value),colour="black",linetype="dashed")+
  facet_wrap(~variable,scales="free")+
  xlim(c(1970,2100))



# Arrange some outputs to take out to excel (all in numbers)

cbind(out[,"time"],
      1000*out[,"Total"],
      1000*(out[,"Total_DS"]),  # Prev (DS)
      1000*(out[,"Total_MDR"]), # Prev (MDR)
      1000*(out[,"Cases_neg"]), # Inc (neg)         
      1000*(out[,"Cases_pos"]), # Inc (pos)  
      1000*(out[,"Cases_ART"]), # Inc (art)  
      1000*(out[,"TB_deaths"]), # Mort (all)
      1000*(out[,"Total_L"]))   # LTBI (all)
      
# distribution CD4 no ART
temp <- rbind(out[,"time"],1000*out[,"CD4500"],1000*out[,"CD4350_500"],1000*out[,"CD4250_349"],1000*out[,"CD4200_249"],
                   1000*out[,"CD4100_199"],1000*out[,"CD450_99"],1000*out[,"CD450"])
write.table(temp,file="CD4_no_ART.txt",sep=" ")

# distribution CD4 with ART
temp <- rbind(out[,"time"],1000*out[,"ART500"],1000*out[,"ART350_500"],1000*out[,"ART250_349"],1000*out[,"ART200_249"],
              1000*out[,"ART100_199"],1000*out[,"ART50_99"],1000*out[,"ART50"])
write.table(temp,file="CD4_ART.txt",sep=" ")

# HIV prevalence 15+
cbind(out[,"time"],100*(1000*out[,"CD4500"]+1000*out[,"CD4350_500"]+1000*out[,"CD4250_349"]+1000*out[,"CD4200_249"]+
      1000*out[,"CD4100_199"]+1000*out[,"CD450_99"]+1000*out[,"CD450"]+1000*out[,"ART500"]+1000*out[,"ART350_500"]+
      1000*out[,"ART250_349"]+1000*out[,"ART200_249"]+1000*out[,"ART100_199"]+1000*out[,"ART50_99"]+1000*out[,"ART50"])/(1000*rowSums(tot[,4:17])))

# Number of HIV positives on and off ART - we currently ignore childhood HIV so equivalent to 15+ in TIME
cbind(out[,"time"],1000*out[,"CD4500"]+1000*out[,"CD4350_500"]+1000*out[,"CD4250_349"]+1000*out[,"CD4200_249"]+
                          1000*out[,"CD4100_199"]+1000*out[,"CD450_99"]+1000*out[,"CD450"],1000*out[,"ART500"]+1000*out[,"ART350_500"]+
                          1000*out[,"ART250_349"]+1000*out[,"ART200_249"]+1000*out[,"ART100_199"]+1000*out[,"ART50_99"]+1000*out[,"ART50"])

cbind(out[,"time"],1000*out[,"ART_new"])
