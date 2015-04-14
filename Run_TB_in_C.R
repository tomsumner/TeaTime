## Script to run the TB model as compiled C code

## Set working directory
setwd("C:/Users/TOM SUMMER/Documents/TIME_research/TeaTime")

## Load packages
library(deSolve)
library(reshape2)
library(ggplot2)

## Compile and load the C code
system("R CMD SHLIB TB_model_v2.c") # Compile
dyn.load("TB_model_v2.dll") # Load
dyn.unload("TB_model_v2.dll") # Unload - need to do this before recompiling

##############################################################################################################################

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

# Combine forcing functions into a list
force <- list(birth_rate,s_birth,s5,s10,s15,s20,s25,s30,s35,s40,s45,s50,s55,s60,s65,s70,s75,s80)

# run model for different values of beta

beta = seq(5,25,by=0.5)
prev <- mat.or.vec(101,length(beta))

system.time(for (run_no in 1:length(beta)){

##############################################################################################################################
# Model initialisation
# run the model from 1970 pop with 1970 birth/death rates and care and control parameters for 400 years to get stable age structure
# Then rerun with 100 TB cases at time 0, no MDR or HIV. Run for 400 years to get stable disease state

# Times to run model for
times <- seq(0,400 , by=1)

# Parameters - proportion primary, proportion smear pos and mortality rates depend on age as in TIME
fit_cost=0.7
g = fit_cost/(1+fit_cost) # superinfections 

parms <- c(age1 = 1/5, age2 = 1/21, beta = beta[run_no], 
           a_a = 0.14, a0 = 0.26432, a5 = 0.14056, a10 = 0.056,  
           p = 0.65, v = 0.001, 
           sig_a = 0.45, sig0 = 0.0684, sig5 = 0.0414, sig10 = 0.0846, rel_inf = 0.25, theta = 0.02, r = 0.25, 
           mu_N = 0.25, mu_N0 = 0.426, mu_I = 0.35, mu_I0 = 0.59, fit_cost = 0.7, e = 0, g=g, k = 0.3, l_s = 0.83, l_m = 0.0, d = 0.8, tau_s = 0.76, tau_m = 0.0,
           eff_n = 0.0, eff_p = 0.0, dst_n = 0.0, dst_p = 0.0) 

# Initial conditions - all susceptible
UN_pop_age <- as.data.frame(read.table("SA_pop_age.txt",header=TRUE)) # Load UN Population data
# add total to data
UN_pop_age_t <- cbind(UN_pop_age,rowSums(UN_pop_age[,2:18]))
colnames(UN_pop_age_t) <- c(colnames(UN_pop_age),"Total")

temp <- c()
for (i in 1:num_ages){temp[i]<-UN_pop_age[21,i+1]}
xstart <- c(S=c(temp),
            Lsn=rep(0,num_ages),Lsp=rep(0,num_ages),Lmn=rep(0,num_ages),Lmp=rep(0,num_ages),
            Nsn=rep(0,num_ages),Nsp=rep(0,num_ages),Nmn=rep(0,num_ages),Nmp=rep(0,num_ages),
            Isn=rep(0,num_ages),Isp=rep(0,num_ages),Imn=rep(0,num_ages),Imp=rep(0,num_ages))

# Run the model

system.time(out_pop <- ode(y=xstart, times, func = "derivsc",
            parms = parms, dllname = "TB_model_v2",initforc = "forcc",
            forcings=force, initfunc = "parmsc", nout = 15,
            outnames = c("Total","Total_S","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                         "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM"), method = rkMethod("rk34f")))

# Update initial conditions based on end of last run and add 100 TB cases
temp <- c()
for (i in 1:num_ages){temp[i]<-out_pop[401,i+1]}
xstart <- c(S=c(temp),
            Lsn=rep(0,num_ages),Lsp=rep(0,num_ages),Lmn=rep(0,num_ages),Lmp=rep(0,num_ages),
            Nsn=rep(0,num_ages),Nsp=rep(0,num_ages),Nmn=rep(0,num_ages),Nmp=rep(0,num_ages),
            Isn=c(rep(0,5),0.1,rep(0,11)),Isp=rep(0,num_ages),Imn=rep(0,num_ages),Imp=rep(0,num_ages))

# Run the model

system.time(out_TB <- ode(y=xstart, times, func = "derivsc",
                       parms = parms, dllname = "TB_model_v2",initforc = "forcc",
                       forcings=force, initfunc = "parmsc", nout = 15,
                       outnames = c("Total","Total_S","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                                    "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM"), method = rkMethod("rk34f")))


# Is TB in equilibrium? # Looks like it's okay but is way higher than TIME

#plot(100*out_TB[,"Total_DS"]/out_TB[,"Total"]) # prevalence (%)

# Is age structure in equilibrium and in line with data? # Looks like it's okay
#pop_sum <- mat.or.vec(dim(out_TB)[1],num_ages)
#for (i in 1:num_ages){
  
#  temp <- seq((i+1),(12*17)+(i+1),17)
#  for (j in 1:length(temp)) pop_sum[,i] <- pop_sum[,i]+out_TB[,temp[j]]
  
#}
#total_pop <- rowSums(pop_sum)
#pop_prop <- pop_sum/total_pop

#par(mfrow=c(4,5))
#for (i in 1:17){
#  plot(pop_prop[,i],ylim=c(0,0.2))
#  abline(h=UN_pop_age_t[21,i+1]/UN_pop_age_t[21,19],col="red")
#}

# Adjust pop down to 1970 values and reassign initial conditions - model can now be run from 1970
temp <- out_TB[dim(out_TB)[1],2:222]
temp <- temp/(sum(temp)/22502) # 22502 is total pop from UN estimates in 1970
xstart <- temp

##############################################################################################################################

# Set times to run for
times <- seq(1970,2070 , by=1)
# Run the model
time_TB<-system.time(out <- ode(y=xstart, times, func = "derivsc",
           parms = parms, dllname = "TB_model_v2",initforc = "forcc",
           forcings=force, initfunc = "parmsc", nout = 15,
           outnames = c("Total","Total_S","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                        "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM"), method = rkMethod("rk34f")))

prev[,run_no]<-100*out[,"Total_DS"]/out[,"Total"]

})






# Plot pop against UN data

# add total to data and convert to long format
UN_pop_age_t <- cbind(UN_pop_age,rowSums(UN_pop_age[,2:18]))
colnames(UN_pop_age_t) <- c(colnames(UN_pop_age),"Total")
temp_data <- melt(UN_pop_age_t,id="Year")

# sum up model outputs over age groups and turn into long format
tot<-mat.or.vec(131,17)
for (i in 1:17){
  temp <- seq(i+1,i+205,17)
  for (j in 1:13){
    tot[,i]<-tot[,i]+out[,temp[j]]
  }
}

temp_model1 <- as.data.frame(cbind(seq(1970,2070),tot,out[,"Total"]))
colnames(temp_model1) <- colnames(UN_pop_age_t)
temp_model1 <- melt(temp_model1,id="Year")

# and plot
plot_pop <- ggplot(temp_model,aes(x=Year,y=value))+
  geom_line(colour="red")+
  geom_point(data=temp_data,aes(x=Year,y=value))+
  geom_line(data=temp_model1,aes(x=Year,y=value),colour="green")+
  facet_wrap(~variable,scales="free")+
  xlim(c(1970,2100))

lines(100*out[,"Total_DS"]/out[,"Total"]) # prevalence (%)

plot(100*out[,"Total_DS"]/out[,"Total"]) # prevalence (%)
