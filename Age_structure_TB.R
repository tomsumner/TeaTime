## Code for age structured TB model

## Age structure is based on http://kinglab.eeb.lsa.umich.edu/EEID/eeid/2011_eco/waifw.pdf

## Set working directory
setwd("C:/Users/TOM SUMMER/Documents/TIME_research/TeaTime")

## Load packages
library(deSolve)
library(reshape2)
library(ggplot2)

# Set up age structure
ages <- c(4,9,14,19,24,29,34,39,44,49,54,59,64,69,74,79,100) # upper end of age classes
num_ages <- length(ages) # calculates the number of age classes
da <- diff(c(0,ages)) # calculates the widths of age classes
da [1] <- da[1]+1 # add one year to first age group as start of group is birth (i.e. 0) not 1

## Initial conditions - 1950 population, 25% infected, 1% prevalence
UN_pop_age <- as.data.frame(read.table("SA_pop_age.txt",header=TRUE)) # Load UN Population data
temp <- c()
for (i in 1:num_ages){temp[i]<-UN_pop_age[1,i+1]}
yinit <- c(S=c(0.9999*temp),
           L_s_n=rep(0,num_ages),L_s_p=rep(0,num_ages),L_m_n=rep(0,num_ages),L_m_p=rep(0,num_ages),
           N_s_n=rep(0,num_ages),N_s_p=rep(0,num_ages),N_m_n=rep(0,num_ages),N_m_p=rep(0,num_ages),
           I_s_n=c(0.0001*temp),I_s_p=rep(0,num_ages),I_m_n=rep(0,num_ages),I_m_p=rep(0,num_ages))
s_index <- 1:num_ages
Lsn_index <- (num_ages+1):(2*num_ages)
Lsp_index <- (2*num_ages+1):(3*num_ages)
Lmn_index <- (3*num_ages+1):(4*num_ages)
Lmp_index <- (4*num_ages+1):(5*num_ages)
Nsn_index <- (5*num_ages+1):(6*num_ages)
Nsp_index <- (6*num_ages+1):(7*num_ages)
Nmn_index <- (7*num_ages+1):(8*num_ages)
Nmp_index <- (8*num_ages+1):(9*num_ages)
Isn_index <- (9*num_ages+1):(10*num_ages)
Isp_index <- (10*num_ages+1):(11*num_ages)
Imn_index <- (11*num_ages+1):(12*num_ages)
Imp_index <- (12*num_ages+1):(13*num_ages)

## Fertility
# Uses crude birth rate (per 1000 population) from UN pop for South Africa which has values for 5 year periods
# Puts these values at midpoint of period and interpolates between 
# 2010-2015 (2012.5) value based on medium fertility projection
birth_rate <- approxfun(x=seq(1952.5,2047.5,5),
                        y=c(44,42,41,38,38,36,34,31,27,25,24,22,21,20,18,17,17,16,15,14),rule=2)

## Survival
Survive_age <- as.data.frame(read.table("SA_survival_age.txt",header=TRUE)) # Load survival proportions calculated from life tables
# Proportion surviving from age 0 to 1 - used to determine entry to first age group
s_birth <- approxfun(x=Survive_age$Year,y=Survive_age$X1,rule=2)
# and to other ages
s5 <- approxfun(x=Survive_age$Year,y=Survive_age$X5,rule=2)
s10 <- approxfun(x=Survive_age$Year,y=Survive_age$X10,rule=2)
s15 <- approxfun(x=Survive_age$Year,y=Survive_age$X15,rule=2)
s20 <- approxfun(x=Survive_age$Year,y=Survive_age$X20,rule=2)
s25 <- approxfun(x=Survive_age$Year,y=Survive_age$X25,rule=2)
s30 <- approxfun(x=Survive_age$Year,y=Survive_age$X30,rule=2)
s35 <- approxfun(x=Survive_age$Year,y=Survive_age$X35,rule=2)
s40 <- approxfun(x=Survive_age$Year,y=Survive_age$X40,rule=2)
s45 <- approxfun(x=Survive_age$Year,y=Survive_age$X45,rule=2)
s50 <- approxfun(x=Survive_age$Year,y=Survive_age$X50,rule=2)
s55 <- approxfun(x=Survive_age$Year,y=Survive_age$X55,rule=2)
s60 <- approxfun(x=Survive_age$Year,y=Survive_age$X60,rule=2)
s65 <- approxfun(x=Survive_age$Year,y=Survive_age$X65,rule=2)
s70 <- approxfun(x=Survive_age$Year,y=Survive_age$X70,rule=2)
s75 <- approxfun(x=Survive_age$Year,y=Survive_age$X75,rule=2)
s80 <- approxfun(x=Survive_age$Year,y=Survive_age$X80,rule=2)

## To model the aging process, define a matrix to hold the rates of movement between age classes.
aging <- diag(-1/da) # This is aging out of class
aging[row(aging)-col(aging)==1] <- 1/head(da,-1) # This is movement into next class assuming 100% survival

# Parameter values

# Core natural history
beta = 10       # Contact rate
a = 0.115       # Proportion developing primary TB
p = 0.65        # Protection due to previous infection
v = 0.0001      # Reactivation rate
sig = 0.62      # Proportion smear positive
rel_inf = 0.22  # relative infectiousness of smear negative TB
theta = 0.015   # smear conversion rate
r = 0.2         # selfcure rate
mu_N = 0.2      # smear negative mortality
mu_I = 0.3      # smear positive mortality

# MDR
fit_cost = 0.7            # fitness cost for MDR
e = 0.014                 # MDR acquistion rate
g = fit_cost/(1+fit_cost) # superinfections 

# Care and control parameters

k = 0     # diagnositic rate (assummed to be the same for DS and MDR cases)         
l_s = 0   # linkage to care (DS TB)
l_m = 0   # linkage to care (MDR TB)
d = 0     # relative detection of smear negative cases
tau_s = 0 # treatment success for DS TB
tau_m = 0 # treatment success for MDR TB
eff_n = 0 # efficacy of first line treatment in new MDR cases
eff_p = 0 # efficacy of first line treatment in retreatment MDR cases
dst_n = 0 # DST coverage in new TB cases
dst_p = 0 # DST coverage in retreatment TB cases

## Now we can put the pieces together to write a simulator TB
ja.multistage.model <- function (t, x, ...) {

  # Notation
  # S: susceptible, L: Latently infected, N: Smear negative TB, I: Smear positive TB
  # s: drug, susceptible, m: drug resistant, n: treatment naive, p: previoulsy treated
  
  S <- x[s_index]
  Lsn <- x[Lsn_index]; Lsp <- x[Lsp_index]; Lmn <- x[Lmn_index]; Lmp <- x[Lmp_index]
  Nsn <- x[Nsn_index]; Nsp <- x[Nsp_index]; Nmn <- x[Nmn_index]; Nmp <- x[Nmp_index]
  Isn <- x[Isn_index]; Isp <- x[Isp_index]; Imn <- x[Imn_index]; Imp <- x[Imp_index]
  
  # create a vector of survial at time t
  surv <- c(s5(t),s10(t),s15(t),s20(t),s25(t),s30(t),s35(t),s40(t),s45(t),s50(t),s55(t),s60(t),s65(t),s70(t),s75(t),s80(t))
  aging_temp <- aging
  aging_temp[row(aging_temp)-col(aging_temp)==1] <- surv*1/head(da,-1)
  
  # Expressions for total populations - also want age group totals at some point
  Total_S <- sum(S)
  
  Total_L_DS <- sum(Lsn) + sum(Lsp)
  Total_L_MDR <- sum(Lmn) + sum(Lmp)
  Total_L <- Total_L_DS + Total_L_MDR
  
  Total_N_DS <- sum(Nsn) + sum(Nsp)
  Total_N_MDR <- sum(Nmn) + sum(Nmp)
  Total_N <- Total_N_DS + Total_N_MDR
  
  Total_I_DS <- sum(Isn) + sum(Isp)
  Total_I_MDR <- sum(Imn) + sum(Imp)
  Total_I <- Total_I_DS + Total_I_MDR
 
  Total_DS <- Total_N_DS + Total_I_DS
  Total_MDR <- Total_N_MDR + Total_I_MDR

  Total <- Total_S + Total_L + Total_N + Total_I
  
  # Expressions for force of infection acounting for reduced infectiousness of smear negative cases
  FS <- beta*(Total_N_DS*rel_inf + Total_I_DS)/Total # Susceptible
  FM <- fit_cost*beta*(Total_N_MDR*rel_inf + Total_I_MDR)/Total # Resistant
  
  # Differential equations
  dSdt <- aging_temp%*%S - (FS + FM)*S

  dLsndt <- aging_temp%*%Lsn + FS*((1-a)*S + (1-a)*(1-p)*(1-g)*Lmn) - 
                               (v + FS*a*(1-p) + FM*a*(1-p) + FM*(1-a)*(1-p)*g)*Lsn + 
                               r*(Isn + Nsn) 

  dLspdt <- aging_temp%*%Lsp + FS*(1-a)*(1-p)*(1-g)*Lmp - 
                               (v + FS*a*(1-p) + FM*a*(1-p) + FM*(1-a)*(1-p)*g)*Lsp + 
                               r*(Isp + Nsp) + k*l_s*(1-e)*tau_s*(Isn + Isp) + k*l_s*(1-e)*tau_s*d*(Nsn + Nsp) 

  dLmndt <- aging_temp%*%Lmn + FM*((1-a)*S + (1-a)*(1-p)*g*Lsn) - 
                               (v + FM*a*(1-p) + FS*a*(1-p) + FS*(1-a)*(1-p)*(1-g))*Lmn + 
                               r*(Imn + Nmn)

  dLmpdt <- aging_temp%*%Lmp + FM*(1-a)*(1-p)*g*Lsp - 
                               (v + FM*a*(1-p) + FS*a*(1-p) + FS*(1-a)*(1-p)*(1-g))*Lmp + 
                               r*(Imp + Nmp) + 
                               k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(Imn + d*Nmn) + 
                               k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(Imp + d*Nmp)

  dNsndt <- aging_temp%*%Nsn + (v*(1-sig) + FS*a*(1-p)*(1-sig))*Lsn + 
                               FS*a*(1-sig)*(S + (1-p)*Lmn) - 
                               (theta + r + k*l_s*d + mu_N)*Nsn

  dNspdt <- aging_temp%*%Nsp + (v*(1-sig) + FS*a*(1-p)*(1-sig))*Lsp + 
                               FS*a*(1-sig)*(1-p)*Lmp - 
                               (theta + r + k*l_s*d*(1-e)*tau_s + k*l_s*d*e + mu_N)*Nsp + 
                               k*l_s*d*(1-e)*(1-tau_s)*Nsn

  dNmndt <- aging_temp%*%Nmn + (v*(1-sig) + FM*a*(1-p)*(1-sig))*Lmn + 
                               FM*a*(1-sig)*(S + (1-p)*Lsn) - 
                               (theta + r + k*l_m*d*dst_n + k*l_s*d*(1-dst_n) + mu_N)*Nmn
  

  dNmpdt <- aging_temp%*%Nmp + (v*(1-sig) + FM*a*(1-p)*(1-sig))*Lmp + 
                               FM*a*(1-sig)*(1-p)*Lsp - 
                               (theta + r + k*l_m*d*dst_p*tau_m + k*l_s*d*(1-dst_p)*tau_s*eff_p + mu_N)*Nmp + 
                               k*l_s*d*e*(Nsn + Nsp) + (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*Nmn

  dIsndt <- aging_temp%*%Isn + (v*sig + FS*a*sig*(1-p))*Lsn + 
                               FS*a*sig*S + FS*a*(1-p)*sig*Lmn + 
                               theta*Nsn - (r + k*l_s + mu_I)*Isn

  dIspdt <- aging_temp%*%Isp + (v*sig + FS*a*sig*(1-p))*Lsp + 
                               FS*a*sig*(1-p)*Lmp + theta*Nsp + 
                               k*l_s*(1-e)*(1-tau_s)*Isn - 
                               (r + k*l_s*(1-e)*tau_s + k*l_s*e + mu_I)*Isp

  dImndt <- aging_temp%*%Imn + (v*sig + FM*a*sig*(1-p))*Lmn + 
                               FM*a*sig*S + FM*a*(1-p)*sig*Lsn + 
                               theta*Nmn - (r + k*l_m*dst_n + k*l_s*(1-dst_n) + mu_I)*Imn

  dImpdt <- aging_temp%*%Imp + (v*sig + FM*a*sig*(1-p))*Lmp + 
                               FM*a*sig*(1-p)*Lsp + theta*Nmp + 
                               (k*l_m*dst_n*(1-tau_m) + k*l_s*(1-dst_n)*(1-(tau_s*eff_n)))*Imn + 
                               k*l_s*e*(Isn + Isp) - 
                               (r + k*l_m*dst_p*tau_m + k*l_s*(1-dst_p)*tau_s*eff_p + mu_I)*Imp                                                                   
  
  Births <- s_birth(t)*birth_rate(t)*Total/1000
  
  dSdt[1] <- dSdt[1]+s_birth(t)*birth_rate(t)*Total/1000
  
  list(c(dSdt,
         dLsndt,dLspdt,dLmndt,dLmpdt,
         dNsndt,dNspdt,dNmndt,dNmpdt,
         dIsndt,dIspdt,dImndt,dImpdt),
         FM=FM,FS=FS,Births=Births,
         Total=Total,Total_S=Total_S,Total_L=Total_L,
         Total_I=Total_I,Total_N=Total_N,Total_DS=Total_DS,Total_MDR=Total_MDR)
}

# And run it
system.time(sol <- ode(y=yinit,times=seq(1950,2250,by=1),func=ja.multistage.model))




