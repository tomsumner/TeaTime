## Script to run the TB model as compiled C code

## Set working directory
setwd("C:/Users/TOM SUMMER/Documents/TIME_research/TeaTime")

## Load packages
library(deSolve)
library(reshape2)
library(ggplot2)
library(gdata)

dyn.unload("TB_model_v6.dll") # Unload - need to do this before recompiling
system("R CMD SHLIB TB_model_v6.c") # Compile
dyn.load("TB_model_v6.dll") # Load

##############################################################################################################################

cn <- "SA"

## Load UN population data
UN_pop_age <- as.data.frame(read.table(paste("Demog/",cn,"_pop_single.txt",sep=""),header=FALSE)) # Load UN Population data
#UN_pop_age_low <- as.data.frame(read.table(paste("Demog/",cn,"_pop_age_low.txt",sep=""),header=TRUE)) # Load UN Population data
#UN_pop_age_high <- as.data.frame(read.table(paste("Demog/",cn,"_pop_age_high.txt",sep=""),header=TRUE)) # Load UN Population data
# Load number of births
births <- as.data.frame(read.table(paste("Demog/",cn,"_births_number.txt",sep=""),header=TRUE))
# Load number of deaths
deaths <- as.data.frame(read.table(paste("Demog/",cn,"_deaths.txt",sep=""),header=TRUE))
# add total, births and deaths to data
UN_pop_age_t <- cbind(births,UN_pop_age[,2:82],rowSums(UN_pop_age[,2:82]),deaths[,2])
colnames(UN_pop_age_t) <- c("Year","births",colnames(UN_pop_age[2:82]),"Total","Deaths")
#UN_pop_age_low_t <- cbind(UN_pop_age_low,rowSums(UN_pop_age_low[,2:18]))
#colnames(UN_pop_age_low_t) <- c(colnames(UN_pop_age_low),"Total")
#UN_pop_age_high_t <- cbind(UN_pop_age_high,rowSums(UN_pop_age_high[,2:18]))
#colnames(UN_pop_age_high_t) <- c(colnames(UN_pop_age_high),"Total")

# Set up age structure
ages <- c(seq(1,80),100) # upper end of age classes
num_ages <- length(ages) # calculates the number of age classes

# Set up the forcing functions for birth and death - all from 1970 onwards

## Fertility
# Uses crude birth rate (per 1000 population) from UN pop 
# values post 2010 based on medium fertiliy
temp <- as.data.frame(read.table(paste("Demog/",cn,"_crude_birth.txt",sep=""),header=TRUE)) # Load crude birth rate
birth_rate <- cbind(temp[,1],temp[,2])

# Load mortality rates taken from demproj

mort_age <- as.data.frame(read.table(paste("Demog/",cn,"_mort_age.txt",sep=""),header=FALSE)) 
for (i in 1:81){  
  assign(paste("s",i,sep=""), cbind(seq(1971,2050),mort_age[,i+1]))
}

################### JUST PLAYING WITH MORTALITY PARAMETERS TO TRY AND IMPROVE FIT ############################
################### NEED TO THINK ABOUT WHETHER TO CHANGE THESE ##############################################
################### COULD FIT DEMOG MODEL WITHOUT TB TO GET THESE  ###########################################
#s81[,2] <- s81[,2]/2
#s1[,2] <- s1[,2]*0.5

# HIV Incidence by age and year
temp <- as.data.frame(read.table(paste("HIV/",cn,"_HIV_Inc_age.txt",sep=""),header=TRUE,fill=TRUE)) # Load HIV incidence data taken from AIM   
# Need to re-arrage to get in year vs age format
HIV_Inc_age <- mat.or.vec(81,18)
HIV_Inc_age[,1] <- seq(1970,2050)
for (i in 1:81){
  j <- (i-1)*19+4
  HIV_Inc_age[i,2:18]=as.numeric(levels(temp[j:(j+16),2]))[temp[j:(j+16),2]]
}

# Data from AIM is rate per 1000 
#HIV_Inc_age[,2:18]=HIV_Inc_age[,2:18]*0   ### This line will turn HIV off 
temp <- c(rep(c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),each=5),18)   # used to duplicate the HIV incidence (reported by 5 year age bins) for single years
for (i in 0:80){
  assign(paste("h",i,sep=""),cbind(seq(1970,2050),HIV_Inc_age[,temp[i+1]]/1000))
}

# ART coverage - based on AIM, use CD4 eligibility threshold and % of those in need on ART
ART_data <- as.data.frame(read.table(paste("HIV/",cn,"_ART_data.txt",sep=""),header=TRUE)) # Load data
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

# Case detection rate by HIV (neg/pos)
# generalised logistic function - NEED TO FIND A NEAT WAY TO DO THIS
base=50
target=90
growth=0.3
shape=1
ybase=1990
ytarget=2020
y=mean(c(ybase,ytarget))
  
kneg <- cbind(seq(1970,2050),(base + ((target-base)/((1+exp(-growth*(seq(1970,2050)-y)))^(1/shape))))/100)
kpos <- cbind(seq(1970,2050),(base + ((target-base)/((1+exp(-growth*(seq(1970,2050)-y)))^(1/shape))))/100)

# Relative detection smear neg
rel_d <- cbind(seq(1970,2050),0.8)

# DST coverage among new and previously treated cases
base=0
target=95
growth=1
shape=2
ybase=1975
ytarget=2010
y=mean(c(ybase,ytarget))

dstneg_n <- cbind(seq(1970,2050),(base + ((target-base)/((1+exp(-growth*(seq(1970,2050)-y)))^(1/shape))))/100)
dstneg_p <- cbind(seq(1970,2050),(base + ((target-base)/((1+exp(-growth*(seq(1970,2050)-y)))^(1/shape))))/100)
dstpos_n <- cbind(seq(1970,2050),(base + ((target-base)/((1+exp(-growth*(seq(1970,2050)-y)))^(1/shape))))/100)
dstpos_p <- cbind(seq(1970,2050),(base + ((target-base)/((1+exp(-growth*(seq(1970,2050)-y)))^(1/shape))))/100) 

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
              kneg,kpos,
              rel_d,
              dstneg_n,dstneg_p,dstpos_n,dstpos_p,
              l_s,l_m,
              tneg_s,tpos_s,tART_s,tneg_m,tpos_m,tART_m)

# Set up TB parameters

# Fitness of MDR, used to calculate parameter for superinfections
# Both of these are passed into "parms" together with e, the MDR acquisition rate (we set this to zero to exclude MDR in equilibirum phase)
fit_cost=0.72
g = fit_cost/(1+fit_cost) # superinfections 
e = 0.01

# proportion primary (a), proportion smear pos (sig) and mortality rates (muN and mu_I) take different values for 
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

##############################################################################################################################
# Run the model 

# First run the model from 1970 pop with 1970 birth/death rates and care and control parameters with 100 TB cases (no MDR or HIV) 
# for 100 years to get stable age structure and disease state

# Times to run model for
times <- seq(0,60, by=1)

# Initial conditions - all susceptible
temp <- c()
for (i in 1:num_ages){temp[i]<-UN_pop_age[21,i+1]}
xstart <- c(S=c(temp),
            Lsn=rep(0,num_ages),Lsp=rep(0,num_ages),Lmn=rep(0,num_ages),Lmp=rep(0,num_ages),
            Nsn=rep(0,num_ages),Nsp=rep(0,num_ages),Nmn=rep(0,num_ages),Nmp=rep(0,num_ages),
            Isn=c(rep(0,25),100,rep(0,55)),Isp=rep(0,num_ages),Imn=rep(0,num_ages),Imp=rep(0,num_ages),
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
            parms = parms, dllname = "TB_model_v6",initforc = "forcc",
            forcings=force, initfunc = "parmsc", nout = 46,
            outnames = c("Total","Total_S","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                         "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM",
                         "CD4500","CD4350_500","CD4250_349","CD4200_249","CD4100_199","CD450_99","CD450",
                         "ART500","ART350_500","ART250_349","ART200_249","ART100_199","ART50_99","ART50",
                         "v1","v2","v3","v4","v5","v6","v7","ART_tot","ART_need","ART_new","ART_on","TB_deaths",
                         "Cases_neg","Cases_pos","Cases_ART",
                         "births","deaths"), 
            events = list(func="event",time=seq(0,60)),
            method = rkMethod("rk4")))

# Adjust pop down to 1970 values (by age) and reassign initial conditions - model can now be run from 1970 with TB and HIV

temp <- out_eq[dim(out_eq)[1],2:30538]

for(i in 1:81){ 
  temp[seq(i,30537,81)] <- temp[seq(i,30537,81)]/(sum(temp[seq(i,30537,81)])/UN_pop_age_t[UN_pop_age_t$Year==1970,i+2])
}

xstart <- temp

# Reset e to allow MDR
parms["e"]=e

# Set times to run for
times <- seq(1970,2050 , by=0.5) # run with 6 month time step using a fixed time step solver - this is faster than adaptive methds but seems to give good accuracy
# Run the model
time_run <-system.time(out <- ode(y=xstart, times, func = "derivsc",
           parms = parms, dllname = "TB_model_v6",initforc = "forcc",
           forcings=force, initfunc = "parmsc", nout = 46,
           outnames = c("Total","Total_S","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                        "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM",
                        "CD4500","CD4350_500","CD4250_349","CD4200_249","CD4100_199","CD450_99","CD450",
                        "ART500","ART350_500","ART250_349","ART200_249","ART100_199","ART50_99","ART50",
                        "v1","v2","v3","v4","v5","v6","v7","ART_tot","ART_need","ART_new","ART_on","TB_deaths",
                        "Cases_neg","Cases_pos","Cases_ART",
                        "births","deaths"), 
           events = list(func="event",time=seq(1970,2050)),
           method = rkMethod("rk4")))

# Just keep every other output now we are running with 6 month time step
out <- out[seq(1,length(times),2),]

######## Some plots for testing things against spectrum and other data ##############################################

#### Compare to TIME TB output

# Load the TIME results (Incidence, Prevalence, Mortality)
TIME_out <- as.data.frame(read.table(paste("TB/",cn,"_TIME_TB.txt",sep=""),header=TRUE,fill=TRUE))
TIME_out <- cbind(TIME_out,"TIME")
colnames(TIME_out) <- c("Year","Prevalence","Incidence","Mortality","Model")

# Arrange model output
R_out <- as.data.frame(cbind(out[,"time"],
                             100000*(out[,"Total_DS"]+out[,"Total_MDR"])/out[,"Total"],  # Prev
                             100000*(out[,"Cases_neg"]+out[,"Cases_pos"]+out[,"Cases_ART"])/out[,"Total"], # Inc 
                             100000*out[,"TB_deaths"]/out[,"Total"])) # Mort (all)
R_out <- cbind(R_out,"R")
colnames(R_out) <- c("Year","Prevalence","Incidence","Mortality","Model")
# combine and melt
Models_out <- rbind(TIME_out,R_out)
Models_out <- melt(Models_out,id=c("Year","Model"))

# Calculate % diff between models and add to df
temp=100*(Models_out[Models_out$Model=="TIME",4]-Models_out[Models_out$Model=="R",4])/Models_out[Models_out$Model=="TIME",4]
Models_out<-cbind(temp,Models_out)

plot_models <- ggplot(Models_out[Models_out$Model!="Data",],aes(x=Year,y=value,colour=Model,linetype=variable))+
  geom_line()+
  xlim(c(1970,2035))+
  ylim(c(0,1100))

plot_diff <- ggplot(Models_out[Models_out$Model=="TIME",],aes(x=Year,y=temp,linetype=variable))+
  geom_line()+
  xlim(c(1970,2050))+
  ylim(c(-20,20))


head(Models_out)






##### NEED TO UPDATE ALL OF THESE TO SUM UP SINGLE YEARS INTO CORRESPONDING 5 YEAR BINS

###################### POPULATION ###################################################################################

## Load and manipulate UN numbers 
## Load UN population data
UN_pop_age <- as.data.frame(read.table(paste("Demog/",cn,"_pop_age.txt",sep=""),header=TRUE)) # Load UN Population data
# Load number of births
births <- as.data.frame(read.table(paste("Demog/",cn,"_births_number.txt",sep=""),header=TRUE))
# Load number of deaths
deaths <- as.data.frame(read.table(paste("Demog/",cn,"_deaths.txt",sep=""),header=TRUE))
# add total, births and deaths to data
UN_pop_age_t <- cbind(births,UN_pop_age[,2:18],rowSums(UN_pop_age[,2:18]),deaths[,2])
colnames(UN_pop_age_t) <- c("Year","births",colnames(UN_pop_age[2:18]),"Total","Deaths")
temp_data <- melt(UN_pop_age_t,id="Year")

# sum up model outputs over age groups and turn into long format
tot <- mapply(function(x,y) sum(out[x,seq(y+1,30538,81)]),rep(seq(1,81),each=81),seq(1,81))
dim(tot) <- c(81,81)
tot <- t(tot)

model_temp <- mat.or.vec(81,21)
model_temp[,1] <- out[,"time"]
model_temp[,2] <- out[,"births"]
model_temp[,20] <- out[,"Total"]
model_temp[,21] <- out[,"deaths"]
model_temp[,19] <- tot[,81]

# Then aggregate into 5 year bins used in TIME
for (i in 1:16){
  t1 <- (i)+(i-1)*4
  t2 <- t1+4
  model_temp[,i+2] <- rowSums(tot[,t1:t2])
}
temp_model <- as.data.frame(model_temp)
colnames(temp_model) <- colnames(UN_pop_age_t)
temp_model_m <- melt(temp_model,id="Year")

# TIME values
temp <- as.data.frame(read.table(paste("Demog/",cn,"_TIME_pop_age.txt",sep=""),header=TRUE,fill=TRUE)) 
# Need to re-arrage to get in year vs age format
TIME_pop <- mat.or.vec(81,19)
TIME_pop[,1] <- seq(1970,2050)
for (i in 1:81){
  j <- (i-1)*19+2
  TIME_pop[i,2:19]=temp[j:(j+17),2]
}
TIME_births <- as.data.frame(read.table(paste("Demog/",cn,"_TIME_births.txt",sep=""),header=TRUE))
TIME_deaths <- as.data.frame(read.table(paste("Demog/",cn,"_TIME_deaths.txt",sep=""),header=TRUE))
TIME_pop <- cbind(TIME_births,TIME_pop[,2:19],TIME_deaths[,2])
colnames(TIME_pop) <- colnames(UN_pop_age_t)
TIME_pop_m <- melt(TIME_pop,id="Year")

# and plot
plot_pop <- ggplot(temp_model_m,aes(x=Year,y=value))+
  geom_line(colour="red")+
  geom_line(data=temp_data,aes(x=Year,y=value),colour="black")+
  #geom_line(data=temp_data_l,aes(x=Year,y=value),colour="black",linetype="dashed")+
  #geom_line(data=temp_data_h,aes(x=Year,y=value),colour="black",linetype="dashed")+
  geom_line(data=TIME_pop_m,aes(x=Year,y=value/1000),colour="green",linetype="dashed")+
  facet_wrap(~variable,scales="free")+
  ggtitle("Population, births and deaths")+
  xlim(c(1970,2050))

# Compare population age structure (i.e % of total pop)
temp_data_s <- as.data.frame(cbind(UN_pop_age_t[,1],100*UN_pop_age_t[,3:20]/UN_pop_age_t[,20]))
colnames(temp_data_s)<-c("Year","x4","X9","X14","X19","X24","X29","X34","X39","X44","X49","X54","X59","X64","X69","X74","X79","X100","Total")
temp_data_s <- melt(temp_data_s,id="Year")

temp_model_s <- as.data.frame(cbind(seq(1970,2050),100*temp_model[,3:20]/temp_model[,20]))
colnames(temp_model_s)<-c("Year","x4","X9","X14","X19","X24","X29","X34","X39","X44","X49","X54","X59","X64","X69","X74","X79","X100","Total")
temp_model_s <- melt(temp_model_s,id="Year")

temp_TIME_s <- as.data.frame(cbind(seq(1970,2050),100*TIME_pop[,3:20]/TIME_pop[,20]))
colnames(temp_TIME_s)<-c("Year","x4","X9","X14","X19","X24","X29","X34","X39","X44","X49","X54","X59","X64","X69","X74","X79","X100","Total")
temp_TIME_s <- melt(temp_TIME_s,id="Year")

# and plot
plot_pop_s <- ggplot(temp_model_s,aes(x=Year,y=value))+
  geom_line(colour="red")+
  geom_line(data=temp_data_s,aes(x=Year,y=value),colour="black")+
  geom_line(data=temp_TIME_s,aes(x=Year,y=value),colour="green",linetype="dashed")+
  facet_wrap(~variable,scales="free")+
  ggtitle("age structure of population (% in age group)")+
  xlim(c(1970,2050))

####### HIV ######################################################################################################

# Load TIME values
temp <- as.data.frame(read.table(paste("HIV/",cn,"_HIV_numbers_age.txt",sep=""),header=TRUE,fill=TRUE)) # Load HIV numbers (output in TIME)  
# Need to re-arrage to get in year vs age format
HIV_number_age <- mat.or.vec(81,18)
HIV_number_age[,1] <- seq(1970,2050)
for (i in 1:81){
  j <- (i-1)*19+2
  HIV_number_age[i,2:18]=temp[j:(j+16),2]
}
HIV_number_age <- as.data.frame(HIV_number_age)
colnames(HIV_number_age) <- colnames(UN_pop_age)
HIV_TIME <- melt(HIV_number_age,id="Year")


# Sum up model outputs for HIV+ categories by age
tot_HIV <- mapply(function(x,y) sum(out[x,seq(y+1054,30538,81)]),rep(seq(1,81),each=81),seq(1,81))
dim(tot_HIV) <- c(81,81)
tot_HIV <- t(tot_HIV)

HIV_model <- mat.or.vec(81,18)
HIV_model[,1] <- seq(1970,2050)
HIV_model[,18] <- tot_HIV[,81]
for (i in 1:16){
  t1 <- (i)+(i-1)*4
  t2 <- t1+4
  HIV_model[,i+1] <- rowSums(tot_HIV[,t1:t2])
}

HIV_model <- as.data.frame(HIV_model )
colnames(HIV_model) <- colnames(UN_pop_age)
HIV_model_m <- melt(HIV_model,id="Year")

# and plot
plot_HIV <- ggplot(HIV_model_m,aes(x=Year,y=value))+
  geom_line(colour="red")+
  geom_line(data=HIV_TIME,aes(x=Year,y=value/1000),colour="black",linetype="dashed")+
  facet_wrap(~variable,scales="free")+
  ggtitle("HIV_numbers")+
  xlim(c(1970,2050))

### Compare HIV prevalence by age group

# TIME
H_p_TIME <- as.data.frame(cbind(HIV_model[,1],100*HIV_number_age[,2:18]/TIME_pop[,3:19]))
colnames(H_p_TIME) <- colnames(UN_pop_age)
H_p_TIME_m <- melt(H_p_TIME,id="Year")
# Model
H_p_model <- as.data.frame(cbind(HIV_model[,1],100*HIV_model[,2:18]/model_temp[,3:19]))
colnames(H_p_model) <- colnames(UN_pop_age)
H_p_model_m <- melt(H_p_model,id="Year")

plot_HIV_p <- ggplot(H_p_model_m,aes(x=Year,y=value))+
  geom_line(colour="red")+
  geom_line(data=H_p_TIME_m,aes(x=Year,y=value),colour="black",linetype="dashed")+
  facet_wrap(~variable,scales="free")+
  ggtitle("HIV_prevalence (age)")+
  xlim(c(1970,2050))




########  SAVE PLOTS TO A PDF FILE ##################################################################################
pdf_name <- paste("C:/Users/TOM SUMMER/Filr/My Files/sync/TIME/TIME Research/",cn,"_TB.pdf",sep="")
pdf(pdf_name,width=10,height=7)
print(plot_pop)
print(plot_pop_s)
print(plot_deaths)
print(plot_death_p)
print(plot_HIV)
print(plot_HIV_p)
dev.off()

#####################################################################################################################


### EXTRA STUFF TO USE WHEN CHECKING CD4 OUTPUTS


# Plot of CD4 distribution against Spectrum - this is normalised by total population to account for differences in population size #####################

# Load TIME outputs
TIME_no_ART <- as.data.frame(read.table(paste("HIV/",cn,"_TIME_no_ART.txt",sep=""),header=TRUE,fill=TRUE))
TIME_no_ART <- cbind(TIME_no_ART[,1],100*TIME_no_ART[,2:8]/TIME_pop[,"Total"],"No_ART","TIME")
colnames(TIME_no_ART) <- c("Year",colnames(TIME_no_ART[2:8]),"Type","Model")

TIME_ART <- as.data.frame(read.table(paste("HIV/",cn,"_TIME_on_ART.txt",sep=""),header=TRUE,fill=TRUE))
TIME_ART <- cbind(TIME_ART[,1],100*TIME_ART[,2:8]/TIME_pop[,"Total"],"ART","TIME")
colnames(TIME_ART) <- c("Year",colnames(TIME_ART[2:8]),"Type","Model")

temp1 <- as.data.frame(cbind(out[,"time"],100*out[,30554:30560]/out[,"Total"]))
temp1 <- cbind(temp1,"No_ART","R")
colnames(temp1) <- c("Year",colnames(temp1[,2:8]),"Type","Model")

temp2 <- as.data.frame(cbind(seq(1970,2050),100*out[,30561:30567]/out[,"Total"]))
temp2 <- cbind(temp2,"ART","R")
colnames(temp2) <- c("Year",colnames(temp1[,2:8]),"Type","Model")

temp_CD4 <- rbind(temp1,temp2,TIME_ART,TIME_no_ART)
temp_CD4 <- melt(temp_CD4,id=c("Year","Type","Model"))

plot_CD4 <- ggplot(temp_CD4[temp_CD4$Model=="R",],aes(x=as.numeric(as.character(Year)),y=as.numeric(as.character(value)),colour=variable))+
  geom_line() +
  geom_line(data=temp_CD4[temp_CD4$Model=="TIME",],aes(x=as.numeric(as.character(Year)),y=as.numeric(as.character(value)),colour=variable),linetype="dashed")+
  facet_wrap(~Type)+
  xlim(c(2000,2020))
















# Arrange some outputs to take out to excel (all in numbers)
cbind(out[,"time"],
      1000*(out[,"Total_DS"]),  # Prev (DS)
      1000*(out[,"Total_MDR"]), # Prev (MDR)
      1000*(out[,"Cases_neg"]), # Inc (neg)         
      1000*(out[,"Cases_pos"]), # Inc (pos)  
      1000*(out[,"Cases_ART"]), # Inc (art)  
      1000*(out[,"TB_deaths"]), # Mort (all)
      1000*(out[,"Total_L"]),   # LTBI (all)
      1000*out[,"Total"])


# distribution CD4 no ART
temp <- rbind(out[,"time"],1000*out[,"CD4500"],1000*out[,"CD4350_500"],1000*out[,"CD4250_349"],1000*out[,"CD4200_249"],
                   1000*out[,"CD4100_199"],1000*out[,"CD450_99"],1000*out[,"CD450"])
write.table(temp,file="CD4_no_ART.txt",sep=" ")

# distribution CD4 with ART
temp <- rbind(out[,"time"],1000*out[,"ART500"],1000*out[,"ART350_500"],1000*out[,"ART250_349"],1000*out[,"ART200_249"],
              1000*out[,"ART100_199"],1000*out[,"ART50_99"],1000*out[,"ART50"])
write.table(temp,file="CD4_ART.txt",sep=" ")

# Population, Number of HIV positives on and off ART (we currently ignore childhood HIV so equivalent to 15+ in TIME), HIV prevalence 15+, new ART
cbind(out[,"time"],1000*out[,"Total"],1000*out[,"CD4500"]+1000*out[,"CD4350_500"]+1000*out[,"CD4250_349"]+1000*out[,"CD4200_249"]+
        1000*out[,"CD4100_199"]+1000*out[,"CD450_99"]+1000*out[,"CD450"],1000*out[,"ART500"]+1000*out[,"ART350_500"]+
        1000*out[,"ART250_349"]+1000*out[,"ART200_249"]+1000*out[,"ART100_199"]+1000*out[,"ART50_99"]+1000*out[,"ART50"],100*(1000*out[,"CD4500"]+1000*out[,"CD4350_500"]+1000*out[,"CD4250_349"]+1000*out[,"CD4200_249"]+
      1000*out[,"CD4100_199"]+1000*out[,"CD450_99"]+1000*out[,"CD450"]+1000*out[,"ART500"]+1000*out[,"ART350_500"]+
      1000*out[,"ART250_349"]+1000*out[,"ART200_249"]+1000*out[,"ART100_199"]+1000*out[,"ART50_99"]+1000*out[,"ART50"])/(1000*rowSums(tot[,4:17])),1000*out[,"ART_new"])



