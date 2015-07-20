## Script to run the demog model with single year age groups
# Adapted to use events to determine births and aging - all births occur at the start of a year and all move to next age group at start of year

## Set working directory
setwd("C:/Users/TOM SUMMER/Documents/TIME_research/TeaTime")

## Load packages
library(deSolve)
library(reshape2)
library(ggplot2)
library(gdata)

dyn.unload("Demog_single.dll") # Unload - need to do this before recompiling
system("R CMD SHLIB Demog_single.c") # Compile
dyn.load("Demog_single.dll") # Load

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

# Combine forcing functions into a list
force <- list(birth_rate,
              s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,
              s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,s32,s33,s34,s35,s36,s37,s38,s39,s40,
              s41,s42,s43,s44,s45,s46,s47,s48,s49,s50,s51,s52,s53,s54,s55,s56,s57,s58,s59,s60,
              s61,s62,s63,s64,s65,s66,s67,s68,s69,s70,s71,s72,s73,s74,s75,s76,s77,s78,s79,s80,
              s81)

# the age1 and age2 parameters govern the width of the age bins (5 years for all except the final bin which we take to be 10 years (80-90))
parms <- c(age_1 = 1/1, age_2 = 1/20)

##############################################################################################################################
# Run the model 

# Times to run model for
times <- seq(1970,2050, by=1)

# Initial conditions - all susceptible
temp <- c()
for (i in 1:num_ages){temp[i]<-UN_pop_age[21,i+1]}
xstart <- c(S=c(temp))

# Run the model
time_eq <- system.time(out <- ode(y=xstart, times, func = "derivsc",
            parms = parms, dllname = "Demog_single",initforc = "forcc",
            forcings=force, initfunc = "parmsc", nout = 3,
            outnames = c("Total","Births","Deaths"), 
            events = list(func="event",time=seq(1971,2050)),
            method = rkMethod("rk34f")))

######## Some plots for testing things against data ##############################################

###################### POPULATION ###################################################################################

# convert UN data to long format
#temp_data <- melt(UN_pop_age_t,id="Year")
#temp_data_l <- melt(UN_pop_age_low_t,id="Year")
#temp_data_h <- melt(UN_pop_age_high_t,id="Year")

#temp_model <- as.data.frame(cbind(seq(1970,2050),out[,"Births"],out[,2:82],out[,"Total"],out[,"Deaths"]))
#colnames(temp_model) <- colnames(UN_pop_age_t)
#temp_model_m <- melt(temp_model,id="Year")

# and plot by single age bins
#plot_pop <- ggplot(temp_model_m,aes(x=Year,y=value))+
#  geom_line(colour="red")+
#  geom_line(data=temp_data,aes(x=Year,y=value),colour="black")+
#  facet_wrap(~variable,scales="free")+
#  ggtitle("Population, births and deaths")+
#  xlim(c(1970,2050))

##### Now sum up into 5 year bins used in TIME and compare ##########################################################

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

model_temp <- mat.or.vec(81,21)

model_temp[,1] <- times
model_temp[,2] <- out[,"Births"]
model_temp[,20] <- out[,"Total"]
model_temp[,21] <- out[,"Deaths"]
model_temp[,19] <- out[,"S81"]

for (i in 1:16){
  
  t1 <- (i+1)+(i-1)*4
  t2 <- t1+4
  model_temp[,i+2] <- rowSums(out[,t1:t2])
  
}
model_temp <- as.data.frame(model_temp)
colnames(model_temp) <- colnames(UN_pop_age_t)
temp_model_m <- melt(model_temp,id="Year")

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

############# DEATHS ################################################################################################

# TIME values
temp <- as.data.frame(read.table(paste("Demog/",cn,"_TIME_deaths_age.txt",sep=""),header=TRUE,fill=TRUE)) # Load HIV numbers (output in TIME)  
# Need to re-arrage to get in year vs age format
deaths_number_age <- mat.or.vec(81,18)
deaths_number_age[,1] <- seq(1970,2050)
for (i in 1:81){
  j <- (i-1)*19+2
  deaths_number_age[i,2:18]=temp[j:(j+16),2]
}
deaths_number_age <- as.data.frame(deaths_number_age)
colnames(deaths_number_age) <- colnames(UN_pop_age)
deaths_TIME <- melt(deaths_number_age,id="Year")

# Arrange model output
model_deaths <- as.data.frame(cbind(out[,1],out[,6455:6471]))
colnames(model_deaths) <- colnames(UN_pop_age)
model_deaths_m <- melt(model_deaths,id="Year")

# and plot
plot_deaths <- ggplot(model_deaths_m,aes(x=Year,y=value))+
  geom_line(colour="red")+
  geom_line(data=deaths_TIME,aes(x=Year,y=value/1000),colour="black",linetype="dashed")+
  facet_wrap(~variable,scales="free")+
  ggtitle("Deaths")+
  xlim(c(1970,2050))

# and also look at rates
# TIME
m_TIME <- as.data.frame(cbind(seq(1970,2050),100*deaths_number_age[,2:18]/TIME_pop[,3:19]))
colnames(m_TIME) <- colnames(UN_pop_age)
m_TIME_m <- melt(m_TIME,id="Year")
# Model
m_model <- as.data.frame(cbind(seq(1970,2050),100*model_deaths[,2:18]/tot))
colnames(m_model) <- colnames(UN_pop_age)
m_model_m <- melt(m_model,id="Year")

plot_death_p <- ggplot(m_model_m,aes(x=Year,y=value))+
  geom_line(colour="red")+
  geom_line(data=m_TIME_m,aes(x=Year,y=value),colour="black",linetype="dashed")+
  facet_wrap(~variable,scales="free")+
  ggtitle("Death_rate")+
  xlim(c(1970,2050))

####### HIV by ######################################################################################################

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
tot_HIV<-mat.or.vec(81,17)
for(i in 1:17){
  tot_HIV[,i] <- apply(out,1,function(x) sum(x[seq(i+222,6410,17)]))
}

HIV_model <- as.data.frame(cbind(seq(1970,2050),tot_HIV))
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
H_p_model <- as.data.frame(cbind(HIV_model[,1],100*tot_HIV/tot))
colnames(H_p_model) <- colnames(UN_pop_age)
H_p_model_m <- melt(H_p_model,id="Year")

plot_HIV_p <- ggplot(H_p_model_m,aes(x=Year,y=value))+
  geom_line(colour="red")+
  geom_line(data=H_p_TIME_m,aes(x=Year,y=value),colour="black",linetype="dashed")+
  facet_wrap(~variable,scales="free")+
  ggtitle("HIV_prevalence (age)")+
  xlim(c(1970,2050))

########  SAVE PLOTS TO A PDF FILE ##################################################################################
pdf_name <- paste("C:/Users/TOM SUMMER/Filr/My Files/sync/TIME/TIME Research/",cn,"_s5_mort_lowered.pdf",sep="")
pdf(pdf_name,width=10,height=7)
print(plot_pop)
print(plot_pop_s)
print(plot_deaths)
print(plot_death_p)
print(plot_HIV)
print(plot_HIV_p)
dev.off()

#####################################################################################################################


### EXTRA STUFF TO USE WHEN CHECKING TB AND CD4 OUTPUTS


# Plot of CD4 distribution against Spectrum - this is normalised by total population to account for differences in population size #####################

# Load TIME outputs
TIME_no_ART <- as.data.frame(read.table("TIME_no_ART_out.txt",header=TRUE))
TIME_no_ART <- cbind(TIME_no_ART,"No_ART","TIME")
colnames(TIME_no_ART) <- c(colnames(TIME_no_ART[1:8]),"Type","Model")

TIME_ART <- as.data.frame(read.table("TIME_ART_out.txt",header=TRUE))
TIME_ART <- cbind(TIME_ART,"ART","TIME")
colnames(TIME_ART) <- c(colnames(TIME_ART[1:8]),"Type","Model")

temp1 <- as.data.frame(cbind(out[,"time"],100*out[,6426:6432]/out[,"Total"]))
temp1 <- cbind(temp1,"No_ART","R")
colnames(temp1) <- c("Year",colnames(temp1[,2:8]),"Type","Model")

temp2 <- as.data.frame(cbind(seq(1970,2050),100*out[,6433:6439]/out[,"Total"]))
temp2 <- cbind(temp2,"ART","R")
colnames(temp2) <- c("Year",colnames(temp1[,2:8]),"Type","Model")

temp_CD4 <- rbind(temp1,temp2,TIME_no_ART,TIME_ART)
temp_CD4 <- melt(temp_CD4,id=c("Year","Type","Model"))

plot_CD4 <- ggplot(temp_CD4[temp_CD4$Model=="R",],aes(x=as.numeric(as.character(Year)),y=as.numeric(as.character(value)),colour=variable))+
  geom_line() +
  geom_line(data=temp_CD4[temp_CD4$Model=="TIME",],aes(x=as.numeric(as.character(Year)),y=as.numeric(as.character(value)),colour=variable),linetype="dashed")+
  facet_wrap(~Type)+
  xlim(c(1970,2050))



#### Compare to TIME output

# Load the TIME results (Incidence, Prevalence, Mortality)
TIME_out <- as.data.frame(read.table("TIME_out.txt",header=TRUE))
TIME_out <- cbind(TIME_out,"TIME")
colnames(TIME_out) <- c("Year","Prevalence","Incidence","Mortality","Model")
# Load the WHO data 
Data_out <- as.data.frame(read.table("Data_out.txt",header=TRUE))
Data_out <- cbind(Data_out,"Data")
colnames(Data_out) <- c("Year","Prevalence","Incidence","Mortality","Model")

# Arrange model output
R_out <- as.data.frame(cbind(out[,"time"],
      100000*(out[,"Total_DS"]+out[,"Total_MDR"])/out[,"Total"],  # Prev
      100000*(out[,"Cases_neg"]+out[,"Cases_pos"]+out[,"Cases_ART"])/out[,"Total"], # Inc 
      100000*out[,"TB_deaths"]/out[,"Total"])) # Mort (all)
R_out <- cbind(R_out,"R")
colnames(R_out) <- c("Year","Prevalence","Incidence","Mortality","Model")
# combine and melt
Models_out <- rbind(TIME_out,R_out,Data_out)
Models_out <- melt(Models_out,id=c("Year","Model"))

plot_models <- ggplot(Models_out[Models_out$Model!="Data",],aes(x=Year,y=value,colour=Model))+
  facet_wrap(~variable,scales="free")+
  geom_line()+
  geom_point(data=Models_out[Models_out$Model=="Data",],aes(x=Year,y=value),colour="black")+
  xlim(c(1970,2050))














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



