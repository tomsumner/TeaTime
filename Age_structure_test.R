## Code to test approach for age structured model
## and initialisation of age structure
## 
## Start model with 1970 pop and rates and run for 400 years to get equilibrium pop structure
## Adjust pop down to 1970 size 
## Run forward to 2050 and check fit

## Based on http://kinglab.eeb.lsa.umich.edu/EEID/eeid/2011_eco/waifw.pdf

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

## Initial conditions - 1970 population
UN_pop_age <- as.data.frame(read.table("SA_pop_age.txt",header=TRUE)) # Load UN Population data
temp <- c()
for (i in 1:num_ages){temp[i]<-UN_pop_age[21,i+1]}
yinit <- c(S=c(temp))
sindex <- 1:num_ages

## Fertility
# Uses crude birth rate (per 1000 population) from UN pop for South Africa which has values for 5 year periods
# Puts these values at midpoint of period and interpolates between 
# 2010-2015 (2012.5) value based on medium fertility projection
birth_rate <- approxfun(x=seq(1972.5,2047.5,5),
                        y=c(38,36,34,31,27,25,24,22,21,20,18,17,17,16,15,14),rule=2)

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

## Now we can put the pieces together to write a simulator for the demogrpahics
ja.multistage.model <- function (t, x, ...) {

  s <- x[sindex] # susceptibles
  
  # create a vector of survial at time t
  surv <- c(s5(t),s10(t),s15(t),s20(t),s25(t),s30(t),s35(t),s40(t),s45(t),s50(t),s55(t),s60(t),s65(t),s70(t),s75(t),s80(t))
  aging_temp <- aging
  aging_temp[row(aging_temp)-col(aging_temp)==1] <- surv*1/head(da,-1)
  
  dsdt <- aging_temp%*%s
  
  Total <- sum(s)
    
  Births <- s_birth(t)*birth_rate(t)*Total/1000
  
  dsdt[1] <- dsdt[1]+s_birth(t)*birth_rate(t)*Total/1000
  
  list(c(dsdt),Total=Total,Births=Births)

}

# Run it for 400 years
sol <- ode(y=yinit,times=seq(0,400,by=1),func=ja.multistage.model)

# Adjust pop down to 1970 values and reassign initial conditions
temp <- sol[dim(sol)[1],2:18]
temp <- temp/(sum(temp)/sol[1,19])

temp2 <- c()
for (i in 1:num_ages){temp2[i]<-temp[i]}
yinit <- c(S=c(temp2))
sindex <- 1:num_ages

# Run from 1970 onwards
sol <- ode(y=yinit,times=seq(1970,2050,by=1),func=ja.multistage.model)


## Plot populations in each age group over time

# add total to data and convert to long format
UN_pop_age_t <- cbind(UN_pop_age,rowSums(UN_pop_age[,2:18]))
colnames(UN_pop_age_t) <- c(colnames(UN_pop_age),"Total")
temp_data <- melt(UN_pop_age_t,id="Year")

# then turn the model output into long format
temp_model <- as.data.frame(sol[,1:19])
colnames(temp_model) <- colnames(UN_pop_age_t)
temp_model <- melt(temp_model,id="Year")

temp2 <- sol[,2:18]/sol[,19]
temp2 <- as.data.frame(cbind(sol[,1],temp2))
colnames(temp2) <- colnames(UN_pop_age_t[1:18])
temp_2 <- melt(temp2,id="Year")

# and plot
plot_temp <- ggplot(temp_model,aes(x=Year,y=value))+
                    geom_line(colour="red")+
                    geom_point(data=temp_data,aes(x=Year,y=value))+
                    facet_wrap(~variable,scales="free")+
                    xlim(c(1950,2050))

plot_temp <- ggplot(temp_model,aes(x=Year,y=value/sum(value)))+
  geom_line(colour="red")+
  #geom_point(data=temp_data,aes(x=Year,y=value))+
  facet_wrap(~variable,scales="free")+
  xlim(c(0,400))+
  ylim(c(0,0.0005))

