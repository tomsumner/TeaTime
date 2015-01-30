## Code to test approach for age structured model - currently limited to trying to get the demography to work

## Based on http://kinglab.eeb.lsa.umich.edu/EEID/eeid/2011_eco/waifw.pdf

## Set working directory
setwd("C:/Users/TOM SUMMER/Documents/TIME_research/TeaTime")

## Load packages
library(deSolve)
library(reshape2)
library(ggplot2)

# Set up age structure
ages <- c(seq(4,79,by=5),100) # upper end of age classes: 0-4 then five year bions to 79 then 80+
num_ages <- length(ages) # calculates the number of age classes
da <- diff(c(0,ages)) # calculates the widths of age classes

## Initial conditions - 1950 population 
yinit <- c(S=c(2088,1706,1485,1340,1183,1040,960,835,725,589,502,415,324,228,148,77,37))
sindex <- 1:num_ages

## To capture the aging process, define a matrix to hold the rates of movement between age classes.
aging <- diag(-1/da)
aging[row(aging)-col(aging)==1] <- 1/head(da,-1)
aging[num_ages,num_ages]=0   # set rate out of last box = 0 - will account for flow out of here via death rate

## Mortality rates
# Uses age specific LE tables from UN pop for South Africa, currently just the average over time, but could be made to vary with time
LE<-c(55.00,50.60,46.02,41.64,37.51,33.66,29.94,26.44,23.01,19.73,16.65,13.75,11.22,9.04,7.20,5.62,3.09)
death <- diag(1/LE)

## Fertility
# Uses crude birth rate (per 1000 population) from UN pop for South Africa, currently takes values at lower time band and linearly interpolates 
births <- approxfun(x=c(1950,1955,1960,1965,1970,1975,1980,1985,1990,1995,2000,2005),y=c(44,42,41,38,38,36,34,31,27,25,24,22),rule=2)

## Now we can put the pieces together to write a simulator for the demogrpahics
ja.multistage.model <- function (t, x, ...) {
  s <- x[sindex] # susceptibles

  Total <- sum(s)
  
  dsdt <- aging%*%s - death%*%s
  
  dsdt[1] <- dsdt[1]+births(t)*Total/1000
  
  list(c(dsdt),Total=Total)
}

# And run it
sol <- ode(y=yinit,times=seq(1950,2010,by=1),func=ja.multistage.model)

## Plot populations in each age group over time

# First load the UN pop estimates and convert to long format
UN_pop_age <- as.data.frame(read.table("SA_pop_age.txt",header=TRUE))
UN_pop_age <- cbind(UN_pop_age,rowSums(UN_pop_age[,2:18]))
colnames(UN_pop_age) <- c("time","0-4","5-9","10-14","15-19","20-24","25-29","30-34",
                          "35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80+","Total")
temp_data <- melt(UN_pop_age,id="time")

# then turn the model output into long format
temp <- rowSums(sol[,2:18])
temp_model <- as.data.frame(cbind(sol,temp))
colnames(temp_model) <- colnames(UN_pop_age)
temp_model <- melt(temp_model,id="time")

# and plot
plot_temp <- ggplot(temp_model,aes(x=time,y=value))+
                    geom_line(colour="red")+
                    geom_point(data=temp_data,aes(x=time,y=value))+
                    facet_wrap(~variable,scales="free")



