## Code to test approach for age structured model - currently limited to trying to get the demography to work

## Based on http://kinglab.eeb.lsa.umich.edu/EEID/eeid/2011_eco/waifw.pdf

## Set working directory
setwd("C:/Users/TOM SUMMER/Documents/TIME_research/TeaTime")

## Load packages
library(deSolve)
library(reshape2)
library(ggplot2)

# Set up age structure
ages <- seq(1,101,by=1) # upper end of age classes: 1 year age bins to 100
num_ages <- length(ages) # calculates the number of age classes
da <- diff(c(0,ages)) # calculates the widths of age classes

## Initial conditions - 1950 population
UN_pop_age <- as.data.frame(read.table("SA_pop_age.txt",header=TRUE)) # Load UN Population data
temp <- c()
for (i in 1:101){temp[i]<-UN_pop_age[1,i+1]}
yinit <- c(S=c(temp))
sindex <- 1:num_ages





## To capture the aging process, define a matrix to hold the rates of movement between age classes.
aging <- diag(-1/da)
aging[row(aging)-col(aging)==1] <- 1/head(da,-1)
aging[num_ages,num_ages]=0   # set rate out of last box = 0 - will account for flow out of here via death rate

## Mortality rates
# Uses age specific LE tables from UN pop for South Africa, currently just the average over time, but could be made to vary with time
#LE<-c(55.00)
#death <- diag(1/LE)

## Fertility
# Uses crude birth rate (per 1000 population) from UN pop for South Africa, currently takes values at lower time band and linearly interpolates 
births <- approxfun(x=c(1952.5,1957.5,1962.5,1967.5,1972.7,1977.5,1982.5,1987.5,1992.5,1997.5,2002.5,2007.5),y=c(44,42,41,38,38,36,34,31,27,25,24,22),rule=2)

## Now we can put the pieces together to write a simulator for the demogrpahics
ja.multistage.model <- function (t, x, ...) {
  s <- x[sindex] # susceptibles

  Total <- sum(s)
  
  dsdt <- aging%*%s #- death%*%s
  
  #dsdt[1] <- dsdt[1]+births(t)*Total/1000
  
  list(c(dsdt),Total=Total)
}

# And run it
sol <- ode(y=yinit,times=seq(1950,2010,by=1),func=ja.multistage.model)

## Plot populations in each age group over time

# add total to data and convert to long format
UN_pop_age <- cbind(UN_pop_age,rowSums(UN_pop_age[,2:102]))
colnames(UN_pop_age) <- c("time",seq(0,100,by=1),"Total")
temp_data <- melt(UN_pop_age,id="time")

# then turn the model output into long format
temp_model <- as.data.frame(sol)
colnames(temp_model) <- colnames(UN_pop_age)
temp_model <- melt(temp_model,id="time")

# and plot
plot_temp <- ggplot(temp_model,aes(x=time,y=value))+
                    geom_line(colour="red")+
                    geom_point(data=temp_data,aes(x=time,y=value))+
                    facet_wrap(~variable,scales="free")



