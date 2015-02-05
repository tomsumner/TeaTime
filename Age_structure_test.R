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
LE <-c(54.18,rep(57.24,4),rep(55.00,5),rep(50.60,5),rep(46.02,5),rep(41.64,5),rep(37.51,5),
       rep(33.66,5),rep(29.94,5),rep(26.44,5),rep(23.01,5),rep(19.73,5),rep(16.65,5),rep(13.75,5),
       rep(11.22,5),rep(9.04,5),rep(7.20,5),rep(5.62,5),rep(4.36,5),rep(3.36,5),rep(2.58,5),2.07)
death <- diag(1/LE)

mort <-c(0.1105,rep(0.0125,5),rep(0.00427,5),rep(0.003346,5),rep(0.005892,5),rep(0.0083542,5),rep(0.0096787,5),
         rep(0.0111635,5),rep(0.0132887,5),rep(0.0162737,5),rep(0.020528,5),rep(0.0274226,5),rep(0.0355126,5),
         rep(0.048119,5),rep(0.063236,5),rep(0.0813139,5),rep(0.1011768,5),rep(0.1219189,5),rep(0.066666667,15))
death <- diag(mort)

## Fertility
# Uses crude birth rate (per 1000 population) from UN pop for South Africa, currently takes values at lower time band and linearly interpolates 
birth_rate <- approxfun(x=c(1952.5,1957.5,1962.5,1967.5,1972.7,1977.5,1982.5,1987.5,1992.5,1997.5,2002.5,2007.5),y=c(44,42,41,38,38,36,34,31,27,25,24,22),rule=2)

## Now we can put the pieces together to write a simulator for the demogrpahics
ja.multistage.model <- function (t, x, ...) {

  s <- x[sindex] # susceptibles

  dsdt <- aging%*%s #- death%*%s
  
  Total <- sum(s)
  Births <- birth_rate(t)*Total/1000
  
  dsdt[1] <- dsdt[1]+birth_rate(t)*Total/1000
  
  list(c(dsdt),Total=Total,Births=Births)
}

# And run it
sol <- ode(y=yinit,times=seq(1950,2010,by=1),func=ja.multistage.model)

## Plot populations in each age group over time

# add total to data and convert to long format
UN_pop_age_t <- cbind(UN_pop_age,rowSums(UN_pop_age[,2:102]))
colnames(UN_pop_age_t) <- c("time",seq(0,100,by=1),"Total")
temp_data <- melt(UN_pop_age_t,id="time")

# then turn the model output into long format
temp_model <- as.data.frame(sol[,1:103])
colnames(temp_model) <- colnames(UN_pop_age_t)
temp_model <- melt(temp_model,id="time")

# and plot
plot_temp <- ggplot(temp_model,aes(x=time,y=value))+
                    geom_line(colour="red")+
                    geom_point(data=temp_data,aes(x=time,y=value))+
                    facet_wrap(~variable,scales="free")+
                    xlim(c(1950,2010))


