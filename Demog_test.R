## Demog_test.R
#
## Attempt to write demographic model as difference equtaions type thing

## Load packages
library(reshape2)
library(ggplot2)

## Set working directory
setwd("C:/Users/TOM SUMMER/Documents/TIME_research/TeaTime")

ages <- seq(0,100,by=1) # upper end of age classes: 1 year age bins to 100
num_ages <- length(ages) # calculates the number of age classes

## Load UN_pop data
UN_pop_age <- read.table("SA_pop_age.txt",header=TRUE) # Load UN Population data

## Create array to store model and populate 1950 values - 60 years (1950-2010), 101 ages (0 to 100)
pop <- mat.or.vec(61,num_ages)
for (i in 1:101){
  pop[1,i] <- UN_pop_age[1,i+1] 
}

## Fertility
# Uses crude birth rate (per 1000 population) from UN pop for South Africa, currently takes values at lower time band and linearly interpolates 
birth_rate <- approxfun(x=c(1952.5,1957.5,1962.5,1967.5,1972.7,1977.5,1982.5,1987.5,1992.5,1997.5,2002.5,2007.5),y=c(44,42,41,38,38,36,34,31,27,25,24,22),rule=2)

# mortality - based on life table for SA - currently fixed over time
#mort <-c(0.1105,rep(0.0125,5),rep(0.00427,5),rep(0.003346,5),rep(0.005892,5),rep(0.0083542,5),rep(0.0096787,5),
#         rep(0.0111635,5),rep(0.0132887,5),rep(0.0162737,5),rep(0.020528,5),rep(0.0274226,5),rep(0.0355126,5),
#         rep(0.048119,5),rep(0.063236,5),rep(0.0813139,5),rep(0.1011768,5),rep(0.1219189,5),rep(0.066666667,15))

mort <- read.table("SA_survival_age.txt",header=FALSE) # Load UN Population data
mort_col <- c("V2",rep("V3",5),rep("V4",5),rep("V5",5),rep("V6",5),rep("V7",5),rep("V8",5),rep("V9",5),rep("V10",5),
              rep("V11",5),rep("V12",5),rep("V13",5),rep("V14",5),rep("V15",5),rep("V16",5),rep("V17",5),rep("V18",5),
              rep("V19",5),rep("V20",15))

for (i in 1:60){
  
  death <- approxfun(x=mort[,"V1"],y=mort[,mort_col[1]],rule=2)
  pop[i+1,1] <- (1-death(i+1950))*birth_rate(i+1950)*sum(pop[i,])/1000 
  
  for (j in 2:101){
 
    death <- approxfun(x=mort[,"V1"],y=mort[,mort_col[j]],rule=2)
    pop[i+1,j] <- pop[i,j-1]*(1-death(i+1950))
  
  }
}
## Plot populations in each age group over time

# add total to data and convert to long format
UN_pop_age_t <- cbind(UN_pop_age,rowSums(UN_pop_age[,2:102]))
colnames(UN_pop_age_t) <- c("time",seq(0,100,by=1),"Total")
temp_data <- melt(UN_pop_age_t,id="time")

# then turn the model output into long format
temp_model <- cbind(seq(1950,2010,1),pop,rowSums(pop))
temp_model <- as.data.frame(temp_model)
colnames(temp_model) <- colnames(UN_pop_age_t)
temp_model <- melt(temp_model,id="time")

# and plot
plot_temp <- ggplot(temp_model,aes(x=time,y=value))+
  geom_line(colour="red")+
  geom_point(data=temp_data,aes(x=time,y=value))+
  facet_wrap(~variable,scales="free")+
  xlim(c(1950,2010))

jpeg(filename="demog_jpg",width=10,height=7,units="in",res=500)
print(plot_temp) 
dev.off() 

