## Code to test approach for age structured model 

## Based on http://kinglab.eeb.lsa.umich.edu/EEID/eeid/2011_eco/waifw.pdf

library(deSolve)

ages <- seq(1,100,by=1) # upper end of age classes - 100 1 year age bins
num_ages <- length(ages) # number of age classes
da <- diff(c(0,ages)) # widths of age classes

## Initial conditions
                                             
yinit <- c(S=c(rep(100,50),90,rep(100,49)),
           I=c(rep(0,50),10,rep(0,49)),
           R=c(rep(0,100)))

sindex <- 1:num_ages
iindex <- (num_ages+1):(2*num_ages)
rindex <- (2*num_ages+1):(3*num_ages)

## To capture the aging process, define a matrix to hold the rates of movement between age classes.
aging <- diag(-1/da)
aging[row(aging)-col(aging)==1] <- 1/head(da,-1)

# Parameters
beta <- 0.005
gamma <- 10
births <- 100

## Now we can put the pieces together to write a simulator for the age-structured SIR dynamics
ja.multistage.model <- function (t, x, ...) {
  s <- x[sindex] # susceptibles
  i <- x[iindex] # infecteds
  r <- x[rindex] # recovereds
  lambda <- beta*i # force of infection
  dsdt <- -lambda*s+aging%*%s
  didt <- lambda*s+aging%*%i-gamma*i
  drdt <- aging%*%r+gamma*i
  dsdt[1] <- dsdt[1]+births
  list(c(dsdt,didt,drdt))
}

# We can plug this into ode just as we did the simpler models to simulate an epidemic

sol <- ode(y=yinit,times=seq(0,100,by=0.1),func=ja.multistage.model)


time <- sol[,1]
infects <- sol[,1+iindex]
plot(time,apply(infects,1,sum),type='l')

