## Function to set up logistic curves for forcing fucntions ############################################
## "base" and "target" are the initial and final values
## "shape" determines the symmetry of the curve
## with "shape=1" the curve is symmetric around the midpoint of "base_year" and "target_year"
## "growth" determines the steepness of the curve
logcurve <- function(base,target,base_year,target_year,growth,shape){
  x <- seq(1970,2050)
  year <- mean(c(base_year,target_year))
  y <- (base + ((target-base)/((1+exp(-growth*(x-year)))^(1/shape))))/100
  z <- cbind(x,y)
  return(z)
}