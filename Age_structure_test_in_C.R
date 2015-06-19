## Code to test approach for age structured model as C code called from R
## Testing use of age specific fertility rates to see if it matches TIME and UN pop data better

## Based on http://kinglab.eeb.lsa.umich.edu/EEID/eeid/2011_eco/waifw.pdf

## Load packages
library(deSolve)
library(reshape2)
library(ggplot2)
library(gdata)

## Set working directory for the model files
setwd("C:/Users/TOM SUMMER/Documents/TIME_research/TeaTime")

# Compile and load C model
dyn.unload("Demog_model.dll")
system("R CMD SHLIB Demog_model.c") # Compile
dyn.load("Demog_model.dll")

cn_list <- c("SouthAfrica","Vietnam","Bangladesh","India")

plotList <- vector('list', length(cn_list))

for (cn in 1:length(cn_list)){

  # Enter the country name
  cntry <- cn_list[cn]

  # Load the demographic parameters

  # TIME population output
  setwd("C:/Users/TOM SUMMER/Filr/My Files/sync/TIME/TIME Research/Demographics_files")
  fname = paste(cntry,".xlsx",sep="")
  temp <- read.xls(fname, sheet="TIME_age_population", verbose=FALSE, "1970",header=FALSE,stringsAsFactors=FALSE)
  pop <- mat.or.vec(81,19)
  pop[,1] <- seq(1970,2050)
  pop_m <- pop
  pop_f <- pop

  for (i in 1:81){
  
    j <- (i-1)*19+2
  
    pop[i,2:19]=as.numeric(temp[j:(j+17),2])
    pop_m[i,2:19]=as.numeric(temp[j:(j+17),3])
    pop_f[i,2:19]=as.numeric(temp[j:(j+17),4])
  
  }
  pop <- as.data.frame(pop)
  pop_m <- as.data.frame(pop_m)
  pop_f <- as.data.frame(pop_f)

  # UN population estimates
  temp <- read.xls(fname, sheet="UN_pop", verbose=FALSE, "1950",header=FALSE,stringsAsFactors=FALSE)
  tempf <- read.xls(fname, sheet="UN_pop_f", verbose=FALSE, "1950",header=FALSE,stringsAsFactors=FALSE)
  tempm <- read.xls(fname, sheet="UN_pop_m", verbose=FALSE, "1950",header=FALSE,stringsAsFactors=FALSE)

  UN_pop_age <- mat.or.vec(101,19)
  UN_pop_age_f <- UN_pop_age
  UN_pop_age_m <- UN_pop_age

  for (i in 1:101){
    UN_pop_age[i,1:19] <- as.numeric(temp[i,1:19])
    UN_pop_age_f[i,1:19] <- as.numeric(tempf[i,1:19])
    UN_pop_age_m[i,1:19] <- as.numeric(tempm[i,1:19])
  }
  colnames(UN_pop_age) <- c("Year",
                          "0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49",
                          "50-54","55-59","60-64","65-69","70-74","75-79","80+","Total")

  colnames(UN_pop_age_f) <- colnames(UN_pop_age)
  colnames(UN_pop_age_m) <- colnames(UN_pop_age)

  # TFR - taken from DemProj
  temp <- read.xls(fname, sheet="TFR", verbose=FALSE, "TFR",header=FALSE,stringsAsFactors=FALSE)
  TFR <- cbind(seq(1970,2050),t(temp[2:82]))

  # ASFR (% of births in age group) - taken from Demproj
  temp <- read.xls(fname, sheet="ASFR", verbose=FALSE, "Age",header=FALSE,stringsAsFactors=FALSE) 
  b20 <- cbind(as.numeric(t(temp[1,2:82])),as.numeric(t(temp[2,2:82])))
  b25 <- cbind(as.numeric(t(temp[1,2:82])),as.numeric(t(temp[3,2:82])))
  b30 <- cbind(as.numeric(t(temp[1,2:82])),as.numeric(t(temp[4,2:82])))
  b35 <- cbind(as.numeric(t(temp[1,2:82])),as.numeric(t(temp[5,2:82])))
  b40 <- cbind(as.numeric(t(temp[1,2:82])),as.numeric(t(temp[6,2:82])))
  b45 <- cbind(as.numeric(t(temp[1,2:82])),as.numeric(t(temp[7,2:82])))
  b50 <- cbind(as.numeric(t(temp[1,2:82])),as.numeric(t(temp[8,2:82])))

  # Male to female birth ratio
  temp <- read.xls(fname, sheet="M_F_births", verbose=FALSE, "Birth",header=FALSE,stringsAsFactors=FALSE)
  m_to_f <- cbind(seq(1970,2050),as.numeric(t(temp[2:82]))/100)

  ## Survival - by age and sex based on UN life tables
  temp <- read.xls(fname, sheet="Life_table_m", verbose=FALSE,header=FALSE,stringsAsFactors=FALSE) 
  s5m <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,4])/as.numeric(temp[,2]))
  s10m <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,5])/as.numeric(temp[,4]))
  s15m <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,6])/as.numeric(temp[,5]))
  s20m <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,7])/as.numeric(temp[,6]))
  s25m <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,8])/as.numeric(temp[,7]))
  s30m <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,9])/as.numeric(temp[,8]))
  s35m <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,10])/as.numeric(temp[,9]))
  s40m <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,11])/as.numeric(temp[,10]))
  s45m <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,12])/as.numeric(temp[,11]))
  s50m <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,13])/as.numeric(temp[,12]))
  s55m <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,14])/as.numeric(temp[,13]))
  s60m <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,15])/as.numeric(temp[,14]))
  s65m <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,16])/as.numeric(temp[,15]))
  s70m <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,17])/as.numeric(temp[,16]))
  s75m <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,18])/as.numeric(temp[,17]))
  s80m <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,19])/as.numeric(temp[,18]))
  s100m <- cbind(seq(1952.5,2097.5,5),0)
  temp <- read.xls(fname, sheet="Life_table_f", verbose=FALSE,header=FALSE,stringsAsFactors=FALSE) 
  s5f <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,4])/as.numeric(temp[,2]))
  s10f <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,5])/as.numeric(temp[,4]))
  s15f <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,6])/as.numeric(temp[,5]))
  s20f <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,7])/as.numeric(temp[,6]))
  s25f <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,8])/as.numeric(temp[,7]))
  s30f <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,9])/as.numeric(temp[,8]))
  s35f <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,10])/as.numeric(temp[,9]))
  s40f <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,11])/as.numeric(temp[,10]))
  s45f <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,12])/as.numeric(temp[,11]))
  s50f <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,13])/as.numeric(temp[,12]))
  s55f <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,14])/as.numeric(temp[,13]))
  s60f <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,15])/as.numeric(temp[,14]))
  s65f <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,16])/as.numeric(temp[,15]))
  s70f <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,17])/as.numeric(temp[,16]))
  s75f <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,18])/as.numeric(temp[,17]))
  s80f <- cbind(seq(1952.5,2097.5,5),as.numeric(temp[,19])/as.numeric(temp[,18]))
  s100f <- cbind(seq(1952.5,2097.5,5),0)

#########################################################################
# Now run the model

  # Number of age groups in model
  num_ages = 17

  ## Initial conditions - extract 1970 population from UN data
  temp <- c()
  for (i in 1:num_ages){temp[i]<-UN_pop_age_f[21,i+1]}
  for (i in (num_ages+1):(2*num_ages)){temp[i]<-UN_pop_age_m[21,(i+1)-17]}
  xstart <- c(S=c(temp))

  # Parameters - width of age bins (5 years for all but the final group which is 21 years)
  parms <- c(d = 0.2, d2 = 0.1)
  ## vector of timesteps
  times <- seq(1970,2050 , by=1)
  # combine all forcing functions
  force <- list(b20,b25,b30,b35,b40,b45,b50,
                s5f,s10f,s15f,s20f,s25f,s30f,s35f,s40f,s45f,s50f,s55f,s60f,s65f,s70f,s75f,s80f,s100f,
                s5m,s10m,s15m,s20m,s25m,s30m,s35m,s40m,s45m,s50m,s55m,s60m,s65m,s70m,s75m,s80m,s100m,
                TFR,m_to_f)
  # Run the model
  out <- ode(y=xstart, times, func = "derivsc",
                              parms = parms, dllname = "Demog_model",initforc = "forcc",
                              forcings=force, initfunc = "parmsc", nout = 4,
                              outnames = c("Total_f","Total_m","Total","Births"), method = rkMethod("rk34f"))

  ## Plot populations in each age group over time

  # add total to UN data and convert to long format
  temp_data_f <- melt(as.data.frame(UN_pop_age_f),id="Year")
  temp_data_m <- melt(as.data.frame(UN_pop_age_m),id="Year")
  temp_data <- melt(as.data.frame(UN_pop_age),id="Year")

  # then turn the model output into long format
  temp_model_f <- as.data.frame(out[,c(seq(1,18),36)])
  colnames(temp_model_f) <- colnames(UN_pop_age_f)
  temp_model_f <- melt(temp_model_f,id="Year")

  temp_model_m <- as.data.frame(out[,c(1,seq(19,35),37)])
  colnames(temp_model_m) <- colnames(UN_pop_age_m)
  temp_model_m <- melt(temp_model_m,id="Year")

  temp_model <- as.data.frame(cbind(out[,1],out[,2:18]+out[,19:35],out[,38]))
  colnames(temp_model) <- colnames(UN_pop_age)
  temp_model <- melt(temp_model,id="Year")

  # and do the same for TIME estimates
  colnames(pop) <- colnames(UN_pop_age)
  TIME_pop <- melt(pop,id="Year")
  colnames(pop_m) <- colnames(UN_pop_age)
  TIME_pop_m <- melt(pop_m,id="Year")
  colnames(pop_f) <- colnames(UN_pop_age)
  TIME_pop_f <- melt(pop_f,id="Year")

  # and plot
  plotList[[cn]] <- ggplot(temp_model_f,aes(x=Year,y=value))+
    geom_line(colour="red")+
    geom_line(data=temp_model_m,aes(x=Year,y=value),colour="blue")+
    geom_line(data=temp_model,aes(x=Year,y=value),colour="black")+
    geom_point(data=temp_data_f,aes(x=Year,y=value),colour="red")+
    geom_point(data=temp_data_m,aes(x=Year,y=value),colour="blue")+
    geom_point(data=temp_data,aes(x=Year,y=value),colour="black")+
    geom_line(data=TIME_pop_f,aes(x=Year,y=value/1000),colour="red",linetype="dashed")+
    geom_line(data=TIME_pop_m,aes(x=Year,y=value/1000),colour="blue",linetype="dashed")+
    geom_line(data=TIME_pop,aes(x=Year,y=value/1000),colour="black",linetype="dashed")+
    facet_wrap(~variable,scales="free")+
    xlim(c(1970,2050))

}
