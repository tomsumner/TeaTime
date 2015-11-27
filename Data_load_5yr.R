### Loads all external data sources and where appropriate creates forcing functions

## This version is for the 5 yr age bin model

## UN population data ###########################################################################################################

# single year age groups for 1970 - used as initial population
UN_pop_start <- as.data.frame(read.table(paste("Demog/Initial_Pop_single_age.txt",sep=""),header=FALSE)) # Load UN Population data
# Load number of births
# add total to pop size data
UN_pop_start_t1 <- cbind(UN_pop_start[,1:83],rowSums(UN_pop_start[,3:83]))
# Then pull out country of interest
UN_pop_start_t1 <- UN_pop_start_t1[UN_pop_start_t1[,1]==cn,]
# Then sum into 5 year age bins
UN_pop_start_t<-mat.or.vec(1,19)
for (i in 1:16){
  j <- (3*i)+2*(i-1)
  UN_pop_start_t[i+2]=sum(UN_pop_start_t1[j:(j+4)])
}
UN_pop_start_t[1]<-as.character(droplevels(UN_pop_start_t1[,1]))
UN_pop_start_t[2]<-UN_pop_start_t1[,2]
UN_pop_start_t[19]<-UN_pop_start_t1[,83]
colnames(UN_pop_start_t) <- c("Country","Year",colnames(UN_pop_start[3:19]))

## UN indicators - Crude birth rate (CBR), births, deaths #######################################################################
UN_ind <- as.data.frame(read.table("Demog/UN_indicators_all_test.txt",header=TRUE))
UN_ind <- UN_ind[UN_ind$Country==cn,]
# Convert CBR into forcing function to be used in C code
birth_rate <- cbind(seq(1970,2050),approx(UN_ind$Year,UN_ind$CBR,seq(1970,2050),rule=2)$y)
# births and deaths will be used in plots later

## Mortality ####################################################################################################################

###### THIS NEEDS TO BE FIXED TO BE 5 YEAR MORTALITY RATES ################

# Age specific mortality taken from demproj outputs
mort_age <- as.data.frame(read.table(paste("Demog/",cn,"/",cn,"_mort_5yr.txt",sep=""),header=FALSE)) 
mort_age[,2:17] <- (1-mort_age[,2:17])/5
mort_age[,18] <- mort_age[,17]/2
# Convert into forcing functions to be used in C code
temp <- c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,100)
for (i in 1:17){  
  assign(paste("s",temp[i],sep=""), cbind(seq(1970,2050),approx(mort_age[,1],mort_age[,i+1],seq(1970,2050),rule=2)$y))
  #assign(paste("s",temp[i],sep=""), cbind(seq(1970,2050),approx(mort_age[,1],rep(0,30),seq(1970,2050),rule=2)$y))
}

## MIGRATION ####################################################################################################################

# Migration data taken from DemProj - duplicate for single year age bins
mig <- as.data.frame(read.table(paste("Demog/",cn,"/",cn,"_migrants_age.txt",sep=""),header=FALSE))
temp <- c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80)
for (i in 1:17){
  assign(paste("mig",temp[i],sep=""),cbind(seq(1970,2050),approx(mig[,1],mig[,i+1]/1000,seq(1970,2050),rule=2)$y))
}

## MORTALITY ADJUSTMENT #########################################################################################################

# Pop adjust - to turn off population adjustment for TB/HIV deaths from 2015 onwards
pop_ad <- cbind(seq(1970,2050),c(rep(1,45),rep(0,35),0))


## HIV ##########################################################################################################################

# HIV Incidence by age and year - taken from AIM (rate per 1000)
temp <- as.data.frame(read.table(paste("HIV/",cn,"/",cn,"_HIV_Inc_age.txt",sep=""),header=TRUE,fill=TRUE))
# Re-arrage to get in year vs age format
HIV_Inc_age <- mat.or.vec(81,18)
HIV_Inc_age[,1] <- seq(1970,2050)
for (i in 1:81){
  j <- (i-1)*19+4
  HIV_Inc_age[i,2:18]=as.numeric(levels(temp[j:(j+16),2]))[temp[j:(j+16),2]]
}
# Now duplicate 5 year age bin values to give values for single year bins and convert into forcing functions
temp <- c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80)
for (i in 1:17){
  assign(paste("h",temp[i],sep=""),cbind(seq(1970,2050),approx(HIV_Inc_age[,1],HIV_Inc_age[,i+1]/1000,seq(1970,2050),rule=2)$y))
}

## ART ##########################################################################################################################

# Take "need for" and "number on" by age and year from AIM, convert to % of eligible on ART 

# ART need (by age and year taken from AIM)
temp <- as.data.frame(read.table(paste("HIV/",cn,"/",cn,"_ART_need_age.txt",sep=""),header=TRUE,fill=TRUE))
# Re-arrage to get in year vs age format
ART_need_age <- mat.or.vec(81,18)
ART_need_age[,1] <- seq(1970,2050)
for (i in 1:81){
  j <- (i-1)*19+4
  ART_need_age[i,2:18]=as.numeric(levels(temp[j:(j+16),2]))[temp[j:(j+16),2]]
}

# ART numbers (by age and year taken from AIM)
temp <- as.data.frame(read.table(paste("HIV/",cn,"/",cn,"_ART_numbers_age.txt",sep=""),header=TRUE,fill=TRUE))
# Re-arrage to get in year vs age format
ART_number_age <- mat.or.vec(81,18)
ART_number_age[,1] <- seq(1970,2050)
for (i in 1:81){
  j <- (i-1)*19+4
  ART_number_age[i,2:18]=as.numeric(levels(temp[j:(j+16),2]))[temp[j:(j+16),2]]
}

# Calculate % coverage
ART_percent <- cbind(seq(1970,2050),ART_number_age[,2:18]/ART_need_age[,2:18])
# Correct for "nan" if "need" is zero
ART_percent[is.nan(ART_percent)] <- 0

# Now duplicate 5 year age bin values to give values for single year bins and convert into forcing functions
temp <- c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80)
for (i in 1:17){
  assign(paste("A",temp[i],sep=""),cbind(seq(1970,2050),approx(seq(1970,2050),ART_percent[,i+1],seq(1970,2050),rule=2)$y))
}

# We still need the ART eligibility threshold (in terms of model CD4 categories) - adapt the code to calculate this from numeric CD4 threshold
ART_data <- as.data.frame(read.table(paste("HIV/",cn,"/",cn,"_ART_data.txt",sep=""),header=TRUE,fill=TRUE))
# ART eligibility threshold
Athresh <- cbind(ART_data[,"Year"],ART_data[,"CD4_cat"])


