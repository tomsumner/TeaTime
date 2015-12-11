### Loads all external data sources and where appropriate creates forcing functions

# Used to pick the country from cn
c_list <- c("South_Africa","Vietnam","Bangladesh")
cn <- c_list[cn]

# define number of TB, HIV and ART states
n_HIV <- 7
n_ART <- 3
n_dis <- 15

## UN population data ###########################################################################################################

# single year age groups for 1970 - used as initial population
UN_pop_start <- as.data.frame(read.table(paste("Demog/Initial_Pop_single_age.txt",sep=""),header=FALSE)) # Load UN Population data
# Load number of births
# add total, births and deaths to pop size data
UN_pop_start_t <- cbind(UN_pop_start[,1:83],rowSums(UN_pop_start[,2:82]))
UN_pop_start_t <- UN_pop_start_t[UN_pop_start_t[,1]==cn,]
colnames(UN_pop_start_t) <- c("Country","Year",colnames(UN_pop_start[3:83]),"Total")

## UN indicators - Crude birth rate (CBR), births, deaths #######################################################################
UN_ind <- as.data.frame(read.table("Demog/UN_indicators_all.txt",header=TRUE))
UN_ind <- UN_ind[UN_ind$Country==cn,]
# Convert CBR into forcing function to be used in C code
birth_rate <- cbind(seq(1970,2050),approx(UN_ind$Year,UN_ind$CBR,seq(1970,2050),rule=2)$y)
# births and deaths will be used in plots later

## Mortality ####################################################################################################################

# Age specific mortality taken from demproj outputs
mort_age <- as.data.frame(read.table(paste("Demog/",cn,"/",cn,"_mort_age.txt",sep=""),header=FALSE)) 

mort_age[,82] <-mort_age[,82]/2 

# Convert into forcing functions to be used in C code
for (i in 1:81){  
  assign(paste("s",i,sep=""), cbind(seq(1970,2050),approx(mort_age[,1],mort_age[,i+1],seq(1970,2050),rule=2)$y))
}

## MIGRATION ####################################################################################################################

# Migration data taken from DemProj - duplicate for single year age bins
mig <- as.data.frame(read.table(paste("Demog/",cn,"/",cn,"_migrants_age.txt",sep=""),header=FALSE))
temp <- c(rep(c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),each=5),18)
for (i in 0:80){
  assign(paste("mig",i,sep=""),cbind(seq(1970,2050),approx(mig[,1],mig[,temp[i+1]]/1000,seq(1970,2050),rule=2)$y))
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
temp <- c(rep(c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),each=5),18)
for (i in 0:80){
  assign(paste("h",i,sep=""),cbind(seq(1970,2050),approx(HIV_Inc_age[,1],HIV_Inc_age[,temp[i+1]]/1000,seq(1970,2050),rule=2)$y))
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
temp <- c(rep(c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),each=5),18)
for (i in 0:80){
  assign(paste("A",i,sep=""),cbind(seq(1970,2050),approx(seq(1970,2050),ART_percent[,temp[i+1]],seq(1970,2050),rule=2)$y))
}

# We still need the ART eligibility threshold (in terms of model CD4 categories) - adapt the code to calculate this from numeric CD4 threshold
ART_data <- as.data.frame(read.table(paste("HIV/",cn,"/",cn,"_ART_data.txt",sep=""),header=TRUE,fill=TRUE))
# ART eligibility threshold
Athresh <- cbind(ART_data[,"Year"],ART_data[,"CD4_cat"])


