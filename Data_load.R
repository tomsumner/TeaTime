### Loads all external data sources and where appropriate creates forcing functions

## UN population data ###########################################################################################################

# single year age groups for 1970 - used as initial population
UN_pop_start <- as.data.frame(read.table(paste("Demog/Initial_Pop_single_age.txt",sep=""),header=FALSE)) # Load UN Population data
# Load number of births
# add total, births and deaths to pop size data
UN_pop_start_t <- cbind(UN_pop_start[,1:83],rowSums(UN_pop_start[,2:82]))
UN_pop_start_t <- UN_pop_start_t[UN_pop_start_t[,1]==cn,]
colnames(UN_pop_start_t) <- c("Country","Year",colnames(UN_pop_start[3:83]),"Total")

# UN indicators - Crude birth rate (CBR), births, deaths ########################################################################
UN_ind <- as.data.frame(read.table("Demog/UN_indicators_all_test.txt",header=TRUE))
UN_ind <- UN_ind[UN_ind$Country==cn,]
# Convert CBR into forcing function to be used in C code
birth_rate <- cbind(seq(1970,2050,by=ss),approx(UN_ind$Year,UN_ind$CBR,seq(1970,2050,by=ss),rule=2)$y)
# births and deaths will be used in plots later

# Mortality #####################################################################################################################

# Age specific mortality taken from demproj outputs
mort_age <- as.data.frame(read.table(paste("Demog/",cn,"/",cn,"_mort_age.txt",sep=""),header=FALSE)) 

mort_age[,82] <-mort_age[,82]/2 

# Convert into forcing functions to be used in C code
for (i in 1:81){  
  assign(paste("s",i,sep=""), cbind(seq(1970,2050,by=ss),approx(mort_age[,1],mort_age[,i+1],seq(1970,2050,by=ss),rule=2)$y))
}


## HIV and ART ##################################################################################################################

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
  assign(paste("h",i,sep=""),cbind(seq(1970,2050,by=ss),approx(HIV_Inc_age[,1],HIV_Inc_age[,temp[i+1]]/1000,seq(1970,2050,by=ss),rule=2)$y))
}

# ART coverage - based on AIM, use CD4 eligibility threshold and % of those in need on ART
ART_data <- as.data.frame(read.table(paste("HIV/",cn,"/",cn,"_ART_data.txt",sep=""),header=TRUE)) # Load data

#ART_data[,"Percent"] <- 0

# Create forcing function of threshold category
Athresh <- cbind(ART_data[,"Year"],ART_data[,"CD4_cat"])
# Create forcing functions which account for threshold and coverage
A50 <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>50,1,0)/100)
A99 <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>99,1,0)/100)
A199 <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>199,1,0)/100)
A249 <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>249,1,0)/100)
A349 <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>349,1,0)/100)
A500 <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>499,1,0)/100)
Ahigh <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>500,1,0)/100)

# Migration data taken from DemProj - duplicate for single year age bins
mig <- as.data.frame(read.table(paste("Demog/",cn,"/",cn,"_migrants_age.txt",sep=""),header=FALSE))
temp <- c(rep(c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),each=5),18)
for (i in 0:80){
  assign(paste("mig",i,sep=""),cbind(seq(1970,2050,by=ss),approx(mig[,1],mig[,temp[i+1]]/1000,seq(1970,2050,by=ss),rule=2)$y))
}

# Pop adjust - to turn off population adjustment for TB/HIV deaths from 2015 onwards
pop_ad <- cbind(seq(1970,2050,by=ss),c(rep(1,45/ss),rep(0,35/ss),0))

