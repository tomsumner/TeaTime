######## Some plots for testing things against spectrum and other data ###########################################

# TB outputs #####################################################################################################

# Load the TIME results (Incidence, Prevalence, Mortality)
TIME_out <- as.data.frame(read.table(paste("TB/",cn,"_TIME_TB.txt",sep=""),header=TRUE,fill=TRUE))
TIME_out <- cbind(TIME_out,"TIME")
colnames(TIME_out) <- c("Year","Prevalence","Incidence","Mortality","Model")

# Arrange model output
R_out <- as.data.frame(cbind(out[,"time"],
                             100000*(out[,"Total_DS"]+out[,"Total_MDR"])/out[,"Total"],  # Prev
                             100000*(out[,"Cases_neg"]+out[,"Cases_pos"]+out[,"Cases_ART"])/out[,"Total"], # Inc 
                             100000*out[,"TB_deaths"]/out[,"Total"])) # Mort (all)
R_out <- cbind(R_out,"R")
colnames(R_out) <- c("Year","Prevalence","Incidence","Mortality","Model")
# combine and melt
Models_out <- rbind(TIME_out,R_out)
Models_out <- melt(Models_out,id=c("Year","Model"))
Models_out <- Models_out[Models_out$Year<=2035,]

# Calculate % diff between models and add to df
temp=100*(Models_out[Models_out$Model=="TIME",4]-Models_out[Models_out$Model=="R",4])/Models_out[Models_out$Model=="TIME",4]
Models_out<-cbind(temp,Models_out)

plot_models <- ggplot(Models_out[Models_out$Model!="Data",],aes(x=Year,y=value,colour=Model))+
  geom_line()+
  facet_wrap(~variable,scales="free_y")+
  xlim(c(1970,2035))+
  ylab("Rate /100,000")+
  scale_y_continuous(expand = c(0, 0))+
  expand_limits(y = 0)
plot_models

plot_diff <- ggplot(Models_out[Models_out$Model=="TIME",],aes(x=Year,y=temp))+
  geom_line()+
  facet_wrap(~variable,scales="free_y")+
  xlim(c(1970,2035))+
  ylab("% difference bewtween R and TIME")
plot_diff

# POPULATION #####################################################################################################

## Load and manipulate UN numbers 
## Load UN population data
UN_pop_age <- as.data.frame(read.table(paste("Demog/",cn,"_pop_age.txt",sep=""),header=TRUE)) # Load UN Population data
# Load number of births
births <- as.data.frame(read.table(paste("Demog/",cn,"_births_number.txt",sep=""),header=TRUE))
# Load number of deaths
deaths <- as.data.frame(read.table(paste("Demog/",cn,"_deaths.txt",sep=""),header=TRUE))
# add total, births and deaths to data
UN_pop_age_t <- cbind(births,UN_pop_age[,2:18],rowSums(UN_pop_age[,2:18]),deaths[,2])
colnames(UN_pop_age_t) <- c("Year","births",colnames(UN_pop_age[2:18]),"Total","Deaths")
temp_data <- melt(UN_pop_age_t,id="Year")

temp_data <- cbind(temp_data,"UNPP")
colnames(temp_data) <- c(colnames(temp_data)[1:3],"Model")

# sum up model outputs over age groups and turn into long format
tot <- mapply(function(x,y) sum(out[x,seq(y+1,30538,81)]),rep(seq(1,81),each=81),seq(1,81))
dim(tot) <- c(81,81)
tot <- t(tot)

model_temp <- mat.or.vec(81,21)
model_temp[,1] <- out[,"time"]
model_temp[,2] <- out[,"births"]
model_temp[,20] <- out[,"Total"]
model_temp[,21] <- out[,"deaths"]
model_temp[,19] <- tot[,81]

# Then aggregate into 5 year bins used in TIME
for (i in 1:16){
  t1 <- (i)+(i-1)*4
  t2 <- t1+4
  model_temp[,i+2] <- rowSums(tot[,t1:t2])
}
temp_model <- as.data.frame(model_temp)
colnames(temp_model) <- colnames(UN_pop_age_t)
temp_model_m <- melt(temp_model,id="Year")
temp_model_m <- cbind(temp_model_m,"R")
colnames(temp_model_m) <- c(colnames(temp_model_m)[1:3],"Model")

# TIME values
temp <- as.data.frame(read.table(paste("Demog/",cn,"_TIME_pop_age.txt",sep=""),header=TRUE,fill=TRUE)) 
# Need to re-arrage to get in year vs age format
TIME_pop <- mat.or.vec(81,19)
TIME_pop[,1] <- seq(1970,2050)
for (i in 1:81){
  j <- (i-1)*19+2
  TIME_pop[i,2:19]=temp[j:(j+17),2]
}
TIME_births <- as.data.frame(read.table(paste("Demog/",cn,"_TIME_births.txt",sep=""),header=TRUE))
TIME_deaths <- as.data.frame(read.table(paste("Demog/",cn,"_TIME_deaths.txt",sep=""),header=TRUE))
TIME_pop <- cbind(TIME_births/1000,TIME_pop[,2:19]/1000,TIME_deaths[,2]/1000)
TIME_pop[,1] <- TIME_pop[,1]*1000
colnames(TIME_pop) <- colnames(UN_pop_age_t)
TIME_pop_m <- melt(TIME_pop,id="Year")
TIME_pop_m <- cbind(TIME_pop_m,"TIME")
colnames(TIME_pop_m) <- c(colnames(TIME_pop_m)[1:3],"Model")

# demproj values
temp <- as.data.frame(read.table(paste("Demog/",cn,"_dem_pop_age.txt",sep=""),header=TRUE,fill=TRUE)) 
# Need to re-arrage to get in year vs age format
dem_pop <- mat.or.vec(81,19)
dem_pop[,1] <- seq(1970,2050)
for (i in 1:81){
  j <- (i-1)*19+2
  dem_pop[i,2:19]=temp[j:(j+17),2]
}
TIME_births <- as.data.frame(read.table(paste("Demog/",cn,"_TIME_births.txt",sep=""),header=TRUE))
TIME_deaths <- as.data.frame(read.table(paste("Demog/",cn,"_TIME_deaths.txt",sep=""),header=TRUE))
dem_pop <- cbind(TIME_births/1000,dem_pop[,2:19]/1000,TIME_deaths[,2]/1000)
dem_pop[,1] <- dem_pop[,1]*1000
colnames(dem_pop) <- colnames(UN_pop_age_t)
dem_pop_m <- melt(dem_pop,id="Year")

dem_pop_m <- cbind(dem_pop_m,"DemProj")
colnames(dem_pop_m) <- c(colnames(dem_pop_m)[1:3],"Model")

dat_to_plot <- rbind(temp_model_m,dem_pop_m,TIME_pop_m,temp_data)

plot_pop <- ggplot(dat_to_plot,aes(x=Year,y=value,colour=Model))+
  geom_line()+
  facet_wrap(~variable,scales="free")+
  ggtitle("Population (by 5 year age bins), total births and total deaths")+
  xlim(c(1970,2050))+
  ylab("Thousands")
plot_pop


# Compare population age structure (i.e % of total pop)
temp_data_s <- as.data.frame(cbind(UN_pop_age_t[,1],100*UN_pop_age_t[,3:20]/UN_pop_age_t[,20]))
colnames(temp_data_s)<-c("Year","x4","X9","X14","X19","X24","X29","X34","X39","X44","X49","X54","X59","X64","X69","X74","X79","X100","Total")
temp_data_s <- melt(temp_data_s,id="Year")
temp_data_s <- cbind(temp_data_s,"UNPP")
colnames(temp_data_s) <- c(colnames(temp_data_s)[1:3],"Model")

temp_model_s <- as.data.frame(cbind(seq(1970,2050),100*temp_model[,3:20]/temp_model[,20]))
colnames(temp_model_s)<-c("Year","x4","X9","X14","X19","X24","X29","X34","X39","X44","X49","X54","X59","X64","X69","X74","X79","X100","Total")
temp_model_s <- melt(temp_model_s,id="Year")
temp_model_s <- cbind(temp_model_s,"R")
colnames(temp_model_s) <- c(colnames(temp_model_s)[1:3],"Model")

temp_TIME_s <- as.data.frame(cbind(seq(1970,2050),100*TIME_pop[,3:20]/TIME_pop[,20]))
colnames(temp_TIME_s)<-c("Year","x4","X9","X14","X19","X24","X29","X34","X39","X44","X49","X54","X59","X64","X69","X74","X79","X100","Total")
temp_TIME_s <- melt(temp_TIME_s,id="Year")
temp_TIME_s <- cbind(temp_TIME_s,"TIME")
colnames(temp_TIME_s) <- c(colnames(temp_TIME_s)[1:3],"Model")

dat_to_plot <- rbind(temp_data_s,temp_model_s,temp_TIME_s)

# and plot
plot_pop_s <- ggplot(dat_to_plot,aes(x=Year,y=value,color=Model))+
  geom_line()+
  facet_wrap(~variable,scales="free")+
  ggtitle("Age structure of population by 5 year age bins")+
  ylab("% of total population")+
  xlim(c(1970,2050))
plot_pop_s

# HIV ############################################################################################################

# Load TIME values
temp <- as.data.frame(read.table(paste("HIV/",cn,"_HIV_numbers_age.txt",sep=""),header=TRUE,fill=TRUE)) # Load HIV numbers (output in TIME)  
# Need to re-arrage to get in year vs age format
HIV_number_age <- mat.or.vec(81,18)
HIV_number_age[,1] <- seq(1970,2050)
for (i in 1:81){
  j <- (i-1)*19+2
  #HIV_number_age[i,2:18]=as.numeric(levels(temp[j:(j+16),2]))[temp[j:(j+16),2]]
  HIV_number_age[i,2:18]=temp[j:(j+16),2]
}
HIV_number_age <- as.data.frame(HIV_number_age)
colnames(HIV_number_age) <- colnames(UN_pop_age)
HIV_TIME <- melt(HIV_number_age,id="Year")
HIV_TIME <- cbind(HIV_TIME,"TIME")
colnames(HIV_TIME) <- c(colnames(HIV_TIME)[1:3],"Model")
HIV_TIME$value <- HIV_TIME$value/1000

# Sum up model outputs for HIV+ categories by age
tot_HIV <- mapply(function(x,y) sum(out[x,seq(y+1054,30538,81)]),rep(seq(1,81),each=81),seq(1,81))
dim(tot_HIV) <- c(81,81)
tot_HIV <- t(tot_HIV)

HIV_model <- mat.or.vec(81,18)
HIV_model[,1] <- seq(1970,2050)
HIV_model[,18] <- tot_HIV[,81]
for (i in 1:16){
  t1 <- (i)+(i-1)*4
  t2 <- t1+4
  HIV_model[,i+1] <- rowSums(tot_HIV[,t1:t2])
}

HIV_model <- as.data.frame(HIV_model )
colnames(HIV_model) <- colnames(UN_pop_age)
HIV_model_m <- melt(HIV_model,id="Year")
HIV_model_m <- cbind(HIV_model_m,"R")
colnames(HIV_model_m) <- c(colnames(HIV_model_m)[1:3],"Model")

dat_to_plot <- rbind(HIV_TIME,HIV_model_m)

# and plot
plot_HIV <- ggplot(dat_to_plot,aes(x=Year,y=value,color=Model))+
  geom_line()+
  facet_wrap(~variable,scales="free")+
  ggtitle("HIV_numbers")+
  xlim(c(1970,2050))
plot_HIV

### Compare HIV prevalence by age group

# TIME
H_p_TIME <- as.data.frame(cbind(HIV_model[,1],100*HIV_number_age[,2:18]/(1000*TIME_pop[,3:19])))
colnames(H_p_TIME) <- colnames(UN_pop_age)
H_p_TIME_m <- melt(H_p_TIME,id="Year")
H_p_TIME_m <- cbind(H_p_TIME_m,"TIME")
colnames(H_p_TIME_m) <- c(colnames(H_p_TIME_m)[1:3],"Model")
# Model
H_p_model <- as.data.frame(cbind(HIV_model[,1],100*HIV_model[,2:18]/model_temp[,3:19]))
colnames(H_p_model) <- colnames(UN_pop_age)
H_p_model_m <- melt(H_p_model,id="Year")
H_p_model_m <- cbind(H_p_model_m,"R")
colnames(H_p_model_m) <- c(colnames(H_p_model_m)[1:3],"Model")

dat_to_plot <- rbind(H_p_TIME_m,H_p_model_m)

plot_HIV_p <- ggplot(dat_to_plot,aes(x=Year,y=value,color=Model))+
  geom_line()+
  facet_wrap(~variable,scales="free")+
  ggtitle("HIV_prevalence (age)")+
  xlim(c(1970,2050))
plot_HIV_p

# CD4 distributions ##############################################################################################

# Plot of CD4 distribution against Spectrum
# this is normalised by total population to account for differences in population size 

# Load TIME outputs
TIME_no_ART <- as.data.frame(read.table(paste("HIV/",cn,"_TIME_no_ART.txt",sep=""),header=TRUE,fill=TRUE))
TIME_no_ART <- cbind(TIME_no_ART[,1],100*TIME_no_ART[,2:8]/(1000*TIME_pop[,"Total"]),"No_ART","TIME")
colnames(TIME_no_ART) <- c("Year",colnames(TIME_no_ART[2:8]),"Type","Model")

TIME_ART <- as.data.frame(read.table(paste("HIV/",cn,"_TIME_on_ART.txt",sep=""),header=TRUE,fill=TRUE))
TIME_ART <- cbind(TIME_ART[,1],100*TIME_ART[,2:8]/(1000*TIME_pop[,"Total"]),"ART","TIME")
colnames(TIME_ART) <- c("Year",colnames(TIME_ART[2:8]),"Type","Model")

temp1 <- as.data.frame(cbind(out[,"time"],100*out[,30554:30560]/out[,"Total"]))
temp1 <- cbind(temp1,"No_ART","R")
colnames(temp1) <- c("Year",colnames(temp1[,2:8]),"Type","Model")

temp2 <- as.data.frame(cbind(seq(1970,2050),100*out[,30561:30567]/out[,"Total"]))
temp2 <- cbind(temp2,"ART","R")
colnames(temp2) <- c("Year",colnames(temp1[,2:8]),"Type","Model")

temp_CD4 <- rbind(temp1,temp2,TIME_ART,TIME_no_ART)
temp_CD4 <- melt(temp_CD4,id=c("Year","Type","Model"))

plot_CD4 <- ggplot(temp_CD4,aes(x=as.numeric(as.character(Year)),y=as.numeric(as.character(value)),colour=variable,linetype=Model))+
  geom_line() +
  facet_wrap(~Type)+
  xlim(c(2000,2020))
plot_CD4

# TB prevalence by age and HIV ###################################################################################

###### HIV_neg ######

# sum up TB cases over age groups and turn into long format
TB_neg <- mapply(function(x,y) sum(out[x,seq(y+406,1054,81)]),rep(seq(1,81),each=81),seq(1,81))
dim(TB_neg) <- c(81,81)
TB_neg <- t(TB_neg)

model_temp <- mat.or.vec(81,18)
model_temp[,1] <- out[,"time"]
model_temp[,18] <- TB_neg[,81]

# Then aggregate into 5 year bins used in TIME
for (i in 1:16){
  t1 <- (i)+(i-1)*4
  t2 <- t1+4
  model_temp[,i+1] <- rowSums(TB_neg[,t1:t2])
}

colnames(model_temp) <- c("Year",colnames(UN_pop_age_t)[3:19])
temp_model_m <- melt(as.data.frame(model_temp),id="Year")

# Load TIME values
temp <- as.data.frame(read.table(paste("TB/",cn,"_TB_prev_numbers_neg_age.txt",sep=""),header=TRUE,fill=TRUE))   
# Need to re-arrage to get in year vs age format
TB_number_age <- mat.or.vec(81,18)
TB_number_age[,1] <- seq(1970,2050)
for (i in 1:81){
  j <- (i-1)*19+2
  TB_number_age[i,2:18]=temp[j:(j+16),2]
}
TB_number_age <- as.data.frame(TB_number_age)
colnames(TB_number_age) <- colnames(UN_pop_age)
TB_TIME <- melt(TB_number_age,id="Year")

plot_TB_prev_neg <- ggplot(temp_model_m,aes(x=Year,y=value))+
  geom_line(colour="red")+
  geom_line(data=TB_TIME,aes(x=Year,y=value/1000),colour="black",linetype="dashed")+
  facet_wrap(~variable,scales="free")+
  ggtitle("HIV- TB prev (by age)")+
  xlim(c(1970,2050))
plot_TB_prev_neg  

###### HIV+ no ART ######
  
# sum up TB cases over age groups and turn into long format
TB_neg <- mapply(function(x,y) sum(out[x,seq(y+3889,8425,81)]),rep(seq(1,81),each=81),seq(1,81))
dim(TB_neg) <- c(81,81)
TB_neg <- t(TB_neg)

model_temp <- mat.or.vec(81,18)
model_temp[,1] <- out[,"time"]
model_temp[,18] <- TB_neg[,81]

# Then aggregate into 5 year bins used in TIME
for (i in 1:16){
  t1 <- (i)+(i-1)*4
  t2 <- t1+4
  model_temp[,i+1] <- rowSums(TB_neg[,t1:t2])
}

colnames(model_temp) <- c("Year",colnames(UN_pop_age_t)[3:19])
temp_model_m <- melt(as.data.frame(model_temp),id="Year")

# Load TIME values
temp <- as.data.frame(read.table(paste("TB/",cn,"_TB_prev_numbers_noART_age.txt",sep=""),header=TRUE,fill=TRUE))   
# Need to re-arrage to get in year vs age format
TB_number_age <- mat.or.vec(81,18)
TB_number_age[,1] <- seq(1970,2050)
for (i in 1:81){
  j <- (i-1)*19+2
  TB_number_age[i,2:18]=temp[j:(j+16),2]
}
TB_number_age <- as.data.frame(TB_number_age)
colnames(TB_number_age) <- colnames(UN_pop_age)
TB_TIME <- melt(TB_number_age,id="Year")

plot_TB_prev_noART <- ggplot(temp_model_m,aes(x=Year,y=value))+
  geom_line(colour="red")+
  geom_line(data=TB_TIME,aes(x=Year,y=value/1000),colour="black",linetype="dashed")+
  facet_wrap(~variable,scales="free")+
  ggtitle("HIV+ no ART TB prev (by age)")+
  xlim(c(1970,2050))
plot_TB_prev_noART
###### HIV+ ART ######

# sum up TB cases over age groups and turn into long format
TB_neg <- mapply(function(x,y) sum(out[x,seq(y+16930,30538,81)]),rep(seq(1,81),each=81),seq(1,81))
dim(TB_neg) <- c(81,81)
TB_neg <- t(TB_neg)

model_temp <- mat.or.vec(81,18)
model_temp[,1] <- out[,"time"]
model_temp[,18] <- TB_neg[,81]

# Then aggregate into 5 year bins used in TIME
for (i in 1:16){
  t1 <- (i)+(i-1)*4
  t2 <- t1+4
  model_temp[,i+1] <- rowSums(TB_neg[,t1:t2])
}

colnames(model_temp) <- c("Year",colnames(UN_pop_age_t)[3:19])
temp_model_m <- melt(as.data.frame(model_temp),id="Year")

# Load TIME values
temp <- as.data.frame(read.table(paste("TB/",cn,"_TB_prev_numbers_ART_age.txt",sep=""),header=TRUE,fill=TRUE))   
# Need to re-arrage to get in year vs age format
TB_number_age <- mat.or.vec(81,18)
TB_number_age[,1] <- seq(1970,2050)
for (i in 1:81){
  j <- (i-1)*19+2
  TB_number_age[i,2:18]=temp[j:(j+16),2]
}
TB_number_age <- as.data.frame(TB_number_age)
colnames(TB_number_age) <- colnames(UN_pop_age)
TB_TIME <- melt(TB_number_age,id="Year")

plot_TB_prev_ART <- ggplot(temp_model_m,aes(x=Year,y=value))+
  geom_line(colour="red")+
  geom_line(data=TB_TIME,aes(x=Year,y=value/1000),colour="black",linetype="dashed")+
  facet_wrap(~variable,scales="free")+
  ggtitle("HIV+, ART TB prev (by age)")+
  xlim(c(1970,2050))
plot_TB_prev_ART
## Deaths by age - compare R vs TIME vs demproj ##################################################################

# sum up deaths over age groups and turn into long format
deaths_age <- out[,30586:30666]

model_temp <- mat.or.vec(81,18)
model_temp[,1] <- out[,"time"]
model_temp[,18] <- deaths_age[,81]

# Then aggregate into 5 year bins used in TIME
for (i in 1:16){
  t1 <- (i)+(i-1)*4
  t2 <- t1+4
  model_temp[,i+1] <- rowSums(deaths_age[,t1:t2])
}

colnames(model_temp) <- c("Year",colnames(UN_pop_age_t)[3:19])
temp_model_m <- melt(as.data.frame(model_temp),id="Year")

# Load TIME values
temp <- as.data.frame(read.table(paste("Demog/",cn,"_TIME_deaths_age.txt",sep=""),header=TRUE,fill=TRUE))   
# Need to re-arrage to get in year vs age format
TIME_deaths_age <- mat.or.vec(81,18)
TIME_deaths_age[,1] <- seq(1970,2050)
for (i in 1:81){
  j <- (i-1)*19+2
  TIME_deaths_age[i,2:18]=temp[j:(j+16),2]
}
TIME_deaths_age <- as.data.frame(TIME_deaths_age)
colnames(TIME_deaths_age) <- colnames(UN_pop_age)
deaths_TIME <- melt(TIME_deaths_age,id="Year")

# load demproj values

# Load TIME values
temp <- as.data.frame(read.table(paste("Demog/",cn,"_dem_deaths_age.txt",sep=""),header=TRUE,fill=TRUE))   
# Need to re-arrage to get in year vs age format
dem_deaths_age <- mat.or.vec(81,18)
dem_deaths_age[,1] <- seq(1970,2050)
for (i in 1:81){
  j <- (i-1)*19+2
  dem_deaths_age[i,2:18]=temp[j:(j+16),2]
}
dem_deaths_age <- as.data.frame(dem_deaths_age)
colnames(dem_deaths_age) <- colnames(UN_pop_age)
deaths_dem <- melt(dem_deaths_age,id="Year")


plot_age_deaths <- ggplot(temp_model_m,aes(x=Year,y=value))+
  geom_line(colour="red")+
  geom_line(data=deaths_TIME,aes(x=Year,y=value/1000),colour="black",linetype="dashed")+
  geom_line(data=deaths_dem,aes(x=Year,y=value/1000),colour="green",linetype="dashed")+
  facet_wrap(~variable,scales="free")+
  xlim(c(1970,2050))
plot_age_deaths
