# Test stuff 

## Plot age specific prevalence from eq_run together with model run

# Main run

# sum up model outputs over age groups and turn into long format to get total pop by age
tot <- mapply(function(x,y) sum(out[x,seq(y+1,35236,81)]),rep(seq(1,81),each=81),seq(1,81))
dim(tot) <- c(81,81)
tot <- t(tot)

pop <- mat.or.vec(81,17)
pop[,17] <- tot[,81]

# Then aggregate into 5 year bins used in TIME
for (i in 1:16){
  t1 <- (i)+(i-1)*4
  t2 <- t1+4
  pop[,i] <- rowSums(tot[,t1:t2])
}

# sum up TB cases over age groups and turn into long format
TB_neg <- mapply(function(x,y) sum(out[x,seq(y+406,1054,81)]),rep(seq(1,81),each=81),seq(1,81))
dim(TB_neg) <- c(81,81)
TB_neg <- t(TB_neg)

prev <- mat.or.vec(81,18)
prev[,1] <- seq(1970,2050)
prev[,18] <- 100000*TB_neg[,81]/pop[,17]

# Then aggregate into 5 year bins used in TIME
for (i in 1:16){
  t1 <- (i)+(i-1)*4
  t2 <- t1+4
  prev[,i+1] <- 100000*rowSums(TB_neg[,t1:t2])/pop[,i]
}

colnames(prev) <- c("Year",colnames(UN_pop_age_t)[3:19])
prev_run <- melt(as.data.frame(prev),id="Year")

# Then check if it matches TIME if we shift cases and pop by one year before we sum up
# Remember, in TIME TB model has 5 year age bins so need to shift cases from across age bin

# Work out what proportion of each 5yr age bin is in the oldest age group - this is the proportion of prevalent cases that should be moved
test<-tot[,seq(5,81,by=5)]
p <- test[,1:16]/pop[,1:16]

# Now shift the age groups, add in births and sum up model outputs over age groups again
tot <- cbind(out[,"births"],tot)

pop <- mat.or.vec(81,17)
pop[,17] <- tot[,81]

# Then aggregate into 5 year bins used in TIME
for (i in 1:16){
  t1 <- (i)+(i-1)*4
  t2 <- t1+4
  pop[,i] <- rowSums(tot[,t1:t2])
}

temp <- mat.or.vec(81,18)
temp[,1] <- seq(1970,2050)
temp[,18] <- TB_neg[,81]

# Now shift TB cases
# First sum into 5 year bins used in TIME
for (i in 1:16){
  t1 <- (i)+(i-1)*4
  t2 <- t1+4
  temp[,i+1] <- rowSums(TB_neg[,t1:t2])
}

prev <- mat.or.vec(81,18)
prev[,1] <- seq(1970,2050)
prev[,18] <- 100000*temp[,81]/pop[,17]
prev[,2] <- 100000*(temp[,2]*(1-p[,1]))/pop[,1]

for (i in 3:17){
  prev[,i] <- 100000*(temp[,i]*(1-p[,i-1]) + temp[,i-1]*p[,i-2])/pop[,i-1]
}

colnames(prev) <- c("Year",colnames(UN_pop_age_t)[3:19])
prev_run_shift <- melt(as.data.frame(prev),id="Year")

# Add TIME age specific prevalence

# Load TIME values
temp <- as.data.frame(read.table(paste("TB/",cn,"/",cn,"_TB_prev_numbers_neg_age.txt",sep=""),header=TRUE,fill=TRUE))   
# Need to re-arrage to get in year vs age format
TB_number_age <- mat.or.vec(81,18)
TB_number_age[,1] <- seq(1970,2050)
for (i in 1:81){
  j <- (i-1)*19+2
  TB_number_age[i,2:18]=temp[j:(j+16),2]/1000
}

prev_TIME <- cbind(seq(1970,2050),100000*TB_number_age[,2:18]/TIME_pop[,3:19])
#prev_TIME <- cbind(seq(1970,2050),TB_number_age[,2:18])
colnames(prev_TIME) <- c("Year",colnames(UN_pop_age_t)[3:19])
prev_TIME <- melt(as.data.frame(prev_TIME),id.vars="Year")

plot_TB_prev_eq <- ggplot(prev_run,aes(x=Year,y=value))+
  geom_line(color="blue")+
  geom_line(data=prev_TIME,aes(x=Year,y=value),color="green")+
  geom_line(data=prev_run_shift,aes(x=Year,y=value),color="blue",lty=2)+
  facet_wrap(~variable,scales="free")+
  ggtitle("HIV- TB prev (by age)")+
  ylab("Prevalence/100,000")+
  xlim(c(1950,2050))












#################################################################################################
## Plot age specific incidence from eq_run together with model run

# Main run

# sum up model outputs over age groups and turn into long format to get total pop by age
tot <- mapply(function(x,y) sum(out[x,seq(y+1,30538,81)]),rep(seq(1,81),each=81),seq(1,81))
dim(tot) <- c(81,81)
tot <- t(tot)

pop <- mat.or.vec(81,17)
pop[,17] <- tot[,81]

# Then aggregate into 5 year bins used in TIME
for (i in 1:16){
  t1 <- (i)+(i-1)*4
  t2 <- t1+4
  pop[,i] <- rowSums(tot[,t1:t2])
}

Inc <- mat.or.vec(81,18)
Inc[,1] <- seq(1970,2050)
Inc[,18] <- 100000*out[,"128"]/pop[,17]

# Then aggregate into 5 year bins used in TIME
for (i in 1:16){
  t1 <- (i-1)*5+30586
  t2 <- t1+4
  Inc[,i+1] <- 100000*rowSums(out[,t1:t2])/pop[,i]
}

colnames(Inc) <- c("Year",colnames(UN_pop_age_t)[3:19])
Inc_run <- melt(as.data.frame(Inc),id="Year")

# Eq_run

# sum up model outputs over age groups and turn into long format to get total pop by age
tot <- mapply(function(x,y) sum(out_eq[x,seq(y+1,30538,81)]),rep(seq(1,101),each=81),seq(1,81))
dim(tot) <- c(81,101)
tot <- t(tot)

pop <- mat.or.vec(101,17)
pop[,17] <- tot[,81]

# Then aggregate into 5 year bins used in TIME
for (i in 1:16){
  t1 <- (i)+(i-1)*4
  t2 <- t1+4
  pop[,i] <- rowSums(tot[,t1:t2])
}

Inc <- mat.or.vec(101,18)
Inc[,1] <- seq(1870,1970)
Inc[,18] <- 100000*out_eq[,"128"]/pop[,17]

# Then aggregate into 5 year bins used in TIME
for (i in 1:16){
  t1 <- (i-1)*5+30586
  t2 <- t1+4
  Inc[,i+1] <- 100000*rowSums(out_eq[,t1:t2])/pop[,i]
}

colnames(Inc) <- c("Year",colnames(UN_pop_age_t)[3:19])
Inc_eq <- melt(as.data.frame(Inc),id="Year")
Inc_eq <- Inc_eq[Inc_eq$Year>1950,]

# Add TIME age specific prevalence

# Load TIME values
temp <- as.data.frame(read.table(paste("TB/",cn,"/",cn,"_TB_Inc_age.txt",sep=""),header=TRUE,fill=TRUE))   
# Need to re-arrage to get in year vs age format
TB_number_age <- mat.or.vec(81,18)
TB_number_age[,1] <- seq(1970,2050)
for (i in 1:81){
  j <- (i-1)*19+2
  TB_number_age[i,2:18]=temp[j:(j+16),2]/1000
}

Inc_TIME <- cbind(out[,"time"],100000*TB_number_age[,2:18]/TIME_pop[,3:19])
colnames(Inc_TIME) <- c("Year",colnames(UN_pop_age_t)[3:19])
Inc_TIME <- melt(Inc_TIME,id.vars="Year")

plot_TB_Inc_eq <- ggplot(Inc_eq,aes(x=Year,y=value))+
  geom_line()+
  geom_line(data=Inc_run,aes(x=Year,y=value),color="blue")+
  geom_line(data=Inc_TIME,aes(x=Year,y=value),color="green")+
  #geom_line(data=prev_eq_a,aes(x=Year,y=value),color="red")+
  #geom_line(data=prev_run_a,aes(x=Year,y=value),color="red")+
  facet_wrap(~variable,scales="free")+
  ggtitle("TB Inc (by age)")+
  ylab("Incidence/100,000")+
  xlim(c(1950,2050))

# plot ratio of prev to inc (approx duration of disease)
# if TIME and R match it suggests that prevalence is out because incidence is out
# if they don't it suggests that something is different in the models  
dur_R<-prev_run$value/Inc_run$value
dur_R_shift<-prev_run_shift$value/Inc_run$value
dur_TIME <-prev_TIME$value/Inc_TIME$value


plot(dur_R,ylim=c(0,4))
points(dur_TIME,col="red")
points(dur_R_shift,col="green")
legend("bottom",c("R","TIME","shifted R"),col=c("black","red","green"),lwd=2)


############################ 
# CHECK POP IS IN EQ IN EQ RUN

# POPULATION #####################################################################################################

## Load and manipulate UN numbers 
## Load UN population data
UN_pop_age <- as.data.frame(read.table(paste("Demog/",cn,"/",cn,"_pop_age.txt",sep=""),header=TRUE)) # Load UN Population data
# add total, births and deaths to data
UN_pop_age_t <- cbind(UN_ind$Year,UN_ind$Births,UN_pop_age[,2:18],rowSums(UN_pop_age[,2:18]),UN_ind$Deaths)
colnames(UN_pop_age_t) <- c("Year","births",colnames(UN_pop_age[2:18]),"Total","Deaths")
temp_data <- melt(UN_pop_age_t,id="Year")

temp_data <- cbind(temp_data,"UNPP")
colnames(temp_data) <- c(colnames(temp_data)[1:3],"Model")

# sum up model outputs over age groups and turn into long format
tot <- mapply(function(x,y) sum(out_eq[x,seq(y+1,30538,81)]),rep(seq(1,101),each=81),seq(1,81))
dim(tot) <- c(81,101)
tot <- t(tot)

model_temp <- mat.or.vec(101,21)
model_temp[,1] <- out_eq[,"time"]
model_temp[,2] <- out_eq[,"births"]
model_temp[,20] <- out_eq[,"Total"]
model_temp[,21] <- out_eq[,"deaths"]
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

# Calculate population age structure
temp_model_s <- as.data.frame(cbind(seq(0,100),100*temp_model[,3:20]/temp_model[,20]))
colnames(temp_model_s)<-c("Year","x4","X9","X14","X19","X24","X29","X34","X39","X44","X49","X54","X59","X64","X69","X74","X79","X100","Total")
temp_model_s <- melt(temp_model_s,id="Year")
temp_model_s <- cbind(temp_model_s,"R")
colnames(temp_model_s) <- c(colnames(temp_model_s)[1:3],"Model")

dat_to_plot <- rbind(temp_model_s,temp_model_noTB)

# and plot
plot_pop_s <- ggplot(dat_to_plot,aes(x=Year,y=value,color=Model))+
  geom_line()+
  facet_wrap(~variable,scales="free")+
  ggtitle("Age structure of population by 5 year age bins")+
  ylab("% of total population")





######### Checking new ART code

# Number who should be on ART by age
ART_on <- out[,35441:35521]

# Number eligible for ART by age
ART_el <- out[,35360:35440]

# Number who need to start ART by age
ART_new <- out[,35279:35359]


## Add up adults and kids who should be on - are these the same as AIM?

on_c <- rowSums(ART_on[,1:15])*1000
on_a <- rowSums(ART_on[,16:81])*1000
temp <- cbind(on_a,on_c)


## Add up adults and kids who are eligible - are these the same as AIM?

el_c <- rowSums(ART_el[,1:15])*1000
el_a <- rowSums(ART_el[,16:81])*1000
temp <- cbind(el_a,el_c)

#### eligilbe has to be bigger than eligible or will break the model 

temp <- ART_el-ART_new




