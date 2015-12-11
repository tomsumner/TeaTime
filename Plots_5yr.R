# CORE PLOTS - DEMOGRPAHY AND TB #################################################################################

# Calculate number of model states

n_state <- n_dis*n_age*(1+n_HIV+n_HIV*n_ART)

# Demography #####################################################################################################

## Load and manipulate UN numbers 
## Load UN population data
UN_pop_age <- as.data.frame(read.table(paste("Demog/",cn,"/",cn,"_pop_age.txt",sep=""),header=TRUE)) # Load UN Population data
# add total, births and deaths to data
UN_pop_age_t <- cbind(UN_ind$Year,UN_pop_age[,2:18],rowSums(UN_pop_age[,2:18]),UN_ind$Births,UN_ind$Deaths)
colnames(UN_pop_age_t) <- c("Year","0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49",
                            "50-54","55-59","60-64","65_69","70-74","75-79","80+","Total","Births","Deaths")
temp_data <- melt(UN_pop_age_t,id="Year")

temp_data <- cbind(temp_data,"UNPP")
colnames(temp_data) <- c(colnames(temp_data)[1:3],"Model")

# sum up model outputs over age groups and turn into long format
tot <- mapply(function(x,y) sum(out[x,seq(y+1,n_state+1,17)]),rep(seq(1,81),each=17),seq(1,17))
dim(tot) <- c(17,81)
tot <- t(tot)

model_temp <- mat.or.vec(81,21)
model_temp[,1] <- out[,"time"]
model_temp[,20] <- out[,"Births"]
model_temp[,19] <- out[,"Total"]
model_temp[,21] <- out[,"Deaths"]
model_temp[,2:18] <- tot[,1:17]

temp_model <- as.data.frame(model_temp)
colnames(temp_model) <- colnames(UN_pop_age_t)
temp_model_m <- melt(temp_model,id="Year")
temp_model_m <- cbind(temp_model_m,"R")
colnames(temp_model_m) <- c(colnames(temp_model_m)[1:3],"Model")

dat_to_plot <- rbind(temp_data,temp_model_m)

plot_pop_5yr <- ggplot(dat_to_plot,aes(x=Year,y=value,colour=Model))+
  geom_line()+
  facet_wrap(~variable,scales="free")+
  ggtitle("Population (by 5 year age bins), total births and total deaths")+
  xlim(c(1970,2050))+
  ylab("Thousands")+theme_bw()

# TB outputs #####################################################################################################

# Mortality (HIV- and HIV+ lines)
# Prevalence (All)
# Incidence (All, HIV+, notified)

# All show model run and WHO envelope 

# Load the WHO estimates
WHO_out <- as.data.frame(read.table(paste("TB/",cn,"/",cn,"_WHO_TB.txt",sep=""),header=TRUE,fill=TRUE))

# Arrange model output
R_out <- as.data.frame(cbind(out[,"time"],
                             100000*out[,"TB_deaths"]/out[,"Total"], # Mort (all)
                             100000*(out[,"Total_DS"]+out[,"Total_MDR"])/out[,"Total"],  # Prev
                             100000*(out[,"Cases_neg"]+out[,"Cases_pos"]+out[,"Cases_ART"])/out[,"Total"], # Inc, all 
                             100000*(out[,"Cases_pos"]+out[,"Cases_ART"])/out[,"Total"], # Inc, HIV+ 
                             100000*(out[,"DS_correct"]+out[,"DS_incorrect"]+out[,"MDR_correct"]+out[,"MDR_incorrect"]+out[,"FP"])/out[,"Total"])) # Total notifcations

colnames(R_out) <- c("Year","V1","V2","V3","V4","V5")
R_out <- R_out[R_out$Year<=2050,]

# combine and melt
Models_out <- R_out
Models_out <- melt(Models_out,id=c("Year"))
Models_out <- cbind(Models_out,c(rep("Mortality",81),rep("Prevalence",81),rep("Incidence",3*81)),c(rep("HIV-",81),rep("All",2*81),rep("HIV+",81),rep("Notifications",81)))
colnames(Models_out) <- c("Year","variable","value","type","group")

plot_TB <- ggplot(Models_out)+
  geom_line(aes(x=Year,y=value,colour=group),size=1)+
  geom_line(data=WHO_out,aes(x=Year,y=mid,colour=group),linetype="dashed")+
  geom_line(data=WHO_out,aes(x=Year,y=lo,colour=group))+
  geom_line(data=WHO_out,aes(x=Year,y=hi,colour=group))+
  geom_ribbon(data=WHO_out,aes(x=Year,ymin=lo,ymax=hi,fill=group),alpha=0.2)+
  facet_wrap(~type,scales="free_y",nrow=3)+
  xlim(c(1990,2050))+
  ylab("")+
  scale_y_continuous(expand = c(0, 0))+
  expand_limits(y = 0)+theme_bw()+theme(legend.title = element_blank()

