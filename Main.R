## Script to run the TB model as compiled C code
## It uses the doParallell and foreach packages to run the model on multiple cores on a single computer

## Set working directory - input folders (Demog, HIV, TB) need to be subfolders of this directory
setwd("C:/Users/TOM SUMMER/Documents/TIME_research/TeaTime")

## Load packages
library(deSolve)
library(reshape2)
library(ggplot2)
library(foreach)
library(doParallel)

## Compile the dll of the model (which is coded in C) ############################################################
system("R CMD SHLIB TB_model_v10.c") # Compile

## Define Country (currently South_Africa, Vietnam, Bangladesh) ##################################################
cn <- "South_Africa"

# Load external data sources and create additional forcing functions where necessary #############################
# This takes ~ 0.3 sec so don't want to do it on every iteration of the model - assumes we will ad no uncertianty to these external things (Demographics, HIV, ART) 
source("Data_load.R",local=TRUE)

## each core needs to be able to see all variables set up in data load ###########################################
## Create list of variables to export to each core
to_export <- c("birth_rate",
               "s1","s2","s3","s4","s5","s6","s7","s8","s9","s10","s11","s12","s13","s14","s15","s16","s17","s18","s19","s20",
               "s21","s22","s23","s24","s25","s26","s27","s28","s29","s30","s31","s32","s33","s34","s35","s36","s37","s38","s39","s40",
               "s41","s42","s43","s44","s45","s46","s47","s48","s49","s50","s51","s52","s53","s54","s55","s56","s57","s58","s59","s60",
               "s61","s62","s63","s64","s65","s66","s67","s68","s69","s70","s71","s72","s73","s74","s75","s76","s77","s78","s79","s80","s81",
               "h0","h1","h2","h3","h4","h5","h6","h7","h8","h9","h10","h11","h12","h13","h14","h15","h16","h17","h18","h19","h20",
               "h21","h22","h23","h24","h25","h26","h27","h28","h29","h30","h31","h32","h33","h34","h35","h36","h37","h38","h39","h40",
               "h41","h42","h43","h44","h45","h46","h47","h48","h49","h50","h51","h52","h53","h54","h55","h56","h57","h58","h59","h60",
               "h61","h62","h63","h64","h65","h66","h67","h68","h69","h70","h71","h72","h73","h74","h75","h76","h77","h78","h79","h80",
               "pop_ad",
               "mig0","mig1","mig2","mig3","mig4","mig5","mig6","mig7","mig8","mig9","mig10","mig11","mig12","mig13","mig14","mig15","mig16","mig17","mig18","mig19","mig20",
               "mig21","mig22","mig23","mig24","mig25","mig26","mig27","mig28","mig29","mig30","mig31","mig32","mig33","mig34","mig35","mig36","mig37","mig38","mig39","mig40",
               "mig41","mig42","mig43","mig44","mig45","mig46","mig47","mig48","mig49","mig50","mig51","mig52","mig53","mig54","mig55","mig56","mig57","mig58","mig59","mig60",
               "mig61","mig62","mig63","mig64","mig65","mig66","mig67","mig68","mig69","mig70","mig71","mig72","mig73","mig74","mig75","mig76","mig77","mig78","mig79","mig80",
               "A0","A1","A2","A3","A4","A5","A6","A7","A8","A9","A10","A11","A12","A13","A14","A15","A16","A17","A18","A19","A20",
               "A21","A22","A23","A24","A25","A26","A27","A28","A29","A30","A31","A32","A33","A34","A35","A36","A37","A38","A39","A40",
               "A41","A42","A43","A44","A45","A46","A47","A48","A49","A50","A51","A52","A53","A54","A55","A56","A57","A58","A59","A60",
               "A61","A62","A63","A64","A65","A66","A67","A68","A69","A70","A71","A72","A73","A74","A75","A76","A77","A78","A79","A80",
               "Athresh","UN_pop_start_t")

## Setup parallel backend to use n-1 processors/cores ############################################################
cl<-makeCluster(detectCores()-1)
registerDoParallel(cl)

## The number of model iterations to do ##########################################################################
iters<-500

# start time
strt<-Sys.time()

## The loop ######################################################################################################

# It uses foreach instead of for and %dopar% to launch in parallel (use %do% to run sequentially on single core)
# Outputs are stored in ls as a list by default (returns "TB_out" and "params") 
ls <- foreach(i=1:iters,.packages=c("deSolve"),.export=to_export)%dopar%{

  source("logcurve.R",local=TRUE)
  source(paste("Para_",cn,".R",sep=""),local=TRUE)
  dyn.load("TB_model_v10.dll") # Load
  source("Run_model.R",local=TRUE) 
  list(parms_samp,TB_out)
}

#end time
print(Sys.time()-strt)

stopCluster(cl)


## ls[[i]] gives all of run i
## ls[[i]][[j]] gives j=1 parsm, and j=2 outputs for run i
## ls[[i]][[1]]["name"] gives parameter name for run i
## ls[[i]][[2]][,"name"] gives output name for run i

## Set this up to sample parameters in Para_.... 
## Need to consider what to include - natural history, care and control?
## Need to return all things being varied to ls

## As a first pass, try fitting to incidence, prevalence, notifications, mortality

## Load the WHO estimates
WHO_out <- as.data.frame(read.table(paste("TB/",cn,"/",cn,"_WHO_TB.txt",sep=""),header=TRUE,fill=TRUE))

# Make Notifications range +/- 20% of reported
WHO_out[WHO_out$variable=="Notifications",4] <- WHO_out[WHO_out$variable=="Notifications",3]*0.8
WHO_out[WHO_out$variable=="Notifications",5] <- WHO_out[WHO_out$variable=="Notifications",3]*1.2
# Variance in WHO data asuming range is 95% CI of normal distribution 
WHO_out <- cbind(WHO_out,((((log(WHO_out$High)-log(WHO_out$Low))/2)/1.96)^2))
colnames(WHO_out)[6] <- "Var"

TB_out <- c()
pars_out <- c()

## Calculate likelihood, combine outputs and parameters
L <- rep(0,iters)
for (i in 1:iters){
  
  temp <- melt(ls[[i]][[2]],id.vars="Year")
  temp <- temp[temp$Year>1989 & temp$Year<2014,]
  L[i] <- prod(((2*pi*WHO_out$Var)^(-1/2))*exp(-(1/(2*WHO_out$Var))*((log(temp$value)-log(WHO_out$Mid))^2)))
  
  TB_out <- rbind(TB_out,cbind(ls[[i]][[2]],as.character(i)))
  pars_out <- cbind(pars_out,ls[[i]][[1]])
}

## Resample runs based on likelihood
N_resamp<-100000
t<-sample(seq(1:iters),N_resamp,replace=TRUE,prob=L/sum(L))
unique_t<-unique(t)


# convert runs into plot-able format
colnames(TB_out)[6] <-"Run"
TB_out <- melt(TB_out,id=c("Year","Run"))

# Plots all runs, best in red and WHO envelope in grey
plot_models <- ggplot(TB_out)+
  #geom_line(aes(x=Year,y=value,colour=Run))+
  geom_line(data=WHO_out,aes(x=Year,y=Mid),colour="black",linetype="dashed")+
  geom_line(data=WHO_out,aes(x=Year,y=Low),colour="black")+
  geom_line(data=WHO_out,aes(x=Year,y=High),colour="black")+
  geom_ribbon(data=WHO_out,aes(x=Year,ymin=Low,ymax=High),alpha=0.2)+
  geom_line(data=TB_out[TB_out$Run==which.max(L),],aes(x=Year,y=value),colour="red",size=1.5)+
  facet_wrap(~variable,scales="free_y")+
  xlim(c(1970,2050))+
  ylab("")+
  scale_y_continuous(expand = c(0, 0))+
  expand_limits(y = 0)+theme_bw()






