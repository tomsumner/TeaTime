## This script runs the model with the inputs defined elsewhere
## Runs for 100 years to get to equilibrium then rescales pop and runs from 1970 to 2050

# Times to run model for
times <- seq(0,60, by=1)

# Initial conditions - all susceptible
temp <- c()
for (i in 1:num_ages){temp[i]<-UN_pop_age[21,i+1]}
xstart <- c(S=c(temp),
            Lsn=rep(0,num_ages),Lsp=rep(0,num_ages),Lmn=rep(0,num_ages),Lmp=rep(0,num_ages),
            Nsn=rep(0,num_ages),Nsp=rep(0,num_ages),Nmn=rep(0,num_ages),Nmp=rep(0,num_ages),
            Isn=c(rep(0,25),100,rep(0,55)),Isp=rep(0,num_ages),Imn=rep(0,num_ages),Imp=rep(0,num_ages),
            S_H=rep(0,num_ages*7),
            Lsn_H=rep(0,num_ages*7),Lsp_H=rep(0,num_ages*7),Lmn_H=rep(0,num_ages*7),Lmp_H=rep(0,num_ages*7),
            Nsn_H=rep(0,num_ages*7),Nsp_H=rep(0,num_ages*7),Nmn_H=rep(0,num_ages*7),Nmp_H=rep(0,num_ages*7),
            Isn_H=rep(0,num_ages*7),Isp_H=rep(0,num_ages*7),Imn_H=rep(0,num_ages*7),Imp_H=rep(0,num_ages*7),
            S_A=rep(0,num_ages*7*3),
            Lsn_A=rep(0,num_ages*7*3),Lsp_A=rep(0,num_ages*7*3),Lmn_A=rep(0,num_ages*7*3),Lmp_A=rep(0,num_ages*7*3),
            Nsn_A=rep(0,num_ages*7*3),Nsp_A=rep(0,num_ages*7*3),Nmn_A=rep(0,num_ages*7*3),Nmp_A=rep(0,num_ages*7*3),
            Isn_A=rep(0,num_ages*7*3),Isp_A=rep(0,num_ages*7*3),Imn_A=rep(0,num_ages*7*3),Imp_A=rep(0,num_ages*7*3))

# For initialisation run turn off MDR by setting e = 0
parms["e"]=0

# Run the model
time_eq <- system.time(out_eq <- ode(y=xstart, times, func = "derivsc",
                                     parms = parms, dllname = "TB_model_v6",initforc = "forcc",
                                     forcings=force, initfunc = "parmsc", nout = 46,
                                     outnames = c("Total","Total_S","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                                                  "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM",
                                                  "CD4500","CD4350_500","CD4250_349","CD4200_249","CD4100_199","CD450_99","CD450",
                                                  "ART500","ART350_500","ART250_349","ART200_249","ART100_199","ART50_99","ART50",
                                                  "v1","v2","v3","v4","v5","v6","v7","ART_tot","ART_need","ART_new","ART_on","TB_deaths",
                                                  "Cases_neg","Cases_pos","Cases_ART",
                                                  "births","deaths"), 
                                     events = list(func="event",time=seq(0,60)),
                                     method = rkMethod("rk4")))

# Adjust pop down to 1970 values (by age) and reassign initial conditions - model can now be run from 1970 with TB and HIV

temp <- out_eq[dim(out_eq)[1],2:30538]

for(i in 1:81){ 
  temp[seq(i,30537,81)] <- temp[seq(i,30537,81)]/(sum(temp[seq(i,30537,81)])/UN_pop_age_t[UN_pop_age_t$Year==1970,i+2])
}

xstart <- temp

# Reset e to allow MDR
parms["e"]=e

# Set times to run for
times <- seq(1970,2050 , by=0.5) # run with 6 month time step using a fixed time step solver - this is faster than adaptive methds but seems to give good accuracy
# Run the model
time_run <-system.time(out <- ode(y=xstart, times, func = "derivsc",
                                  parms = parms, dllname = "TB_model_v6",initforc = "forcc",
                                  forcings=force, initfunc = "parmsc", nout = 46,
                                  outnames = c("Total","Total_S","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                                               "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM",
                                               "CD4500","CD4350_500","CD4250_349","CD4200_249","CD4100_199","CD450_99","CD450",
                                               "ART500","ART350_500","ART250_349","ART200_249","ART100_199","ART50_99","ART50",
                                               "v1","v2","v3","v4","v5","v6","v7","ART_tot","ART_need","ART_new","ART_on","TB_deaths",
                                               "Cases_neg","Cases_pos","Cases_ART",
                                               "births","deaths"), 
                                  events = list(func="event",time=seq(1970,2050)),
                                  method = rkMethod("rk4")))

# Just keep every other output now we are running with 6 month time step
out <- out[seq(1,length(times),2),]