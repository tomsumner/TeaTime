## This script runs the model with the inputs defined elsewhere
## Runs for 100 years to get to equilibrium then rescales pop and runs from 1970 to 2050

# Set up age structure ###########################################################################################
ages <- c(seq(1,80),100) # upper end of age classes
num_ages <- length(ages) # calculates the number of age classes

# Times to run model for #########################################################################################
times <- seq(0,100, by=1)

## Now put together all the forcing fucntions in a list to be passed to the C code ###############################

force <- list(birth_rate,
              s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,
              s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,s32,s33,s34,s35,s36,s37,s38,s39,s40,
              s41,s42,s43,s44,s45,s46,s47,s48,s49,s50,s51,s52,s53,s54,s55,s56,s57,s58,s59,s60,
              s61,s62,s63,s64,s65,s66,s67,s68,s69,s70,s71,s72,s73,s74,s75,s76,s77,s78,s79,s80,s81,
              h0,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h14,h15,h16,h17,h18,h19,h20,
              h21,h22,h23,h24,h25,h26,h27,h28,h29,h30,h31,h32,h33,h34,h35,h36,h37,h38,h39,h40,
              h41,h42,h43,h44,h45,h46,h47,h48,h49,h50,h51,h52,h53,h54,h55,h56,h57,h58,h59,h60,
              h61,h62,h63,h64,h65,h66,h67,h68,h69,h70,h71,h72,h73,h74,h75,h76,h77,h78,h79,h80,
              Ahigh,A500,A349,A249,A199,A99,A50,Athresh,
              BCG_cov,pop_ad,
              kneg,kpos,rel_d,dstneg_n,dstneg_p,dstpos_n,dstpos_p,l_s,l_m,tneg_s,tpos_s,tART_s,tneg_m,tpos_m,tART_m,
              mig0,mig1,mig2,mig3,mig4,mig5,mig6,mig7,mig8,mig9,mig10,mig11,mig12,mig13,mig14,mig15,mig16,mig17,mig18,mig19,mig20,
              mig21,mig22,mig23,mig24,mig25,mig26,mig27,mig28,mig29,mig30,mig31,mig32,mig33,mig34,mig35,mig36,mig37,mig38,mig39,mig40,
              mig41,mig42,mig43,mig44,mig45,mig46,mig47,mig48,mig49,mig50,mig51,mig52,mig53,mig54,mig55,mig56,mig57,mig58,mig59,mig60,
              mig61,mig62,mig63,mig64,mig65,mig66,mig67,mig68,mig69,mig70,mig71,mig72,mig73,mig74,mig75,mig76,mig77,mig78,mig79,mig80,
              se_I_neg,se_N_neg,se_m_neg,sp_I_neg,sp_N_neg,sp_m_neg,
              se_I_pos,se_N_pos,se_m_pos,sp_I_pos,sp_N_pos,sp_m_pos,
              health)

# EQUILIBRIUM RUN ################################################################################################

# Initial conditions - all susceptible
temp <- c()
for (i in 1:num_ages){temp[i]<-as.numeric(UN_pop_start_t[i+2])}
xstart <- c(S=c(temp),
            Lsn=rep(0,num_ages),Lsp=rep(0,num_ages),Lmn=rep(0,num_ages),Lmp=rep(0,num_ages),
            Nsn=rep(0,num_ages),Nsp=rep(0,num_ages),Nmn=rep(0,num_ages),Nmp=rep(0,num_ages),
            Isn=c(rep(0,25),0.1,rep(0,55)),Isp=rep(0,num_ages),Imn=rep(0,num_ages),Imp=rep(0,num_ages),
            PTn=rep(0,num_ages),PTp=rep(0,num_ages),
            S_H=rep(0,num_ages*7),
            Lsn_H=rep(0,num_ages*7),Lsp_H=rep(0,num_ages*7),Lmn_H=rep(0,num_ages*7),Lmp_H=rep(0,num_ages*7),
            Nsn_H=rep(0,num_ages*7),Nsp_H=rep(0,num_ages*7),Nmn_H=rep(0,num_ages*7),Nmp_H=rep(0,num_ages*7),
            Isn_H=rep(0,num_ages*7),Isp_H=rep(0,num_ages*7),Imn_H=rep(0,num_ages*7),Imp_H=rep(0,num_ages*7),
            PTn_H=rep(0,num_ages*7),PTp_H=rep(0,num_ages*7),
            S_A=rep(0,num_ages*7*3),
            Lsn_A=rep(0,num_ages*7*3),Lsp_A=rep(0,num_ages*7*3),Lmn_A=rep(0,num_ages*7*3),Lmp_A=rep(0,num_ages*7*3),
            Nsn_A=rep(0,num_ages*7*3),Nsp_A=rep(0,num_ages*7*3),Nmn_A=rep(0,num_ages*7*3),Nmp_A=rep(0,num_ages*7*3),
            Isn_A=rep(0,num_ages*7*3),Isp_A=rep(0,num_ages*7*3),Imn_A=rep(0,num_ages*7*3),Imp_A=rep(0,num_ages*7*3),
            PTn_A=rep(0,num_ages*7*3),PTp_A=rep(0,num_ages*7*3))

# For initialisation run turn off MDR by setting e = 0
parms["e"]=0
# For initialisation run no HIV so don't run those parameters
parms["HIV_run"]=0

# Run the model
time_eq <- system.time(out_eq <- ode(y=xstart, times, func = "derivsc",
                                     parms = parms, dllname = "TB_model_v8",initforc = "forcc",
                                     forcings=force, initfunc = "parmsc", nout = 53,
                                     outnames = c("Total","Total_S","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                                                  "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM",
                                                  "CD4500","CD4350_500","CD4250_349","CD4200_249","CD4100_199","CD450_99","CD450",
                                                  "ART500","ART350_500","ART250_349","ART200_249","ART100_199","ART50_99","ART50",
                                                  "v1","v2","v3","v4","v5","v6","v7","ART_tot","ART_need","ART_new","ART_on","TB_deaths",
                                                  "Cases_neg","Cases_pos","Cases_ART",
                                                  "births","deaths","ART_deaths","Total_PT",
                                                  "DS_correct","DS_incorrect","MDR_correct","MDR_incorrect","FP"), 
                                     events = list(func="event",time=seq(0,100)),
                                     method = rkMethod("rk45dp7",hmax=1)))
                                     
# PROJECTION RUN #################################################################################################

# Adjust pop down to 1970 values (by age) and reassign initial conditions - model can now be run from 1970 with TB and HIV

temp <- out_eq[dim(out_eq)[1],2:35236]

for(i in 1:81){ 
  temp[seq(i,35235,81)] <- temp[seq(i,35235,81)]/(sum(temp[seq(i,35235,81)])/as.numeric(UN_pop_start_t[i+2]))
}

xstart <- temp

# Reset e to allow MDR
parms["e"]=e
# Reset to model HIV
parms["HIV_run"]=1

# Set times to run for
times <- seq(1970,2050 , by=ss) # run with 6 month time step using a fixed time step solver - this is faster than adaptive methds but seems to give good accuracy
# Run the model
time_run <-system.time(out <- ode(y=xstart, times, func = "derivsc",
                                  parms = parms, dllname = "TB_model_v8",initforc = "forcc",
                                  forcings=force, initfunc = "parmsc", nout = 53,
                                  outnames = c("Total","Total_S","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                                               "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM",
                                               "CD4500","CD4350_500","CD4250_349","CD4200_249","CD4100_199","CD450_99","CD450",
                                               "ART500","ART350_500","ART250_349","ART200_249","ART100_199","ART50_99","ART50",
                                               "v1","v2","v3","v4","v5","v6","v7","ART_tot","ART_need","ART_new","ART_on","TB_deaths",
                                               "Cases_neg","Cases_pos","Cases_ART",
                                               "births","deaths","ART_deaths","Total_PT",
                                               "DS_correct","DS_incorrect","MDR_correct","MDR_incorrect","FP"), 
                                  events = list(func="event",time=seq(1970,2050)),
                                  method = rkMethod("rk45dp7",hmax=1)))
                                  
# Just keep every other output now we are running with 6 month time step
#out <- out[seq(1,length(times),1/ss),]
