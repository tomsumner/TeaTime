# Create forcing functions - these are time varying inputs used in the C code ####################################

# They need to consist of x,y pairs where x is the year and y is the value in that year
# They can be entered as numeric values or using the "logcurve" function to simulate a logistic curve:

# logcurve(base, target, base_year, target_year, growth, shape)
# "base" and "target" are the initial and final values
# "shape" determines the symmetry of the curve
# with "shape=1" the curve is symmetric around the midpoint of "base_year" and "target_year"
# "growth" determines the steepness of the curve
# "ss" determines how often to output values

# BCG coverage - currently assume 100% at all times
BCG_cov <- cbind(seq(0,2050),1) # FIXED

# Case detection rate by HIV (neg/pos) - sample parameters of this, need to think about data to constrain the values here
kbase <- runif(1,30,60)
ktarget <- runif(1,90,120)
kbase_year <- runif(1,1970,1990)
ktarget_year <- runif(1,2000,2010)   
kgrowth <- runif(1,0.1,1)

kneg <- logcurve(kbase,ktarget,kbase_year,ktarget_year,kgrowth,1)
kpos <- kneg
#kpos <- logcurve(50,100,1985,2005,0.8,2)  # Old values

# Relative detection smear neg - SAMPLE
rel_d_samp <- runif(1,0.2,0.8)
rel_d <- cbind(seq(1970,2050),rel_d_samp) 

# Relative presentation healthy 
health <- cbind(seq(1970,2050),0.015) # FIX? - TO FIT THIS WOULD NEED NUMBERS SCREENED, NOT SURE THIS DATA IS AVAILABLE FOR ALL COUNTRIES?

# DST coverage among new and previously treated cases 
dstneg_n <- logcurve(0,39,1975,2010,1,2) # FIX BASED ON DATA - as in Targets work
dstneg_p <- logcurve(0,39,1975,2010,1,2) # FIX BASED ON DATA
dstpos_n <- logcurve(0,39,1975,2010,1,2) # FIX BASED ON DATA
dstpos_p <- logcurve(0,39,1975,2010,1,2) # FIX BASED ON DATA

# Sens (se) and spec (sp) of algorithms for sm+ (I), sm- (N) and DST testing (m) - TRY AND BASE THESE ON KNOWN ALGORITHMS?
se_I_neg <- logcurve(100,100,1970,2010,1,2)  # NOT SURE THAT SENS IS IDENTIFIABLE IF DIAGNOSTIC RATES NOT KNOWN? 
se_N_neg <- logcurve(100,100,1970,2010,1,2)
se_m_neg <- logcurve(100,100,1970,2010,1,2)
sp_I_neg <- logcurve(95,95,1970,2010,1,2) # TO FIT SPEC WOULD NEED EITHER NUMBERS SCREENED OR PPV, DON'T THINK IT IS IDENTIFIABLE
sp_N_neg <- logcurve(95,95,1970,2010,1,2)
sp_m_neg <- logcurve(95,95,1970,2010,1,2)

se_I_pos <- logcurve(100,100,1970,2010,1,2)
se_N_pos <- logcurve(100,100,1970,2010,1,2)
se_m_pos <- logcurve(100,100,1970,2010,1,2)
sp_I_pos <- logcurve(95,95,1970,2010,1,2)
sp_N_pos <- logcurve(95,95,1970,2010,1,2)
sp_m_pos <- logcurve(95,95,1970,2010,1,2)

# Linkage to care
l_s <- cbind(seq(1970,2050),0.8) # VARY
l_m <- cbind(seq(1970,2050),0.8) # VARY

# Treatment success by HIV (neg, pos no ART, pos on ART) and susceptibility
tneg_s <- cbind(seq(1970,2050),0.76) # FIX BASED ON OUTCOMES DATA
tpos_s <- cbind(seq(1970,2050),0.76) # FIX BASED ON OUTCOMES DATA
tART_s <- cbind(seq(1970,2050),0.76) # FIX BASED ON OUTCOMES DATA
tneg_m <- cbind(seq(1970,2050),0.5) # FIX BASED ON OUTCOMES DATA
tpos_m <- cbind(seq(1970,2050),0.5) # FIX BASED ON OUTCOMES DATA
tART_m <- cbind(seq(1970,2050),0.5) # FIX BASED ON OUTCOMES DATA

# Set up TB parameters ###########################################################################################

# Fitness of MDR (fit_cost), used to calculate parameter for superinfections (g)
# Both of these are passed into "parms" together with the MDR acquisition rate (e)
fit_cost=0.73 # VARY
g = fit_cost/(1+fit_cost) # superinfections 
e = 0.014 # VARY

# beta = contact rate; a = proportion developing primary TB; p = protection due to previous infection; v = reactivation rate; 
# sig = proportion of cases developing smear positive TB; rel_inf = relative infectiousness of smear negative cases;
# theta = rate of conversion from smear negative to smear positive; r = self cure rate; 
# mu_N = mortality for smear negative cases; mu_I = mortality for smear positive cases; 
# eff_n/p = relative efficacy of first line treatment in new/previously treated MDR cases; 

# _H indicates values for HIV+ (mu_N, mu_I, sig, r and rel_inf)
# other natural history parameters for HIV+ (a,v,p) are adjusted using the rate ratio parameters RR1a etc
# ART modifies a,v,p, muN_H and muI_H by ART_TB1 etc 

# proportion primary (a), proportion smear pos (sig) and mortality rates (muN and muI) take different values for 
# adults (>15) (_a), 0-4 (_0), 5-9 (_5) and 10-14 (_10)

parms <- c(beta = runif(1,10,25), # VARY
           a_a = runif(1,0.08,0.15), a0 = 0.2551, a5 = 0.1357, a10 = 0.0541, # VARY a_a, BUT FIX RRs FOR KIDS?
           v = runif(1,0.0001,0.0025), # VARY
           p = runif(1,0.37,0.9), # VARY
           sig_a = runif(1,0.4,0.5), sig0 = 0.0804, sig5 = 0.0486, sig10 = 0.0994, rel_inf = runif(1,0.1,0.37), theta = runif(1,0.007,0.03), r = 0.2,
           mu_N = runif(1,0.18,0.25), mu_N0 = 0.4205, mu_I = runif(1,0.20,0.41), mu_I0 = 0.6007, fit_cost = fit_cost, e = e, g=g,
           eff_n = 0.61, eff_p = 0.45, 
           muN_H = 0.42, muI_H = 0.6, RR1a = 2.6, RR2a = 1.36, RR1v = 2.6, RR2v = 1.36, RR1p = 0.8, RR2p = 1.3,
           ART_TB1 = 0.204, ART_TB2 = 0.554, ART_TB3 = 0.70, ART_mort1 = 0.232, ART_mort2 = 0.629, ART_mort3 = 0.795,
           BCG_eff = 0.56,
           sig_H = runif(1,0.219,0.425),r_H=0.1,rel_inf_H=runif(1,0.1,0.37),theta_H=runif(1,0.007,0.03))


#### Parameters to return to loop - this includes those above and all others sampled for forcings

parms_samp <- c(parms,
                "kbase"=kbase,"ktarget"=ktarget,"kbase_year"=kbase_year,"ktarget_year"=ktarget_year,"kgrowth"=kgrowth,
                "rel_d"=rel_d_samp) 





