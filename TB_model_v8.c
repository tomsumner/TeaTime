/* TB model in C code to call from R */

/* Version 8 - converted to single year age bins and aging done via events- this ensures individuals remain in the correct age cohorts */
/*             matches TIME within reasonable bounds with the exception of ART                                                         */
/*             this version updated to include true and false positive diagnosis - HIV- eqns done need to add parameters               */
/*             and post PT compartment added - just for false positives at present - no PT modelled yet                                */

/* Can be compiled within R with system("R CMD SHLIB TB_model_v8.c") */
/* This creates a dynamic linked library (.dll) which can be loaded (dyn.load(TB_model_v4.dll)) into R and used as the model fucntion in a call to desolve */

/* C libraries needed */
#include <R.h>
#include <math.h>

/* You need to define number of parameters and forcing functions passed to the model here */
/* These must match number in intializer functions below */
static double parms[43];
static double forc[282];

/* ###### A TRICK TO KEEP UP WITH THE PARAMETERS AND FORCINGS ###### */

/* !!!!!!! NOTE IN C INDEXING STARTS AT 0 AND NOT 1 (i.e. parms[0] is the first value in parms) !!!!!! */

/* Parameters are constant over time */
#define beta parms[0]       /* effective contact rate */
#define a_a parms[1]        /* proportion developing primary disease in adults (>15) */
#define a0 parms[2]         /* proportion developing primary disease in 0-4 year olds */ 
#define a5 parms[3]         /* proportion developing primary disease in 5-9 year olds */
#define a10 parms[4]        /* proportion developing primary disease in 10-14 year olds */
#define p parms[5]          /* protection against disease due to prior infection */
#define v parms[6]          /* reactivation rate */
#define sig_a parms[7]      /* proportion developing smear positive disease in adults (>15) */
#define sig0 parms[8]       /* proportion developing smear positive disease in 0-4 year olds */
#define sig5 parms[9]       /* proportion developing smear positive disease in 5-9 year olds */
#define sig10 parms[10]     /* proportion developing smear positive disease in 10-14 year olds */
#define rel_inf parms[11]   /* relative infectiousness of smear negative TB */
#define theta parms[12]     /* rate of conversion from smear negative to smear positive TB */
#define r parms[13]         /* rate of self cure from TB */
#define mu_N parms[14]      /* mortality rate from smear negative TB (>5 years old) */
#define mu_N0 parms[15]     /* mortality rate from smear negative TB in 0-4 year olds */
#define mu_I parms[16]      /* mortality rate from smear positive TB (>5 years old) */
#define mu_I0 parms[17]     /* mortality rate from smear positive TB in 0-4 year olds */
#define fit_cost parms[18]  /* relative fitness of MDR strains (transmission only) */
#define e parms[19]         /* rate of acquisition of MDR */
#define g parms[20]         /* superinfections */
#define eff_n parms[21]     /* efficacy of first line drugs in treating new MDR cases */
#define eff_p parms[22]     /* efficacy of first line drugs in treating retreatment MDR cases */
#define muN_H parms[23]     /* mortaltiy rate from smear negative TB in HIV+ */
#define muI_H parms[24]     /* mortaltiy rate from smear positive TB in HIV+ */
#define RR1a parms[25]      /* Relative risk for primary disease following HIV infection */
#define RR2a parms[26]      /* Relative risk for primary disease by CD4 */
#define RR1v parms[27]      /* Relative risk for reactivation following HIV infection */
#define RR2v parms[28]      /* Relative risk for reactivation by CD4 */
#define RR1p parms[29]      /* Relative risk of protection against reinfection following HIV infection */
#define RR2p parms[30]      /* Relative risk of protection against reinfection by CD4 */
#define ART_TB1 parms[31]   /* Reduction in TB parameters on ART (<6m) */
#define ART_TB2 parms[32]   /* Reduction in TB parameters on ART (7-12m) */
#define ART_TB3 parms[33]   /* Reduction in TB parameters on ART (>12m) */
#define ART_mort1 parms[34] /* Reduction in TB mortality on ART (<6m) */
#define ART_mort2 parms[35] /* Reduction in TB mortality on ART (7-12m) */
#define ART_mort3 parms[36] /* Reduction in TB mortality on ART (>12m) */
#define BCG_eff parms[37]   /* Efficacy of BCG (reduces primary and reactivation risks) */
#define sig_H parms[38]     /* proportion developing smear positive disease in HIV+ */
#define r_H parms[39]       /* rate of self cure from TB in HIV+ */
#define rel_inf_H parms[40] /* relative infectiousness of smear negative cases in HIV+ */
#define theta_H parms[41]   /* rate of conversion from smear negative to smear positive TB in HIV+ */ 

#define HIV_run parms[42]   /* used to skip HIV equations in equilibrium run */

/* Forcings are time dependant functions e.g. mortality rates */
#define birth_rate forc[0] /* Birth rate */
#define s1 forc[1]         /* Survival to age x */
#define s2 forc[2] 
#define s3 forc[3] 
#define s4 forc[4] 
#define s5 forc[5]
#define s6 forc[6]
#define s7 forc[7]
#define s8 forc[8]
#define s9 forc[9]
#define s10 forc[10]
#define s11 forc[11]
#define s12 forc[12]
#define s13 forc[13]
#define s14 forc[14]
#define s15 forc[15]
#define s16 forc[16]
#define s17 forc[17]
#define s18 forc[18]       
#define s19 forc[19] 
#define s20 forc[20] 
#define s21 forc[21] 
#define s22 forc[22]
#define s23 forc[23]
#define s24 forc[24]
#define s25 forc[25]
#define s26 forc[26]
#define s27 forc[27]
#define s28 forc[28]
#define s29 forc[29]
#define s30 forc[30]
#define s31 forc[31]
#define s32 forc[32]
#define s33 forc[33]
#define s34 forc[34]
#define s35 forc[35]
#define s36 forc[36]
#define s37 forc[37]
#define s38 forc[38]         
#define s39 forc[39] 
#define s40 forc[40] 
#define s41 forc[41] 
#define s42 forc[42]
#define s43 forc[43]
#define s44 forc[44]
#define s45 forc[45]
#define s46 forc[46]
#define s47 forc[47]
#define s48 forc[48]
#define s49 forc[49]
#define s50 forc[50]
#define s51 forc[51]
#define s52 forc[52]
#define s53 forc[53]
#define s54 forc[54]
#define s55 forc[55]
#define s56 forc[56]
#define s57 forc[57]
#define s58 forc[58]        
#define s59 forc[59] 
#define s60 forc[60] 
#define s61 forc[61] 
#define s62 forc[62]
#define s63 forc[63]
#define s64 forc[64]
#define s65 forc[65]
#define s66 forc[66]
#define s67 forc[67]
#define s68 forc[68]
#define s69 forc[69]
#define s70 forc[70]
#define s71 forc[71] 
#define s72 forc[72]
#define s73 forc[73]
#define s74 forc[74]
#define s75 forc[75]
#define s76 forc[76]
#define s77 forc[77]
#define s78 forc[78]
#define s79 forc[79]
#define s80 forc[80]
#define s81 forc[81]

#define h0 forc[82]  /* HIV incidence at age x */
#define h1 forc[83]
#define h2 forc[84]
#define h3 forc[85]
#define h4 forc[86]
#define h5 forc[87]
#define h6 forc[88]
#define h7 forc[89]
#define h8 forc[90]
#define h9 forc[91]
#define h10 forc[92]
#define h11 forc[93]
#define h12 forc[94]
#define h13 forc[95]
#define h14 forc[96]
#define h15 forc[97]
#define h16 forc[98]
#define h17 forc[99]
#define h18 forc[100]       
#define h19 forc[101] 
#define h20 forc[102] 
#define h21 forc[103] 
#define h22 forc[104]
#define h23 forc[105]
#define h24 forc[106]
#define h25 forc[107]
#define h26 forc[108]
#define h27 forc[109]
#define h28 forc[110]
#define h29 forc[111]
#define h30 forc[112]
#define h31 forc[113]
#define h32 forc[114]
#define h33 forc[115]
#define h34 forc[116]
#define h35 forc[117]
#define h36 forc[118]
#define h37 forc[119]
#define h38 forc[120]         
#define h39 forc[121] 
#define h40 forc[122] 
#define h41 forc[123] 
#define h42 forc[124]
#define h43 forc[125]
#define h44 forc[126]
#define h45 forc[127]
#define h46 forc[128]
#define h47 forc[129]
#define h48 forc[130]
#define h49 forc[131]
#define h50 forc[132]
#define h51 forc[133]
#define h52 forc[134]
#define h53 forc[135]
#define h54 forc[136]
#define h55 forc[137]
#define h56 forc[138]
#define h57 forc[139]
#define h58 forc[140]        
#define h59 forc[141] 
#define h60 forc[142] 
#define h61 forc[143] 
#define h62 forc[144]
#define h63 forc[145]
#define h64 forc[146]
#define h65 forc[147]
#define h66 forc[148]
#define h67 forc[149]
#define h68 forc[150]
#define h69 forc[151]
#define h70 forc[152]
#define h71 forc[153] 
#define h72 forc[154]
#define h73 forc[155]
#define h74 forc[156]
#define h75 forc[157]
#define h76 forc[158]
#define h77 forc[159]
#define h78 forc[160]
#define h79 forc[161]
#define h80 forc[162]

#define Ahigh forc[163] /* ART coverage by CD4 - depends on overall coverage and eligibility threshold */
#define A500 forc[164]
#define A349 forc[165]
#define A249 forc[166]
#define A199 forc[167]
#define A99 forc[168]
#define A50 forc[169]
#define Athresh forc[170] /* ART eligibility threshold in terms of model CD4 categories i.e. 1 means anyone <500 */

#define BCG_cov forc[171] /* BCG coverage */

#define pop_ad forc[172]  /* used to turn on/off adjustment of population to account for disease induced mortality - idea is to turn it off from 2015 onwards*/

#define kneg forc[173]  /* detection rate */
#define kpos forc[174]

#define rel_d forc[175] /* relative detection of smear negatives */

#define dstneg_n forc[176]   /* dst coverage in new cases */
#define dstneg_p forc[177]   /* dst coverage in previoulsy treated cases */
#define dstpos_n forc[178]  
#define dstpos_p forc[179] 

#define l_s forc[180]  /* linkage to care */
#define l_m forc[181]

#define tneg_s forc[182]  /* treatment success */
#define tpos_s forc[183]
#define tART_s forc[184]
#define tneg_m forc[185]
#define tpos_m forc[186]
#define tART_m forc[187]

#define mig0 forc[188]  /* migration by age x */
#define mig1 forc[189]
#define mig2 forc[190]
#define mig3 forc[191]
#define mig4 forc[192]
#define mig5 forc[193]
#define mig6 forc[194]
#define mig7 forc[195]
#define mig8 forc[196]
#define mig9 forc[197]
#define mig10 forc[198]
#define mig11 forc[199]
#define mig12 forc[200]
#define mig13 forc[201]
#define mig14 forc[202]
#define mig15 forc[203]
#define mig16 forc[204]
#define mig17 forc[205]
#define mig18 forc[206]       
#define mig19 forc[207] 
#define mig20 forc[208] 
#define mig21 forc[209] 
#define mig22 forc[210]
#define mig23 forc[211]
#define mig24 forc[212]
#define mig25 forc[213]
#define mig26 forc[214]
#define mig27 forc[215]
#define mig28 forc[216]
#define mig29 forc[217]
#define mig30 forc[218]
#define mig31 forc[219]
#define mig32 forc[220]
#define mig33 forc[221]
#define mig34 forc[222]
#define mig35 forc[223]
#define mig36 forc[224]
#define mig37 forc[225]
#define mig38 forc[226]         
#define mig39 forc[227] 
#define mig40 forc[228] 
#define mig41 forc[229] 
#define mig42 forc[230]
#define mig43 forc[231]
#define mig44 forc[232]
#define mig45 forc[233]
#define mig46 forc[234]
#define mig47 forc[245]
#define mig48 forc[236]
#define mig49 forc[237]
#define mig50 forc[238]
#define mig51 forc[239]
#define mig52 forc[240]
#define mig53 forc[241]
#define mig54 forc[242]
#define mig55 forc[243]
#define mig56 forc[244]
#define mig57 forc[245]
#define mig58 forc[246]        
#define mig59 forc[247] 
#define mig60 forc[248] 
#define mig61 forc[249] 
#define mig62 forc[250]
#define mig63 forc[251]
#define mig64 forc[252]
#define mig65 forc[253]
#define mig66 forc[254]
#define mig67 forc[255]
#define mig68 forc[256]
#define mig69 forc[257]
#define mig70 forc[258]
#define mig71 forc[259] 
#define mig72 forc[260]
#define mig73 forc[261]
#define mig74 forc[262]
#define mig75 forc[263]
#define mig76 forc[264]
#define mig77 forc[265]
#define mig78 forc[266]
#define mig79 forc[267]
#define mig80 forc[268]

/* sens and spec for tests */
#define se_I_neg forc[269]
#define se_N_neg forc[270]
#define se_m_neg forc[271]

#define sp_I_neg forc[272]
#define sp_N_neg forc[273]
#define sp_m_neg forc[274]

#define se_I_pos forc[275]
#define se_N_pos forc[276]
#define se_m_pos forc[277]

#define sp_I_pos forc[278]
#define sp_N_pos forc[279]
#define sp_m_pos forc[280]

/* RR for presentation for healthy individuals */
#define health forc[281]

/* ###### FUNCTION TO SUM ARRAY FROM ELEMENT i_start TO i_end ###### */
double sumsum(double ar[], int i_start, int i_end)
{
   int i=0;
   double sum=0;
   for (i=i_start; i<=i_end; i++)
   {
    sum = sum + ar[i];
   }
   return(sum);
}

/* ###### FUNCTION TO INITIALIZE PARAMETERS PASSED FROM R - if the number of parameters is changed you must update N here ###### */
void parmsc(void (* odeparms)(int *, double *))
{
    int N=43;
    odeparms(&N, parms);
}

/* ###### FUNCTION TO INITIALIZE FORCINGS PASSED FROM R - if the number of parameters is changed you must update N here ###### */
void forcc(void (* odeforcs)(int *, double *))
{
    int N=282;
    odeforcs(&N, forc);
}

/* ##### EVENTS ARE USED TO ADD BIRTHS AND SHIFT THE POPULATION BY AGE - EQUIVALENT TO THE METHOD OF SCHENZLE ###### */

void event(int *n, double *t, double *y) 
{
  int i;
  
  /* Store current population in temp and shift every age group forward one */
  double temp[35235];
  temp[0] = y[0];
  for (i=1; i<35235; i++){
    temp[i] = y[i];
    y[i] = temp[i-1];
  }
  /* Set every 0 age group to zero and every >80 age group to the previous age group plus those still surviving  */
  for (i=0; i<35235; i+=81) {
    y[i] = 0;
    y[i+80] = temp[i+80] + temp[i+79];  
  }
  /* Then add births into group 0 - only susceptibles get born */ 
  y[0] = birth_rate*sumsum(temp,0,35234)/1000;
}

/* ###### DERIVATIVE FUNCTIONS - THIS IS THE MODEL ITSELF ###### */

void derivsc(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if (ip[0] <2) error("nout should be at least 2");
    
    /* Expand state variables so can use more meaningful names than y and ydot (the variables and rates of change are vectors y and ydot here) */
    /* There are 81 age groups, 7 HIV positive categories and 3 times on ART */
    /* _H = HIV+; _A = HIV+ on ART */
    
    /* These are the variables */
    double S[81]={0};   double S_H[81][7]={{0}};   double S_A[81][7][3]={{{0}}};      /* Susceptible */
    double Lsn[81]={0}; double Lsn_H[81][7]={{0}}; double Lsn_A[81][7][3]={{{0}}};    /* Latent, DS, new */
    double Lsp[81]={0}; double Lsp_H[81][7]={{0}}; double Lsp_A[81][7][3]={{{0}}};    /* Latent, DS, previous */
    double Lmn[81]={0}; double Lmn_H[81][7]={{0}}; double Lmn_A[81][7][3]={{{0}}};    /* Latent, DR, new */
    double Lmp[81]={0}; double Lmp_H[81][7]={{0}}; double Lmp_A[81][7][3]={{{0}}};    /* Latent, DR, previous */
    double Nsn[81]={0}; double Nsn_H[81][7]={{0}}; double Nsn_A[81][7][3]={{{0}}};    /* Smear negative, DS, new */
    double Nsp[81]={0}; double Nsp_H[81][7]={{0}}; double Nsp_A[81][7][3]={{{0}}};    /* Smear negative, DS, previous */
    double Nmn[81]={0}; double Nmn_H[81][7]={{0}}; double Nmn_A[81][7][3]={{{0}}};    /* Smear negative, DR, new */
    double Nmp[81]={0}; double Nmp_H[81][7]={{0}}; double Nmp_A[81][7][3]={{{0}}};    /* Smear negative, DR, previous */
    double Isn[81]={0}; double Isn_H[81][7]={{0}}; double Isn_A[81][7][3]={{{0}}};    /* Smear positive, DS, new */
    double Isp[81]={0}; double Isp_H[81][7]={{0}}; double Isp_A[81][7][3]={{{0}}};    /* Smear positive, DS, previous */
    double Imn[81]={0}; double Imn_H[81][7]={{0}}; double Imn_A[81][7][3]={{{0}}};    /* Smear positive, DR, new */
    double Imp[81]={0}; double Imp_H[81][7]={{0}}; double Imp_A[81][7][3]={{{0}}};    /* Smear positive, DR, previous */
    double PTn[81]={0}; double PTn_H[81][7]={{0}}; double PTn_A[81][7][3]={{{0}}};     /* Post PT, new - also move people here if they are false positive for TB and receive Rx */
    double PTp[81]={0}; double PTp_H[81][7]={{0}}; double PTp_A[81][7][3]={{{0}}};     /* Post PT, previous - also move people here if they are false positive for TB and receive Rx */

    /* These are the rates of change (same names but prefixed with d) */
    double dS[81]={0};   double dS_H[81][7]={{0}};   double dS_A[81][7][3]={{{0}}};
    double dLsn[81]={0}; double dLsn_H[81][7]={{0}}; double dLsn_A[81][7][3]={{{0}}};
    double dLsp[81]={0}; double dLsp_H[81][7]={{0}}; double dLsp_A[81][7][3]={{{0}}};
    double dLmn[81]={0}; double dLmn_H[81][7]={{0}}; double dLmn_A[81][7][3]={{{0}}};
    double dLmp[81]={0}; double dLmp_H[81][7]={{0}}; double dLmp_A[81][7][3]={{{0}}};
    double dNsn[81]={0}; double dNsn_H[81][7]={{0}}; double dNsn_A[81][7][3]={{{0}}};
    double dNsp[81]={0}; double dNsp_H[81][7]={{0}}; double dNsp_A[81][7][3]={{{0}}};
    double dNmn[81]={0}; double dNmn_H[81][7]={{0}}; double dNmn_A[81][7][3]={{{0}}};
    double dNmp[81]={0}; double dNmp_H[81][7]={{0}}; double dNmp_A[81][7][3]={{{0}}};
    double dIsn[81]={0}; double dIsn_H[81][7]={{0}}; double dIsn_A[81][7][3]={{{0}}};
    double dIsp[81]={0}; double dIsp_H[81][7]={{0}}; double dIsp_A[81][7][3]={{{0}}};
    double dImn[81]={0}; double dImn_H[81][7]={{0}}; double dImn_A[81][7][3]={{{0}}};
    double dImp[81]={0}; double dImp_H[81][7]={{0}}; double dImp_A[81][7][3]={{{0}}};
    double dPTn[81]={0}; double dPTn_H[81][7]={{0}}; double dPTn_A[81][7][3]={{{0}}}; 
    double dPTp[81]={0}; double dPTp_H[81][7]={{0}}; double dPTp_A[81][7][3]={{{0}}}; 

    /* intergers to use as counters */ 
    int i,j,l,ij;
     
    /* Then need to assign the variables to the correct bit of y (we do the same for the rates of change after solving the equations) */
     
    /* HIV- */ 
    
    int n_age = 81;     /* Number of age groups */
    int n_HIV = 7;      /* Number of HIV pos groups */
    int n_ART = 3;      /* Number of ART groups */
    int n_disease = 15; /* Number of disease states */
    
    for (i=0; i<n_age; i++) S[i] = y[i];             
    for (i=n_age; i<n_age*2; i++) Lsn[i-n_age] = y[i];       
    for (i=n_age*2; i<n_age*3; i++) Lsp[i-n_age*2] = y[i];       
    for (i=n_age*3; i<n_age*4; i++) Lmn[i-n_age*3] = y[i];    
    for (i=n_age*4; i<n_age*5; i++) Lmp[i-n_age*4] = y[i];     
    for (i=n_age*5; i<n_age*6; i++) Nsn[i-n_age*5] = y[i];      
    for (i=n_age*6; i<n_age*7; i++) Nsp[i-n_age*6] = y[i];    
    for (i=n_age*7; i<n_age*8; i++) Nmn[i-n_age*7] = y[i];    
    for (i=n_age*8; i<n_age*9; i++) Nmp[i-n_age*8] = y[i];    
    for (i=n_age*9; i<n_age*10; i++) Isn[i-n_age*9] = y[i];    
    for (i=n_age*10; i<n_age*11; i++) Isp[i-n_age*10] = y[i];    
    for (i=n_age*11; i<n_age*12; i++) Imn[i-n_age*11] = y[i];    
    for (i=n_age*12; i<n_age*13; i++) Imp[i-n_age*12] = y[i];  
    for (i=n_age*13; i<n_age*14; i++) PTn[i-n_age*13] = y[i];
    for (i=n_age*14; i<n_age*15; i++) PTp[i-n_age*14] = y[i];
    /* HIV+ */
    ij = n_age*n_disease;  
    for (j=0; j<n_HIV; j++){
      for (i=0; i<n_age; i++){  
        S_H[i][j] = y[ij];
        Lsn_H[i][j] = y[ij+(n_age*n_HIV)];
        Lsp_H[i][j] = y[ij+(2*n_age*n_HIV)];
        Lmn_H[i][j] = y[ij+(3*n_age*n_HIV)];
        Lmp_H[i][j] = y[ij+(4*n_age*n_HIV)];  
        Nsn_H[i][j] = y[ij+(5*n_age*n_HIV)];
        Nsp_H[i][j] = y[ij+(6*n_age*n_HIV)];
        Nmn_H[i][j] = y[ij+(7*n_age*n_HIV)];
        Nmp_H[i][j] = y[ij+(8*n_age*n_HIV)];
        Isn_H[i][j] = y[ij+(9*n_age*n_HIV)];
        Isp_H[i][j] = y[ij+(10*n_age*n_HIV)];
        Imn_H[i][j] = y[ij+(11*n_age*n_HIV)];
        Imp_H[i][j] = y[ij+(12*n_age*n_HIV)];
        PTn_H[i][j] = y[ij+(13*n_age*n_HIV)];
        PTp_H[i][j] = y[ij+(14*n_age*n_HIV)];
        ij = ij+1;
      }
    }
    /* HIV+, on ART */
    ij = (n_HIV+1)*n_age*n_disease;   
    for (l=0; l<n_ART; l++){                          
      for (j=0; j<n_HIV; j++){
        for (i=0; i<n_age; i++){
          S_A[i][j][l] = y[ij];
          Lsn_A[i][j][l] = y[ij+(n_age*n_ART*n_HIV)];
          Lsp_A[i][j][l] = y[ij+(2*n_age*n_ART*n_HIV)];
          Lmn_A[i][j][l] = y[ij+(3*n_age*n_ART*n_HIV)];
          Lmp_A[i][j][l] = y[ij+(4*n_age*n_ART*n_HIV)];  
          Nsn_A[i][j][l] = y[ij+(5*n_age*n_ART*n_HIV)];
          Nsp_A[i][j][l] = y[ij+(6*n_age*n_ART*n_HIV)];
          Nmn_A[i][j][l] = y[ij+(7*n_age*n_ART*n_HIV)];
          Nmp_A[i][j][l] = y[ij+(8*n_age*n_ART*n_HIV)];
          Isn_A[i][j][l] = y[ij+(9*n_age*n_ART*n_HIV)];
          Isp_A[i][j][l] = y[ij+(10*n_age*n_ART*n_HIV)];
          Imn_A[i][j][l] = y[ij+(11*n_age*n_ART*n_HIV)];
          Imp_A[i][j][l] = y[ij+(12*n_age*n_ART*n_HIV)];
          PTn_A[i][j][l] = y[ij+(13*n_age*n_ART*n_HIV)];
          PTp_A[i][j][l] = y[ij+(14*n_age*n_ART*n_HIV)];
          ij = ij+1;
        }
      }
    }
      
    /* Adjust TB model parameters for age, HIV and ART */

    /* Create vectors of disease parameters (by age - now single year bins) to use in derivatives - includes BCG effect on risk of primary disease */
    double a_age[81];
    double sig_age[81];
    double v_age[81];
    double muN_age[81];
    double muI_age[81];
    double bcg = (BCG_cov*(1-BCG_eff)+(1-BCG_cov));  /* those on bcg (BCG_COV) have RR of (1-BCG_eff) */

    for (i=0; i<5; i++){
      a_age[i] = a0*bcg;
      a_age[i+5] = a5*bcg;
      a_age[i+10] = a10*bcg;
      
      sig_age[i] = sig0;
      sig_age[i+5] = sig5;
      sig_age[i+10] = sig10;
      
      v_age[i] = v;
      v_age[i+5] = v;
      v_age[i+10] = v;
      
      muN_age[i] = mu_N0;
      muN_age[i+5] = mu_N;
      muN_age[i+10] = mu_N;
      
      muI_age[i] = mu_I0;
      muI_age[i+5] = mu_I;
      muI_age[i+10] = mu_I;
    }
    for (i=15; i<81; i++){
      a_age[i] = a_a;
      sig_age[i] = sig_a;
      v_age[i] = v; 
      muN_age[i] = mu_N;
      muI_age[i] = mu_I;
    }
    
    /* Now adjust parameters for HIV and ART */
    /* HIV mortality rates are passed in directly and ART mortality reduction is implemented in the derivatives */

    double mid_CD4[7] = {500,425,300,225,150,75,25};       /* mid points of CD4 categories */
    double ART_TB[3] = {ART_TB1,ART_TB2,ART_TB3};          /* vector of ART relative risks for TB disease */
    double ART_mort[3] = {ART_mort1,ART_mort2,ART_mort3};  /* vector of ART relative risks for TB mortlaity */
    /* then adjust parameters */
    double a_age_H[81][7];
    double p_H[7];
    double v_age_H[81][7];
    double a_age_A[81][7][3];
    double p_A[7][3];
    double v_age_A[81][7][3];
    double muN_H_A[3];
    double muI_H_A[3];
    for (j=0; j<n_HIV; j++){
      p_H[j] = p*RR1p*pow(RR2p,-1*(500-mid_CD4[j])/100);
      for (i=0; i<n_age; i++){
        a_age_H[i][j] = fmin(a_age[i]*RR1a*pow(RR2a,(500-mid_CD4[j])/100),1); /* fmin ensures proportion developing active disease can't go above 1 - assume this cap is also applied in TIME */
        v_age_H[i][j] = v_age[i]*RR1v*pow(RR2v,(500-mid_CD4[j])/100);
        for (l=0; l<n_ART; l++){
          a_age_A[i][j][l] = fmax(a_age_H[i][j]*(1-ART_TB[l]),a_age[i]);   /* fmax or fmin ensures being on ART can't be better than being HIV- */
          v_age_A[i][j][l] = fmax(v_age_H[i][j]*(1-ART_TB[l]),v_age[i]);
          p_A[j][l] = fmin(p_H[j]*(1+ART_TB[l]),p);                        /* protection term gets higher as ART is taken */
        
          muN_H_A[l] = fmax(muN_H*(1-ART_mort[l]),muN_age[i]); /* make sure mortality can't go lower than HIV- */ 
          muI_H_A[l] = fmax(muI_H*(1-ART_mort[l]),muI_age[i]);
        
        }
      }
    }
  
    /* Set up parameters for HIV model - these are taken from AIM */

    double H_CD4[7][81] = {{0},{0},{0},{0},{0},{0},{0}}; /* Distribution of new HIV infections (age,CD4) - assume distribution of new child infections mirrors adults */
    for (i=0; i<25; i++) {
      H_CD4[0][i] = 0.643; H_CD4[1][i] = 0.357;
    }
    for (i=25; i<35; i++){
      H_CD4[0][i] = 0.607; H_CD4[1][i] = 0.393;
    }
    for (i=35; i<45; i++){
      H_CD4[0][i] = 0.585; H_CD4[1][i] = 0.415;
    }
    for (i=45; i<81; i++){
      H_CD4[0][i] = 0.552; H_CD4[1][i] = 0.448;
    }

    /* Have updated these values based on the durations in the AIM manual (rate = 1/duration) as they are different from rates in AIM editor in software */
    /* those are actually risks (i.e. 1-exp(-rate)) */
    double H_prog[8][81] = {{0},{0},{0},{0},{0},{0},{0},{0}};  /* Progression through CD4 categories (age, CD4) - has extra row to avoid progression in/out of first/last groups */
    for (i=0; i<15; i++){
      H_prog[1][i] = 0.298; H_prog[2][i] = 0.239; H_prog[3][i] = 0.183; H_prog[4][i] = 0.183; H_prog[5][i] = 0.130; H_prog[6][i] = 0.130;
    }
    for (i=15; i<25; i++){
      H_prog[1][i] = 0.117; H_prog[2][i] = 0.223; H_prog[3][i] = 0.294; H_prog[4][i] = 0.508; H_prog[5][i] = 0.214; H_prog[6][i] = 0.348;
    }
    for (i=25; i<35; i++){
      H_prog[1][i] = 0.147; H_prog[2][i] = 0.240; H_prog[3][i] = 0.452; H_prog[4][i] = 1.087; H_prog[5][i] = 0.637; H_prog[6][i] = 1.449;
    }
    for (i=35; i<45; i++){
      H_prog[1][i] = 0.183; H_prog[2][i] = 0.355; H_prog[3][i] = 0.581; H_prog[4][i] = 1.250; H_prog[5][i] = 0.676; H_prog[6][i] = 1.449;
    }
    for (i=45; i<81; i++){
      H_prog[1][i] = 0.213; H_prog[2][i] = 0.535; H_prog[3][i] = 0.855; H_prog[4][i] = 1.818; H_prog[5][i] = 0.952; H_prog[6][i] = 2.000;
    }

    double H_mort[7][81]; /* Mortality due to HIV (no ART) (age, CD4) */
    for (i=0; i<5; i++){
      H_mort[0][i] = 0.312; H_mort[1][i] = 0.382; H_mort[2][i] = 0.466; H_mort[3][i] = 0.466; H_mort[4][i] = 0.569; H_mort[5][i] = 0.569; H_mort[6][i] = 0.569;
    }
    for (i=5; i<15; i++){
      H_mort[0][i] = 0.039; H_mort[1][i] = 0.048; H_mort[2][i] = 0.058; H_mort[3][i] = 0.058; H_mort[4][i] = 0.071; H_mort[5][i] = 0.071; H_mort[6][i] = 0.071;
    }
    for (i=15; i<25; i++){
      H_mort[0][i] = 0.005; H_mort[1][i] = 0.011; H_mort[2][i] = 0.026; H_mort[3][i] = 0.061; H_mort[4][i] = 0.139; H_mort[5][i] = 0.321; H_mort[6][i] = 0.737;
    }    
    for (i=25; i<35; i++){
      H_mort[0][i] = 0.004; H_mort[1][i] = 0.010; H_mort[2][i] = 0.026; H_mort[3][i] = 0.069; H_mort[4][i] = 0.185; H_mort[5][i] = 0.499; H_mort[6][i] = 1.342;
    } 
    for (i=35; i<45; i++){
      H_mort[0][i] = 0.005; H_mort[1][i] = 0.013; H_mort[2][i] = 0.036; H_mort[3][i] = 0.096; H_mort[4][i] = 0.258; H_mort[5][i] = 0.691; H_mort[6][i] = 1.851;
    } 
    for (i=45; i<81; i++){
      H_mort[0][i] = 0.005; H_mort[1][i] = 0.013; H_mort[2][i] = 0.032; H_mort[3][i] = 0.080; H_mort[4][i] = 0.203; H_mort[5][i] = 0.513; H_mort[6][i] = 1.295;
    } 

    double A_mort[3][7][81]; /* On ART mortality (age,starting CD4, time on ART)  - this is an average of male and female values weigthed by sex of those on ART  */
    
    for (i=0; i<5; i++){
      A_mort[0][0][i] = 0.062; A_mort[0][1][i] = 0.282;  A_mort[0][2][i] = 0.214; A_mort[0][3][i] = 0.214; A_mort[0][4][i] = 0.749;  A_mort[0][5][i] = 0.749;  A_mort[0][6][i] = 0.749; 
      A_mort[1][0][i] = 0.147; A_mort[1][1][i] = 0.210;  A_mort[1][2][i] = 0.199; A_mort[1][3][i] = 0.199; A_mort[1][4][i] = 0.376;  A_mort[1][5][i] = 0.376;  A_mort[1][6][i] = 0.376; 
      A_mort[2][0][i] = 0.062; A_mort[2][1][i] = 0.089;  A_mort[2][2][i] = 0.084; A_mort[2][3][i] = 0.084; A_mort[2][4][i] = 0.159;  A_mort[2][5][i] = 0.159;  A_mort[2][6][i] = 0.159; 
    }
    for (i=5; i<10; i++){
      A_mort[0][0][i] = 0.008; A_mort[0][1][i] = 0.035;  A_mort[0][2][i] = 0.027; A_mort[0][3][i] = 0.027; A_mort[0][4][i] = 0.094;  A_mort[0][5][i] = 0.094;  A_mort[0][6][i] = 0.094; 
      A_mort[1][0][i] = 0.018; A_mort[1][1][i] = 0.026;  A_mort[1][2][i] = 0.025; A_mort[1][3][i] = 0.025; A_mort[1][4][i] = 0.047;  A_mort[1][5][i] = 0.047;  A_mort[1][6][i] = 0.047; 
      A_mort[2][0][i] = 0.008; A_mort[2][1][i] = 0.011;  A_mort[2][2][i] = 0.011; A_mort[2][3][i] = 0.011; A_mort[2][4][i] = 0.020;  A_mort[2][5][i] = 0.020;  A_mort[2][6][i] = 0.020; 
    }
    for (i=10; i<15; i++){
      A_mort[0][0][i] = 0.007; A_mort[0][1][i] = 0.033;  A_mort[0][2][i] = 0.025; A_mort[0][3][i] = 0.025; A_mort[0][4][i] = 0.088;  A_mort[0][5][i] = 0.088;  A_mort[0][6][i] = 0.088; 
      A_mort[1][0][i] = 0.009; A_mort[1][1][i] = 0.012;  A_mort[1][2][i] = 0.012; A_mort[1][3][i] = 0.012; A_mort[1][4][i] = 0.022;  A_mort[1][5][i] = 0.022;  A_mort[1][6][i] = 0.022; 
      A_mort[2][0][i] = 0.004; A_mort[2][1][i] = 0.005;  A_mort[2][2][i] = 0.005; A_mort[2][3][i] = 0.005; A_mort[2][4][i] = 0.009;  A_mort[2][5][i] = 0.009;  A_mort[2][6][i] = 0.009; 
    }
    for (i=15; i<25; i++){
      A_mort[0][0][i] = 0.0050; A_mort[0][1][i] = 0.0115;  A_mort[0][2][i] = 0.0264; A_mort[0][3][i] = 0.0557; A_mort[0][4][i] = 0.0953;  A_mort[0][5][i] = 0.1553;  A_mort[0][6][i] = 0.3418; 
      A_mort[1][0][i] = 0.0050; A_mort[1][1][i] = 0.0115;  A_mort[1][2][i] = 0.0204; A_mort[1][3][i] = 0.0216; A_mort[1][4][i] = 0.0272;  A_mort[1][5][i] = 0.0343;  A_mort[1][6][i] = 0.0496; 
      A_mort[2][0][i] = 0.0050; A_mort[2][1][i] = 0.0068;  A_mort[2][2][i] = 0.0073; A_mort[2][3][i] = 0.0079; A_mort[2][4][i] = 0.0105;  A_mort[2][5][i] = 0.0139;  A_mort[2][6][i] = 0.0211; 
    }
    for (i=25; i<35; i++){
      A_mort[0][0][i] = 0.0035; A_mort[0][1][i] = 0.0095;  A_mort[0][2][i] = 0.0256; A_mort[0][3][i] = 0.0477; A_mort[0][4][i] = 0.0792;  A_mort[0][5][i] = 0.1305;  A_mort[0][6][i] = 0.2897; 
      A_mort[1][0][i] = 0.0035; A_mort[1][1][i] = 0.0095;  A_mort[1][2][i] = 0.0242; A_mort[1][3][i] = 0.0278; A_mort[1][4][i] = 0.0351;  A_mort[1][5][i] = 0.0442;  A_mort[1][6][i] = 0.0640; 
      A_mort[2][0][i] = 0.0035; A_mort[2][1][i] = 0.0081;  A_mort[2][2][i] = 0.0093; A_mort[2][3][i] = 0.0100; A_mort[2][4][i] = 0.0134;  A_mort[2][5][i] = 0.0177;  A_mort[2][6][i] = 0.0271; 
    }
    for (i=35; i<45; i++){
      A_mort[0][0][i] = 0.0050; A_mort[0][1][i] = 0.0134;  A_mort[0][2][i] = 0.0359; A_mort[0][3][i] = 0.0495; A_mort[0][4][i] = 0.0833;  A_mort[0][5][i] = 0.1384;  A_mort[0][6][i] = 0.3094; 
      A_mort[1][0][i] = 0.0050; A_mort[1][1][i] = 0.0134;  A_mort[1][2][i] = 0.0267; A_mort[1][3][i] = 0.0283; A_mort[1][4][i] = 0.0362;  A_mort[1][5][i] = 0.0461;  A_mort[1][6][i] = 0.0675; 
      A_mort[2][0][i] = 0.0050; A_mort[2][1][i] = 0.0076;  A_mort[2][2][i] = 0.0083; A_mort[2][3][i] = 0.0091; A_mort[2][4][i] = 0.0128;  A_mort[2][5][i] = 0.0174;  A_mort[2][6][i] = 0.0275; 
    }
    for (i=45; i<81; i++){
      A_mort[0][0][i] = 0.0050; A_mort[0][1][i] = 0.0126;  A_mort[0][2][i] = 0.0319; A_mort[0][3][i] = 0.0489; A_mort[0][4][i] = 0.0872;  A_mort[0][5][i] = 0.1496;  A_mort[0][6][i] = 0.3431; 
      A_mort[1][0][i] = 0.0050; A_mort[1][1][i] = 0.0126;  A_mort[1][2][i] = 0.0306; A_mort[1][3][i] = 0.0366; A_mort[1][4][i] = 0.0481;  A_mort[1][5][i] = 0.0625;  A_mort[1][6][i] = 0.0935; 
      A_mort[2][0][i] = 0.0045; A_mort[2][1][i] = 0.0066;  A_mort[2][2][i] = 0.0076; A_mort[2][3][i] = 0.0087; A_mort[2][4][i] = 0.0141;  A_mort[2][5][i] = 0.0209;  A_mort[2][6][i] = 0.0355; 
    }
    
    double A_prog[4] = {0,2,2,0}; /* Progression through time on ART, 6 monthly time blocks - 0 ensure no progression into first catergory and no progression out of last category*/
    double A_start[3] = {1,0,0};  /* Used to make sure ART initiations are only added to the fist time on ART box */ 
    
    /* sum up various totals */

    /* Use sumsum function to add up HIV- */
    double Total_S = sumsum(S,0,80);                        /* Total susceptible */
    double Total_Ls = sumsum(Lsn,0,80)+sumsum(Lsp,0,80);    /* Total LTBI with drug susceptible (DS) strain */
    double Total_Lm = sumsum(Lmn,0,80)+sumsum(Lmp,0,80);    /* Total LTBI with drug resistant (DR) strain */
    double Total_Ns_N = sumsum(Nsn,0,80)+sumsum(Nsp,0,80);    /* Total DS smear negative TB */
    double Total_Nm_N = sumsum(Nmn,0,80)+sumsum(Nmp,0,80);    /* Total DR smear negative TB */
    double Total_Is_N = sumsum(Isn,0,80)+sumsum(Isp,0,80);    /* Total DS smear positive TB */
    double Total_Im_N = sumsum(Imn,0,80)+sumsum(Imp,0,80);    /* Total DR smear positive TB */
    double Total_PT = sumsum(PTn,0,80)+sumsum(PTp,0,80);      /* Post PT */
    
    /* Now loop through HIV and ART and add them in */
    double Total_Ns_H =0; double Total_Nm_H = 0; double Total_Is_H = 0; double Total_Im_H =0;
    
    for (j=0; j<n_HIV; j++){
      for (i=0; i<n_age; i++){
        Total_S = Total_S + S_H[i][j];
        Total_Ls = Total_Ls + Lsn_H[i][j]+Lsp_H[i][j];
        Total_Lm = Total_Lm + Lmn_H[i][j]+Lmp_H[i][j];
        Total_Ns_H = Total_Ns_H + Nsn_H[i][j]+Nsp_H[i][j];
        Total_Nm_H = Total_Nm_H + Nmn_H[i][j]+Nmp_H[i][j];
        Total_Is_H = Total_Is_H + Isn_H[i][j]+Isp_H[i][j];
        Total_Im_H = Total_Im_H + Imn_H[i][j]+Imp_H[i][j];
        Total_PT = Total_PT + PTn_H[i][j] + PTp_H[i][j];
        for (l=0; l<n_ART; l++){
          Total_S = Total_S + S_A[i][j][l];
          Total_Ls = Total_Ls + Lsn_A[i][j][l]+Lsp_A[i][j][l];
          Total_Lm = Total_Lm + Lmn_A[i][j][l]+Lmp_A[i][j][l];
          Total_Ns_H = Total_Ns_H + Nsn_A[i][j][l]+Nsp_A[i][j][l];
          Total_Nm_H = Total_Nm_H + Nmn_A[i][j][l]+Nmp_A[i][j][l];
          Total_Is_H = Total_Is_H + Isn_A[i][j][l]+Isp_A[i][j][l];
          Total_Im_H = Total_Im_H + Imn_A[i][j][l]+Imp_A[i][j][l];
          Total_PT = Total_PT + PTn_A[i][j][l] + PTp_A[i][j][l];
        }
      }
    }
    double Total_L = Total_Ls + Total_Lm;           /* Total LTBI */
    double Total_N_N = Total_Ns_N + Total_Nm_N;     /* Total smear negative TB (HIV-) */
    double Total_N_H = Total_Ns_H + Total_Nm_H;     /* Total smear negative TB (HIV+) */
    double Total_N = Total_N_N + Total_N_H;         /* Total smear negative TB */
    double Total_I_N = Total_Is_N + Total_Im_N;     /* Total smear positive TB (HIV-) */
    double Total_I_H = Total_Is_H + Total_Im_H;     /* Total smear positive TB (HIV+) */
    double Total_I = Total_I_N + Total_I_H;         /* Total smear positive TB */
    
    double Total_DS = Total_Ns_N + Total_Ns_H + Total_Is_N + Total_Is_H;    /* Total DS TB */
    double Total_MDR = Total_Nm_N + Total_Nm_H + Total_Im_N + Total_Im_H;   /* Total DR TB */
    double Total = Total_S+Total_L+Total_N+Total_I+Total_PT; /* Total */
    
    /* Mortality calculations and adjustments */
    /* HIV mortality rates include TB deaths */
    /* Background mortality rates include HIV and TB deaths */
    /* Need to make an adjustment to both HIV and background rates to avoid double counting */
    
    /* Calculate deaths due to disease (TB and HIV) by age, CD4 and ART */
    /* work out disease induced mortality rate if pre 2015 - if after 2015 just use 2015 value */
    /* and adjust background mortality rate accordingly */
    /* Calculate total population in the same loop and prevalence of TB in HIV- */
    double TB_deaths_neg[81];
    double TB_deaths_HIV[81][7];
    double TB_deaths_ART[81][7][3];
    double TB_deaths_HIV_age[81] = {0};
    double TB_deaths_ART_age[81] = {0};
    double TB_deaths[81];
    
    double HIV_deaths_HIV[81] = {0}; ;
    double HIV_deaths_ART[81] = {0};

    double up_H_mort[7][81];
    double up_A_mort[3][7][81];
    double m_b[81];
    double rate_dis_death[81];
    
    double tot_age[81] = {0};
    double tot_age_HIV[81][7];
    double tot_age_ART[81][7][3];
    
    /*double Tot_deaths = 0;*/
    double Tot_deaths_age[81];
    double ART_deaths = 0;
    double Tot_deaths=0;
    double tot_age_neg[81] = {0};
    
 
    for (i=0; i<n_age; i++) {
      
      /* Calculate HIV- TB deaths */
      TB_deaths_neg[i] = (Nsn[i]+Nsp[i]+Nmn[i]+Nmp[i])*muN_age[i] + (Isn[i]+Isp[i]+Imn[i]+Imp[i])*muI_age[i];
      TB_deaths[i] = TB_deaths_neg[i];
      /* Calculate size of age group */
      tot_age_neg[i] = S[i]+Lsn[i]+Lsp[i]+Lmn[i]+Lmp[i]+Nsn[i]+Nsp[i]+Nmn[i]+Nmp[i]+Isn[i]+Isp[i]+Imn[i]+Imp[i]+PTn[i]+PTp[i];
      tot_age[i] = tot_age[i] + tot_age_neg[i];
      
      for(j=0; j<n_HIV; j++){
        
        /* Calculate HIV+ TB deaths */
        TB_deaths_HIV[i][j] = (Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j])*muN_H + (Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j])*muI_H;
        TB_deaths_HIV_age[i] = TB_deaths_HIV_age[i] + TB_deaths_HIV[i][j];
        TB_deaths[i] = TB_deaths[i] + TB_deaths_HIV[i][j];
        /* Calculate size of HIV+ age group */                              
        tot_age_HIV[i][j] = S_H[i][j]+Lsn_H[i][j]+Lsp_H[i][j]+Lmn_H[i][j]+Lmp_H[i][j]+Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j]+Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j]+PTn_H[i][j]+PTp_H[i][j]; 
        /* Update size of age group */
        tot_age[i] = tot_age[i] + tot_age_HIV[i][j];
        /* Adjust HIV mortality probability to remove TB deaths (only if there is any HIV yet (otherwise we get a divide by 0 error)) */
        if (tot_age_HIV[i][j]>0.0){
          up_H_mort[j][i] = fmax(0,(H_mort[j][i] - (TB_deaths_HIV[i][j]/tot_age_HIV[i][j])));
        }
        else{
          up_H_mort[j][i] = H_mort[j][i];
        }
        
        /* Calculate HIV deaths using the updated probabilities */  
        HIV_deaths_HIV[i] = HIV_deaths_HIV[i] + up_H_mort[j][i]*(S_H[i][j]+Lsn_H[i][j]+Lsp_H[i][j]+Lmn_H[i][j]+Lmp_H[i][j]+
                                                                 Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j]+
                                                                 Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j]+
                                                                 PTn_H[i][j]+PTp_H[i][j]);
      
        for (l=0; l<n_ART; l++){
          
          /* Calculate ART TB deaths */
          TB_deaths_ART[i][j][l] = (Nsn_A[i][j][l]+Nsp_A[i][j][l]+Nmn_A[i][j][l]+Nmp_A[i][j][l])*muN_H_A[l] + 
                                   (Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l])*muI_H_A[l]; 
          TB_deaths_ART_age[i] = TB_deaths_ART_age[i] + TB_deaths_ART[i][j][l];
          TB_deaths[i] = TB_deaths[i] + TB_deaths_ART[i][j][l];
          /* Calculate size of ART age group */
          tot_age_ART[i][j][l] = S_A[i][j][l]+Lsn_A[i][j][l]+Lsp_A[i][j][l]+Lmn_A[i][j][l]+Lmp_A[i][j][l]+Nsn_A[i][j][l]+Nsp_A[i][j][l]+
                                 Nmn_A[i][j][l]+Nmp_A[i][j][l]+Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l]+PTn_A[i][j][l]+PTp_A[i][j][l];
          /* Update size of age group */
          tot_age[i] = tot_age[i] + tot_age_ART[i][j][l];
          /* Adjust ART mortality probability to remove TB deaths (only if there is any HIV yet (otherwise we get a divide by 0 error)) */
          if (tot_age_ART[i][j][l]>0.0){
            up_A_mort[l][j][i] = fmax(0,A_mort[l][j][i] - (TB_deaths_ART[i][j][l]/tot_age_ART[i][j][l]));
          }
          else{
            up_A_mort[l][j][i] = A_mort[l][j][i];
          }
          
          /* Calculate ART deaths using the updated probabilites */
          HIV_deaths_ART[i] = HIV_deaths_ART[i] + up_A_mort[l][j][i]*(S_A[i][j][l]+Lsn_A[i][j][l]+Lsp_A[i][j][l]+Lmn_A[i][j][l]+Lmp_A[i][j][l]+
                                                                      Nsn_A[i][j][l]+Nsp_A[i][j][l]+Nmn_A[i][j][l]+Nmp_A[i][j][l]+ 
                                                                      Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l]+
                                                                      PTn_A[i][j][l]+PTp_A[i][j][l]);    
                                                                      
          /* Add up all deaths on ART - used to put new people on ART */
          ART_deaths = ART_deaths + TB_deaths_ART[i][j][l] + up_A_mort[l][j][i]*(S_A[i][j][l]+Lsn_A[i][j][l]+Lsp_A[i][j][l]+Lmn_A[i][j][l]+Lmp_A[i][j][l]+
                                                                                 Nsn_A[i][j][l]+Nsp_A[i][j][l]+Nmn_A[i][j][l]+Nmp_A[i][j][l]+ 
                                                                                 Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l]+
                                                                                 PTn_A[i][j][l]+PTp_A[i][j][l]);
          
        }                               
      }

      if (pop_ad>0) rate_dis_death[i] = (TB_deaths_neg[i]+TB_deaths_HIV_age[i]+TB_deaths_ART_age[i]+HIV_deaths_HIV[i]+HIV_deaths_ART[i])/tot_age[i];
      m_b[i] = fmax(0,forc[i+1]-rate_dis_death[i]);
      
      Tot_deaths_age[i] = m_b[i]*tot_age[i] + TB_deaths[i] + HIV_deaths_HIV[i] + HIV_deaths_ART[i];  
      Tot_deaths = Tot_deaths + Tot_deaths_age[i];
      
      /* Add in background deaths on ART - used to calculate new people to put on ART */
      for(j=0; j<n_HIV; j++){
        for (l=0; l<n_ART; l++){
          ART_deaths = ART_deaths + m_b[i]*(S_A[i][j][l]+Lsn_A[i][j][l]+Lsp_A[i][j][l]+Lmn_A[i][j][l]+Lmp_A[i][j][l]+Nsn_A[i][j][l]+Nsp_A[i][j][l]+
                                               Nmn_A[i][j][l]+Nmp_A[i][j][l]+Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l]+PTn_A[i][j][l]+PTp_A[i][j][l]);
        }
      }
      
    } 
    double TB_deaths_tot = sumsum(TB_deaths,0,80);
    
    /* Sum up populations over CD4 categories, with and without ART and calculate rates of ART initiation */
    double CD4_dist[7] = {0,0,0,0,0,0,0};     /* Not on ART by CD4 */
    double CD4_dist_ART[7] = {0,0,0,0,0,0,0}; /* On ART by CD4 */
    double CD4_deaths[7] = {0,0,0,0,0,0,0};   /* Deaths by CD4 (no ART) */
    double Tot_ART = 0;                       /* Total on ART */ 
    double ART_need = 0;                      /* Number needing ART - this is all those below the eligibility threshold*/
    double ART_on = 0;                        /* Number who should be on ART based on reported % */
    double ART_new = 0;                       /* Number who need to start ART */
    double ART_el = 0;                        /* Number who are eligible but not on ART */
    double ART_el_deaths = 0;                 /* Number eligible who will die */
    double ART_prop[7] = {0,0,0,0,0,0,0};     /* Proportion of CD4 category who should start ART */

    for (j=0; j<n_HIV; j++){
      for (i=0; i<n_age; i++){
        CD4_dist[j] = CD4_dist[j]+S_H[i][j]+Lsn_H[i][j]+Lsp_H[i][j]+Lmn_H[i][j]+Lmp_H[i][j]+
                      Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j]+ 
                      Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j]+
                      PTn_H[i][j]+PTp_H[i][j];
        CD4_deaths[j] = CD4_deaths[j]+
                        up_H_mort[j][i]*(S_H[i][j]+Lsn_H[i][j]+Lsp_H[i][j]+Lmn_H[i][j]+Lmp_H[i][j]+
                        Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j]+ 
                        Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j]+
                        PTn_H[i][j]+PTp_H[i][j]);
        for (l=0; l<n_ART; l++){
          CD4_dist_ART[j] = CD4_dist_ART[j]+S_A[i][j][l]+Lsn_A[i][j][l]+Lsp_A[i][j][l]+Lmn_A[i][j][l]+Lmp_A[i][j][l]+
                            Nsn_A[i][j][l]+Nsp_A[i][j][l]+Nmn_A[i][j][l]+Nmp_A[i][j][l]+ 
                            Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l]+
                            PTn_A[i][j][l]+PTp_A[i][j][l];         
                      
        }            
      } 
      Tot_ART = Tot_ART + CD4_dist_ART[j];  /* sum up number currently on ART */
      ART_on = ART_on + (CD4_dist[j] + CD4_dist_ART[j])*forc[j+163]; /* number who should be on ART - HIV+ population times coverage (by CD4) */
    }
    ART_new = fmax(0,ART_on - (Tot_ART - ART_deaths));   /* number who need to start is number who should be on minus those already on plus those on ART who will die in current time */ 
        
    /* Then work out where these should go by CD4 - based on proportion of eligible population in CD4 group and proportion of deaths which occuring in CD4 group */
    for (j=Athresh; j<n_HIV; j++) {
      ART_el = ART_el + CD4_dist[j];
      ART_el_deaths = ART_el_deaths + CD4_deaths[j];
      ART_need = ART_need + CD4_dist[j] + CD4_dist_ART[j];
    }
    if (ART_el > 0){
      for (j=Athresh; j<n_HIV; j++) {
        if (CD4_dist[j] > 0) {
          ART_prop[j] = (((CD4_dist[j]/ART_el)+(CD4_deaths[j]/ART_el_deaths))/2)*(ART_new/CD4_dist[j]); /* applies weighting and size of CD4 group to work out % of CD4 group that should move */
           
        }
      }
    }
    
    /* Force of infection */
    double FS = beta*(Total_Ns_N*rel_inf + Total_Ns_H*rel_inf_H + Total_Is_N + Total_Is_H)/Total; 
    double FM = fit_cost*beta*(Total_Nm_N*rel_inf + Total_Nm_H*rel_inf_H + Total_Im_N + Total_Im_H)/Total; 
    
    /* Variables to store numbers of new cases */
    double TB_cases_age[81] = {0};
    double TB_cases_neg_age[81];
    double TB_cases_neg = 0;
    double TB_cases_pos_age[81][7];
    double TB_cases_pos = 0;
    double TB_cases_ART_age[81][7][3];
    double TB_cases_ART = 0;
        
    /* Derivatives */ 
 
    /* HIV-: loop through ages*/ 
    
    double births = birth_rate*Total/1000;

    for (i=0; i<n_age; i++){
      
      /* Calculate the disease flows here and use these in the derivatives - intention is to make the model more flexible/easier to understand */
      
      double S_to_Lsn = FS*(1-a_age[i])*S[i];                                     /* Susceptible to latent DS infection (no disease history) */
      double S_to_Nsn = FS*a_age[i]*(1-sig_age[i])*S[i];                          /* Susceptible to primary DS smear negative disease (no disease history) */
      double S_to_Isn = FS*a_age[i]*sig_age[i]*S[i];                              /* Susceptible to primary DS smear positive disease (no disease history) */
      double S_to_Lmn = FM*(1-a_age[i])*S[i];                                     /* Susceptible to latent DR infection (no disease history) */
      double S_to_Nmn = FM*a_age[i]*(1-sig_age[i])*S[i];                          /* Susceptible to primary DR smear negative disease (no disease history) */
      double S_to_Imn = FM*a_age[i]*sig_age[i]*S[i];                              /* Susceptible to primary DR smear positive disease (no disease history) */
      
      double Lsn_to_Nsn = (v_age[i] + FS*a_age[i]*(1-p))*(1-sig_age[i])*Lsn[i];   /* Latent DS to smear negative DS disease (no disease history) - reactivation and reinfection */ 
      double Lsn_to_Isn = (v_age[i] + FS*a_age[i]*(1-p))*sig_age[i]*Lsn[i];       /* Latent DS to smear positive DS disease (no disease history) - reactivation and reinfection */ 
      double Lsn_to_Nmn = FM*a_age[i]*(1-p)*(1-sig_age[i])*Lsn[i];                /* Latent DS to smear negative DR disease (no disease history) - co-infection */ 
      double Lsn_to_Imn = FM*a_age[i]*(1-p)*sig_age[i]*Lsn[i];                    /* Latent DS to smear positive DR disease (no disease history) - co-infection */
      double Lsn_to_Lmn = FM*(1-a_age[i])*(1-p)*g*Lsn[i];                         /* Latent DS to latent DR (no disease history) */
      
      double Lmn_to_Nmn = (v_age[i] + FM*a_age[i]*(1-p))*(1-sig_age[i])*Lmn[i];   /* Latent DR to smear negative DR disease (no disease history) - reactivation and reinfection */ 
      double Lmn_to_Imn = (v_age[i] + FM*a_age[i]*(1-p))*sig_age[i]*Lmn[i];       /* Latent DR to smear positive DR disease (no disease history) - reactivation and reinfection */ 
      double Lmn_to_Nsn = FS*a_age[i]*(1-sig_age[i])*(1-p)*Lmn[i];                /* Latent DR to smear negative DS disease (no disease history) - co-infection */
      double Lmn_to_Isn = FS*a_age[i]*sig_age[i]*(1-p)*Lmn[i];                    /* Latent DR to smear positive DS disease (no disease history) - co-infection */
      double Lmn_to_Lsn = FS*(1-a_age[i])*(1-p)*(1-g)*Lmn[i];                     /* Latent DR to latent DS (no disease history) */

      double Lsp_to_Nsp = (v_age[i] + FS*a_age[i]*(1-p))*(1-sig_age[i])*Lsp[i];   /* Latent DS to smear negative DS disease (prior Rx) - reactivation and reinfection */ 
      double Lsp_to_Isp = (v_age[i] + FS*a_age[i]*(1-p))*sig_age[i]*Lsp[i];       /* Latent DS to smear positive DS disease (prior Rx) - reactivation and reinfection */     
      double Lsp_to_Nmp = FM*a_age[i]*(1-p)*(1-sig_age[i])*Lsp[i];                /* Latent DS to smear negative DR disease (prior_rx) - reactivation and reinfection */ 
      double Lsp_to_Imp = FM*a_age[i]*(1-p)*sig_age[i]*Lsp[i];                    /* Latent DS to smear positive DR disease (prior_Rx) - reactivation and reinfection */
      double Lsp_to_Lmp = FM*(1-a_age[i])*(1-p)*g*Lsp[i];                         /* Latent DS to latent DR (prior Rx) */
      
      double Lmp_to_Nmp = (v_age[i] + FM*a_age[i]*(1-p))*(1-sig_age[i])*Lmp[i];   /* Latent DR to smear negative DR disease (prior Rx) - reactivation and reinfection */ 
      double Lmp_to_Imp = (v_age[i] + FM*a_age[i]*(1-p))*sig_age[i]*Lmp[i];       /* Latent DR to smear positive DR disease (prior Rx) - reactivation and reinfection */     
      double Lmp_to_Nsp = FS*a_age[i]*(1-p)*(1-sig_age[i])*Lmp[i];                /* Latent DR to smear negative DS disease (prior_rx) - reactivation and reinfection */ 
      double Lmp_to_Isp = FS*a_age[i]*(1-p)*sig_age[i]*Lmp[i];                    /* Latent DR to smear positive DS disease (prior_Rx) - reactivation and reinfection */
      double Lmp_to_Lsp = FS*(1-a_age[i])*(1-p)*(1-g)*Lmp[i];                     /* Latent DR to latent DS (prior Rx) */
      
      double PTn_to_Lsn = FS*(1-a_age[i])*(1-p)*PTn[i];                           /* Post PT to latent DS (no disease history) */
      double PTn_to_Nsn = FS*a_age[i]*(1-p)*(1-sig_age[i])*PTn[i];                /* Post PT to smear negative DS disease (no disease history) */
      double PTn_to_Isn = FS*a_age[i]*(1-p)*sig_age[i]*PTn[i];                    /* Post PT to smear positive DS disease (no disease history) */
      double PTn_to_Lmn = FM*(1-a_age[i])*(1-p)*g*PTn[i];                         /* Post PT to latent DR (no disease history) */
      double PTn_to_Nmn = FM*a_age[i]*(1-p)*(1-sig_age[i])*PTn[i];                /* Post PT to smear negative DR disease (no disease history) */
      double PTn_to_Imn = FM*a_age[i]*(1-p)*sig_age[i]*PTn[i];                    /* Post PT to smear positive DR disease (no disease history) */
      
      double PTp_to_Lsp = FS*(1-a_age[i])*(1-p)*PTp[i];                           /* Post PT to latent DS (prior Rx) */
      double PTp_to_Nsp = FS*a_age[i]*(1-p)*(1-sig_age[i])*PTp[i];                /* Post PT to smear negative DS disease (prior Rx) */
      double PTp_to_Isp = FS*a_age[i]*(1-p)*sig_age[i]*PTp[i];                    /* Post PT to smear positive DS disease (prior Rx) */
      double PTp_to_Lmp = FM*(1-a_age[i])*(1-p)*g*PTp[i];                         /* Post PT to latent DR (prior Rx) */
      double PTp_to_Nmp = FM*a_age[i]*(1-p)*(1-sig_age[i])*PTp[i];                /* Post PT to smear negative DR disease (prior Rx) */
      double PTp_to_Imp = FM*a_age[i]*(1-p)*sig_age[i]*PTp[i];                    /* Post PT to smear positive DR disease (prior Rx) */
          
      /* Susceptible - NOTE BIRTHS ARE ADDED TO HERE IN THE EVENTS FUNCTION*/
      dS[i] = - (FS + FM)*S[i] - forc[i+82]*S[i] - m_b[i]*S[i] + (S[i]/tot_age[i])*(forc[i+188]/5);
       
      /* Latent, ds, naive */
      dLsn[i] = - m_b[i]*Lsn[i] + S_to_Lsn + Lmn_to_Lsn + r*(Isn[i] + Nsn[i]) - Lsn_to_Lmn - Lsn_to_Nsn - Lsn_to_Isn  - Lsn_to_Nmn - Lsn_to_Imn -
                forc[i+82]*Lsn[i] + (Lsn[i]/tot_age[i])*(forc[i+188]/5) + PTn_to_Lsn -
                health*kneg*(1-sp_I_neg*sp_N_neg)*l_s*tneg_s*Lsn[i];                   
      
      /* Latent, ds, prev */
      dLsp[i] = - m_b[i]*Lsp[i] + Lmp_to_Lsp + r*(Isp[i] + Nsp[i]) - Lsp_to_Lmp - Lsp_to_Nsp - Lsp_to_Isp  - Lsp_to_Nmp - Lsp_to_Imp +     
                kneg*(l_s*(1-e)*tneg_s*((dstneg_p*sp_m_neg)+(1-dstneg_p)) + l_m*tneg_m*dstneg_p*(1-sp_m_neg))*(se_I_neg*Isp[i] + se_N_neg*rel_d*Nsp[i]) + 
                kneg*(l_s*(1-e)*tneg_s*((dstneg_n*sp_m_neg)+(1-dstneg_n)) + l_m*tneg_m*dstneg_n*(1-sp_m_neg))*(se_I_neg*Isn[i] + se_N_neg*rel_d*Nsn[i]) - /* Care */
                forc[i+82]*Lsp[i] + (Lsp[i]/tot_age[i])*(forc[i+188]/5) + PTp_to_Lsp -
                health*kneg*(1-sp_I_neg*sp_N_neg)*l_s*tneg_s*Lsp[i];         

      /* Latent, mdr, naive */ 
      dLmn[i] = - m_b[i]*Lmn[i] + S_to_Lmn + Lsn_to_Lmn + r*(Imn[i] + Nmn[i]) - Lmn_to_Lsn - Lmn_to_Nsn - Lmn_to_Isn - Lmn_to_Nmn - Lmn_to_Imn  -                
                forc[i+82]*Lmn[i] + (Lmn[i]/tot_age[i])*(forc[i+188]/5) + PTn_to_Lmn;              
      
      /* Latent, mdr, prev */
      dLmp[i] = - m_b[i]*Lmp[i] + Lsp_to_Lmp + r*(Imp[i] + Nmp[i]) - Lmp_to_Lsp - Lmp_to_Nsp - Lmp_to_Isp - Lmp_to_Nmp - Lmp_to_Imp +
                kneg*(dstneg_p*se_m_neg*l_m*tneg_m + l_s*tneg_s*eff_p*((1-dstneg_p)+dstneg_p*(1-se_m_neg)))*(se_I_neg*Imp[i] + se_N_neg*rel_d*Nmp[i]) +
                kneg*(dstneg_n*se_m_neg*l_m*tneg_m + l_s*tneg_s*eff_n*((1-dstneg_n)+dstneg_n*(1-se_m_neg)))*(se_I_neg*Imn[i] + se_N_neg*rel_d*Nmn[i]) - /* Care */
                forc[i+82]*Lmp[i] + (Lmp[i]/tot_age[i])*(forc[i+188]/5) + PTp_to_Lmp;  
                
      /* Smear neg, ds, new */
      dNsn[i] = - m_b[i]*Nsn[i] + S_to_Nsn + Lsn_to_Nsn + Lmn_to_Nsn - (theta + r + muN_age[i])*Nsn[i] -
                kneg*se_N_neg*rel_d*(l_s*(dstneg_n*sp_m_neg + (1-dstneg_n)) + dstneg_n*(1-sp_m_neg)*l_m)*Nsn[i] - /* care */ 
                forc[i+82]*Nsn[i] + (Nsn[i]/tot_age[i])*(forc[i+188]/5) + PTn_to_Nsn;  
      
      /* Smear neg, ds, prev */                             
      dNsp[i] = - m_b[i]*Nsp[i] + Lsp_to_Nsp + Lmp_to_Nsp - (theta + r + muN_age[i])*Nsp[i] -
                kneg*rel_d*se_N_neg*((1-e)*tneg_s*l_s*(sp_m_neg*dstneg_p + (1-dstneg_p)) + dstneg_p*(1-sp_m_neg)*l_m*tneg_m + l_s*e*(sp_m_neg*dstneg_p +(1-dstneg_p)))*Nsp[i] + 
                kneg*rel_d*se_N_neg*(l_s*(1-e)*(1-tneg_s)*(dstneg_n*sp_m_neg + (1-dstneg_n)) + dstneg_n*(1-sp_m_neg)*(1-tneg_m)*l_m)*Nsn[i] - /* Care */
                forc[i+82]*Nsp[i] + (Nsp[i]/tot_age[i])*(forc[i+188]/5) + PTp_to_Nsp;  

      /* Smear neg, mdr, new */
      dNmn[i] = - m_b[i]*Nmn[i] + S_to_Nmn + Lsn_to_Nmn + Lmn_to_Nmn - (theta + r + muN_age[i])*Nmn[i] -
                kneg*rel_d*se_N_neg*(l_m*dstneg_n*se_m_neg + l_s*((1-dstneg_n)+dstneg_n*(1-se_m_neg)))*Nmn[i] - /* Care */
                forc[i+82]*Nmn[i] + (Nmn[i]/tot_age[i])*(forc[i+188]/5) + PTn_to_Nmn;   
                   
      /* Smear neg, mdr, prev */
      dNmp[i] = - m_b[i]*Nmp[i] + Lsp_to_Nmp + Lmp_to_Nmp - (theta + r + muN_age[i])*Nmp[i] -
                kneg*rel_d*se_N_neg*(dstneg_p*se_m_neg*l_m*tneg_m + l_s*tneg_s*eff_p*((1-dstneg_p)+dstneg_p*(1-se_m_neg)))*Nmp[i] + 
                kneg*l_s*rel_d*e*se_N_neg*((dstneg_n*sp_m_neg+(1-dstneg_n))*Nsn[i]+(dstneg_p*sp_m_neg+(1-dstneg_p))*Nsp[i]) + 
                kneg*se_N_neg*rel_d*(dstneg_n*l_m*se_m_neg*(1-tneg_m) + l_s*(1-(tneg_s*eff_n))*((1-dstneg_n)+dstneg_n*(1-se_m_neg)))*Nmn[i] - /* Care */
                forc[i+82]*Nmp[i] + (Nmp[i]/tot_age[i])*(forc[i+188]/5) + PTp_to_Nmp;  

      /* Smear pos, ds, new */
      dIsn[i] = - m_b[i]*Isn[i] + S_to_Isn + Lsn_to_Isn + Lmn_to_Isn + theta*Nsn[i] - (r + muI_age[i])*Isn[i] - 
                kneg*se_I_neg*((dstneg_n*sp_m_neg + (1-dstneg_n))*l_s +(dstneg_n*(1-sp_m_neg)*l_m))*Isn[i] - /* Care */  
                forc[i+82]*Isn[i] + (Isn[i]/tot_age[i])*(forc[i+188]/5) + PTn_to_Isn;  
      
      /* Smear pos, ds, prev */
      dIsp[i] = - m_b[i]*Isp[i] + Lsp_to_Isp + Lmp_to_Isp + theta*Nsp[i] - (r + muI_age[i])*Isp[i] -   
                kneg*se_I_neg*(l_s*(1-e)*tneg_s*(sp_m_neg*dstneg_p + (1-dstneg_p)) + dstneg_p*(1-sp_m_neg)*l_m*tneg_m + l_s*e*(sp_m_neg*dstneg_p + (1-dstneg_p)))*Isp[i] + 
                kneg*se_I_neg*(l_s*(1-e)*(1-tneg_s)*(dstneg_n*sp_m_neg + (1-dstneg_n)) + (dstneg_n*(1-sp_m_neg)*l_m*(1-tneg_m)))*Isn[i] - /* Care */
                forc[i+82]*Isp[i] + (Isp[i]/tot_age[i])*(forc[i+188]/5) + PTp_to_Isp; 

      /* Smear pos, mdr, new */
      dImn[i] = - m_b[i]*Imn[i] + S_to_Imn + Lsn_to_Imn + Lmn_to_Imn + theta*Nmn[i] - (r + muI_age[i])*Imn[i] -
                kneg*se_I_neg*(l_m*dstneg_n*se_m_neg + (1-dstneg_n)*l_s + dstneg_n*(1-se_m_neg)*l_s)*Imn[i] - /* Care */ 
                forc[i+82]*Imn[i] + (Imn[i]/tot_age[i])*(forc[i+188]/5) + PTn_to_Imn;  
                
                
      /* Smear pos, mdr, prev */
      dImp[i] = - m_b[i]*Imp[i] + Lsp_to_Imp + Lmp_to_Imp + theta*Nmp[i] - (r + muI_age[i])*Imp[i] -
                kneg*se_I_neg*(l_m*dstneg_p*tneg_m*se_m_neg + l_s*tneg_s*eff_p*((1-dstneg_p)+dstneg_p*(1-se_m_neg)))*Imp[i] + 
                kneg*se_I_neg*(se_m_neg*l_m*dstneg_n*(1-tneg_m) + l_s*(1-(tneg_s*eff_n))*((1-dstneg_n)+dstneg_n*(1-se_m_neg)))*Imn[i] + 
                kneg*se_I_neg*l_s*e*((dstneg_n*sp_m_neg + (1-dstneg_n))*Isn[i]+(dstneg_p*sp_m_neg+ (1-dstneg_p))*Isp[i]) -
                forc[i+82]*Imp[i] + (Imp[i]/tot_age[i])*(forc[i+188]/5) + PTp_to_Imp;  
                
      /* Post PT, ds, new */          
      dPTn[i] = - m_b[i]*PTn[i] - PTn_to_Lsn - PTn_to_Nsn - PTn_to_Isn - PTn_to_Lmn - PTn_to_Nmn - PTn_to_Imn +
                health*kneg*(1-sp_I_neg*sp_N_neg)*l_s*tneg_s*Lsn[i] -
                forc[i+82]*PTn[i] + (PTn[i]/tot_age[i])*(forc[i+188]/5);
      
      /* Post PT, ds, prev */ 
      dPTp[i] = - m_b[i]*PTp[i] - PTp_to_Lsp - PTp_to_Nsp - PTp_to_Isp - PTp_to_Lmp - PTp_to_Nmp - PTp_to_Imp +
                health*kneg*(1-sp_I_neg*sp_N_neg)*l_s*tneg_s*Lsp[i] -
                forc[i+82]*PTp[i] + (PTp[i]/tot_age[i])*(forc[i+188]/5);
                        
      /* sum up new HIV- cases */           
      TB_cases_neg_age[i] = (v_age[i]*(1-sig_age[i]) + FS*a_age[i]*(1-p)*(1-sig_age[i]))*Lsn[i] + FS*a_age[i]*(1-sig_age[i])*(S[i] + (1-p)*(Lmn[i] + PTn[i])) + /*sneg,sus,new*/
                            (v_age[i]*(1-sig_age[i]) + FS*a_age[i]*(1-p)*(1-sig_age[i]))*Lsp[i] + FS*a_age[i]*(1-sig_age[i])*(1-p)*(Lmp[i] + PTp[i]) +          /*sneg,sus,prev*/
                            (v_age[i]*(1-sig_age[i]) + FM*a_age[i]*(1-p)*(1-sig_age[i]))*Lmn[i] + FM*a_age[i]*(1-sig_age[i])*(S[i] + (1-p)*(Lsn[i] + PTn[i])) + /*sneg,mdr,new*/
                            (v_age[i]*(1-sig_age[i]) + FM*a_age[i]*(1-p)*(1-sig_age[i]))*Lmp[i] + FM*a_age[i]*(1-sig_age[i])*(1-p)*(Lsp[i] + PTp[i]) +          /*sneg,mdr,prev*/
                            (v_age[i]*sig_age[i] + FS*a_age[i]*sig_age[i]*(1-p))*Lsn[i] + FS*a_age[i]*sig_age[i]*(S[i] + (1-p)*(Lmn[i] + PTn[i])) +             /*spos,sus,new*/
                            (v_age[i]*sig_age[i] + FS*a_age[i]*sig_age[i]*(1-p))*Lsp[i] + FS*a_age[i]*sig_age[i]*(1-p)*(Lmp[i] + PTp[i]) +                      /*spos,sus,prev*/
                            (v_age[i]*sig_age[i] + FM*a_age[i]*sig_age[i]*(1-p))*Lmn[i] + FM*a_age[i]*sig_age[i]*(S[i] + (1-p)*(Lsn[i] + PTn[i]))  +            /*spos,mdr,new*/
                            (v_age[i]*sig_age[i] + FM*a_age[i]*sig_age[i]*(1-p))*Lmp[i] + FM*a_age[i]*sig_age[i]*(1-p)*(Lsp[i] + PTp[i]);                       /*spos,mdr,prev*/
                            
      TB_cases_neg = TB_cases_neg + TB_cases_neg_age[i];
    
      TB_cases_age[i] = TB_cases_age[i] + TB_cases_neg_age[i];
    
        /* HIV+: Loop through CD4 categories */
    
        if (HIV_run>0.0){   /* don't run these for equilibrium as no HIV - save time */
    
        for (j=0; j<n_HIV; j++){      /* CD4 */


          double PTnH_to_LsnH = FS*(1-a_age_H[i][j])*(1-p_H[j])*PTn_H[i][j];            /* Post PT to latent DS (no disease history) */
          double PTnH_to_NsnH = FS*a_age_H[i][j]*(1-p_H[j])*(1-sig_H)*PTn_H[i][j];      /* Post PT to smear negative DS disease (no disease history) */
          double PTnH_to_IsnH = FS*a_age_H[i][j]*(1-p_H[j])*sig_H*PTn_H[i][j];          /* Post PT to smear positive DS disease (no disease history) */
          double PTnH_to_LmnH = FM*(1-a_age_H[i][j])*(1-p_H[j])*g*PTn_H[i][j];          /* Post PT to latent DR (no disease history) */
          double PTnH_to_NmnH = FM*a_age_H[i][j]*(1-p_H[j])*(1-sig_H)*PTn_H[i][j];      /* Post PT to smear negative DR disease (no disease history) */
          double PTnH_to_ImnH = FM*a_age_H[i][j]*(1-p)*sig_H*PTn_H[i][j];               /* Post PT to smear positive DR disease (no disease history) */
      
          double PTpH_to_LspH = FS*(1-a_age_H[i][j])*(1-p_H[j])*PTp_H[i][j];            /* Post PT to latent DS (prior Rx) */
          double PTpH_to_NspH = FS*a_age_H[i][j]*(1-p_H[j])*(1-sig_H)*PTp_H[i][j];      /* Post PT to smear negative DS disease (prior Rx) */
          double PTpH_to_IspH = FS*a_age_H[i][j]*(1-p_H[j])*sig_H*PTp_H[i][j];          /* Post PT to smear positive DS disease (prior Rx) */
          double PTpH_to_LmpH = FM*(1-a_age_H[i][j])*(1-p_H[j])*g*PTp_H[i][j];          /* Post PT to latent DR (prior Rx) */
          double PTpH_to_NmpH = FM*a_age_H[i][j]*(1-p_H[j])*(1-sig_H)*PTp_H[i][j];      /* Post PT to smear negative DR disease (prior Rx) */
          double PTpH_to_ImpH = FM*a_age_H[i][j]*(1-p_H[j])*sig_H*PTp_H[i][j];          /* Post PT to smear positive DR disease (prior Rx) */

          dS_H[i][j] = - m_b[i]*S_H[i][j] - 
                      (FS + FM)*S_H[i][j] +
                      forc[i+82]*H_CD4[j][i]*S[i] - H_prog[j+1][i]*S_H[i][j] + H_prog[j][i]*S_H[i][j-1] - 
                      up_H_mort[j][i]*S_H[i][j] - ART_prop[j]*S_H[i][j] + (S_H[i][j]/tot_age[i])*(forc[i+188]/5);

          dLsn_H[i][j] = - m_b[i]*Lsn_H[i][j] +
                        FS*((1-a_age_H[i][j])*S_H[i][j] + (1-a_age_H[i][j])*(1-p_H[j])*(1-g)*Lmn_H[i][j]) -
                        (v_age_H[i][j] + FS*a_age_H[i][j]*(1-p_H[j]) + FM*a_age_H[i][j]*(1-p_H[j]) + FM*(1-a_age_H[i][j])*(1-p_H[j])*g)*Lsn_H[i][j] +
                        r*(Isn_H[i][j] + Nsn_H[i][j]) +
                        forc[i+82]*H_CD4[j][i]*Lsn[i] - H_prog[j+1][i]*Lsn_H[i][j] + H_prog[j][i]*Lsn_H[i][j-1] - 
                        up_H_mort[j][i]*Lsn_H[i][j] - ART_prop[j]*Lsn_H[i][j] + (Lsn_H[i][j]/tot_age[i])*(forc[i+188]/5)-
                        health*kpos*(1-sp_I_pos*sp_N_pos)*l_s*tpos_s*Lsn_H[i][j] + PTnH_to_LsnH;

          dLsp_H[i][j] = - m_b[i]*Lsp_H[i][j] +
                        FS*(1-a_age_H[i][j])*(1-p_H[j])*(1-g)*Lmp_H[i][j] -
                        (v_age_H[i][j] + FS*a_age_H[i][j]*(1-p_H[j]) + FM*a_age_H[i][j]*(1-p_H[j]) + FM*(1-a_age_H[i][j])*(1-p_H[j])*g)*Lsp_H[i][j] +
                        r_H*(Isp_H[i][j]+Nsp_H[i][j]) +
                        kpos*(l_s*(1-e)*tpos_s*((dstpos_p*sp_m_pos)+(1-dstpos_p)) + l_m*tpos_m*dstpos_p*(1-sp_m_pos))*(se_I_pos*Isp_H[i][j] + se_N_pos*rel_d*Nsp_H[i][j]) + 
                        kpos*(l_s*(1-e)*tpos_s*((dstpos_n*sp_m_pos)+(1-dstpos_n)) + l_m*tpos_m*dstpos_n*(1-sp_m_pos))*(se_I_pos*Isn_H[i][j] + se_N_pos*rel_d*Nsn_H[i][j]) +
                        forc[i+82]*H_CD4[j][i]*Lsp[i] - H_prog[j+1][i]*Lsp_H[i][j] + H_prog[j][i]*Lsp_H[i][j-1] - 
                        up_H_mort[j][i]*Lsp_H[i][j] - ART_prop[j]*Lsp_H[i][j] + (Lsp_H[i][j]/tot_age[i])*(forc[i+188]/5)-
                        health*kpos*(1-sp_I_pos*sp_N_pos)*l_s*tpos_s*Lsp_H[i][j] + PTpH_to_LspH;

          dLmn_H[i][j] = - m_b[i]*Lmn_H[i][j] +
                        FM*((1-a_age_H[i][j])*S_H[i][j] + (1-a_age_H[i][j])*(1-p_H[j])*g*Lsn_H[i][j]) -
                        (v_age_H[i][j] + FM*a_age_H[i][j]*(1-p_H[j]) + FS*a_age_H[i][j]*(1-p_H[j]) + FS*(1-a_age_H[i][j])*(1-p_H[j])*(1-g))*Lmn_H[i][j] +
                        r_H*(Imn_H[i][j]+Nmn_H[i][j]) +
                        forc[i+82]*H_CD4[j][i]*Lmn[i] - H_prog[j+1][i]*Lmn_H[i][j] + H_prog[j][i]*Lmn_H[i][j-1] - 
                        up_H_mort[j][i]*Lmn_H[i][j] - ART_prop[j]*Lmn_H[i][j] + (Lmn_H[i][j]/tot_age[i])*(forc[i+188]/5) + PTnH_to_LmnH;
                     
          dLmp_H[i][j] = - m_b[i]*Lmp_H[i][j] + 
                        FM*(1-a_age_H[i][j])*(1-p_H[j])*g*Lsp_H[i][j] -
                        (v_age_H[i][j] + FM*a_age_H[i][j]*(1-p_H[j]) + FS*a_age_H[i][j]*(1-p_H[j]) + FS*(1-a_age_H[i][j])*(1-p_H[j])*(1-g))*Lmp_H[i][j] +
                        r_H*(Imp_H[i][j]+Nmp_H[i][j]) +
                        kpos*(dstpos_p*se_m_pos*l_m*tpos_m + l_s*tpos_s*eff_p*((1-dstpos_p)+dstpos_p*(1-se_m_pos)))*(se_I_pos*Imp_H[i][j] + se_N_pos*rel_d*Nmp_H[i][j]) +
                        kpos*(dstpos_n*se_m_pos*l_m*tpos_m + l_s*tpos_s*eff_n*((1-dstpos_n)+dstpos_n*(1-se_m_pos)))*(se_I_pos*Imn_H[i][j] + se_N_pos*rel_d*Nmn_H[i][j]) + 
                        forc[i+82]*H_CD4[j][i]*Lmp[i] - H_prog[j+1][i]*Lmp_H[i][j] + H_prog[j][i]*Lmp_H[i][j-1] - 
                        up_H_mort[j][i]*Lmp_H[i][j] - ART_prop[j]*Lmp_H[i][j] + (Lmp_H[i][j]/tot_age[i])*(forc[i+188]/5) + PTpH_to_LmpH;

          dNsn_H[i][j] = - m_b[i]*Nsn_H[i][j] +
                        (v_age_H[i][j]*(1-sig_H) + FS*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lsn_H[i][j] +
                        FS*a_age_H[i][j]*(1-sig_H)*(S_H[i][j] + (1-p_H[j])*Lmn_H[i][j]) -
                        (theta_H + r_H + muN_H)*Nsn_H[i][j] -
                        kpos*se_N_pos*rel_d*(l_s*(dstpos_n*sp_m_pos + (1-dstpos_n)) + dstpos_n*(1-sp_m_pos)*l_m)*Nsn_H[i][j] + 
                        forc[i+82]*H_CD4[j][i]*Nsn[i] - H_prog[j+1][i]*Nsn_H[i][j] + H_prog[j][i]*Nsn_H[i][j-1] - 
                        up_H_mort[j][i]*Nsn_H[i][j] - ART_prop[j]*Nsn_H[i][j] + (Nsn_H[i][j]/tot_age[i])*(forc[i+188]/5) + PTnH_to_NsnH;
                                                
          dNsp_H[i][j] = - m_b[i]*Nsp_H[i][j] + 
                        (v_age_H[i][j]*(1-sig_H) + FS*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lsp_H[i][j] +
                        FS*a_age_H[i][j]*(1-sig_H)*(1-p_H[j])*Lmp_H[i][j] -
                        (theta_H + r_H + muN_H)*Nsp_H[i][j] -
                        kpos*rel_d*se_N_pos*((1-e)*tpos_s*l_s*(sp_m_pos*dstpos_p + (1-dstpos_p)) + dstpos_p*(1-sp_m_pos)*l_m*tpos_m + l_s*e*(sp_m_pos*dstpos_p +(1-dstpos_p)))*Nsp_H[i][j] + 
                        kpos*rel_d*se_N_pos*(l_s*(1-e)*(1-tpos_s)*(dstpos_n*sp_m_pos + (1-dstpos_n)) + dstpos_n*(1-sp_m_pos)*(1-tpos_m)*l_m)*Nsn_H[i][j] +                       
                        forc[i+82]*H_CD4[j][i]*Nsp[i] - H_prog[j+1][i]*Nsp_H[i][j] + H_prog[j][i]*Nsp_H[i][j-1] - 
                        up_H_mort[j][i]*Nsp_H[i][j] - ART_prop[j]*Nsp_H[i][j] + (Nsp_H[i][j]/tot_age[i])*(forc[i+188]/5) + PTpH_to_NspH;

          dNmn_H[i][j] = - m_b[i]*Nmn_H[i][j] +
                        (v_age_H[i][j]*(1-sig_H) + FM*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lmn_H[i][j] +
                        FM*a_age_H[i][j]*(1-sig_H)*(S_H[i][j] + (1-p_H[j])*Lsn_H[i][j]) -
                        (theta_H + r_H + muN_H)*Nmn_H[i][j] -
                        kpos*rel_d*se_N_pos*(l_m*dstpos_n*se_m_pos + l_s*((1-dstpos_n)+dstpos_n*(1-se_m_pos)))*Nmn_H[i][j] +
                        forc[i+82]*H_CD4[j][i]*Nmn[i] - H_prog[j+1][i]*Nmn_H[i][j] + H_prog[j][i]*Nmn_H[i][j-1] - 
                        up_H_mort[j][i]*Nmn_H[i][j] - ART_prop[j]*Nmn_H[i][j] + (Nmn_H[i][j]/tot_age[i])*(forc[i+188]/5) + PTnH_to_NmnH;
                       
          dNmp_H[i][j] = - m_b[i]*Nmp_H[i][j] +
                        (v_age_H[i][j]*(1-sig_H) + FM*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lmp_H[i][j] +
                        FM*a_age_H[i][j]*(1-sig_H)*(1-p_H[j])*Lsp_H[i][j] -
                        (theta_H + r_H + muN_H)*Nmp_H[i][j] -
                        kpos*rel_d*se_N_pos*(dstpos_p*se_m_pos*l_m*tpos_m + l_s*tpos_s*eff_p*((1-dstpos_p)+dstpos_p*(1-se_m_pos)))*Nmp_H[i][j] + 
                        kpos*l_s*rel_d*e*se_N_pos*((dstpos_n*sp_m_pos+(1-dstpos_n))*Nsn_H[i][j]+(dstpos_p*sp_m_pos+(1-dstpos_p))*Nsp_H[i][j]) + 
                        kpos*se_N_pos*rel_d*(dstpos_n*l_m*se_m_pos*(1-tpos_m) + l_s*(1-(tpos_s*eff_n))*((1-dstpos_n)+dstpos_n*(1-se_m_pos)))*Nmn_H[i][j] +
                        forc[i+82]*H_CD4[j][i]*Nmp[i] - H_prog[j+1][i]*Nmp_H[i][j] + H_prog[j][i]*Nmp_H[i][j-1] - 
                        up_H_mort[j][i]*Nmp_H[i][j] - ART_prop[j]*Nmp_H[i][j] + (Nmp_H[i][j]/tot_age[i])*(forc[i+188]/5) + PTpH_to_NmpH;

          dIsn_H[i][j] = - m_b[i]*Isn_H[i][j] +
                        (v_age_H[i][j]*sig_H + FS*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lsn_H[i][j] +
                        FS*a_age_H[i][j]*sig_H*S_H[i][j] +
                        FS*a_age_H[i][j]*(1-p_H[j])*sig_H*Lmn_H[i][j] +
                        theta_H*Nsn_H[i][j] -
                        (r_H + muI_H)*Isn_H[i][j] -
                        kpos*se_I_pos*((dstpos_n*sp_m_pos + (1-dstpos_n))*l_s +(dstpos_n*(1-sp_m_pos)*l_m))*Isn_H[i][j] +
                        forc[i+82]*H_CD4[j][i]*Isn[i] - H_prog[j+1][i]*Isn_H[i][j] + H_prog[j][i]*Isn_H[i][j-1] - 
                        up_H_mort[j][i]*Isn_H[i][j] - ART_prop[j]*Isn_H[i][j] + (Isn_H[i][j]/tot_age[i])*(forc[i+188]/5) + PTnH_to_IsnH;

          dIsp_H[i][j] = - m_b[i]*Isp_H[i][j] +
                        (v_age_H[i][j]*sig_H + FS*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lsp_H[i][j] +
                        FS*a_age_H[i][j]*sig_H*(1-p_H[j])*Lmp_H[i][j] +
                        theta_H*Nsp_H[i][j] +
                        kpos*se_I_pos*(l_s*(1-e)*(1-tpos_s)*(dstpos_n*sp_m_pos + (1-dstpos_n)) + (dstpos_n*(1-sp_m_pos)*l_m*(1-tpos_m)))*Isn_H[i][j] - 
                        (r_H + muI_H)*Isp_H[i][j] -
                        kpos*se_I_pos*(l_s*(1-e)*tpos_s*(sp_m_pos*dstpos_p + (1-dstpos_p)) + dstpos_p*(1-sp_m_pos)*l_m*tpos_m + l_s*e*(sp_m_pos*dstpos_p + (1-dstpos_p)))*Isp_H[i][j] +
                        forc[i+82]*H_CD4[j][i]*Isp[i] - H_prog[j+1][i]*Isp_H[i][j] + H_prog[j][i]*Isp_H[i][j-1] - 
                        up_H_mort[j][i]*Isp_H[i][j] - ART_prop[j]*Isp_H[i][j] + (Isp_H[i][j]/tot_age[i])*(forc[i+188]/5) + PTpH_to_IspH;       
     
          dImn_H[i][j] = - m_b[i]*Imn_H[i][j] +
                        (v_age_H[i][j]*sig_H + FM*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lmn_H[i][j] +
                        FM*a_age_H[i][j]*sig_H*S_H[i][j] +
                        FM*a_age_H[i][j]*(1-p_H[j])*sig_H*Lsn_H[i][j] +
                        theta_H*Nmn_H[i][j] -
                        kpos*se_I_pos*(l_m*dstpos_n*se_m_pos + (1-dstpos_n)*l_s + dstpos_n*(1-se_m_pos)*l_s)*Imn_H[i][j] -
                        (r_H + muI_H)*Imn_H[i][j] +
                        forc[i+82]*H_CD4[j][i]*Imn[i] - H_prog[j+1][i]*Imn_H[i][j] + H_prog[j][i]*Imn_H[i][j-1] - 
                        up_H_mort[j][i]*Imn_H[i][j] - ART_prop[j]*Imn_H[i][j] + (Imn_H[i][j]/tot_age[i])*(forc[i+188]/5) + PTnH_to_ImnH;

          dImp_H[i][j] = - m_b[i]*Imp_H[i][j] +
                        (v_age_H[i][j]*sig_H + FM*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lmp_H[i][j] +
                        FM*a_age_H[i][j]*sig_H*(1-p_H[j])*Lsp_H[i][j] +
                        theta_H*Nmp_H[i][j] + 
                        kpos*se_I_pos*(se_m_pos*l_m*dstpos_n*(1-tpos_m) + l_s*(1-(tpos_s*eff_n))*((1-dstpos_n)+dstpos_n*(1-se_m_pos)))*Imn_H[i][j] + 
                        kpos*se_I_pos*l_s*e*((dstpos_n*sp_m_pos + (1-dstpos_n))*Isn_H[i][j]+(dstpos_p*sp_m_pos+ (1-dstpos_p))*Isp_H[i][j]) -
                        (r_H + muI_H)*Imp_H[i][j] -
                        kpos*se_I_pos*(l_m*dstpos_p*tpos_m*se_m_pos + l_s*tpos_s*eff_p*((1-dstpos_p)+dstpos_p*(1-se_m_pos)))*Imp_H[i][j]+
                        forc[i+82]*H_CD4[j][i]*Imp[i] - H_prog[j+1][i]*Imp_H[i][j] + H_prog[j][i]*Imp_H[i][j-1] - 
                        up_H_mort[j][i]*Imp_H[i][j] - ART_prop[j]*Imp_H[i][j] + (Imp_H[i][j]/tot_age[i])*(forc[i+188]/5) + PTpH_to_ImpH;       
                                                      
          dPTn_H[i][j] = -m_b[i]*PTn_H[i][j] - PTnH_to_LsnH - PTnH_to_NsnH - PTnH_to_IsnH - PTnH_to_LmnH - PTnH_to_NmnH - PTnH_to_ImnH +
                         health*kpos*(1-sp_I_pos*sp_N_pos)*l_s*tpos_s*Lsn_H[i][j] +
                         forc[i+82]*H_CD4[j][i]*PTn[i] - H_prog[j+1][i]*PTn_H[i][j] + H_prog[j][i]*PTn_H[i][j-1] - 
                         up_H_mort[j][i]*PTn_H[i][j] - ART_prop[j]*PTn_H[i][j] + (PTn_H[i][j]/tot_age[i])*(forc[i+188]/5);
          
          dPTp_H[i][j] = - m_b[i]*PTp_H[i][j] - PTpH_to_LspH - PTpH_to_NspH - PTpH_to_IspH - PTpH_to_LmpH - PTpH_to_NmpH - PTpH_to_ImpH +
                         health*kpos*(1-sp_I_pos*sp_N_pos)*l_s*tpos_s*Lsp_H[i][j] +
                         forc[i+82]*H_CD4[j][i]*PTp[i] - H_prog[j+1][i]*PTp_H[i][j] + H_prog[j][i]*PTp_H[i][j-1] - 
                         up_H_mort[j][i]*PTp_H[i][j] - ART_prop[j]*PTp_H[i][j] + (PTp_H[i][j]/tot_age[i])*(forc[i+188]/5);
                   
          TB_cases_pos_age[i][j] =(v_age_H[i][j]*(1-sig_H) + FS*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lsn_H[i][j] + FS*a_age_H[i][j]*(1-sig_H)*(S_H[i][j] + (1-p_H[j])*(Lmn_H[i][j] + PTn_H[i][j])) + 
                                  (v_age_H[i][j]*(1-sig_H) + FS*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lsp_H[i][j] + FS*a_age_H[i][j]*(1-sig_H)*(1-p_H[j])*(Lmp_H[i][j] + PTp_H[i][j]) +
                                  (v_age_H[i][j]*(1-sig_H) + FM*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lmn_H[i][j] + FM*a_age_H[i][j]*(1-sig_H)*(S_H[i][j] + (1-p_H[j])*(Lsn_H[i][j] + PTn_H[i][j])) +         
                                  (v_age_H[i][j]*(1-sig_H) + FM*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lmp_H[i][j] + FM*a_age_H[i][j]*(1-sig_H)*(1-p_H[j])*(Lsp_H[i][j] + PTp_H[i][j]) +
                                  (v_age_H[i][j]*sig_H + FS*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lsn_H[i][j] + FS*a_age_H[i][j]*sig_H*(S_H[i][j] + (1-p_H[j])*(Lmn_H[i][j] + PTn_H[i][j]))+
                                  (v_age_H[i][j]*sig_H + FS*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lsp_H[i][j] + FS*a_age_H[i][j]*sig_H*(1-p_H[j])*(Lmp_H[i][j] + PTp_H[i][j]) +
                                  (v_age_H[i][j]*sig_H + FM*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lmn_H[i][j] + FM*a_age_H[i][j]*sig_H*(S_H[i][j] + (1-p_H[j])*(Lsn_H[i][j] + PTn_H[i][j]))+
                                  (v_age_H[i][j]*sig_H + FM*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lmp_H[i][j] + FM*a_age_H[i][j]*sig_H*(1-p_H[j])*(Lsp_H[i][j] + PTp_H[i][j]);

          TB_cases_pos = TB_cases_pos + TB_cases_pos_age[i][j];
          
          TB_cases_age[i] = TB_cases_age[i] + TB_cases_pos_age[i][j];
       
          if(ART_on>0.0){ /* don't run ART equations if no ART yet - save some time */
       /* HIV+ on ART: loop through time on ART, CD4 at initiation, age */
          for (l=0; l<n_ART; l++){
        
            dS_A[i][j][l] = - m_b[i]*S_A[i][j][l] - (FS + FM)*S_A[i][j][l] +
                            ART_prop[j]*A_start[l]*S_H[i][j] + A_prog[l]*S_A[i][j][l-1] - A_prog[l+1]*S_A[i][j][l] - up_A_mort[l][j][i]*S_A[i][j][l] + 
                            (S_A[i][j][l]/tot_age[i])*(forc[i+188]/5);

            dLsn_A[i][j][l] = - m_b[i]*Lsn_A[i][j][l] +
                            FS*((1-a_age_A[i][j][l])*S_A[i][j][l] + (1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g)*Lmn_A[i][j][l]) -
                            (v_age_A[i][j][l] + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g)*Lsn_A[i][j][l] +
                            r_H*(Isn_A[i][j][l] + Nsn_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lsn_H[i][j] + A_prog[l]*Lsn_A[i][j][l-1] - A_prog[l+1]*Lsn_A[i][j][l] - up_A_mort[l][j][i]*Lsn_A[i][j][l] +
                            (Lsn_A[i][j][l]/tot_age[i])*(forc[i+188]/5);
        
            dLsp_A[i][j][l] = - m_b[i]*Lsp_A[i][j][l] +
                            FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g)*Lmp_A[i][j][l] -
                            (v_age_A[i][j][l] + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g)*Lsp_A[i][j][l] +
                            r_H*(Isp_A[i][j][l]+Nsp_A[i][j][l]) +
                            kpos*l_s*(1-e)*tART_s*(Isn_A[i][j][l]+Isp_A[i][j][l]) + 
                            kpos*l_s*(1-e)*tART_s*rel_d*(Nsn_A[i][j][l]+Nsp_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lsp_H[i][j] + A_prog[l]*Lsp_A[i][j][l-1] - A_prog[l+1]*Lsp_A[i][j][l] - up_A_mort[l][j][i]*Lsp_A[i][j][l] +
                            (Lsp_A[i][j][l]/tot_age[i])*(forc[i+188]/5);
                             
            dLmn_A[i][j][l] = - m_b[i]*Lmn_A[i][j][l] +
                            FM*((1-a_age_A[i][j][l])*S_A[i][j][l] + (1-a_age_A[i][j][l])*(1-p_A[j][l])*g*Lsn_A[i][j][l]) -
                            (v_age_A[i][j][l] + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g))*Lmn_A[i][j][l] +
                            r_H*(Imn_A[i][j][l]+Nmn_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lmn_H[i][j] + A_prog[l]*Lmn_A[i][j][l-1] - A_prog[l+1]*Lmn_A[i][j][l] - up_A_mort[l][j][i]*Lmn_A[i][j][l] +
                            (Lmn_A[i][j][l]/tot_age[i])*(forc[i+188]/5);
                     
            dLmp_A[i][j][l] = - m_b[i]*Lmp_A[i][j][l] + 
                            FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g*Lsp_A[i][j][l] -
                            (v_age_A[i][j][l] + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g))*Lmp_A[i][j][l] +
                            r_H*(Imp_A[i][j][l]+Nmp_A[i][j][l]) +
                            kpos*(l_m*dstpos_n*tART_m + l_s*(1-dstpos_n)*tART_s*eff_n)*(Imn_A[i][j][l]+rel_d*Nmn_A[i][j][l]) +
                            kpos*(l_m*dstpos_p*tART_m + l_s*(1-dstpos_p)*tART_s*eff_p)*(Imp_A[i][j][l]+rel_d*Nmp_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lmp_H[i][j] + A_prog[l]*Lmp_A[i][j][l-1] - A_prog[l+1]*Lmp_A[i][j][l] - up_A_mort[l][j][i]*Lmp_A[i][j][l] +
                            (Lmp_A[i][j][l]/tot_age[i])*(forc[i+188]/5);
           
            dNsn_A[i][j][l] = - m_b[i]*Nsn_A[i][j][l]+
                            (v_age_A[i][j][l]*(1-sig_H) + FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lsn_A[i][j][l] +
                            FS*a_age_A[i][j][l]*(1-sig_H)*(S_A[i][j][l] + (1-p_A[j][l])*Lmn_A[i][j][l]) -
                            (theta_H + r_H + kpos*l_s*rel_d + muN_H_A[l])*Nsn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Nsn_H[i][j] + A_prog[l]*Nsn_A[i][j][l-1] - A_prog[l+1]*Nsn_A[i][j][l] - up_A_mort[l][j][i]*Nsn_A[i][j][l] +
                            (Nsn_A[i][j][l]/tot_age[i])*(forc[i+188]/5);
                                                
            dNsp_A[i][j][l] = - m_b[i]*Nsp_A[i][j][l] + 
                            (v_age_A[i][j][l]*(1-sig_H) + FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lsp_A[i][j][l] +
                            FS*a_age_A[i][j][l]*(1-sig_H)*(1-p_A[j][l])*Lmp_A[i][j][l] -
                            (theta_H + r_H + kpos*l_s*rel_d*(1-e)*tART_s + kpos*l_s*rel_d*e + muN_H_A[l])*Nsp_A[i][j][l] +
                            kpos*l_s*rel_d*(1-e)*(1-tART_s)*Nsn_A[i][j][l] +      
                            ART_prop[j]*A_start[l]*Nsp_H[i][j] + A_prog[l]*Nsp_A[i][j][l-1] - A_prog[l+1]*Nsp_A[i][j][l] - up_A_mort[l][j][i]*Nsp_A[i][j][l] +
                            (Nsp_A[i][j][l]/tot_age[i])*(forc[i+188]/5);
                       
            dNmn_A[i][j][l] = - m_b[i]*Nmn_A[i][j][l] +
                            (v_age_A[i][j][l]*(1-sig_H) + FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lmn_A[i][j][l] +
                            FM*a_age_A[i][j][l]*(1-sig_H)*(S_A[i][j][l] + (1-p_A[j][l])*Lsn_A[i][j][l]) -
                            (theta_H + r_H + kpos*l_m*rel_d*dstpos_n + kpos*l_s*rel_d*(1-dstpos_n) + muN_H_A[l])*Nmn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Nmn_H[i][j] + A_prog[l]*Nmn_A[i][j][l-1] - A_prog[l+1]*Nmn_A[i][j][l] - up_A_mort[l][j][i]*Nmn_A[i][j][l] +
                            (Nmn_A[i][j][l]/tot_age[i])*(forc[i+188]/5);
                       
            dNmp_A[i][j][l] = - m_b[i]*Nmp_A[i][j][l] +
                            (v_age_A[i][j][l]*(1-sig_H) + FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lmp_A[i][j][l] +
                            FM*a_age_A[i][j][l]*(1-sig_H)*(1-p_A[j][l])*Lsp_A[i][j][l] -
                            (theta_H + r_H + kpos*l_m*rel_d*dstpos_p*tART_m + kpos*l_s*rel_d*(1-dstpos_p)*tART_s*eff_p + muN_H_A[l])*Nmp_A[i][j][l] +
                            kpos*l_s*rel_d*e*(Nsn_A[i][j][l]+Nsp_A[i][j][l]) +
                            (kpos*l_m*rel_d*dstpos_n*(1-tART_m) + kpos*l_s*rel_d*(1-dstpos_n)*(1-(tART_s*eff_n)))*Nmn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Nmp_H[i][j] + A_prog[l]*Nmp_A[i][j][l-1] - A_prog[l+1]*Nmp_A[i][j][l] - up_A_mort[l][j][i]*Nmp_A[i][j][l] +
                            (Nmp_A[i][j][l]/tot_age[i])*(forc[i+188]/5);
           
            dIsn_A[i][j][l] = - m_b[i]*Isn_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_H + FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lsn_A[i][j][l] +
                            FS*a_age_A[i][j][l]*sig_H*S_A[i][j][l] +
                            FS*a_age_A[i][j][l]*(1-p_A[j][l])*sig_H*Lmn_A[i][j][l] +
                            theta_H*Nsn_A[i][j][l] -
                            (r_H + kpos*l_s + muI_H_A[l])*Isn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Isn_H[i][j] + A_prog[l]*Isn_A[i][j][l-1] - A_prog[l+1]*Isn_A[i][j][l] - up_A_mort[l][j][i]*Isn_A[i][j][l] +
                            (Isn_A[i][j][l]/tot_age[i])*(forc[i+188]/5);

            dIsp_A[i][j][l] = - m_b[i]*Isp_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_H + FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lsp_A[i][j][l] +
                            FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l])*Lmp_A[i][j][l] +
                            theta_H*Nsp_A[i][j][l] +
                            kpos*l_s*(1-e)*(1-tART_s)*Isn_A[i][j][l] -
                            (r_H + kpos*l_s*(1-e)*tART_s + kpos*l_s*e + muI_H_A[l])*Isp_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Isp_H[i][j] + A_prog[l]*Isp_A[i][j][l-1] - A_prog[l+1]*Isp_A[i][j][l] - up_A_mort[l][j][i]*Isp_A[i][j][l] +
                            (Isp_A[i][j][l]/tot_age[i])*(forc[i+188]/5);
                                   
            dImn_A[i][j][l] = - m_b[i]*Imn_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_H + FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lmn_A[i][j][l] +
                            FM*a_age_A[i][j][l]*sig_H*S_A[i][j][l] +
                            FM*a_age_A[i][j][l]*(1-p_A[j][l])*sig_H*Lsn_A[i][j][l] +
                            theta_H*Nmn_A[i][j][l] -
                            (r_H + kpos*l_m*dstpos_n + kpos*l_s*(1-dstpos_n) + muI_H_A[l])*Imn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Imn_H[i][j] + A_prog[l]*Imn_A[i][j][l-1] - A_prog[l+1]*Imn_A[i][j][l] - up_A_mort[l][j][i]*Imn_A[i][j][l] +
                            (Imn_A[i][j][l]/tot_age[i])*(forc[i+188]/5);

            dImp_A[i][j][l] = - m_b[i]*Imp_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_H + FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lmp_A[i][j][l] +
                            FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l])*Lsp_A[i][j][l] +
                            theta_H*Nmp_A[i][j][l] +
                            (kpos*l_m*dstpos_n*(1-tART_m) + kpos*l_s*(1-dstpos_n)*(1-(tART_s*eff_n)))*Imn_A[i][j][l] +
                            kpos*l_s*e*(Isn_A[i][j][l]+Isp_A[i][j][l]) -
                            (r_H + kpos*l_m*dstpos_p*tART_m + kpos*l_s*(1-dstpos_p)*tART_s*eff_p + muI_H_A[l])*Imp_A[i][j][l] +                     
                            ART_prop[j]*A_start[l]*Imp_H[i][j] + A_prog[l]*Imp_A[i][j][l-1] - A_prog[l+1]*Imp_A[i][j][l] - up_A_mort[l][j][i]*Imp_A[i][j][l] +
                            (Imp_A[i][j][l]/tot_age[i])*(forc[i+188]/5);             
        
                                   
            dPTn_A[i][j][l] = 0;
            dPTp_A[i][j][l] = 0;   
        
            TB_cases_ART_age[i][j][l] = (v_age_A[i][j][l]*(1-sig_H) + FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lsn_A[i][j][l] + FS*a_age_A[i][j][l]*(1-sig_H)*(S_A[i][j][l] + (1-p_A[j][l])*Lmn_A[i][j][l]) +
                                      (v_age_A[i][j][l]*(1-sig_H) + FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lsp_A[i][j][l] + FS*a_age_A[i][j][l]*(1-sig_H)*(1-p_A[j][l])*Lmp_A[i][j][l] +
                                      (v_age_A[i][j][l]*(1-sig_H) + FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lmn_A[i][j][l] + FM*a_age_A[i][j][l]*(1-sig_H)*(S_A[i][j][l] + (1-p_A[j][l])*Lsn_A[i][j][l]) +
                                      (v_age_A[i][j][l]*(1-sig_H) + FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lmp_A[i][j][l] + FM*a_age_A[i][j][l]*(1-sig_H)*(1-p_A[j][l])*Lsp_A[i][j][l] +
                                      (v_age_A[i][j][l]*sig_H + FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lsn_A[i][j][l] + FS*a_age_A[i][j][l]*sig_H*S_A[i][j][l] + FS*a_age_A[i][j][l]*(1-p_A[j][l])*sig_H*Lmn_A[i][j][l] +
                                      (v_age_A[i][j][l]*sig_H + FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lsp_A[i][j][l] + FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l])*Lmp_A[i][j][l] +
                                      (v_age_A[i][j][l]*sig_H + FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lmn_A[i][j][l] + FM*a_age_A[i][j][l]*sig_H*S_A[i][j][l] + FM*a_age_A[i][j][l]*(1-p_A[j][l])*sig_H*Lsn_A[i][j][l] +
                                      (v_age_A[i][j][l]*sig_H + FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lmp_A[i][j][l] + FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l])*Lsp_A[i][j][l];
                            
                            
                            
                            
                            
                            
                            
            TB_cases_ART = TB_cases_ART + TB_cases_ART_age[i][j][l];
        
            TB_cases_age[i] = TB_cases_age[i] + TB_cases_ART_age[i][j][l];
        
          }  /* end loop on ART */
          }

        }    /* end loop on HIV */
        }  

    }        /* end loop on age */

    /* Put our calculated rates of change back into ydot */
    /*  AGAIN, MUST BE A NEATER WAY TO DO THIS */
    /* HIV+ */
    for (i=0; i<n_age; i++) ydot[i] = dS[i];             
    for (i=n_age; i<n_age*2; i++) ydot[i] = dLsn[i-n_age];      
    for (i=n_age*2; i<n_age*3; i++) ydot[i] = dLsp[i-n_age*2];    /* Lsp: 34-50 */
    for (i=n_age*3; i<n_age*4; i++) ydot[i] = dLmn[i-n_age*3];    /* Lmn: 51-67 */
    for (i=n_age*4; i<n_age*5; i++) ydot[i] = dLmp[i-n_age*4];    /* Lmp: 68-84 */
    for (i=n_age*5; i<n_age*6; i++) ydot[i] = dNsn[i-n_age*5];    /* Nsn: 85-101 */
    for (i=n_age*6; i<n_age*7; i++) ydot[i] = dNsp[i-n_age*6];    /* Nsp: 102-118 */
    for (i=n_age*7; i<n_age*8; i++) ydot[i] = dNmn[i-n_age*7];    /* Nmn: 119-135 */
    for (i=n_age*8; i<n_age*9; i++) ydot[i] = dNmp[i-n_age*8];    /* Nmp: 136-152 */ 
    for (i=n_age*9; i<n_age*10; i++) ydot[i] = dIsn[i-n_age*9];   /* Isn: 153-169 */
    for (i=n_age*10; i<n_age*11; i++) ydot[i] = dIsp[i-n_age*10]; /* Isp: 170-186 */
    for (i=n_age*11; i<n_age*12; i++) ydot[i] = dImn[i-n_age*11]; /* Imn: 187-203 */
    for (i=n_age*12; i<n_age*13; i++) ydot[i] = dImp[i-n_age*12]; /* Imp: 204-220 */
    for (i=n_age*13; i<n_age*14; i++) ydot[i] = dPTn[i-n_age*13]; /* Imn: 187-203 */
    for (i=n_age*14; i<n_age*15; i++) ydot[i] = dPTp[i-n_age*14]; /* Imp: 204-220 */
    
    /* HIV+ */
    ij = n_age*n_disease;                                     
    for (j=0; j<n_HIV; j++){
      for (i=0; i<n_age; i++){  
        ydot[ij] = dS_H[i][j];
        ydot[ij+(n_age*n_HIV)] = dLsn_H[i][j];
        ydot[ij+(2*n_age*n_HIV)] = dLsp_H[i][j]; 
        ydot[ij+(3*n_age*n_HIV)] = dLmn_H[i][j];
        ydot[ij+(4*n_age*n_HIV)] = dLmp_H[i][j];  
        ydot[ij+(5*n_age*n_HIV)] = dNsn_H[i][j];
        ydot[ij+(6*n_age*n_HIV)] = dNsp_H[i][j];
        ydot[ij+(7*n_age*n_HIV)] = dNmn_H[i][j];
        ydot[ij+(8*n_age*n_HIV)] = dNmp_H[i][j];
        ydot[ij+(9*n_age*n_HIV)] = dIsn_H[i][j];
        ydot[ij+(10*n_age*n_HIV)] = dIsp_H[i][j];
        ydot[ij+(11*n_age*n_HIV)] = dImn_H[i][j];
        ydot[ij+(12*n_age*n_HIV)] = dImp_H[i][j];
        ydot[ij+(13*n_age*n_HIV)] = dPTn_H[i][j];
        ydot[ij+(14*n_age*n_HIV)] = dPTp_H[i][j];
        ij = ij+1;
      }
    }
    /* HIV+, on ART */
    ij = (n_HIV+1)*n_disease*n_age;   
    for (l=0; l<n_ART; l++){                          
      for (j=0; j<n_HIV; j++){
        for (i=0; i<n_age; i++){
          ydot[ij] = dS_A[i][j][l];
          ydot[ij+(n_age*n_ART*n_HIV)] = dLsn_A[i][j][l];
          ydot[ij+(2*n_age*n_ART*n_HIV)] = dLsp_A[i][j][l];
          ydot[ij+(3*n_age*n_ART*n_HIV)] = dLmn_A[i][j][l];
          ydot[ij+(4*n_age*n_ART*n_HIV)] = dLmp_A[i][j][l];  
          ydot[ij+(5*n_age*n_ART*n_HIV)] = dNsn_A[i][j][l];
          ydot[ij+(6*n_age*n_ART*n_HIV)] = dNsp_A[i][j][l];
          ydot[ij+(7*n_age*n_ART*n_HIV)] = dNmn_A[i][j][l];
          ydot[ij+(8*n_age*n_ART*n_HIV)] = dNmp_A[i][j][l];
          ydot[ij+(9*n_age*n_ART*n_HIV)] = dIsn_A[i][j][l];
          ydot[ij+(10*n_age*n_ART*n_HIV)] = dIsp_A[i][j][l];
          ydot[ij+(11*n_age*n_ART*n_HIV)] = dImn_A[i][j][l];
          ydot[ij+(12*n_age*n_ART*n_HIV)] = dImp_A[i][j][l];
          ydot[ij+(13*n_age*n_ART*n_HIV)] = dPTn_A[i][j][l];
          ydot[ij+(14*n_age*n_ART*n_HIV)] = dPTp_A[i][j][l];
          ij = ij+1;
        }
      }
    }
    
    /* Calculate notifications, treatments etc */
    double DS_correct = 0;
    double DS_incorrect = 0; 
    double MDR_correct = 0; 
    double MDR_incorrect = 0;
    double FP = 0; 
    
    for (i=0; i<n_age; i++){
    
      /* ## DS notifications ## */
      
      /* Correct */
      DS_correct =  DS_correct + rel_d*kneg*se_N_neg*(dstneg_n*sp_m_neg + (1-dstneg_n))*l_s*Nsn[i] + /* s-, ds, new */
                                 rel_d*kneg*se_N_neg*(dstneg_p*sp_m_neg + (1-dstneg_p))*l_s*Nsp[i] + /* s-, ds, prev */                                        
                                 kneg*se_I_neg*(dstneg_n*sp_m_neg + (1-dstneg_n))*l_s*Isn[i] +       /* s-, ds, new */
                                 kneg*se_I_neg*(dstneg_p*sp_m_neg + (1-dstneg_p))*l_s*Isp[i];        /* s-, ds, prev */  
             
      /* Incorrect - should be MDR */    
      DS_incorrect = DS_incorrect + rel_d*kneg*se_N_neg*(dstneg_n*(1-se_m_neg) + (1-dstneg_n))*l_s*Nmn[i] + /* s-, dr, new */
                                    rel_d*kneg*se_N_neg*(dstneg_p*(1-se_m_neg) + (1-dstneg_p))*l_s*Nmp[i] + /* s-, dr, prev */                                        
                                    kneg*se_I_neg*(dstneg_n*(1-se_m_neg) + (1-dstneg_n))*l_s*Imn[i] +       /* s-, dr, new */
                                    kneg*se_I_neg*(dstneg_p*(1-se_m_neg) + (1-dstneg_p))*l_s*Imp[i];        /* s-, dr, prev */                            
    
      /* ## MDR notifications ## */
    
      /* Correct */ 
      MDR_correct = MDR_correct + rel_d*kneg*se_N_neg*l_m*(dstneg_n*Nmn[i] + dstneg_p*Nmp[i]) + /* s-, dr */
                                  kneg*se_I_neg*l_m*(dstneg_n*Imn[i] + dstneg_p*Imp[i]);        /* s+, dr */
                          
      /* Incorrect - should be DS */
      MDR_incorrect = MDR_incorrect + rel_d*kneg*se_N_neg*(1-sp_m_neg)*l_m*(dstneg_n*Nsn[i] + dstneg_p*Nsp[i]) + /* s-, ds */
                                      kneg*se_I_neg*(1-sp_m_neg)*l_m*(dstneg_n*Isn[i] + dstneg_p*Isp[i]);        /* s+, ds */ 
                          
      /* ## False positives - don't have TB ## */
      FP = FP + health*kneg*(1-sp_I_neg*sp_N_neg)*l_s*(S[i]+Lsn[i]+Lsp[i]+Lmn[i]+Lmp[i]+PTn[i]+PTp[i]);
               
      for (j=0; j<n_HIV; j++){
                          
          DS_correct = DS_correct + rel_d*kpos*se_N_pos*(dstpos_n*sp_m_pos + (1-dstpos_n))*l_s*Nsn_H[i][j] +
                                    rel_d*kpos*se_N_pos*(dstpos_p*sp_m_pos + (1-dstpos_p))*l_s*Nsp_H[i][j] +                                         
                                    kpos*se_I_pos*(dstpos_n*sp_m_pos + (1-dstpos_n))*l_s*Isn_H[i][j] +       
                                    kpos*se_I_pos*(dstpos_p*sp_m_pos + (1-dstpos_p))*l_s*Isp_H[i][j];
                 
                 
          DS_incorrect = DS_incorrect + rel_d*kpos*se_N_pos*(dstpos_n*(1-se_m_pos) + (1-dstpos_n))*l_s*Nmn_H[i][j] + 
                                        rel_d*kpos*se_N_pos*(dstpos_p*(1-se_m_pos) + (1-dstpos_p))*l_s*Nmp_H[i][j] +                                        
                                        kpos*se_I_pos*(dstpos_n*(1-se_m_pos) + (1-dstpos_n))*l_s*Imn_H[i][j] +       
                                        kpos*se_I_pos*(dstpos_p*(1-se_m_pos) + (1-dstpos_p))*l_s*Imp_H[i][j];  
                 
          MDR_correct = MDR_correct + rel_d*kpos*se_N_pos*l_m*(dstpos_n*Nmn_H[i][j] + dstpos_p*Nmp_H[i][j]) + 
                                      kpos*se_I_pos*l_m*(dstpos_n*Imn_H[i][j] + dstpos_p*Imp_H[i][j]);  
                                      
          MDR_incorrect = MDR_incorrect + rel_d*kpos*se_N_pos*(1-sp_m_pos)*l_m*(dstpos_n*Nsn_H[i][j] + dstpos_p*Nsp_H[i][j]) + 
                                          kpos*se_I_pos*(1-sp_m_pos)*l_m*(dstpos_n*Isn_H[i][j] + dstpos_p*Isp_H[i][j]);  
                 
          FP = FP + health*kpos*(1-sp_I_pos*sp_N_pos)*l_s*(S_H[i][j]+Lsn_H[i][j]+Lsp_H[i][j]+Lmn_H[i][j]+Lmp_H[i][j]+PTn_H[i][j]+PTp_H[i][j]);
      
      }                  
    }                          
              
    /* Finally assign the things we want to use in R (in addition to the state variables) to yout */
    /* The ode call in R needs to define the number of these and give names */
    yout[0] = Total;
    yout[1] = Total_S;
    yout[2] = Total_Ls;
    yout[3] = Total_Lm;
    yout[4] = Total_L;
    yout[5] = Total_Ns_N + Total_Ns_H;
    yout[6] = Total_Nm_N + Total_Nm_H;
    yout[7] = Total_N;
    yout[8] = Total_Is_N + Total_Is_H;
    yout[9] = Total_Im_N + Total_Im_H;
    yout[10] = Total_I;
    yout[11] = Total_DS;
    yout[12] = Total_MDR;
    yout[13] = FS;
    yout[14] = FM;
    yout[15] = CD4_dist[0];
    yout[16] = CD4_dist[1];
    yout[17] = CD4_dist[2];
    yout[18] = CD4_dist[3];
    yout[19] = CD4_dist[4];
    yout[20] = CD4_dist[5];
    yout[21] = CD4_dist[6];
    yout[22] = CD4_dist_ART[0];
    yout[23] = CD4_dist_ART[1];
    yout[24] = CD4_dist_ART[2];
    yout[25] = CD4_dist_ART[3];
    yout[26] = CD4_dist_ART[4];
    yout[27] = CD4_dist_ART[5];
    yout[28] = CD4_dist_ART[6];
    yout[29] = ART_prop[0];
    yout[30] = ART_prop[1];
    yout[31] = ART_prop[2];
    yout[32] = ART_prop[3];
    yout[33] = ART_prop[4];
    yout[34] = ART_prop[5];
    yout[35] = ART_prop[6];
    yout[36] = Tot_ART;
    yout[37] = ART_need;
    yout[38] = ART_new;
    yout[39] = ART_on;
    yout[40] = TB_deaths_tot;
    yout[41] = TB_cases_neg;
    yout[42] = TB_cases_pos;
    yout[43] = TB_cases_ART;
    yout[44] = births;
    yout[45] = Tot_deaths;
    yout[46] = ART_deaths; 
    yout[47] = Total_PT;
    yout[48] = DS_correct;
    yout[49] = DS_incorrect;
    yout[50] = MDR_correct;
    yout[51] = MDR_incorrect;
    yout[52] = FP;
}



