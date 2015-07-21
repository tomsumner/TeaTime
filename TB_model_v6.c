/* TB model in C code to call from R */

/* Version 6 - converted to single year age bins and aging done via events */

/* To add: Update the way background mort and population adjustent is done
           Let ART mortality be an input (to allow different countries to be modelled)
           Turning off population scaling from 2015
           HIV testing and ART inititation
           TB interventions
           Recode equations to pull out new cases terms - so can adjust care flow more easily in future
*/

/* Can be compiled within R with system("R CMD SHLIB TB_model_v6.c") */
/* This creates a dynamic linked library (.dll) which can be loaded (dyn.load(TB_model_v4.dll)) into R and used as the model fucntion in a call to desolve */

/* C libraries needed */
#include <R.h>
#include <math.h>

/* You need to define number of parameters and forcing functions passed to the model here */
/* These must match number in intializer functions below */
static double parms[45];
static double forc[112];

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
#define sig0 parms[8]      /* proportion developing smear positive disease in 0-4 year olds */
#define sig5 parms[9]      /* proportion developing smear positive disease in 5-9 year olds */
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
#define l_s parms[21]       /* linkage to care for diagnosed DS TB */
#define l_m parms[22]       /* linkage to care for diagnosed MDR TB */
#define d parms[23]         /* relative detection rate of smear negative TB */ 
#define tau_s parms[24]     /* treatment success for DS TB */
#define tau_m parms[25]     /* treatment success for MDR TB */
#define eff_n parms[26]     /* efficacy of first line drugs in treating new MDR cases */
#define eff_p parms[27]     /* efficacy of first line drugs in treating retreatment MDR cases */
#define muN_H parms[28]     /* mortaltiy rate from smear negative TB in HIV+ */
#define muI_H parms[29]     /* mortaltiy rate from smear positive TB in HIV+ */
#define RR1a parms[30]      /* Relative risk for primary disease following HIV infection */
#define RR2a parms[31]      /* Relative risk for primary disease by CD4 */
#define RR1v parms[32]      /* Relative risk for reactivation following HIV infection */
#define RR2v parms[33]      /* Relative risk for reactivation by CD4 */
#define RR1p parms[34]      /* Relative risk of protection against reinfection following HIV infection */
#define RR2p parms[35]      /* Relative risk of protection against reinfection by CD4 */
#define ART_TB1 parms[36]   /* Reduction in TB parameters on ART (<6m) */
#define ART_TB2 parms[37]   /* Reduction in TB parameters on ART (7-12m) */
#define ART_TB3 parms[38]   /* Reduction in TB parameters on ART (>12m) */
#define ART_mort1 parms[39] /* Reduction in TB mortality on ART (<6m) */
#define ART_mort2 parms[40] /* Reduction in TB mortality on ART (7-12m) */
#define ART_mort3 parms[41] /* Reduction in TB mortality on ART (>12m) */
#define BCG_eff parms[42]  /* Efficacy of BCG (reduces primary and reactivation risks) */
#define sig_H parms[43]     /* proportion developing smear positive disease in HIV+ */
#define r_H parms[44]       /* rate of self cure from TB in HIV+ */

/* Forcings are time dependant functions */
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
#define h5 forc[83]
#define h10 forc[84]
#define h15 forc[85]
#define h20 forc[86]
#define h25 forc[87]
#define h30 forc[88]
#define h35 forc[89]
#define h40 forc[90]
#define h45 forc[91]
#define h50 forc[92]
#define h55 forc[93]
#define h60 forc[94]
#define h65 forc[95]
#define h70 forc[96]
#define h75 forc[97]
#define h80 forc[98]

#define Ahigh forc[99] /* ART coverage by CD4 - depends on overall coverage and eligibility threshold */
#define A500 forc[100]
#define A349 forc[101]
#define A249 forc[102]
#define A199 forc[103]
#define A99 forc[104]
#define A50 forc[105]
#define Athresh forc[106] /* ART eligibility threshold in terms of model CD4 categories i.e. 1 means anyone <500 */

#define BCG_cov forc[107] /* BCG coverage */
#define pop_ad forc[108]  /* used to turn on/off adjustment of population to account for disease induced mortality - idea is to turn it off from 2015 onwards*/
#define k forc[109]       /* detection rate */
#define dst_n forc[110]   /* dst coverage in new cases */
#define dst_p forc[111]   /* dst coverage in previoulsy treated cases */

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
    int N=45;
    odeparms(&N, parms);
}

/* ###### FUNCTION TO INITIALIZE FORCINGS PASSED FROM R - if the number of parameters is changed you must update N here ###### */
void forcc(void (* odeforcs)(int *, double *))
{
    int N=112;
    odeforcs(&N, forc);
}

/* ##### EVENTS ARE USED TO ADD BIRTHS AND SHIFT THE POPULATION BY AGE - EQUIVALENT TO THE METHOD OF SCHENZLE ###### */

/* ~~~~~~ NEED TO UPDATE THIS TO DEAL WITH ALL DISEASE STATES ~~~~~~~~~~ */
void event(int *n, double *t, double *y) 
{
  int i;
  double temp[81];
  for (i=0; i<81; i++) temp[i] = y[i];
  y[0] = birth_rate*sumsum(y,0,80)/1000;
  for (i=1; i<80; i++) y[i] = temp[i-1];
  y[80] = y[80] + temp[79];
}





/* ###### DERIVATIVE FUNCTION ###### */

/* ~~~~~~~~~~~~~~~~ NEED TO UPDATE THIS TO HAVE 81 AGE CATEGORIES AND REMOVE AGING AND BIRTHS FROM DERIVATIVES ~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~ ALSO NEED TO ENSURE THAT HIV RATES ARE USED CORRECTLY ACROSS AGES ~~~~~~~~~~~~~~~~~~~~~~~~~ */

void derivsc(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if (ip[0] <2) error("nout should be at least 2");
    
    /* Expand state variables so can use more meaningful names than y and ydot (the variables and rates of change are vectors y and ydot here) */
    /* There are 81 age groups, 7 HIV positive categories and 3 times on ART */
    /* _H = HIV+; _A = HIV+ on ART */
    
    /* These are the variables */
    double S[81];          /* Susceptible */
    double Lsn[81];        /* Latent, DS, new */
    double Lsp[81];        /* Latent, DS, previous */
    double Lmn[81];        /* Latent, DR, new */
    double Lmp[81];        /* Latent, DR, previous */
    double Nsn[81];        /* Smear negative, DS, new */
    double Nsp[81];        /* Smear negative, DS, previous */
    double Nmn[81];        /* Smear negative, DR, new */
    double Nmp[81];        /* Smear negative, DR, previous */
    double Isn[81];        /* Smear positive, DS, new */
    double Isp[81];        /* Smear positive, DS, previous */
    double Imn[81];        /* Smear positive, DR, new */
    double Imp[81];        /* Smear positive, DR, previous */
    double S_H[81][7];
    double Lsn_H[81][7];
    double Lsp_H[81][7];
    double Lmn_H[81][7];
    double Lmp_H[81][7];
    double Nsn_H[81][7];
    double Nsp_H[81][7];
    double Nmn_H[81][7];
    double Nmp_H[81][7];
    double Isn_H[81][7];
    double Isp_H[81][7];
    double Imn_H[81][7];
    double Imp_H[81][7];
    double S_A[81][7][3];
    double Lsn_A[81][7][3];
    double Lsp_A[81][7][3];
    double Lmn_A[81][7][3];
    double Lmp_A[81][7][3];
    double Nsn_A[81][7][3];
    double Nsp_A[81][7][3];
    double Nmn_A[81][7][3];
    double Nmp_A[81][7][3];
    double Isn_A[81][7][3];
    double Isp_A[81][7][3];
    double Imn_A[81][7][3];
    double Imp_A[81][7][3];
    /* These are the rates of change (same names but prefixed with d) */
    double dS[81];
    double dLsn[81];
    double dLsp[81];
    double dLmn[81];
    double dLmp[81];
    double dNsn[81];
    double dNsp[81];
    double dNmn[81];
    double dNmp[81];
    double dIsn[81];
    double dIsp[81];
    double dImn[81];
    double dImp[81];
    double dS_H[81][7];
    double dLsn_H[81][7];
    double dLsp_H[81][7];
    double dLmn_H[81][7];
    double dLmp_H[81][7];
    double dNsn_H[81][7];
    double dNsp_H[81][7];
    double dNmn_H[81][7];
    double dNmp_H[81][7];
    double dIsn_H[81][7];
    double dIsp_H[81][7];
    double dImn_H[81][7];
    double dImp_H[81][7];
    double dS_A[81][7][3];
    double dLsn_A[81][7][3];
    double dLsp_A[81][7][3];
    double dLmn_A[81][7][3];
    double dLmp_A[81][7][3];
    double dNsn_A[81][7][3];
    double dNsp_A[81][7][3];
    double dNmn_A[81][7][3];
    double dNmp_A[81][7][3];
    double dIsn_A[81][7][3];
    double dIsp_A[81][7][3];
    double dImn_A[81][7][3];
    double dImp_A[81][7][3]; 
     
    /* intergers to use as counters */ 
    int i;  
    int j;
    int l;
    int ij;
     
    /* Then need to assign the variablesto the correct bit of y (we do the same for the rates of change after solving the equations) */
    /* SEE IF CAN FIND A NEATER WAY TO DO THIS  - MAYBE JUST HAVE IT ALL IN A LOOP WITH A COUNTER TO KEEP TRACK OF WHERE WE ARE ??*/
     
    /* HIV- */ 
    for (i=0; i<17; i++) S[i] = y[i];             /* S: 0-17 */
    for (i=17; i<34; i++) Lsn[i-17] = y[i];       /* Lsn: 17-33 */
    for (i=34; i<51; i++) Lsp[i-34] = y[i];       /* Lsp: 34-50 */
    for (i=51; i<68; i++) Lmn[i-51] = y[i];       /* Lmn: 51-67 */
    for (i=68; i<85; i++) Lmp[i-68] = y[i];       /* Lmp: 68-84 */
    for (i=85; i<102; i++) Nsn[i-85] = y[i];      /* Nsn: 85-101 */
    for (i=102; i<119; i++) Nsp[i-102] = y[i];    /* Nsp: 102-118 */
    for (i=119; i<136; i++) Nmn[i-119] = y[i];    /* Nmn: 119-135 */
    for (i=136; i<153; i++) Nmp[i-136] = y[i];    /* Nmp: 136-152 */
    for (i=153; i<170; i++) Isn[i-153] = y[i];    /* Isn: 153-169 */
    for (i=170; i<187; i++) Isp[i-170] = y[i];    /* Isp: 170-186 */
    for (i=187; i<204; i++) Imn[i-187] = y[i];    /* Imn: 187-203 */
    for (i=204; i<221; i++) Imp[i-204] = y[i];    /* Imp: 204-220 */  
    /* HIV+ */
    ij = 221;                                     
    for (j=0; j<7; j++){
      for (i=0; i<17; i++){  
        S_H[i][j] = y[ij];
        Lsn_H[i][j] = y[ij+119];
        Lsp_H[i][j] = y[ij+(2*119)];
        Lmn_H[i][j] = y[ij+(3*119)];
        Lmp_H[i][j] = y[ij+(4*119)];  
        Nsn_H[i][j] = y[ij+(5*119)];
        Nsp_H[i][j] = y[ij+(6*119)];
        Nmn_H[i][j] = y[ij+(7*119)];
        Nmp_H[i][j] = y[ij+(8*119)];
        Isn_H[i][j] = y[ij+(9*119)];
        Isp_H[i][j] = y[ij+(10*119)];
        Imn_H[i][j] = y[ij+(11*119)];
        Imp_H[i][j] = y[ij+(12*119)];
        ij = ij+1;
      }
    }
    /* HIV+, on ART */
    ij = 1768;   
    for (l=0; l<3; l++){                          
      for (j=0; j<7; j++){
        for (i=0; i<17; i++){
          S_A[i][j][l] = y[ij];
          Lsn_A[i][j][l] = y[ij+357];
          Lsp_A[i][j][l] = y[ij+(2*357)];
          Lmn_A[i][j][l] = y[ij+(3*357)];
          Lmp_A[i][j][l] = y[ij+(4*357)];  
          Nsn_A[i][j][l] = y[ij+(5*357)];
          Nsp_A[i][j][l] = y[ij+(6*357)];
          Nmn_A[i][j][l] = y[ij+(7*357)];
          Nmp_A[i][j][l] = y[ij+(8*357)];
          Isn_A[i][j][l] = y[ij+(9*357)];
          Isp_A[i][j][l] = y[ij+(10*357)];
          Imn_A[i][j][l] = y[ij+(11*357)];
          Imp_A[i][j][l] = y[ij+(12*357)];
          ij = ij+1;
        }
      }
    }
    
    /* Adjust TB model parameters for age, HIV and ART */
    
    /* Create vectors of aging rates to use in derivatives */
    double age_out[17] = {age1,age1,age1,age1,age1,age1,age1,age1,age1,age1,age1,age1,age1,age1,age1,age1,age2};  /* age_out is lower for final age group because it is wider */
    double age_in[17] = {0.0,age1,age1,age1,age1,age1,age1,age1,age1,age1,age1,age1,age1,age1,age1,age1,age1};    /* age_in is set to zero for first age group as new borns only enter the S class */
    
    /* Create vectors of disease parameters (by age) to use in derivatives - includes BCG effect */
    double a_age[17] = {a0*(BCG_cov*(1-BCG_eff)+(1-BCG_cov)),a5*(BCG_cov*(1-BCG_eff)+(1-BCG_cov)),a10*(BCG_cov*(1-BCG_eff)+(1-BCG_cov)),
                       a_a,a_a,a_a,a_a,a_a,a_a,a_a,a_a,a_a,a_a,a_a,a_a,a_a,a_a};
    double v_age[17] = {v*(BCG_cov*(1-BCG_eff)+(1-BCG_cov)),v*(BCG_cov*(1-BCG_eff)+(1-BCG_cov)),v*(BCG_cov*(1-BCG_eff)+(1-BCG_cov)),
                       v,v,v,v,v,v,v,v,v,v,v,v,v,v};
    double sig_age[17] = {sig0,sig5,sig10,sig_a,sig_a,sig_a,sig_a,sig_a,sig_a,sig_a,sig_a,sig_a,sig_a,sig_a,sig_a,sig_a,sig_a};
    double muN_age[17] = {mu_N0,mu_N,mu_N,mu_N,mu_N,mu_N,mu_N,mu_N,mu_N,mu_N,mu_N,mu_N,mu_N,mu_N,mu_N,mu_N,mu_N};
    double muI_age[17] = {mu_I0,mu_I,mu_I,mu_I,mu_I,mu_I,mu_I,mu_I,mu_I,mu_I,mu_I,mu_I,mu_I,mu_I,mu_I,mu_I,mu_I};

    /* Now adjust parameters for HIV and ART */
    /* HIV mortality rates are passed in directly and ART mortality reduction is implemented in the derivatives */

    double mid_CD4[7] = {500,425,300,225,150,75,25};       /* mid points of CD4 categories */
    double ART_TB[3] = {ART_TB1,ART_TB2,ART_TB3};          /* vector of ART relative risks for TB disease */
    double ART_mort[3] = {ART_mort1,ART_mort2,ART_mort3};  /* vector of ART relative risks for TB mortlaity */
    /* then adjust parameters */
    double a_age_H[17][7];
    double p_H[7];
    double v_age_H[17][7];
    double a_age_A[17][7][3];
    double p_A[7][3];
    double v_age_A[17][7][3];
    for (j=0; j<7; j++){
      p_H[j] = p*RR1p*pow(RR2p,-1*(500-mid_CD4[j])/100);
      for (i=0; i<17; i++){
        a_age_H[i][j] = a_age[i]*RR1a*pow(RR2a,(500-mid_CD4[j])/100);
        v_age_H[i][j] = v_age[i]*RR1v*pow(RR2v,(500-mid_CD4[j])/100);
        for (l=0; l<3; l++){
          a_age_A[i][j][l] = fmax(a_age_H[i][j]*ART_TB[l],a_age[i]);   /* fmax or fmin ensures being on ART can't be better than being HIV- */
          v_age_A[i][j][l] = fmax(v_age_H[i][j]*ART_TB[l],v_age[i]);
          p_A[j][l] = fmin(p_H[j]/ART_TB[l],p);
        }
      }
    }
    
    /* Set up parameters for HIV model - these are taken from AIM */
    /* Looks like only ART mortality depends on country - others can remain hard coded but this will need to be an input (and need to remember to take average by sex) */

    double H_CD4[7][17] = { /* Distribution of new HIV infections (age,CD4) - assume distribution of new child infections mirrors adults */
      {0.643,0.643,0.643,0.643,0.643,0.607,0.607,0.585,0.585,0.552,0.552,0.552,0.552,0.552,0.552,0.552,0.552},
      {0.357,0.357,0.357,0.357,0.357,0.393,0.393,0.415,0.415,0.448,0.448,0.448,0.448,0.448,0.448,0.448,0.448},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    };
    
    /* Have updated these values based on the durations in the AIM manual (rate = 1/duration) as they are different from rates in AIM editor in software */
    /* those are actually risks (i.e. 1-exp(-rate)) */
    double H_prog[8][17] = {  /* Progression through CD4 categories (age, CD4) - has extra row to avoid progression in/out of first/last groups */
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0.298,0.298,0.298,0.117,0.117,0.147,0.147,0.183,0.183,0.213,0.212,0.212,0.212,0.212,0.212,0.212,0.212},
      {0.239,0.239,0.239,0.223,0.223,0.240,0.240,0.355,0.355,0.535,0.535,0.535,0.535,0.535,0.535,0.535,0.535},
      {0.183,0.183,0.183,0.294,0.294,0.452,0.452,0.581,0.581,0.855,0.855,0.855,0.855,0.855,0.855,0.855,0.855},
      {0.183,0.183,0.183,0.508,0.508,1.087,1.087,1.250,1.250,1.818,1.818,1.818,1.818,1.818,1.818,1.818,1.818},
      {0.130,0.130,0.130,0.214,0.214,0.637,0.637,0.676,0.676,0.952,0.952,0.952,0.952,0.952,0.952,0.952,0.952},
      {0.130,0.130,0.130,0.348,0.348,1.449,1.449,1.449,1.449,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    };    

    double H_mort[7][17] = { /* Mortality due to HIV (no ART) (age, CD4) */
      {0.312,0.039,0.039,0.005,0.005,0.004,0.004,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005},
      {0.382,0.048,0.048,0.011,0.011,0.010,0.010,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013},
      {0.466,0.058,0.058,0.026,0.026,0.026,0.026,0.036,0.036,0.032,0.032,0.032,0.032,0.032,0.032,0.032,0.032},
      {0.466,0.058,0.058,0.061,0.061,0.069,0.069,0.096,0.096,0.080,0.080,0.080,0.080,0.080,0.080,0.080,0.080},
      {0.569,0.071,0.071,0.139,0.139,0.185,0.185,0.258,0.258,0.203,0.203,0.203,0.203,0.203,0.203,0.203,0.203},
      {0.569,0.071,0.071,0.321,0.321,0.499,0.499,0.691,0.691,0.513,0.513,0.513,0.513,0.513,0.513,0.513,0.513},
      {0.569,0.071,0.071,0.737,0.737,1.342,1.342,1.851,1.851,1.295,1.295,1.295,1.295,1.295,1.295,1.295,1.295}
    };
    double A_mort[3][7][17] = { /* On ART mortality (age,starting CD4, time on ART)  - this is an average of male and female values weigthed by sex of those on ART  */
      {{0.062,0.008,0.007,0.0050,0.0050,0.0035,0.0035,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050},
       {0.282,0.035,0.033,0.0115,0.0115,0.0095,0.0095,0.0134,0.0134,0.0126,0.0126,0.0126,0.0126,0.0126,0.0126,0.0126,0.0126},
       {0.214,0.027,0.025,0.0264,0.0264,0.0256,0.0256,0.0359,0.0359,0.0319,0.0319,0.0319,0.0319,0.0319,0.0319,0.0319,0.0319},
       {0.214,0.027,0.025,0.0557,0.0557,0.0477,0.0477,0.0495,0.0495,0.0489,0.0489,0.0489,0.0489,0.0489,0.0489,0.0489,0.0489},
       {0.749,0.094,0.088,0.0953,0.0953,0.0792,0.0792,0.0833,0.0833,0.0872,0.0872,0.0872,0.0872,0.0872,0.0872,0.0872,0.0872},
       {0.749,0.094,0.088,0.1553,0.1553,0.1305,0.1305,0.1384,0.1384,0.1496,0.1496,0.1496,0.1496,0.1496,0.1496,0.1496,0.1496},
       {0.749,0.094,0.088,0.3418,0.3418,0.2897,0.2897,0.3094,0.3094,0.3431,0.3431,0.3431,0.3431,0.3431,0.3431,0.3431,0.3431}},
      {{0.147,0.018,0.009,0.0050,0.0050,0.0035,0.0035,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050},
       {0.210,0.026,0.012,0.0115,0.0115,0.0095,0.0095,0.0134,0.0134,0.0126,0.0126,0.0126,0.0126,0.0126,0.0126,0.0126,0.0126},
       {0.199,0.025,0.012,0.0204,0.0204,0.0242,0.0242,0.0267,0.0267,0.0306,0.0306,0.0306,0.0306,0.0306,0.0306,0.0306,0.0306},
       {0.199,0.025,0.012,0.0216,0.0216,0.0278,0.0278,0.0283,0.0283,0.0366,0.0366,0.0366,0.0366,0.0366,0.0366,0.0366,0.0366},
       {0.376,0.047,0.022,0.0272,0.0272,0.0351,0.0351,0.0362,0.0362,0.0481,0.0481,0.0481,0.0481,0.0481,0.0481,0.0481,0.0481},
       {0.376,0.047,0.022,0.0343,0.0343,0.0442,0.0442,0.0461,0.0461,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625},
       {0.376,0.047,0.022,0.0496,0.0496,0.0640,0.0640,0.0675,0.0675,0.0935,0.0935,0.0935,0.0935,0.0935,0.0935,0.0935,0.0935}},
      {{0.062,0.008,0.004,0.0050,0.0050,0.0035,0.0035,0.0050,0.0050,0.0045,0.0045,0.0045,0.0045,0.0045,0.0045,0.0045,0.0045},
       {0.089,0.011,0.005,0.0068,0.0068,0.0081,0.0081,0.0076,0.0076,0.0066,0.0066,0.0066,0.0066,0.0066,0.0066,0.0066,0.0066},
       {0.084,0.011,0.005,0.0073,0.0073,0.0093,0.0093,0.0083,0.0083,0.0076,0.0076,0.0076,0.0076,0.0076,0.0076,0.0076,0.0076},
       {0.084,0.011,0.005,0.0079,0.0079,0.0100,0.0100,0.0091,0.0091,0.0087,0.0087,0.0087,0.0087,0.0087,0.0087,0.0087,0.0087},
       {0.159,0.020,0.009,0.0105,0.0105,0.0134,0.0134,0.0128,0.0128,0.0141,0.0141,0.0141,0.0141,0.0141,0.0141,0.0141,0.0141},
       {0.159,0.020,0.009,0.0139,0.0139,0.0177,0.0177,0.0174,0.0174,0.0209,0.0209,0.0209,0.0209,0.0209,0.0209,0.0209,0.0209},
       {0.159,0.020,0.009,0.0211,0.0211,0.0271,0.0271,0.0275,0.0275,0.0355,0.0355,0.0355,0.0355,0.0355,0.0355,0.0355,0.0355}}
    };      
    double A_prog[4] = {0,2,2,0}; /* Progression through time on ART, 6 monthly time blocks - 0 ensure no progression into first catergory and no progression out of last category*/
    double A_start[3] = {1,0,0};  /* Used to make sure ART initiations are only added to the fist time on ART box */ 
    
    /* sum up various totals */

    /* Use sumsum function to add up HIV- */
    double Total_S = sumsum(S,0,16);                        /* Total susceptible */
    double Total_Ls = sumsum(Lsn,0,16)+sumsum(Lsp,0,16);    /* Total LTBI with drug susceptible (DS) strain */
    double Total_Lm = sumsum(Lmn,0,16)+sumsum(Lmp,0,16);    /* Total LTBI with drug resistant (DR) strain */
    double Total_Ns = sumsum(Nsn,0,16)+sumsum(Nsp,0,16);    /* Total DS smear negative TB */
    double Total_Nm = sumsum(Nmn,0,16)+sumsum(Nmp,0,16);    /* Total DR smear negative TB */
    double Total_Is = sumsum(Isn,0,16)+sumsum(Isp,0,16);    /* Total DS smear positive TB */
    double Total_Im = sumsum(Imn,0,16)+sumsum(Imp,0,16);    /* Total DR smear positive TB */
    
    /* Now loop through HIV and ART and add them in */
    for (j=0; j<7; j++){
      for (i=0; i<17; i++){
        Total_S = Total_S + S_H[i][j];
        Total_Ls = Total_Ls + Lsn_H[i][j]+Lsp_H[i][j];
        Total_Lm = Total_Lm + Lmn_H[i][j]+Lmp_H[i][j];
        Total_Ns = Total_Ns + Nsn_H[i][j]+Nsp_H[i][j];
        Total_Nm = Total_Nm + Nmn_H[i][j]+Nmp_H[i][j];
        Total_Is = Total_Is + Isn_H[i][j]+Isp_H[i][j];
        Total_Im = Total_Im + Imn_H[i][j]+Imp_H[i][j];
        for (l=0; l<3; l++){
          Total_S = Total_S + S_A[i][j][l];
          Total_Ls = Total_Ls + Lsn_A[i][j][l]+Lsp_A[i][j][l];
          Total_Lm = Total_Lm + Lmn_A[i][j][l]+Lmp_A[i][j][l];
          Total_Ns = Total_Ns + Nsn_A[i][j][l]+Nsp_A[i][j][l];
          Total_Nm = Total_Nm + Nmn_A[i][j][l]+Nmp_A[i][j][l];
          Total_Is = Total_Is + Isn_A[i][j][l]+Isp_A[i][j][l];
          Total_Im = Total_Im + Imn_A[i][j][l]+Imp_A[i][j][l];
        }
      }
    }
    double Total_L = Total_Ls + Total_Lm;           /* Total LTBI */
    double Total_N = Total_Ns + Total_Nm;           /* Total smear negative TB */ 
    double Total_I = Total_Is + Total_Im;           /* Total smear positive TB */
    double Total_DS = Total_Ns + Total_Is;          /* Total DS TB */
    double Total_MDR = Total_Nm + Total_Im;         /* Total DR TB */
    double Total = Total_S+Total_L+Total_N+Total_I; /* Total */
    
    /* Mortality calculations and adjustments */
    /* Background survival rates used in the demographic model include death due to disease */
    /* Need to make an adjustment to these rates (based on levels of disease induced mortality) to avoid double counting */
    
    /* Calculate deaths due to disease (TB and HIV) by age group */
    double TB_deaths[17];
    double TB_ART_deaths = 0;
    double HIV_deaths[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double ART_deaths[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    for (i=0; i<17; i++) {
      TB_deaths[i] = (Nsn[i]+Nsp[i]+Nmn[i]+Nmp[i])*muN_age[i] + (Isn[i]+Isp[i]+Imn[i]+Imp[i])*muI_age[i];
      for(j=0; j<7; j++){
        TB_deaths[i] = TB_deaths[i] + (Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j])*muN_H + (Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j])*muI_H;
                                      
        HIV_deaths[i] = HIV_deaths[i]+H_mort[j][i]*(S_H[i][j]+Lsn_H[i][j]+Lsp_H[i][j]+Lmn_H[i][j]+Lmp_H[i][j]+
                        Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j]+Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j]);                              
                          
        for (l=0; l<3; l++){
          TB_deaths[i] = TB_deaths[i] + (Nsn_A[i][j][l]+Nsp_A[i][j][l]+Nmn_A[i][j][l]+Nmp_A[i][j][l])*muN_H*ART_mort[l] +
                                        (Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l])*muI_H*ART_mort[l];
                                        
          TB_ART_deaths = TB_ART_deaths + (Nsn_A[i][j][l]+Nsp_A[i][j][l]+Nmn_A[i][j][l]+Nmp_A[i][j][l])*muN_H*ART_mort[l] + 
                                          (Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l])*muI_H*ART_mort[l];  
              
          ART_deaths[i] = ART_deaths[i]+A_mort[l][j][i]*(S_A[i][j][l]+Lsn_A[i][j][l]+Lsp_A[i][j][l]+Lmn_A[i][j][l]+Lmp_A[i][j][l]+
                                                         Nsn_A[i][j][l]+Nsp_A[i][j][l]+Nmn_A[i][j][l]+Nmp_A[i][j][l]+ 
                                                         Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l]);                             
                              
        }                               
      }

    } 
    double TB_deaths_tot = sumsum(TB_deaths,0,16);
    double ART_deaths_tot = sumsum(ART_deaths,0,16) + TB_ART_deaths; 
    
    /* Calculate total population by age */
    double tot_age[17];
    for (i=0; i<17; i++){
      tot_age[i] = S[i]+Lsn[i]+Lsp[i]+Lmn[i]+Lmp[i]+Nsn[i]+Nsp[i]+Nmn[i]+Nmp[i]+Isn[i]+Isp[i]+Imn[i]+Imp[i];
      for (j=0; j<7; j++){
        tot_age[i] = tot_age[i]+S_H[i][j]+
                     Lsn_H[i][j]+Lsp_H[i][j]+Lmn_H[i][j]+Lmp_H[i][j]+
                     Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j]+
                     Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j];
        for (l=0; l<3; l++){
          tot_age[i] = tot_age[i]+S_A[i][j][l]+Lsn_A[i][j][l]+Lsp_A[i][j][l]+Lmn_A[i][j][l]+Lmp_A[i][j][l]+
                            Nsn_A[i][j][l]+Nsp_A[i][j][l]+Nmn_A[i][j][l]+Nmp_A[i][j][l]+ 
                            Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l];
        }
      }
    }
    /* Then calculate % of pop that die of disease by age - this is multiplied by function pop_ad which turns adjustment off from 2015 onwards */
    double prop_dis_death[17];
    if(pop_ad > 0){
      for (i=0; i<17; i++){
        prop_dis_death[i] = (TB_deaths[i]+HIV_deaths[i]+ART_deaths[i])/tot_age[i];
      }
    }

    /* Then calculate the proportion of people who should die per year by age and reduce it by the disease mortality */
    double m_b[17];
    for (i=0; i<17; i++){
      /*m_b[i] = ((1-forc[i+1])/(1/age_out[i]))-prop_dis_death[i];*/
      m_b[i] = forc[i+1]-prop_dis_death[i];
    }
    
    
    /* Calculate total deaths */
    double Tot_deaths = 0;
    double Tot_deaths_age[17];
    for (i=0; i<17; i++){
      Tot_deaths_age[i] = m_b[i]*tot_age[i] + TB_deaths[i] + HIV_deaths[i] + ART_deaths[i];     
      Tot_deaths = Tot_deaths + Tot_deaths_age[i];
    }
    
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

    for (j=0; j<7; j++){
      for (i=0; i<17; i++){
        CD4_dist[j] = CD4_dist[j]+S_H[i][j]+Lsn_H[i][j]+Lsp_H[i][j]+Lmn_H[i][j]+Lmp_H[i][j]+
                      Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j]+ 
                      Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j];
        CD4_deaths[j] = CD4_deaths[j]+
                        H_mort[j][i]*(S_H[i][j]+Lsn_H[i][j]+Lsp_H[i][j]+Lmn_H[i][j]+Lmp_H[i][j]+
                        Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j]+ 
                        Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j]);
        for (l=0; l<3; l++){
          CD4_dist_ART[j] = CD4_dist_ART[j]+S_A[i][j][l]+Lsn_A[i][j][l]+Lsp_A[i][j][l]+Lmn_A[i][j][l]+Lmp_A[i][j][l]+
                            Nsn_A[i][j][l]+Nsp_A[i][j][l]+Nmn_A[i][j][l]+Nmp_A[i][j][l]+ 
                            Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l];         
                      
        }            
      } 
      Tot_ART = Tot_ART + CD4_dist_ART[j];  /* sum up number currently on ART */
      ART_on = ART_on + (CD4_dist[j] + CD4_dist_ART[j])*forc[j+35]; /* number who should be on ART - HIV+ population times coverage (by CD4) */
    }
    ART_new = fmax(0,ART_on - (Tot_ART - ART_deaths_tot));   /* number who need to start is number who should be on minus those already on plus those on ART who will die in current time */ 
    
     /*for (j=0; j<7; j++) ART_prop[j] = 0; /* THIS LINE IS HERE TO ALLOW ART TO BE TURNED OFF */
    
    /* Then work out where these should go by CD4 - based on proportion of eligible population in CD4 group and proportion of deaths which occuring in CD4 group */
    for (j=Athresh; j<7; j++) {
      ART_el = ART_el + CD4_dist[j];
      ART_el_deaths = ART_el_deaths + CD4_deaths[j];
      ART_need = ART_need + CD4_dist[j] + CD4_dist_ART[j];
    }
    if (ART_el > 0){
      for (j=Athresh; j<7; j++) {
        if (CD4_dist[j] > 0) {
          ART_prop[j] = (((CD4_dist[j]/ART_el)+(CD4_deaths[j]/ART_el_deaths))/2)*(ART_new/CD4_dist[j]); /* applies weighting and size of CD4 group to work out % of CD4 group that should move */
           /* ART_prop[j] = (CD4_dist[j]/ART_el)*(ART_new/CD4_dist[j]); */
       }
      }
    }
    
    /* Force of infection */
    double FS = beta*(Total_Ns*rel_inf + Total_Is)/Total; 
    double FM = fit_cost*beta*(Total_Nm*rel_inf + Total_Im)/Total; 
    
    /* Variables to store numbers of new cases */
    double TB_cases_neg_age[17];
    double TB_cases_neg = 0;
    double TB_cases_pos_age[17][7];
    double TB_cases_pos = 0;
    double TB_cases_ART_age[17][7][3];
    double TB_cases_ART = 0;
        
    /* Derivatives */ 
 
    /* HIV-: loop through ages*/ 
    
    double births = birth_rate*Total/1000;
    
    /* Susceptible - note this is done separately to deal with births */
    dS[0] = births - age_out[0]*S[0] - (FS + FM)*S[0] - forc[18]*S[0] - m_b[0]*S[0];
    for (i=1; i<17; i++) dS[i] = age_in[i]*S[i-1] - age_out[i]*S[i] - (FS + FM)*S[i] - forc[i+18]*S[i] - m_b[i]*S[i];
    /* Other TB states */
    for (i=0; i<17; i++){
      
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
      
      /* Latent, ds, naive */
      dLsn[i] = age_in[i]*Lsn[i-1] - age_out[i]*Lsn[i] - m_b[i]*Lsn[i] + S_to_Lsn + Lmn_to_Lsn + r*(Isn[i] + Nsn[i]) - Lsn_to_Lmn - Lsn_to_Nsn - Lsn_to_Isn  - Lsn_to_Nmn - Lsn_to_Imn -
                forc[i+18]*Lsn[i];                   
      /* Latent, ds, prev */
      dLsp[i] = age_in[i]*Lsp[i-1] - age_out[i]*Lsp[i] - m_b[i]*Lsp[i] + Lmp_to_Lsp + r*(Isp[i] + Nsp[i]) - Lsp_to_Lmp - Lsp_to_Nsp - Lsp_to_Isp  - Lsp_to_Nmp - Lsp_to_Imp + 
                k*l_s*(1-e)*tau_s*(Isn[i]+Isp[i]) + k*l_s*(1-e)*tau_s*d*(Nsn[i]+Nsp[i]) - /* Care */
                forc[i+18]*Lsp[i]; /* HIV */          
      /* Latent, mdr, naive */ 
      dLmn[i] = age_in[i]*Lmn[i-1] - age_out[i]*Lmn[i] - m_b[i]*Lmn[i] + S_to_Lmn + Lsn_to_Lmn + r*(Imn[i] + Nmn[i]) - Lmn_to_Lsn - Lmn_to_Nsn - Lmn_to_Isn - Lmn_to_Nmn - Lmn_to_Imn  -                
                forc[i+18]*Lmn[i]; /* HIV */             
      /* Latent, mdr, prev */
      dLmp[i] = age_in[i]*Lmp[i-1] - age_out[i]*Lmp[i] - m_b[i]*Lmp[i] + Lsp_to_Lmp + r*(Imp[i] + Nmp[i]) - Lmp_to_Lsp - Lmp_to_Nsp - Lmp_to_Isp - Lmp_to_Nmp - Lmp_to_Imp +
                k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(Imn[i]+d*Nmn[i]) + k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(Imp[i]+d*Nmp[i]) - /* Care */
                forc[i+18]*Lmp[i]; /* HIV */          
      /* Smear neg, ds, new */
      dNsn[i] = age_in[i]*Nsn[i-1] - age_out[i]*Nsn[i] - m_b[i]*Nsn[i] + S_to_Nsn + Lsn_to_Nsn + Lmn_to_Nsn - (theta + r + muN_age[i])*Nsn[i] -
                k*l_s*d*Nsn[i] - /* care */
                forc[i+18]*Nsn[i]; /* HIV */
      /* Smear neg, ds, prev */                             
      dNsp[i] = age_in[i]*Nsp[i-1] - age_out[i]*Nsp[i] - m_b[i]*Nsp[i] + Lsp_to_Nsp + Lmp_to_Nsp - (theta + r + muN_age[i])*Nsp[i] -
                (k*l_s*d*(1-e)*tau_s + k*l_s*d*e)*Nsp[i] + k*l_s*d*(1-e)*(1-tau_s)*Nsn[i] - /* Care */
                forc[i+18]*Nsp[i]; /* HIV */
      /* Smear neg, mdr, new */
      dNmn[i] = age_in[i]*Nmn[i-1] - age_out[i]*Nmn[i] - m_b[i]*Nmn[i] + S_to_Nmn + Lsn_to_Nmn + Lmn_to_Nmn - (theta + r + muN_age[i])*Nmn[i] -
                k*(l_m*d*dst_n + l_s*d*(1-dst_n))*Nmn[i] - /* Care */
                forc[i+18]*Nmn[i]; /* HIV */        
      /* Smear neg, mdr, prev */
      dNmp[i] = age_in[i]*Nmp[i-1] - age_out[i]*Nmp[i] - m_b[i]*Nmp[i] + Lsp_to_Nmp + Lmp_to_Nmp - (theta + r + muN_age[i])*Nmp[i] -
                k*(l_m*d*dst_p*tau_m + l_s*d*(1-dst_p)*tau_s*eff_p)*Nmp[i] + k*l_s*d*e*(Nsn[i]+Nsp[i]) + (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*Nmn[i] - /* Care */
                forc[i+18]*Nmp[i]; /* HIV */
      /* Smear pos, ds, new */
      dIsn[i] = age_in[i]*Isn[i-1] - age_out[i]*Isn[i] - m_b[i]*Isn[i] + S_to_Isn + Lsn_to_Isn + Lmn_to_Isn + theta*Nsn[i] - (r + muI_age[i])*Isn[i] - 
                k*l_s*Isn[i] - /* Care */
                forc[i+18]*Isn[i]; /* HIV */
      /* Smear pos, ds, prev */
      dIsp[i] = age_in[i]*Isp[i-1] - age_out[i]*Isp[i] - m_b[i]*Isp[i] + Lsp_to_Isp + Lmp_to_Isp + theta*Nsp[i] - (r + muI_age[i])*Isp[i] - 
                (k*l_s*(1-e)*tau_s + k*l_s*e)*Isp[i] + k*l_s*(1-e)*(1-tau_s)*Isn[i] - /* Care */
                forc[i+18]*Isp[i]; /* HIV */
      /* Smear pos, mdr, new */
      dImn[i] = age_in[i]*Imn[i-1] - age_out[i]*Imn[i] - m_b[i]*Imn[i] + S_to_Imn + Lsn_to_Imn + Lmn_to_Imn + theta*Nmn[i] - (r + muI_age[i])*Imn[i] -
                (k*l_m*dst_n + k*l_s*(1-dst_n))*Imn[i] - /* Care */ 
                forc[i+18]*Imn[i]; /* HIV */
      /* Smear pos, mdr, prev */
      dImp[i] = age_in[i]*Imp[i-1] - age_out[i]*Imp[i] - m_b[i]*Imp[i] + Lsp_to_Imp + Lmp_to_Imp + theta*Nmp[i] - (r + muI_age[i])*Imp[i] -
                k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*Imp[i] + (k*l_m*dst_n*(1-tau_m) + k*l_s*(1-dst_n)*(1-(tau_s*eff_n)))*Imn[i] + k*l_s*e*(Isn[i]+Isp[i]) -
                forc[i+18]*Imp[i];
                        
      /* sum up new HIV- cases */           
      TB_cases_neg_age[i] = (v_age[i]*(1-sig_age[i]) + FS*a_age[i]*(1-p)*(1-sig_age[i]))*Lsn[i] + FS*a_age[i]*(1-sig_age[i])*(S[i] + (1-p)*Lmn[i]) +          /*sneg,sus,new*/
                            (v_age[i]*(1-sig_age[i]) + FS*a_age[i]*(1-p)*(1-sig_age[i]))*Lsp[i] + FS*a_age[i]*(1-sig_age[i])*(1-p)*Lmp[i] +                   /*sneg,sus,prev*/
                            (v_age[i]*(1-sig_age[i]) + FM*a_age[i]*(1-p)*(1-sig_age[i]))*Lmn[i] + FM*a_age[i]*(1-sig_age[i])*(S[i] + (1-p)*Lsn[i]) +          /*sneg,mdr,new*/
                            (v_age[i]*(1-sig_age[i]) + FM*a_age[i]*(1-p)*(1-sig_age[i]))*Lmp[i] + FM*a_age[i]*(1-sig_age[i])*(1-p)*Lsp[i] +                   /*sneg,mdr,prev*/
                            (v_age[i]*sig_age[i] + FS*a_age[i]*sig_age[i]*(1-p))*Lsn[i] + FS*a_age[i]*sig_age[i]*S[i] + FS*a_age[i]*(1-p)*sig_age[i]*Lmn[i] + /*spos,sus,new*/
                            (v_age[i]*sig_age[i] + FS*a_age[i]*sig_age[i]*(1-p))*Lsp[i] + FS*a_age[i]*sig_age[i]*(1-p)*Lmp[i] +                               /*spos,sus,prev*/
                            (v_age[i]*sig_age[i] + FM*a_age[i]*sig_age[i]*(1-p))*Lmn[i] + FM*a_age[i]*sig_age[i]*S[i] + FM*a_age[i]*(1-p)*sig_age[i]*Lsn[i] + /*spos,mdr,new*/
                            (v_age[i]*sig_age[i] + FM*a_age[i]*sig_age[i]*(1-p))*Lmp[i] + FM*a_age[i]*sig_age[i]*(1-p)*Lsp[i];                                /*spos,mdr,prev*/
    
      TB_cases_neg = TB_cases_neg + TB_cases_neg_age[i];
    
      /* HIV+: Loop through CD4 categories */
    
      for (j=0; j<7; j++){      /* CD4 */

        dS_H[i][j] = age_in[i]*S_H[i-1][j] - age_out[i]*S_H[i][j] - m_b[i]*S_H[i][j] - 
                     (FS + FM)*S_H[i][j] +
                     forc[i+18]*H_CD4[j][i]*S[i] - H_prog[j+1][i]*S_H[i][j] + H_prog[j][i]*S_H[i][j-1] - 
                     H_mort[j][i]*S_H[i][j] - ART_prop[j]*S_H[i][j];

        dLsn_H[i][j] = age_in[i]*Lsn_H[i-1][j] - age_out[i]*Lsn_H[i][j] - m_b[i]*Lsn_H[i][j] +
                       FS*((1-a_age_H[i][j])*S_H[i][j] + (1-a_age_H[i][j])*(1-p_H[j])*(1-g)*Lmn_H[i][j]) -
                       (v_age_H[i][j] + FS*a_age_H[i][j]*(1-p_H[j]) + FM*a_age_H[i][j]*(1-p_H[j]) + FM*(1-a_age_H[i][j])*(1-p_H[j])*g)*Lsn_H[i][j] +
                       r*(Isn_H[i][j] + Nsn_H[i][j]) +
                       forc[i+18]*H_CD4[j][i]*Lsn[i] - H_prog[j+1][i]*Lsn_H[i][j] + H_prog[j][i]*Lsn_H[i][j-1] - 
                       H_mort[j][i]*Lsn_H[i][j] - ART_prop[j]*Lsn_H[i][j];
        
        dLsp_H[i][j] = age_in[i]*Lsp_H[i-1][j] - age_out[i]*Lsp_H[i][j] - m_b[i]*Lsp_H[i][j] +
                       FS*(1-a_age_H[i][j])*(1-p_H[j])*(1-g)*Lmp_H[i][j] -
                       (v_age_H[i][j] + FS*a_age_H[i][j]*(1-p_H[j]) + FM*a_age_H[i][j]*(1-p_H[j]) + FM*(1-a_age_H[i][j])*(1-p_H[j])*g)*Lsp_H[i][j] +
                       r_H*(Isp_H[i][j]+Nsp_H[i][j]) +
                       k*l_s*(1-e)*tau_s*(Isn_H[i][j]+Isp_H[i][j]) + 
                       k*l_s*(1-e)*tau_s*d*(Nsn_H[i][j]+Nsp_H[i][j]) +
                       forc[i+18]*H_CD4[j][i]*Lsp[i] - H_prog[j+1][i]*Lsp_H[i][j] + H_prog[j][i]*Lsp_H[i][j-1] - 
                       H_mort[j][i]*Lsp_H[i][j] - ART_prop[j]*Lsp_H[i][j];
                             
        dLmn_H[i][j] = age_in[i]*Lmn_H[i-1][j] - age_out[i]*Lmn_H[i][j] - m_b[i]*Lmn_H[i][j] +
                       FM*((1-a_age_H[i][j])*S_H[i][j] + (1-a_age_H[i][j])*(1-p_H[j])*g*Lsn_H[i][j]) -
                       (v_age_H[i][j] + FM*a_age_H[i][j]*(1-p_H[j]) + FS*a_age_H[i][j]*(1-p_H[j]) + FS*(1-a_age_H[i][j])*(1-p_H[j])*(1-g))*Lmn_H[i][j] +
                       r_H*(Imn_H[i][j]+Nmn_H[i][j]) +
                       forc[i+18]*H_CD4[j][i]*Lmn[i] - H_prog[j+1][i]*Lmn_H[i][j] + H_prog[j][i]*Lmn_H[i][j-1] - 
                       H_mort[j][i]*Lmn_H[i][j] - ART_prop[j]*Lmn_H[i][j];
                     
        dLmp_H[i][j] = age_in[i]*Lmp_H[i-1][j] - age_out[i]*Lmp_H[i][j] - m_b[i]*Lmp_H[i][j] + 
                       FM*(1-a_age_H[i][j])*(1-p_H[j])*g*Lsp_H[i][j] -
                       (v_age_H[i][j] + FM*a_age_H[i][j]*(1-p_H[j]) + FS*a_age_H[i][j]*(1-p_H[j]) + FS*(1-a_age_H[i][j])*(1-p_H[j])*(1-g))*Lmp_H[i][j] +
                       r_H*(Imp_H[i][j]+Nmp_H[i][j]) +
                       k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(Imn_H[i][j]+d*Nmn_H[i][j]) +
                       k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(Imp_H[i][j]+d*Nmp_H[i][j]) +
                       forc[i+18]*H_CD4[j][i]*Lmp[i] - H_prog[j+1][i]*Lmp_H[i][j] + H_prog[j][i]*Lmp_H[i][j-1] - 
                       H_mort[j][i]*Lmp_H[i][j] - ART_prop[j]*Lmp_H[i][j];
           
        dNsn_H[i][j] = age_in[i]*Nsn_H[i-1][j] - age_out[i]*Nsn_H[i][j] - m_b[i]*Nsn_H[i][j] +
                       (v_age_H[i][j]*(1-sig_H) + FS*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lsn_H[i][j] +
                       FS*a_age_H[i][j]*(1-sig_H)*(S_H[i][j] + (1-p_H[j])*Lmn_H[i][j]) -
                       (theta + r_H + k*l_s*d + muN_H)*Nsn_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Nsn[i] - H_prog[j+1][i]*Nsn_H[i][j] + H_prog[j][i]*Nsn_H[i][j-1] - 
                       H_mort[j][i]*Nsn_H[i][j] - ART_prop[j]*Nsn_H[i][j];
                                                
        dNsp_H[i][j] = age_in[i]*Nsp_H[i-1][j] - age_out[i]*Nsp_H[i][j] - m_b[i]*Nsp_H[i][j] + 
                       (v_age_H[i][j]*(1-sig_H) + FS*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lsp_H[i][j] +
                       FS*a_age_H[i][j]*(1-sig_H)*(1-p_H[j])*Lmp_H[i][j] -
                       (theta + r_H + k*l_s*d*(1-e)*tau_s + k*l_s*d*e + muN_H)*Nsp_H[i][j] +
                       k*l_s*d*(1-e)*(1-tau_s)*Nsn_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Nsp[i] - H_prog[j+1][i]*Nsp_H[i][j] + H_prog[j][i]*Nsp_H[i][j-1] - 
                       H_mort[j][i]*Nsp_H[i][j] - ART_prop[j]*Nsp_H[i][j];
                       
        dNmn_H[i][j] = age_in[i]*Nmn_H[i-1][j] - age_out[i]*Nmn_H[i][j] - m_b[i]*Nmn_H[i][j] +
                       (v_age_H[i][j]*(1-sig_H) + FM*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lmn_H[i][j] +
                       FM*a_age_H[i][j]*(1-sig_H)*(S_H[i][j] + (1-p_H[j])*Lsn_H[i][j]) -
                       (theta + r_H + k*l_m*d*dst_n + k*l_s*d*(1-dst_n) + muN_H)*Nmn_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Nmn[i] - H_prog[j+1][i]*Nmn_H[i][j] + H_prog[j][i]*Nmn_H[i][j-1] - 
                       H_mort[j][i]*Nmn_H[i][j] - ART_prop[j]*Nmn_H[i][j];
                       
        dNmp_H[i][j] = age_in[i]*Nmp_H[i-1][j] - age_out[i]*Nmp_H[i][j] - m_b[i]*Nmp_H[i][j] +
                       (v_age_H[i][j]*(1-sig_H) + FM*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lmp_H[i][j] +
                       FM*a_age_H[i][j]*(1-sig_H)*(1-p_H[j])*Lsp_H[i][j] -
                       (theta + r_H + k*l_m*d*dst_p*tau_m + k*l_s*d*(1-dst_p)*tau_s*eff_p + muN_H)*Nmp_H[i][j] +
                       k*l_s*d*e*(Nsn_H[i][j]+Nsp_H[i][j]) +
                       (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*Nmn_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Nmp[i] - H_prog[j+1][i]*Nmp_H[i][j] + H_prog[j][i]*Nmp_H[i][j-1] - 
                       H_mort[j][i]*Nmp_H[i][j] - ART_prop[j]*Nmp_H[i][j];
           
        dIsn_H[i][j] = age_in[i]*Isn_H[i-1][j] - age_out[i]*Isn_H[i][j] - m_b[i]*Isn_H[i][j] +
                       (v_age_H[i][j]*sig_H + FS*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lsn_H[i][j] +
                       FS*a_age_H[i][j]*sig_H*S_H[i][j] +
                       FS*a_age_H[i][j]*(1-p_H[j])*sig_H*Lmn_H[i][j] +
                       theta*Nsn_H[i][j] -
                       (r_H + k*l_s + muI_H)*Isn_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Isn[i] - H_prog[j+1][i]*Isn_H[i][j] + H_prog[j][i]*Isn_H[i][j-1] - 
                       H_mort[j][i]*Isn_H[i][j] - ART_prop[j]*Isn_H[i][j];

        dIsp_H[i][j] = age_in[i]*Isp_H[i-1][j] - age_out[i]*Isp_H[i][j] - m_b[i]*Isp_H[i][j] +
                       (v_age_H[i][j]*sig_H + FS*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lsp_H[i][j] +
                       FS*a_age_H[i][j]*sig_H*(1-p_H[j])*Lmp_H[i][j] +
                       theta*Nsp_H[i][j] +
                       k*l_s*(1-e)*(1-tau_s)*Isn_H[i][j] -
                       (r_H + k*l_s*(1-e)*tau_s + k*l_s*e + muI_H)*Isp_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Isp[i] - H_prog[j+1][i]*Isp_H[i][j] + H_prog[j][i]*Isp_H[i][j-1] - 
                       H_mort[j][i]*Isp_H[i][j] - ART_prop[j]*Isp_H[i][j];
                                   
        dImn_H[i][j] = age_in[i]*Imn_H[i-1][j] - age_out[i]*Imn_H[i][j] - m_b[i]*Imn_H[i][j] +
                       (v_age_H[i][j]*sig_H + FM*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lmn_H[i][j] +
                       FM*a_age_H[i][j]*sig_H*S_H[i][j] +
                       FM*a_age_H[i][j]*(1-p_H[j])*sig_H*Lsn_H[i][j] +
                       theta*Nmn_H[i][j] -
                       (r_H + k*l_m*dst_n + k*l_s*(1-dst_n) + muI_H)*Imn_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Imn[i] - H_prog[j+1][i]*Imn_H[i][j] + H_prog[j][i]*Imn_H[i][j-1] - 
                       H_mort[j][i]*Imn_H[i][j] - ART_prop[j]*Imn_H[i][j];

        dImp_H[i][j] = age_in[i]*Imp_H[i-1][j] - age_out[i]*Imp_H[i][j] - m_b[i]*Imp_H[i][j] +
                       (v_age_H[i][j]*sig_H + FM*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lmp_H[i][j] +
                       FM*a_age_H[i][j]*sig_H*(1-p_H[j])*Lsp_H[i][j] +
                       theta*Nmp_H[i][j] +
                       (k*l_m*dst_n*(1-tau_m) + k*l_s*(1-dst_n)*(1-(tau_s*eff_n)))*Imn_H[i][j] +
                       k*l_s*e*(Isn_H[i][j]+Isp_H[i][j]) -
                       (r_H + k*l_m*dst_p*tau_m + k*l_s*(1-dst_p)*tau_s*eff_p + muI_H)*Imp_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Imp[i] - H_prog[j+1][i]*Imp_H[i][j] + H_prog[j][i]*Imp_H[i][j-1] - 
                       H_mort[j][i]*Imp_H[i][j] - ART_prop[j]*Imp_H[i][j];
                       
       TB_cases_pos_age[i][j] = (v_age_H[i][j]*(1-sig_H) + FS*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lsn_H[i][j] + FS*a_age_H[i][j]*(1-sig_H)*(S_H[i][j] + (1-p_H[j])*Lmn_H[i][j]) + 
                                (v_age_H[i][j]*(1-sig_H) + FS*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lsp_H[i][j] + FS*a_age_H[i][j]*(1-sig_H)*(1-p_H[j])*Lmp_H[i][j] +
                                (v_age_H[i][j]*(1-sig_H) + FM*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lmn_H[i][j] + FM*a_age_H[i][j]*(1-sig_H)*(S_H[i][j] + (1-p_H[j])*Lsn_H[i][j]) +         
                                (v_age_H[i][j]*(1-sig_H) + FM*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lmp_H[i][j] + FM*a_age_H[i][j]*(1-sig_H)*(1-p_H[j])*Lsp_H[i][j] +
                                (v_age_H[i][j]*sig_H + FS*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lsn_H[i][j] + FS*a_age_H[i][j]*sig_H*S_H[i][j] + FS*a_age_H[i][j]*(1-p_H[j])*sig_H*Lmn_H[i][j] +
                                (v_age_H[i][j]*sig_H + FS*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lsp_H[i][j] + FS*a_age_H[i][j]*sig_H*(1-p_H[j])*Lmp_H[i][j] +
                                (v_age_H[i][j]*sig_H + FM*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lmn_H[i][j] + FM*a_age_H[i][j]*sig_H*S_H[i][j] + FM*a_age_H[i][j]*(1-p_H[j])*sig_H*Lsn_H[i][j] +
                                (v_age_H[i][j]*sig_H + FM*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lmp_H[i][j] + FM*a_age_H[i][j]*sig_H*(1-p_H[j])*Lsp_H[i][j];
                       
       TB_cases_pos = TB_cases_pos + TB_cases_pos_age[i][j];
       
       /* HIV+ on ART: loop through time on ART, CD4 at initiation, age - NEED TO ADD IN IMPACT OF ART ON TB PARAMETERS */
        for (l=0; l<3; l++){
        
          dS_A[i][j][l] = age_in[i]*S_A[i-1][j][l] - age_out[i]*S_A[i][j][l] - m_b[i]*S_A[i][j][l] - (FS + FM)*S_A[i][j][l] +
                          ART_prop[j]*A_start[l]*S_H[i][j] + A_prog[l]*S_A[i][j][l-1] - A_prog[l+1]*S_A[i][j][l] - A_mort[l][j][i]*S_A[i][j][l];

          dLsn_A[i][j][l] = age_in[i]*Lsn_A[i-1][j][l] - age_out[i]*Lsn_A[i][j][l] - m_b[i]*Lsn_A[i][j][l] +
                            FS*((1-a_age_A[i][j][l])*S_A[i][j][l] + (1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g)*Lmn_A[i][j][l]) -
                            (v_age_A[i][j][l] + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g)*Lsn_A[i][j][l] +
                            r_H*(Isn_A[i][j][l] + Nsn_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lsn_H[i][j] + A_prog[l]*Lsn_A[i][j][l-1] - A_prog[l+1]*Lsn_A[i][j][l] - A_mort[l][j][i]*Lsn_A[i][j][l];
        
          dLsp_A[i][j][l] = age_in[i]*Lsp_A[i-1][j][l] - age_out[i]*Lsp_A[i][j][l] - m_b[i]*Lsp_A[i][j][l] +
                            FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g)*Lmp_A[i][j][l] -
                            (v_age_A[i][j][l] + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g)*Lsp_A[i][j][l] +
                            r_H*(Isp_A[i][j][l]+Nsp_A[i][j][l]) +
                            k*l_s*(1-e)*tau_s*(Isn_A[i][j][l]+Isp_A[i][j][l]) + 
                            k*l_s*(1-e)*tau_s*d*(Nsn_A[i][j][l]+Nsp_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lsp_H[i][j] + A_prog[l]*Lsp_A[i][j][l-1] - A_prog[l+1]*Lsp_A[i][j][l] - A_mort[l][j][i]*Lsp_A[i][j][l];
                             
          dLmn_A[i][j][l] = age_in[i]*Lmn_A[i-1][j][l] - age_out[i]*Lmn_A[i][j][l] - m_b[i]*Lmn_A[i][j][l] +
                            FM*((1-a_age_A[i][j][l])*S_A[i][j][l] + (1-a_age_A[i][j][l])*(1-p_A[j][l])*g*Lsn_A[i][j][l]) -
                            (v_age_A[i][j][l] + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g))*Lmn_A[i][j][l] +
                            r_H*(Imn_A[i][j][l]+Nmn_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lmn_H[i][j] + A_prog[l]*Lmn_A[i][j][l-1] - A_prog[l+1]*Lmn_A[i][j][l] - A_mort[l][j][i]*Lmn_A[i][j][l];
                     
          dLmp_A[i][j][l] = age_in[i]*Lmp_A[i-1][j][l] - age_out[i]*Lmp_A[i][j][l] - m_b[i]*Lmp_A[i][j][l] + 
                            FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g*Lsp_A[i][j][l] -
                            (v_age_A[i][j][l] + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g))*Lmp_A[i][j][l] +
                            r_H*(Imp_A[i][j][l]+Nmp_A[i][j][l]) +
                            k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(Imn_A[i][j][l]+d*Nmn_A[i][j][l]) +
                            k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(Imp_A[i][j][l]+d*Nmp_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lmp_H[i][j] + A_prog[l]*Lmp_A[i][j][l-1] - A_prog[l+1]*Lmp_A[i][j][l] - A_mort[l][j][i]*Lmp_A[i][j][l];
           
          dNsn_A[i][j][l] = age_in[i]*Nsn_A[i-1][j][l]- age_out[i]*Nsn_A[i][j][l] - m_b[i]*Nsn_A[i][j][l]+
                            (v_age_A[i][j][l]*(1-sig_H) + FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lsn_A[i][j][l] +
                            FS*a_age_A[i][j][l]*(1-sig_H)*(S_A[i][j][l] + (1-p_A[j][l])*Lmn_A[i][j][l]) -
                            (theta + r_H + k*l_s*d + muN_H*ART_mort[l])*Nsn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Nsn_H[i][j] + A_prog[l]*Nsn_A[i][j][l-1] - A_prog[l+1]*Nsn_A[i][j][l] - A_mort[l][j][i]*Nsn_A[i][j][l];
                                                
          dNsp_A[i][j][l] = age_in[i]*Nsp_A[i-1][j][l] - age_out[i]*Nsp_A[i][j][l] - m_b[i]*Nsp_A[i][j][l] + 
                            (v_age_A[i][j][l]*(1-sig_H) + FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lsp_A[i][j][l] +
                            FS*a_age_A[i][j][l]*(1-sig_H)*(1-p_A[j][l])*Lmp_A[i][j][l] -
                            (theta + r_H + k*l_s*d*(1-e)*tau_s + k*l_s*d*e + muN_H*ART_mort[l])*Nsp_A[i][j][l] +
                            k*l_s*d*(1-e)*(1-tau_s)*Nsn_A[i][j][l] +      
                            ART_prop[j]*A_start[l]*Nsp_H[i][j] + A_prog[l]*Nsp_A[i][j][l-1] - A_prog[l+1]*Nsp_A[i][j][l] - A_mort[l][j][i]*Nsp_A[i][j][l];
                       
          dNmn_A[i][j][l] = age_in[i]*Nmn_A[i-1][j][l] - age_out[i]*Nmn_A[i][j][l] - m_b[i]*Nmn_A[i][j][l] +
                            (v_age_A[i][j][l]*(1-sig_H) + FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lmn_A[i][j][l] +
                            FM*a_age_A[i][j][l]*(1-sig_H)*(S_A[i][j][l] + (1-p_A[j][l])*Lsn_A[i][j][l]) -
                            (theta + r_H + k*l_m*d*dst_n + k*l_s*d*(1-dst_n) + muN_H*ART_mort[l])*Nmn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Nmn_H[i][j] + A_prog[l]*Nmn_A[i][j][l-1] - A_prog[l+1]*Nmn_A[i][j][l] - A_mort[l][j][i]*Nmn_A[i][j][l];
                       
          dNmp_A[i][j][l] = age_in[i]*Nmp_A[i-1][j][l] - age_out[i]*Nmp_A[i][j][l] - m_b[i]*Nmp_A[i][j][l] +
                            (v_age_A[i][j][l]*(1-sig_H) + FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lmp_A[i][j][l] +
                            FM*a_age_A[i][j][l]*(1-sig_H)*(1-p_A[j][l])*Lsp_A[i][j][l] -
                            (theta + r_H + k*l_m*d*dst_p*tau_m + k*l_s*d*(1-dst_p)*tau_s*eff_p + muN_H*ART_mort[l])*Nmp_A[i][j][l] +
                            k*l_s*d*e*(Nsn_A[i][j][l]+Nsp_A[i][j][l]) +
                            (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*Nmn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Nmp_H[i][j] + A_prog[l]*Nmp_A[i][j][l-1] - A_prog[l+1]*Nmp_A[i][j][l] - A_mort[l][j][i]*Nmp_A[i][j][l];
           
          dIsn_A[i][j][l] = age_in[i]*Isn_A[i-1][j][l] - age_out[i]*Isn_A[i][j][l] - m_b[i]*Isn_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_H + FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lsn_A[i][j][l] +
                            FS*a_age_A[i][j][l]*sig_H*S_A[i][j][l] +
                            FS*a_age_A[i][j][l]*(1-p_A[j][l])*sig_H*Lmn_A[i][j][l] +
                            theta*Nsn_A[i][j][l] -
                            (r_H + k*l_s + muI_H*ART_mort[l])*Isn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Isn_H[i][j] + A_prog[l]*Isn_A[i][j][l-1] - A_prog[l+1]*Isn_A[i][j][l] - A_mort[l][j][i]*Isn_A[i][j][l];

          dIsp_A[i][j][l] = age_in[i]*Isp_A[i-1][j][l] - age_out[i]*Isp_A[i][j][l] - m_b[i]*Isp_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_H + FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lsp_A[i][j][l] +
                            FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l])*Lmp_A[i][j][l] +
                            theta*Nsp_A[i][j][l] +
                            k*l_s*(1-e)*(1-tau_s)*Isn_A[i][j][l] -
                            (r_H + k*l_s*(1-e)*tau_s + k*l_s*e + muI_H*ART_mort[l])*Isp_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Isp_H[i][j] + A_prog[l]*Isp_A[i][j][l-1] - A_prog[l+1]*Isp_A[i][j][l] - A_mort[l][j][i]*Isp_A[i][j][l];
                                   
          dImn_A[i][j][l] = age_in[i]*Imn_A[i-1][j][l] - age_out[i]*Imn_A[i][j][l] - m_b[i]*Imn_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_H + FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lmn_A[i][j][l] +
                            FM*a_age_A[i][j][l]*sig_H*S_A[i][j][l] +
                            FM*a_age_A[i][j][l]*(1-p_A[j][l])*sig_H*Lsn_A[i][j][l] +
                            theta*Nmn_A[i][j][l] -
                            (r_H + k*l_m*dst_n + k*l_s*(1-dst_n) + muI_H*ART_mort[l])*Imn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Imn_H[i][j] + A_prog[l]*Imn_A[i][j][l-1] - A_prog[l+1]*Imn_A[i][j][l] - A_mort[l][j][i]*Imn_A[i][j][l];

          dImp_A[i][j][l] = age_in[i]*Imp_A[i-1][j][l] - age_out[i]*Imp_A[i][j][l] - m_b[i]*Imp_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_H + FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lmp_A[i][j][l] +
                            FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l])*Lsp_A[i][j][l] +
                            theta*Nmp_A[i][j][l] +
                            (k*l_m*dst_n*(1-tau_m) + k*l_s*(1-dst_n)*(1-(tau_s*eff_n)))*Imn_A[i][j][l] +
                            k*l_s*e*(Isn_A[i][j][l]+Isp_A[i][j][l]) -
                            (r_H + k*l_m*dst_p*tau_m + k*l_s*(1-dst_p)*tau_s*eff_p + muI_H*ART_mort[l])*Imp_A[i][j][l] +                     
                            ART_prop[j]*A_start[l]*Imp_H[i][j] + A_prog[l]*Imp_A[i][j][l-1] - A_prog[l+1]*Imp_A[i][j][l] - A_mort[l][j][i]*Imp_A[i][j][l];              
        
          TB_cases_ART_age[i][j][l] = (v_age_A[i][j][l]*(1-sig_H) + FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lsn_A[i][j][l] + FS*a_age_A[i][j][l]*(1-sig_H)*(S_A[i][j][l] + (1-p_A[j][l])*Lmn_A[i][j][l]) +
                                      (v_age_A[i][j][l]*(1-sig_H) + FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lsp_A[i][j][l] + FS*a_age_A[i][j][l]*(1-sig_H)*(1-p_A[j][l])*Lmp_A[i][j][l] +
                                      (v_age_A[i][j][l]*(1-sig_H) + FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lmn_A[i][j][l] + FM*a_age_A[i][j][l]*(1-sig_H)*(S_A[i][j][l] + (1-p_A[j][l])*Lsn_A[i][j][l]) +
                                      (v_age_A[i][j][l]*(1-sig_H) + FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lmp_A[i][j][l] + FM*a_age_A[i][j][l]*(1-sig_H)*(1-p_A[j][l])*Lsp_A[i][j][l] +
                                      (v_age_A[i][j][l]*sig_H + FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lsn_A[i][j][l] + FS*a_age_A[i][j][l]*sig_H*S_A[i][j][l] + FS*a_age_A[i][j][l]*(1-p_A[j][l])*sig_H*Lmn_A[i][j][l] +
                                      (v_age_A[i][j][l]*sig_H + FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lsp_A[i][j][l] + FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l])*Lmp_A[i][j][l] +
                                      (v_age_A[i][j][l]*sig_H + FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lmn_A[i][j][l] + FM*a_age_A[i][j][l]*sig_H*S_A[i][j][l] + FM*a_age_A[i][j][l]*(1-p_A[j][l])*sig_H*Lsn_A[i][j][l] +
                                      (v_age_A[i][j][l]*sig_H + FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lmp_A[i][j][l] + FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l])*Lsp_A[i][j][l];
                            
          TB_cases_ART = TB_cases_ART + TB_cases_ART_age[i][j][l];
        
        }
      }
    }

    /* Put our calculated rates of change back into ydot */
    /*  AGAIN, MUST BE A NEATER WAY TO DO THIS */
    /* HIV+ */
    for (i=0; i<17; i++) ydot[i] = dS[i];             /* S: 0-16 */
    for (i=17; i<34; i++) ydot[i] = dLsn[i-17];       /* Lsn: 17-33 */
    for (i=34; i<51; i++) ydot[i] = dLsp[i-34];       /* Lsp: 34-50 */
    for (i=51; i<68; i++) ydot[i] = dLmn[i-51];       /* Lmn: 51-67 */
    for (i=68; i<85; i++) ydot[i] = dLmp[i-68];       /* Lmp: 68-84 */
    for (i=85; i<102; i++) ydot[i] = dNsn[i-85];      /* Nsn: 85-101 */
    for (i=102; i<119; i++) ydot[i] = dNsp[i-102];    /* Nsp: 102-118 */
    for (i=119; i<136; i++) ydot[i] = dNmn[i-119];    /* Nmn: 119-135 */
    for (i=136; i<153; i++) ydot[i] = dNmp[i-136];    /* Nmp: 136-152 */ 
    for (i=153; i<170; i++) ydot[i] = dIsn[i-153];    /* Isn: 153-169 */
    for (i=170; i<187; i++) ydot[i] = dIsp[i-170];    /* Isp: 170-186 */
    for (i=187; i<204; i++) ydot[i] = dImn[i-187];    /* Imn: 187-203 */
    for (i=204; i<221; i++) ydot[i] = dImp[i-204];    /* Imp: 204-220 */
    /* HIV+ */
    ij = 221;                                     
    for (j=0; j<7; j++){
      for (i=0; i<17; i++){  
        ydot[ij] = dS_H[i][j];
        ydot[ij+119] = dLsn_H[i][j];
        ydot[ij+(2*119)] = dLsp_H[i][j]; 
        ydot[ij+(3*119)] = dLmn_H[i][j];
        ydot[ij+(4*119)] = dLmp_H[i][j];  
        ydot[ij+(5*119)] = dNsn_H[i][j];
        ydot[ij+(6*119)] = dNsp_H[i][j];
        ydot[ij+(7*119)] = dNmn_H[i][j];
        ydot[ij+(8*119)] = dNmp_H[i][j];
        ydot[ij+(9*119)] = dIsn_H[i][j];
        ydot[ij+(10*119)] = dIsp_H[i][j];
        ydot[ij+(11*119)] = dImn_H[i][j];
        ydot[ij+(12*119)] = dImp_H[i][j];
        ij = ij+1;
      }
    }
    /* HIV+, on ART */
    ij = 1768;   
    for (l=0; l<3; l++){                          
      for (j=0; j<7; j++){
        for (i=0; i<17; i++){
          ydot[ij] = dS_A[i][j][l];
          ydot[ij+357] = dLsn_A[i][j][l];
          ydot[ij+(2*357)] = dLsp_A[i][j][l];
          ydot[ij+(3*357)] = dLmn_A[i][j][l];
          ydot[ij+(4*357)] = dLmp_A[i][j][l];  
          ydot[ij+(5*357)] = dNsn_A[i][j][l];
          ydot[ij+(6*357)] = dNsp_A[i][j][l];
          ydot[ij+(7*357)] = dNmn_A[i][j][l];
          ydot[ij+(8*357)] = dNmp_A[i][j][l];
          ydot[ij+(9*357)] = dIsn_A[i][j][l];
          ydot[ij+(10*357)] = dIsp_A[i][j][l];
          ydot[ij+(11*357)] = dImn_A[i][j][l];
          ydot[ij+(12*357)] = dImp_A[i][j][l];
          ij = ij+1;
        }
      }
    }
    
    /* Finally assign the things we want to use in R (in addition to the state variables) to yout */
    /* The ode call in R needs to define the number of these and give names */
    yout[0] = Total;
    yout[1] = Total_S;
    yout[2] = Total_Ls;
    yout[3] = Total_Lm;
    yout[4] = Total_L;
    yout[5] = Total_Ns;
    yout[6] = Total_Nm;
    yout[7] = Total_N;
    yout[8] = Total_Is;
    yout[9] = Total_Im;
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
    yout[44] = Tot_deaths_age[0];
    yout[45] = Tot_deaths_age[1];
    yout[46] = Tot_deaths_age[2];
    yout[47] = Tot_deaths_age[3];
    yout[48] = Tot_deaths_age[4];
    yout[49] = Tot_deaths_age[5];
    yout[50] = Tot_deaths_age[6];
    yout[51] = Tot_deaths_age[7];
    yout[52] = Tot_deaths_age[8];
    yout[53] = Tot_deaths_age[9];
    yout[54] = Tot_deaths_age[10];
    yout[55] = Tot_deaths_age[11];
    yout[56] = Tot_deaths_age[12];
    yout[57] = Tot_deaths_age[13];
    yout[58] = Tot_deaths_age[14];
    yout[59] = Tot_deaths_age[15];
    yout[60] = Tot_deaths_age[16];
    yout[61] = births;
    yout[62] = Tot_deaths;
    
}



