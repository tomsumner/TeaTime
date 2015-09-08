/* TB model in C code to call from R */

/* To add: Update the way background mort and population adjustent is done
           Let ART mortality be an input (to allow different countries to be modelled)
           Turning off population scaling from 2015
           HIV testing and ART inititation
           TB interventions
           Recode equations to pull out new cases terms - so can adjust care flow more easily in future
*/

/* Can be compiled within R with system("R CMD SHLIB TB_model_v4.c") */
/* This creates a dynamic linked library (.dll) which can be loaded (dyn.load(TB_model_v4.dll)) into R and used as the model fucntion in a call to desolve */

/* C libraries needed */
#include <R.h>
#include <math.h>

/* You need to define number of parameters and forcing functions passed to the model here */
/* These must match number in intializer functions below */
static double parms[47];
static double forc[48];

/* ###### A TRICK TO KEEP UP WITH THE PARAMETERS AND FORCINGS ###### */

/* !!!!!!! NOTE IN C INDEXING STARTS AT 0 AND NOT 1 (i.e. parms[0] is the first value in parms) !!!!!! */

/* Parameters are constant over time */
#define age1 parms[0]       /* rate of entry/exit to/from age group = 1/width */
#define age2 parms[1]       /* exit rate is higher from final age group  as bin is wider */
#define beta parms[2]       /* effective contact rate */
#define a_a parms[3]        /* proportion developing primary disease in adults (>15) */
#define a0 parms[4]         /* proportion developing primary disease in 0-4 year olds */ 
#define a5 parms[5]         /* proportion developing primary disease in 5-9 year olds */
#define a10 parms[6]        /* proportion developing primary disease in 10-14 year olds */
#define p parms[7]          /* protection against disease due to prior infection */
#define v parms[8]          /* reactivation rate */
#define sig_a parms[9]      /* proportion developing smear positive disease in adults (>15) */
#define sig0 parms[10]      /* proportion developing smear positive disease in 0-4 year olds */
#define sig5 parms[11]      /* proportion developing smear positive disease in 5-9 year olds */
#define sig10 parms[12]     /* proportion developing smear positive disease in 10-14 year olds */
#define rel_inf parms[13]   /* relative infectiousness of smear negative TB */
#define theta parms[14]     /* rate of conversion from smear negative to smear positive TB */
#define r parms[15]         /* rate of self cure from TB */
#define mu_N parms[16]      /* mortality rate from smear negative TB (>5 years old) */
#define mu_N0 parms[17]     /* mortality rate from smear negative TB in 0-4 year olds */
#define mu_I parms[18]      /* mortality rate from smear positive TB (>5 years old) */
#define mu_I0 parms[19]     /* mortality rate from smear positive TB in 0-4 year olds */
#define fit_cost parms[20]  /* relative fitness of MDR strains (transmission only) */
#define e parms[21]         /* rate of acquisition of MDR */
#define g parms[22]         /* superinfections */
#define l_s parms[23]       /* linkage to care for diagnosed DS TB */
#define l_m parms[24]       /* linkage to care for diagnosed MDR TB */
#define d parms[25]         /* relative detection rate of smear negative TB */ 
#define tau_s parms[26]     /* treatment success for DS TB */
#define tau_m parms[27]     /* treatment success for MDR TB */
#define eff_n parms[28]     /* efficacy of first line drugs in treating new MDR cases */
#define eff_p parms[29]     /* efficacy of first line drugs in treating retreatment MDR cases */
#define muN_H parms[30]     /* mortaltiy rate from smear negative TB in HIV+ */
#define muI_H parms[31]     /* mortaltiy rate from smear positive TB in HIV+ */
#define RR1a parms[32]      /* Relative risk for primary disease following HIV infection */
#define RR2a parms[33]      /* Relative risk for primary disease by CD4 */
#define RR1v parms[34]      /* Relative risk for reactivation following HIV infection */
#define RR2v parms[35]      /* Relative risk for reactivation by CD4 */
#define RR1p parms[36]      /* Relative risk of protection against reinfection following HIV infection */
#define RR2p parms[37]      /* Relative risk of protection against reinfection by CD4 */
#define ART_TB1 parms[38]   /* Reduction in TB parameters on ART (<6m) */
#define ART_TB2 parms[39]   /* Reduction in TB parameters on ART (7-12m) */
#define ART_TB3 parms[40]   /* Reduction in TB parameters on ART (>12m) */
#define ART_mort1 parms[41] /* Reduction in TB mortality on ART (<6m) */
#define ART_mort2 parms[42] /* Reduction in TB mortality on ART (7-12m) */
#define ART_mort3 parms[43] /* Reduction in TB mortality on ART (>12m) */
#define BCG_eff parms[44]   /* Efficacy of BCG (reduces primary and reactivation risks) */
#define sig_H parms[45]     /* proportion developing smear positive disease in HIV+ */
#define r_H parms[46]       /* rate of self cure from TB in HIV+ */

/* Forcings are time dependant functions */
#define birth_rate forc[0] /* Birth rate */
#define s_birth forc[1]    /* Survival at birth */
#define s5 forc[2]         /* Survival to age x */
#define s10 forc[3] 
#define s15 forc[4] 
#define s20 forc[5] 
#define s25 forc[6]
#define s30 forc[7]
#define s35 forc[8]
#define s40 forc[9]
#define s45 forc[10]
#define s50 forc[11]
#define s55 forc[12]
#define s60 forc[13]
#define s65 forc[14]
#define s70 forc[15]
#define s75 forc[16]
#define s80 forc[17]

#define h0 forc[18]  /* HIV incidence at age x */
#define h5 forc[19]
#define h10 forc[20]
#define h15 forc[21]
#define h20 forc[22]
#define h25 forc[23]
#define h30 forc[24]
#define h35 forc[25]
#define h40 forc[26]
#define h45 forc[27]
#define h50 forc[28]
#define h55 forc[29]
#define h60 forc[30]
#define h65 forc[31]
#define h70 forc[32]
#define h75 forc[33]
#define h80 forc[34]

#define Ahigh forc[35] /* ART coverage by CD4 - depends on overall coverage and eligibility threshold */
#define A500 forc[36]
#define A349 forc[37]
#define A249 forc[38]
#define A199 forc[39]
#define A99 forc[40]
#define A50 forc[41]
#define Athresh forc[42] /* ART eligibility threshold in terms of model CD4 categories i.e. 1 means anyone <500 */

#define BCG_cov forc[43] /* BCG coverage */
#define pop_ad forc[44]  /* used to turn on/off adjustment of population to account for disease induced mortality - idea is to turn it off from 2015 onwards*/
#define k forc[45]       /* detection rate */
#define dst_n forc[46]   /* dst coverage in new cases */
#define dst_p forc[47]   /* dst coverage in previoulsy treated cases */

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
    int N=47;
    odeparms(&N, parms);
}

/* ###### FUNCTION TO INITIALIZE FORCINGS PASSED FROM R - if the number of parameters is changed you must update N here ###### */
void forcc(void (* odeforcs)(int *, double *))
{
    int N=48;
    odeforcs(&N, forc);
}

/* ###### DERIVATIVE FUNCTION ###### */

void derivsc(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if (ip[0] <2) error("nout should be at least 2");
    
    /* Expand state variables so can use more meaningful names than y and ydot (the variables and rates of change are vectors y and ydot here) */
    /* There are 17 age groups, 7 HIV positive categories and 3 times on ART */
    /* _H = HIV+; _A = HIV+ on ART */
    
    /* These are the variables */
    double S[17];          /* Susceptible */
    double Lsn[17];        /* Latent, DS, new */
    double Lsp[17];        /* Latent, DS, previous */
    double Lmn[17];        /* Latent, DR, new */
    double Lmp[17];        /* Latent, DR, previous */
    double Nsn[17];        /* Smear negative, DS, new */
    double Nsp[17];        /* Smear negative, DS, previous */
    double Nmn[17];        /* Smear negative, DR, new */
    double Nmp[17];        /* Smear negative, DR, previous */
    double Isn[17];        /* Smear positive, DS, new */
    double Isp[17];        /* Smear positive, DS, previous */
    double Imn[17];        /* Smear positive, DR, new */
    double Imp[17];        /* Smear positive, DR, previous */
    double S_H[17][7];
    double Lsn_H[17][7];
    double Lsp_H[17][7];
    double Lmn_H[17][7];
    double Lmp_H[17][7];
    double Nsn_H[17][7];
    double Nsp_H[17][7];
    double Nmn_H[17][7];
    double Nmp_H[17][7];
    double Isn_H[17][7];
    double Isp_H[17][7];
    double Imn_H[17][7];
    double Imp_H[17][7];
    double S_A[17][7][3];
    double Lsn_A[17][7][3];
    double Lsp_A[17][7][3];
    double Lmn_A[17][7][3];
    double Lmp_A[17][7][3];
    double Nsn_A[17][7][3];
    double Nsp_A[17][7][3];
    double Nmn_A[17][7][3];
    double Nmp_A[17][7][3];
    double Isn_A[17][7][3];
    double Isp_A[17][7][3];
    double Imn_A[17][7][3];
    double Imp_A[17][7][3];
    /* These are the rates of change (same names but prefixed with d) */
    double dS[17];
    double dLsn[17];
    double dLsp[17];
    double dLmn[17];
    double dLmp[17];
    double dNsn[17];
    double dNsp[17];
    double dNmn[17];
    double dNmp[17];
    double dIsn[17];
    double dIsp[17];
    double dImn[17];
    double dImp[17];
    double dS_H[17][7];
    double dLsn_H[17][7];
    double dLsp_H[17][7];
    double dLmn_H[17][7];
    double dLmp_H[17][7];
    double dNsn_H[17][7];
    double dNsp_H[17][7];
    double dNmn_H[17][7];
    double dNmp_H[17][7];
    double dIsn_H[17][7];
    double dIsp_H[17][7];
    double dImn_H[17][7];
    double dImp_H[17][7];
    double dS_A[17][7][3];
    double dLsn_A[17][7][3];
    double dLsp_A[17][7][3];
    double dLmn_A[17][7][3];
    double dLmp_A[17][7][3];
    double dNsn_A[17][7][3];
    double dNsp_A[17][7][3];
    double dNmn_A[17][7][3];
    double dNmp_A[17][7][3];
    double dIsn_A[17][7][3];
    double dIsp_A[17][7][3];
    double dImn_A[17][7][3];
    double dImp_A[17][7][3]; 
     
    /* intergers to use as counters */ 
    int i;  
    int j;
    int l;
    int ij;
     
    /* Then need to assign the variablesto the correct bit of y (we do the same for the rates of change after solving the equations) */
    /* SEE IF CAN FIND A NEATER WAY TO DO THIS */
     
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

    double H_CD4[7][17] = { /* Distribution of new HIV infections (age,CD4) */
      {0,0,0,0.643,0.643,0.607,0.607,0.585,0.585,0.552,0.552,0.552,0.552,0.552,0.552,0.552,0.552},
      {0,0,0,0.357,0.357,0.393,0.393,0.415,0.415,0.448,0.448,0.448,0.448,0.448,0.448,0.448,0.448},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    };
    double H_prog[8][17] = {  /* Progression through CD4 categories (age, CD4) - has extra row to avoid progression in/out of first/last groups */
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0.111,0.111,0.137,0.137,0.168,0.168,0.192,0.192,0.192,0.192,0.192,0.192,0.192,0.192},
      {0,0,0,0.200,0.200,0.213,0.213,0.299,0.299,0.414,0.414,0.414,0.414,0.414,0.414,0.414,0.414},
      {0,0,0,0.255,0.255,0.364,0.364,0.441,0.441,0.575,0.575,0.575,0.575,0.575,0.575,0.575,0.575},
      {0,0,0,0.398,0.398,0.663,0.663,0.713,0.713,0.838,0.838,0.838,0.838,0.838,0.838,0.838,0.838},
      {0,0,0,0.193,0.193,0.471,0.471,0.491,0.491,0.614,0.614,0.614,0.614,0.614,0.614,0.614,0.614},
      {0,0,0,0.294,0.294,0.765,0.765,0.765,0.765,0.865,0.865,0.865,0.865,0.865,0.865,0.865,0.865},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    };    
    double H_mort[7][17] = { /* Mortality due to HIV (no ART) (age, CD4) */
      {0,0,0,0.005,0.005,0.004,0.004,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005},
      {0,0,0,0.011,0.011,0.010,0.010,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013},
      {0,0,0,0.026,0.026,0.026,0.026,0.036,0.036,0.032,0.032,0.032,0.032,0.032,0.032,0.032,0.032},
      {0,0,0,0.061,0.061,0.069,0.069,0.096,0.096,0.080,0.080,0.080,0.080,0.080,0.080,0.080,0.080},
      {0,0,0,0.139,0.139,0.185,0.185,0.258,0.258,0.203,0.203,0.203,0.203,0.203,0.203,0.203,0.203},
      {0,0,0,0.321,0.321,0.499,0.499,0.691,0.691,0.513,0.513,0.513,0.513,0.513,0.513,0.513,0.513},
      {0,0,0,0.737,0.737,1.342,1.342,1.851,1.851,1.295,1.295,1.295,1.295,1.295,1.295,1.295,1.295}
    };
    double A_mort[3][7][17] = { /* On ART mortality (age,starting CD4, time on ART)  - this is an average of male and female values weigthed by sex of those on ART  */
      {{0,0,0,0.0050,0.0050,0.0035,0.0035,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050},
       {0,0,0,0.0115,0.0115,0.0095,0.0095,0.0134,0.0134,0.0126,0.0126,0.0126,0.0126,0.0126,0.0126,0.0126,0.0126},
       {0,0,0,0.0264,0.0264,0.0256,0.0256,0.0359,0.0359,0.0319,0.0319,0.0319,0.0319,0.0319,0.0319,0.0319,0.0319},
       {0,0,0,0.0557,0.0557,0.0477,0.0477,0.0495,0.0495,0.0489,0.0489,0.0489,0.0489,0.0489,0.0489,0.0489,0.0489},
       {0,0,0,0.0953,0.0953,0.0792,0.0792,0.0833,0.0833,0.0872,0.0872,0.0872,0.0872,0.0872,0.0872,0.0872,0.0872},
       {0,0,0,0.1553,0.1553,0.1305,0.1305,0.1384,0.1384,0.1496,0.1496,0.1496,0.1496,0.1496,0.1496,0.1496,0.1496},
       {0,0,0,0.3418,0.3418,0.2897,0.2897,0.3094,0.3094,0.3431,0.3431,0.3431,0.3431,0.3431,0.3431,0.3431,0.3431}},
      {{0,0,0,0.0050,0.0050,0.0035,0.0035,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050},
       {0,0,0,0.0115,0.0115,0.0095,0.0095,0.0134,0.0134,0.0126,0.0126,0.0126,0.0126,0.0126,0.0126,0.0126,0.0126},
       {0,0,0,0.0204,0.0204,0.0242,0.0242,0.0267,0.0267,0.0306,0.0306,0.0306,0.0306,0.0306,0.0306,0.0306,0.0306},
       {0,0,0,0.0216,0.0216,0.0278,0.0278,0.0283,0.0283,0.0366,0.0366,0.0366,0.0366,0.0366,0.0366,0.0366,0.0366},
       {0,0,0,0.0272,0.0272,0.0351,0.0351,0.0362,0.0362,0.0481,0.0481,0.0481,0.0481,0.0481,0.0481,0.0481,0.0481},
       {0,0,0,0.0343,0.0343,0.0442,0.0442,0.0461,0.0461,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625},
       {0,0,0,0.0496,0.0496,0.0640,0.0640,0.0675,0.0675,0.0935,0.0935,0.0935,0.0935,0.0935,0.0935,0.0935,0.0935}},
      {{0,0,0,0.0050,0.0050,0.0035,0.0035,0.0050,0.0050,0.0045,0.0045,0.0045,0.0045,0.0045,0.0045,0.0045,0.0045},
       {0,0,0,0.0068,0.0068,0.0081,0.0081,0.0076,0.0076,0.0066,0.0066,0.0066,0.0066,0.0066,0.0066,0.0066,0.0066},
       {0,0,0,0.0073,0.0073,0.0093,0.0093,0.0083,0.0083,0.0076,0.0076,0.0076,0.0076,0.0076,0.0076,0.0076,0.0076},
       {0,0,0,0.0079,0.0079,0.0100,0.0100,0.0091,0.0091,0.0087,0.0087,0.0087,0.0087,0.0087,0.0087,0.0087,0.0087},
       {0,0,0,0.0105,0.0105,0.0134,0.0134,0.0128,0.0128,0.0141,0.0141,0.0141,0.0141,0.0141,0.0141,0.0141,0.0141},
       {0,0,0,0.0139,0.0139,0.0177,0.0177,0.0174,0.0174,0.0209,0.0209,0.0209,0.0209,0.0209,0.0209,0.0209,0.0209},
       {0,0,0,0.0211,0.0211,0.0271,0.0271,0.0275,0.0275,0.0355,0.0355,0.0355,0.0355,0.0355,0.0355,0.0355,0.0355}}
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
    /*TB_deaths[i] = TB_deaths[i]*pop_ad;*/ /* We may want to turn off this adjustment for projections beyond 2015 */
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
    /* Then calculate % of pop that die of disease by age */
    double prop_dis_death[17];
    for (i=0; i<17; i++){
      prop_dis_death[i] = (TB_deaths[i]+HIV_deaths[i]+ART_deaths[i])/tot_age[i];
    }
    /* Then calculate the proportion of people who should die per year by age and reduce it by the disease mortality */
    double mort_b[17];
    for (i=0; i<17; i++){
      mort_b[i] = ((1-forc[i+1])/(1/age_out[i]))-prop_dis_death[i];
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
    
    /* Susceptible - note this is done separately to deal with births */
    dS[0] = s_birth*birth_rate*Total/1000 - age1*S[0] - (FS + FM)*S[0] - forc[18]*S[0] - mort_b[0]*S[0];
    for (i=1; i<17; i++) dS[i] = age_in[i]*S[i-1] - age_out[i]*S[i] - (FS + FM)*S[i] - forc[i+18]*S[i] - mort_b[i]*S[i];
    /* Other TB states */
    for (i=0; i<17; i++){
      
      /* Latent, ds, naive */
      dLsn[i] = age_in[i]*Lsn[i-1] - age_out[i]*Lsn[i] - mort_b[i]*Lsn[i] +
                FS*((1-a_age[i])*S[i] + (1-a_age[i])*(1-p)*(1-g)*Lmn[i]) - (v_age[i] + FS*a_age[i]*(1-p) + FM*a_age[i]*(1-p) + FM*(1-a_age[i])*(1-p)*g)*Lsn[i] + r*(Isn[i] + Nsn[i]) - 
                forc[i+18]*Lsn[i];
      /* Latent, ds, prev */
      dLsp[i] = age_in[i]*Lsp[i-1] - age_out[i]*Lsp[i] - mort_b[i]*Lsp[i] +
                FS*(1-a_age[i])*(1-p)*(1-g)*Lmp[i] -
                (v_age[i] + FS*a_age[i]*(1-p) + FM*a_age[i]*(1-p) + FM*(1-a_age[i])*(1-p)*g)*Lsp[i] +
                r*(Isp[i]+Nsp[i]) + k*l_s*(1-e)*tau_s*(Isn[i]+Isp[i]) + k*l_s*(1-e)*tau_s*d*(Nsn[i]+Nsp[i]) - 
                forc[i+18]*Lsp[i];
      /* Latent, mdr, naive */ 
      dLmn[i] = age_in[i]*Lmn[i-1] - age_out[i]*Lmn[i] - mort_b[i]*Lmn[i] +
                FM*((1-a_age[i])*S[i] + (1-a_age[i])*(1-p)*g*Lsn[i]) - (v_age[i] + FM*a_age[i]*(1-p) + FS*a_age[i]*(1-p) + FS*(1-a_age[i])*(1-p)*(1-g))*Lmn[i] +
                r*(Imn[i]+Nmn[i]) - 
                forc[i+18]*Lmn[i];
      /* Latent, mdr, prev */
      dLmp[i] = age_in[i]*Lmp[i-1] - age_out[i]*Lmp[i] - mort_b[i]*Lmp[i] + 
                FM*(1-a_age[i])*(1-p)*g*Lsp[i] - (v_age[i] + FM*a_age[i]*(1-p) + FS*a_age[i]*(1-p) + FS*(1-a_age[i])*(1-p)*(1-g))*Lmp[i] +
                r*(Imp[i]+Nmp[i]) + k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(Imn[i]+d*Nmn[i]) + k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(Imp[i]+d*Nmp[i]) - 
                forc[i+18]*Lmp[i];
      /* Smear neg, ds, new */
      dNsn[i] = age_in[i]*Nsn[i-1] - age_out[i]*Nsn[i] - mort_b[i]*Nsn[i] +
                (v_age[i]*(1-sig_age[i]) + FS*a_age[i]*(1-p)*(1-sig_age[i]))*Lsn[i] + FS*a_age[i]*(1-sig_age[i])*(S[i] + (1-p)*Lmn[i]) -
                (theta + r + k*l_s*d + muN_age[i])*Nsn[i] - 
                forc[i+18]*Nsn[i];
      /* Smear neg, ds, prev */                             
      dNsp[i] = age_in[i]*Nsp[i-1] - age_out[i]*Nsp[i] - mort_b[i]*Nsp[i] + 
                (v_age[i]*(1-sig_age[i]) + FS*a_age[i]*(1-p)*(1-sig_age[i]))*Lsp[i] + FS*a_age[i]*(1-sig_age[i])*(1-p)*Lmp[i] -
                (theta + r + k*l_s*d*(1-e)*tau_s + k*l_s*d*e + muN_age[i])*Nsp[i] +
                k*l_s*d*(1-e)*(1-tau_s)*Nsn[i] - 
                forc[i+18]*Nsp[i];
      /* Smear neg, mdr, new */
      dNmn[i] = age_in[i]*Nmn[i-1] - age_out[i]*Nmn[i] - mort_b[i]*Nmn[i] +
                (v_age[i]*(1-sig_age[i]) + FM*a_age[i]*(1-p)*(1-sig_age[i]))*Lmn[i] + FM*a_age[i]*(1-sig_age[i])*(S[i] + (1-p)*Lsn[i]) -
                (theta + r + k*l_m*d*dst_n + k*l_s*d*(1-dst_n) + muN_age[i])*Nmn[i] - 
                forc[i+18]*Nmn[i];
      /* Smear neg, mdr, prev */
      dNmp[i] = age_in[i]*Nmp[i-1] - age_out[i]*Nmp[i] - mort_b[i]*Nmp[i] +
                (v_age[i]*(1-sig_age[i]) + FM*a_age[i]*(1-p)*(1-sig_age[i]))*Lmp[i] + FM*a_age[i]*(1-sig_age[i])*(1-p)*Lsp[i] -
                (theta + r + k*l_m*d*dst_p*tau_m + k*l_s*d*(1-dst_p)*tau_s*eff_p + muN_age[i])*Nmp[i] +
                k*l_s*d*e*(Nsn[i]+Nsp[i]) + (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*Nmn[i] - 
                forc[i+18]*Nmp[i];
      /* Smear pos, ds, new */
      dIsn[i] = age_in[i]*Isn[i-1] - age_out[i]*Isn[i] - mort_b[i]*Isn[i] +
                (v_age[i]*sig_age[i] + FS*a_age[i]*sig_age[i]*(1-p))*Lsn[i] + FS*a_age[i]*sig_age[i]*S[i] + FS*a_age[i]*(1-p)*sig_age[i]*Lmn[i] +
                theta*Nsn[i] - (r + k*l_s + muI_age[i])*Isn[i] - 
                forc[i+18]*Isn[i];
      /* Smear pos, ds, prev */
      dIsp[i] = age_in[i]*Isp[i-1] - age_out[i]*Isp[i] - mort_b[i]*Isp[i] +
                (v_age[i]*sig_age[i] + FS*a_age[i]*sig_age[i]*(1-p))*Lsp[i] + FS*a_age[i]*sig_age[i]*(1-p)*Lmp[i] +
                theta*Nsp[i] + k*l_s*(1-e)*(1-tau_s)*Isn[i] - (r + k*l_s*(1-e)*tau_s + k*l_s*e + muI_age[i])*Isp[i] - 
                forc[i+18]*Isp[i];
      /* Smear pos, mdr, new */
      dImn[i] = age_in[i]*Imn[i-1] - age_out[i]*Imn[i] - mort_b[i]*Imn[i] +
                (v_age[i]*sig_age[i] + FM*a_age[i]*sig_age[i]*(1-p))*Lmn[i] + FM*a_age[i]*sig_age[i]*S[i] + FM*a_age[i]*(1-p)*sig_age[i]*Lsn[i] +
                theta*Nmn[i] - (r + k*l_m*dst_n + k*l_s*(1-dst_n) + muI_age[i])*Imn[i] - 
                forc[i+18]*Imn[i];
      /* Smear pos, mdr, prev */
      dImp[i] = age_in[i]*Imp[i-1] - age_out[i]*Imp[i] - mort_b[i]*Imp[i] +
                (v_age[i]*sig_age[i] + FM*a_age[i]*sig_age[i]*(1-p))*Lmp[i] + FM*a_age[i]*sig_age[i]*(1-p)*Lsp[i] +
                theta*Nmp[i] + (k*l_m*dst_n*(1-tau_m) + k*l_s*(1-dst_n)*(1-(tau_s*eff_n)))*Imn[i] + k*l_s*e*(Isn[i]+Isp[i]) -
                (r + k*l_m*dst_p*tau_m + k*l_s*(1-dst_p)*tau_s*eff_p + muI_age[i])*Imp[i] - 
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

        dS_H[i][j] = age_in[i]*S_H[i-1][j] - age_out[i]*S_H[i][j] - mort_b[i]*S_H[i][j] - 
                     (FS + FM)*S_H[i][j] +
                     forc[i+18]*H_CD4[j][i]*S[i] - H_prog[j+1][i]*S_H[i][j] + H_prog[j][i]*S_H[i][j-1] - 
                     H_mort[j][i]*S_H[i][j] - ART_prop[j]*S_H[i][j];

        dLsn_H[i][j] = age_in[i]*Lsn_H[i-1][j] - age_out[i]*Lsn_H[i][j] - mort_b[i]*Lsn_H[i][j] +
                       FS*((1-a_age_H[i][j])*S_H[i][j] + (1-a_age_H[i][j])*(1-p_H[j])*(1-g)*Lmn_H[i][j]) -
                       (v_age_H[i][j] + FS*a_age_H[i][j]*(1-p_H[j]) + FM*a_age_H[i][j]*(1-p_H[j]) + FM*(1-a_age_H[i][j])*(1-p_H[j])*g)*Lsn_H[i][j] +
                       r*(Isn_H[i][j] + Nsn_H[i][j]) +
                       forc[i+18]*H_CD4[j][i]*Lsn[i] - H_prog[j+1][i]*Lsn_H[i][j] + H_prog[j][i]*Lsn_H[i][j-1] - 
                       H_mort[j][i]*Lsn_H[i][j] - ART_prop[j]*Lsn_H[i][j];
        
        dLsp_H[i][j] = age_in[i]*Lsp_H[i-1][j] - age_out[i]*Lsp_H[i][j] - mort_b[i]*Lsp_H[i][j] +
                       FS*(1-a_age_H[i][j])*(1-p_H[j])*(1-g)*Lmp_H[i][j] -
                       (v_age_H[i][j] + FS*a_age_H[i][j]*(1-p_H[j]) + FM*a_age_H[i][j]*(1-p_H[j]) + FM*(1-a_age_H[i][j])*(1-p_H[j])*g)*Lsp_H[i][j] +
                       r_H*(Isp_H[i][j]+Nsp_H[i][j]) +
                       k*l_s*(1-e)*tau_s*(Isn_H[i][j]+Isp_H[i][j]) + 
                       k*l_s*(1-e)*tau_s*d*(Nsn_H[i][j]+Nsp_H[i][j]) +
                       forc[i+18]*H_CD4[j][i]*Lsp[i] - H_prog[j+1][i]*Lsp_H[i][j] + H_prog[j][i]*Lsp_H[i][j-1] - 
                       H_mort[j][i]*Lsp_H[i][j] - ART_prop[j]*Lsp_H[i][j];
                             
        dLmn_H[i][j] = age_in[i]*Lmn_H[i-1][j] - age_out[i]*Lmn_H[i][j] - mort_b[i]*Lmn_H[i][j] +
                       FM*((1-a_age_H[i][j])*S_H[i][j] + (1-a_age_H[i][j])*(1-p_H[j])*g*Lsn_H[i][j]) -
                       (v_age_H[i][j] + FM*a_age_H[i][j]*(1-p_H[j]) + FS*a_age_H[i][j]*(1-p_H[j]) + FS*(1-a_age_H[i][j])*(1-p_H[j])*(1-g))*Lmn_H[i][j] +
                       r_H*(Imn_H[i][j]+Nmn_H[i][j]) +
                       forc[i+18]*H_CD4[j][i]*Lmn[i] - H_prog[j+1][i]*Lmn_H[i][j] + H_prog[j][i]*Lmn_H[i][j-1] - 
                       H_mort[j][i]*Lmn_H[i][j] - ART_prop[j]*Lmn_H[i][j];
                     
        dLmp_H[i][j] = age_in[i]*Lmp_H[i-1][j] - age_out[i]*Lmp_H[i][j] - mort_b[i]*Lmp_H[i][j] + 
                       FM*(1-a_age_H[i][j])*(1-p_H[j])*g*Lsp_H[i][j] -
                       (v_age_H[i][j] + FM*a_age_H[i][j]*(1-p_H[j]) + FS*a_age_H[i][j]*(1-p_H[j]) + FS*(1-a_age_H[i][j])*(1-p_H[j])*(1-g))*Lmp_H[i][j] +
                       r_H*(Imp_H[i][j]+Nmp_H[i][j]) +
                       k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(Imn_H[i][j]+d*Nmn_H[i][j]) +
                       k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(Imp_H[i][j]+d*Nmp_H[i][j]) +
                       forc[i+18]*H_CD4[j][i]*Lmp[i] - H_prog[j+1][i]*Lmp_H[i][j] + H_prog[j][i]*Lmp_H[i][j-1] - 
                       H_mort[j][i]*Lmp_H[i][j] - ART_prop[j]*Lmp_H[i][j];
           
        dNsn_H[i][j] = age_in[i]*Nsn_H[i-1][j] - age_out[i]*Nsn_H[i][j] - mort_b[i]*Nsn_H[i][j] +
                       (v_age_H[i][j]*(1-sig_H) + FS*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lsn_H[i][j] +
                       FS*a_age_H[i][j]*(1-sig_H)*(S_H[i][j] + (1-p_H[j])*Lmn_H[i][j]) -
                       (theta + r_H + k*l_s*d + muN_H)*Nsn_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Nsn[i] - H_prog[j+1][i]*Nsn_H[i][j] + H_prog[j][i]*Nsn_H[i][j-1] - 
                       H_mort[j][i]*Nsn_H[i][j] - ART_prop[j]*Nsn_H[i][j];
                                                
        dNsp_H[i][j] = age_in[i]*Nsp_H[i-1][j] - age_out[i]*Nsp_H[i][j] - mort_b[i]*Nsp_H[i][j] + 
                       (v_age_H[i][j]*(1-sig_H) + FS*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lsp_H[i][j] +
                       FS*a_age_H[i][j]*(1-sig_H)*(1-p_H[j])*Lmp_H[i][j] -
                       (theta + r_H + k*l_s*d*(1-e)*tau_s + k*l_s*d*e + muN_H)*Nsp_H[i][j] +
                       k*l_s*d*(1-e)*(1-tau_s)*Nsn_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Nsp[i] - H_prog[j+1][i]*Nsp_H[i][j] + H_prog[j][i]*Nsp_H[i][j-1] - 
                       H_mort[j][i]*Nsp_H[i][j] - ART_prop[j]*Nsp_H[i][j];
                       
        dNmn_H[i][j] = age_in[i]*Nmn_H[i-1][j] - age_out[i]*Nmn_H[i][j] - mort_b[i]*Nmn_H[i][j] +
                       (v_age_H[i][j]*(1-sig_H) + FM*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lmn_H[i][j] +
                       FM*a_age_H[i][j]*(1-sig_H)*(S_H[i][j] + (1-p_H[j])*Lsn_H[i][j]) -
                       (theta + r_H + k*l_m*d*dst_n + k*l_s*d*(1-dst_n) + muN_H)*Nmn_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Nmn[i] - H_prog[j+1][i]*Nmn_H[i][j] + H_prog[j][i]*Nmn_H[i][j-1] - 
                       H_mort[j][i]*Nmn_H[i][j] - ART_prop[j]*Nmn_H[i][j];
                       
        dNmp_H[i][j] = age_in[i]*Nmp_H[i-1][j] - age_out[i]*Nmp_H[i][j] - mort_b[i]*Nmp_H[i][j] +
                       (v_age_H[i][j]*(1-sig_H) + FM*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lmp_H[i][j] +
                       FM*a_age_H[i][j]*(1-sig_H)*(1-p_H[j])*Lsp_H[i][j] -
                       (theta + r_H + k*l_m*d*dst_p*tau_m + k*l_s*d*(1-dst_p)*tau_s*eff_p + muN_H)*Nmp_H[i][j] +
                       k*l_s*d*e*(Nsn_H[i][j]+Nsp_H[i][j]) +
                       (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*Nmn_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Nmp[i] - H_prog[j+1][i]*Nmp_H[i][j] + H_prog[j][i]*Nmp_H[i][j-1] - 
                       H_mort[j][i]*Nmp_H[i][j] - ART_prop[j]*Nmp_H[i][j];
           
        dIsn_H[i][j] = age_in[i]*Isn_H[i-1][j] - age_out[i]*Isn_H[i][j] - mort_b[i]*Isn_H[i][j] +
                       (v_age_H[i][j]*sig_H + FS*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lsn_H[i][j] +
                       FS*a_age_H[i][j]*sig_H*S_H[i][j] +
                       FS*a_age_H[i][j]*(1-p_H[j])*sig_H*Lmn_H[i][j] +
                       theta*Nsn_H[i][j] -
                       (r_H + k*l_s + muI_H)*Isn_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Isn[i] - H_prog[j+1][i]*Isn_H[i][j] + H_prog[j][i]*Isn_H[i][j-1] - 
                       H_mort[j][i]*Isn_H[i][j] - ART_prop[j]*Isn_H[i][j];

        dIsp_H[i][j] = age_in[i]*Isp_H[i-1][j] - age_out[i]*Isp_H[i][j] - mort_b[i]*Isp_H[i][j] +
                       (v_age_H[i][j]*sig_H + FS*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lsp_H[i][j] +
                       FS*a_age_H[i][j]*sig_H*(1-p_H[j])*Lmp_H[i][j] +
                       theta*Nsp_H[i][j] +
                       k*l_s*(1-e)*(1-tau_s)*Isn_H[i][j] -
                       (r_H + k*l_s*(1-e)*tau_s + k*l_s*e + muI_H)*Isp_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Isp[i] - H_prog[j+1][i]*Isp_H[i][j] + H_prog[j][i]*Isp_H[i][j-1] - 
                       H_mort[j][i]*Isp_H[i][j] - ART_prop[j]*Isp_H[i][j];
                                   
        dImn_H[i][j] = age_in[i]*Imn_H[i-1][j] - age_out[i]*Imn_H[i][j] - mort_b[i]*Imn_H[i][j] +
                       (v_age_H[i][j]*sig_H + FM*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lmn_H[i][j] +
                       FM*a_age_H[i][j]*sig_H*S_H[i][j] +
                       FM*a_age_H[i][j]*(1-p_H[j])*sig_H*Lsn_H[i][j] +
                       theta*Nmn_H[i][j] -
                       (r_H + k*l_m*dst_n + k*l_s*(1-dst_n) + muI_H)*Imn_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Imn[i] - H_prog[j+1][i]*Imn_H[i][j] + H_prog[j][i]*Imn_H[i][j-1] - 
                       H_mort[j][i]*Imn_H[i][j] - ART_prop[j]*Imn_H[i][j];

        dImp_H[i][j] = age_in[i]*Imp_H[i-1][j] - age_out[i]*Imp_H[i][j] - mort_b[i]*Imp_H[i][j] +
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
        
          dS_A[i][j][l] = age_in[i]*S_A[i-1][j][l] - age_out[i]*S_A[i][j][l] - mort_b[i]*S_A[i][j][l] - (FS + FM)*S_A[i][j][l] +
                          ART_prop[j]*A_start[l]*S_H[i][j] + A_prog[l]*S_A[i][j][l-1] - A_prog[l+1]*S_A[i][j][l] - A_mort[l][j][i]*S_A[i][j][l];

          dLsn_A[i][j][l] = age_in[i]*Lsn_A[i-1][j][l] - age_out[i]*Lsn_A[i][j][l] - mort_b[i]*Lsn_A[i][j][l] +
                            FS*((1-a_age_A[i][j][l])*S_A[i][j][l] + (1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g)*Lmn_A[i][j][l]) -
                            (v_age_A[i][j][l] + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g)*Lsn_A[i][j][l] +
                            r_H*(Isn_A[i][j][l] + Nsn_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lsn_H[i][j] + A_prog[l]*Lsn_A[i][j][l-1] - A_prog[l+1]*Lsn_A[i][j][l] - A_mort[l][j][i]*Lsn_A[i][j][l];
        
          dLsp_A[i][j][l] = age_in[i]*Lsp_A[i-1][j][l] - age_out[i]*Lsp_A[i][j][l] - mort_b[i]*Lsp_A[i][j][l] +
                            FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g)*Lmp_A[i][j][l] -
                            (v_age_A[i][j][l] + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g)*Lsp_A[i][j][l] +
                            r_H*(Isp_A[i][j][l]+Nsp_A[i][j][l]) +
                            k*l_s*(1-e)*tau_s*(Isn_A[i][j][l]+Isp_A[i][j][l]) + 
                            k*l_s*(1-e)*tau_s*d*(Nsn_A[i][j][l]+Nsp_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lsp_H[i][j] + A_prog[l]*Lsp_A[i][j][l-1] - A_prog[l+1]*Lsp_A[i][j][l] - A_mort[l][j][i]*Lsp_A[i][j][l];
                             
          dLmn_A[i][j][l] = age_in[i]*Lmn_A[i-1][j][l] - age_out[i]*Lmn_A[i][j][l] - mort_b[i]*Lmn_A[i][j][l] +
                            FM*((1-a_age_A[i][j][l])*S_A[i][j][l] + (1-a_age_A[i][j][l])*(1-p_A[j][l])*g*Lsn_A[i][j][l]) -
                            (v_age_A[i][j][l] + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g))*Lmn_A[i][j][l] +
                            r_H*(Imn_A[i][j][l]+Nmn_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lmn_H[i][j] + A_prog[l]*Lmn_A[i][j][l-1] - A_prog[l+1]*Lmn_A[i][j][l] - A_mort[l][j][i]*Lmn_A[i][j][l];
                     
          dLmp_A[i][j][l] = age_in[i]*Lmp_A[i-1][j][l] - age_out[i]*Lmp_A[i][j][l] - mort_b[i]*Lmp_A[i][j][l] + 
                            FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g*Lsp_A[i][j][l] -
                            (v_age_A[i][j][l] + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g))*Lmp_A[i][j][l] +
                            r_H*(Imp_A[i][j][l]+Nmp_A[i][j][l]) +
                            k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(Imn_A[i][j][l]+d*Nmn_A[i][j][l]) +
                            k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(Imp_A[i][j][l]+d*Nmp_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lmp_H[i][j] + A_prog[l]*Lmp_A[i][j][l-1] - A_prog[l+1]*Lmp_A[i][j][l] - A_mort[l][j][i]*Lmp_A[i][j][l];
           
          dNsn_A[i][j][l] = age_in[i]*Nsn_A[i-1][j][l]- age_out[i]*Nsn_A[i][j][l] - mort_b[i]*Nsn_A[i][j][l]+
                            (v_age_A[i][j][l]*(1-sig_H) + FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lsn_A[i][j][l] +
                            FS*a_age_A[i][j][l]*(1-sig_H)*(S_A[i][j][l] + (1-p_A[j][l])*Lmn_A[i][j][l]) -
                            (theta + r_H + k*l_s*d + muN_H*ART_mort[l])*Nsn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Nsn_H[i][j] + A_prog[l]*Nsn_A[i][j][l-1] - A_prog[l+1]*Nsn_A[i][j][l] - A_mort[l][j][i]*Nsn_A[i][j][l];
                                                
          dNsp_A[i][j][l] = age_in[i]*Nsp_A[i-1][j][l] - age_out[i]*Nsp_A[i][j][l] - mort_b[i]*Nsp_A[i][j][l] + 
                            (v_age_A[i][j][l]*(1-sig_H) + FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lsp_A[i][j][l] +
                            FS*a_age_A[i][j][l]*(1-sig_H)*(1-p_A[j][l])*Lmp_A[i][j][l] -
                            (theta + r_H + k*l_s*d*(1-e)*tau_s + k*l_s*d*e + muN_H*ART_mort[l])*Nsp_A[i][j][l] +
                            k*l_s*d*(1-e)*(1-tau_s)*Nsn_A[i][j][l] +      
                            ART_prop[j]*A_start[l]*Nsp_H[i][j] + A_prog[l]*Nsp_A[i][j][l-1] - A_prog[l+1]*Nsp_A[i][j][l] - A_mort[l][j][i]*Nsp_A[i][j][l];
                       
          dNmn_A[i][j][l] = age_in[i]*Nmn_A[i-1][j][l] - age_out[i]*Nmn_A[i][j][l] - mort_b[i]*Nmn_A[i][j][l] +
                            (v_age_A[i][j][l]*(1-sig_H) + FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lmn_A[i][j][l] +
                            FM*a_age_A[i][j][l]*(1-sig_H)*(S_A[i][j][l] + (1-p_A[j][l])*Lsn_A[i][j][l]) -
                            (theta + r_H + k*l_m*d*dst_n + k*l_s*d*(1-dst_n) + muN_H*ART_mort[l])*Nmn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Nmn_H[i][j] + A_prog[l]*Nmn_A[i][j][l-1] - A_prog[l+1]*Nmn_A[i][j][l] - A_mort[l][j][i]*Nmn_A[i][j][l];
                       
          dNmp_A[i][j][l] = age_in[i]*Nmp_A[i-1][j][l] - age_out[i]*Nmp_A[i][j][l] - mort_b[i]*Nmp_A[i][j][l] +
                            (v_age_A[i][j][l]*(1-sig_H) + FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lmp_A[i][j][l] +
                            FM*a_age_A[i][j][l]*(1-sig_H)*(1-p_A[j][l])*Lsp_A[i][j][l] -
                            (theta + r_H + k*l_m*d*dst_p*tau_m + k*l_s*d*(1-dst_p)*tau_s*eff_p + muN_H*ART_mort[l])*Nmp_A[i][j][l] +
                            k*l_s*d*e*(Nsn_A[i][j][l]+Nsp_A[i][j][l]) +
                            (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*Nmn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Nmp_H[i][j] + A_prog[l]*Nmp_A[i][j][l-1] - A_prog[l+1]*Nmp_A[i][j][l] - A_mort[l][j][i]*Nmp_A[i][j][l];
           
          dIsn_A[i][j][l] = age_in[i]*Isn_A[i-1][j][l] - age_out[i]*Isn_A[i][j][l] - mort_b[i]*Isn_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_H + FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lsn_A[i][j][l] +
                            FS*a_age_A[i][j][l]*sig_H*S_A[i][j][l] +
                            FS*a_age_A[i][j][l]*(1-p_A[j][l])*sig_H*Lmn_A[i][j][l] +
                            theta*Nsn_A[i][j][l] -
                            (r_H + k*l_s + muI_H*ART_mort[l])*Isn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Isn_H[i][j] + A_prog[l]*Isn_A[i][j][l-1] - A_prog[l+1]*Isn_A[i][j][l] - A_mort[l][j][i]*Isn_A[i][j][l];

          dIsp_A[i][j][l] = age_in[i]*Isp_A[i-1][j][l] - age_out[i]*Isp_A[i][j][l] - mort_b[i]*Isp_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_H + FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lsp_A[i][j][l] +
                            FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l])*Lmp_A[i][j][l] +
                            theta*Nsp_A[i][j][l] +
                            k*l_s*(1-e)*(1-tau_s)*Isn_A[i][j][l] -
                            (r_H + k*l_s*(1-e)*tau_s + k*l_s*e + muI_H*ART_mort[l])*Isp_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Isp_H[i][j] + A_prog[l]*Isp_A[i][j][l-1] - A_prog[l+1]*Isp_A[i][j][l] - A_mort[l][j][i]*Isp_A[i][j][l];
                                   
          dImn_A[i][j][l] = age_in[i]*Imn_A[i-1][j][l] - age_out[i]*Imn_A[i][j][l] - mort_b[i]*Imn_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_H + FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lmn_A[i][j][l] +
                            FM*a_age_A[i][j][l]*sig_H*S_A[i][j][l] +
                            FM*a_age_A[i][j][l]*(1-p_A[j][l])*sig_H*Lsn_A[i][j][l] +
                            theta*Nmn_A[i][j][l] -
                            (r_H + k*l_m*dst_n + k*l_s*(1-dst_n) + muI_H*ART_mort[l])*Imn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Imn_H[i][j] + A_prog[l]*Imn_A[i][j][l-1] - A_prog[l+1]*Imn_A[i][j][l] - A_mort[l][j][i]*Imn_A[i][j][l];

          dImp_A[i][j][l] = age_in[i]*Imp_A[i-1][j][l] - age_out[i]*Imp_A[i][j][l] - mort_b[i]*Imp_A[i][j][l] +
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
    yout[44] = prop_dis_death[0];
    yout[45] = prop_dis_death[1];
    yout[46] = prop_dis_death[2];
    yout[47] = prop_dis_death[3];
    yout[48] = prop_dis_death[4];
    yout[49] = prop_dis_death[5];
    yout[50] = prop_dis_death[6];
    yout[51] = prop_dis_death[7];
    yout[52] = prop_dis_death[8];
    yout[53] = prop_dis_death[9];
    yout[54] = prop_dis_death[10];
    yout[55] = prop_dis_death[11];
    yout[56] = prop_dis_death[12];
    yout[57] = prop_dis_death[13];
    yout[58] = prop_dis_death[14];
    yout[59] = prop_dis_death[15];
    yout[60] = prop_dis_death[16];
}


