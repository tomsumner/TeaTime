/* TB model in C code to call from R */

/* Version 10 - now matches TIME including ART (pretty good fit) */
/* Testing conversion back to 5 year age bins to increase speed  */   

/* Can be compiled within R with system("R CMD SHLIB TB_model_v10.c") */
/* This creates a dynamic linked library (.dll) which can be loaded (dyn.load(TB_model_v10.dll)) into R and used as the model function in a call to desolve */

/* C libraries needed */
#include <R.h>
#include <math.h>

/* You need to define number of parameters and forcing functions passed to the model here */
/* These must match number in intializer functions below */
static double parms[43];
static double forc[100];

/* ###### A TRICK TO KEEP UP WITH THE PARAMETERS AND FORCINGS ###### */

/* !!!!!!! NOTE IN C INDEXING STARTS AT 0 AND NOT 1 (i.e. parms[0] is the first value in parms) !!!!!! */

/* Parameters are constant over time */
#define beta parms[0]       /* effective contact rate */
#define a_a parms[1]        /* proportion developing primary disease in adults (>15) */
#define a0 parms[2]         /* proportion developing primary disease in 0-4 year olds */ 
#define a5 parms[3]         /* proportion developing primary disease in 5-9 year olds */
#define a10 parms[4]        /* proportion developing primary disease in 10-14 year olds */
#define v parms[5]          /* reactivation rate */
#define p parms[6]          /* protection against disease due to prior infection */
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

#define s5 forc[1]         /* Survival to age x */
#define s10 forc[2] 
#define s15 forc[3] 
#define s20 forc[4] 
#define s25 forc[5]
#define s30 forc[6]
#define s35 forc[7]
#define s40 forc[8]
#define s45 forc[9]
#define s50 forc[10]
#define s55 forc[11]
#define s60 forc[12]
#define s65 forc[13]
#define s70 forc[14]
#define s75 forc[15]
#define s80 forc[16]
#define s100 forc[17]

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

#define BCG_cov forc[35] /* BCG coverage */

#define pop_ad forc[36]  /* used to turn on/off adjustment of population to account for disease induced mortality - idea is to turn it off from 2015 onwards*/

#define kneg forc[37]  /* detection rate */
#define kpos forc[38]

#define rel_d forc[39] /* relative detection of smear negatives */

#define dstneg_n forc[40]   /* dst coverage in new cases */
#define dstneg_p forc[41]   /* dst coverage in previoulsy treated cases */
#define dstpos_n forc[42]  
#define dstpos_p forc[43] 

#define l_s forc[44]  /* linkage to care */
#define l_m forc[45]

#define tneg_s forc[46]  /* treatment success */
#define tpos_s forc[47]
#define tART_s forc[48]
#define tneg_m forc[49]
#define tpos_m forc[50]
#define tART_m forc[51]

#define mig0 forc[52]  /* migration by age x */
#define mig5 forc[53]
#define mig10 forc[54]
#define mig15 forc[55]
#define mig20 forc[56]
#define mig25 forc[57]
#define mig30 forc[58]
#define mig35 forc[59]
#define mig40 forc[60]
#define mig45 forc[61]
#define mig50 forc[62]
#define mig55 forc[63]
#define mig60 forc[64]
#define mig65 forc[65]
#define mig70 forc[66]
#define mig75 forc[67]
#define mig80 forc[68]

/* sens and spec for tests */
#define se_I_neg forc[69]
#define se_N_neg forc[70]
#define se_m_neg forc[71]

#define sp_I_neg forc[72]
#define sp_N_neg forc[73]
#define sp_m_neg forc[74]

#define se_I_pos forc[75]
#define se_N_pos forc[76]
#define se_m_pos forc[77]

#define sp_I_pos forc[78]
#define sp_N_pos forc[79]
#define sp_m_pos forc[80]

/* RR for presentation for healthy individuals */
#define health forc[81]

/* ART coverage (% of those eligible by age) */
#define A0 forc[82]
#define A5 forc[83]
#define A10 forc[84]
#define A15 forc[85]
#define A20 forc[86]
#define A25 forc[87]
#define A30 forc[88]
#define A35 forc[89]
#define A40 forc[90]
#define A45 forc[91]
#define A50 forc[92]
#define A55 forc[93]
#define A60 forc[94]
#define A65 forc[95]
#define A70 forc[96]
#define A75 forc[97]
#define A80 forc[98]

#define Athresh forc[99] /* ART eligibility threshold (CD4 category) */

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
    int N=100;
    odeforcs(&N, forc);
}

/* ##### EVENTS ARE USED TO ADD BIRTHS AND SHIFT THE POPULATION BY AGE - EQUIVALENT TO THE METHOD OF SCHENZLE ###### */

void event(int *n, double *t, double *y) 
{
  int i;
  
  /* Store current population in temp and shift 1/5 every age group forward one */
  double temp[7395];
  temp[0] = y[0];
  for (i=1; i<7395; i++){
    temp[i] = y[i];
    y[i] = (temp[i]*4/5) + (temp[i-1]/5);
  }
  /* Set every 0 age group to 4/5 of current value and every >80 age group to 1/5 previous age group plus those still surviving  */
  for (i=0; i<7395; i+=17) {
    y[i] = (temp[i]*4/5);
    y[i+16] = temp[i+16] + temp[i+15]/5;  
  }
  /* Then add births into group 0 - only susceptibles get born */ 
  y[0] = y[0] + birth_rate*sumsum(temp,0,7394)/1000;
/*y[0] = 600;*/
}

/* ###### DERIVATIVE FUNCTIONS - THIS IS THE MODEL ITSELF ###### */

void derivs5(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if (ip[0] <2) error("nout should be at least 2");
    
    /* Expand state variables so can use more meaningful names than y and ydot (the variables and rates of change are vectors y and ydot here) */
    /* There are 17 age groups, 7 HIV positive categories and 3 times on ART */
    /* _H = HIV+; _A = HIV+ on ART */
    
    /* These are the variables */
    double S[17]={0};   double S_H[17][7]={{0}};   double S_A[17][7][3]={{{0}}};      /* Susceptible */
    double Lsn[17]={0}; double Lsn_H[17][7]={{0}}; double Lsn_A[17][7][3]={{{0}}};    /* Latent, DS, new */
    double Lsp[17]={0}; double Lsp_H[17][7]={{0}}; double Lsp_A[17][7][3]={{{0}}};    /* Latent, DS, previous */
    double Lmn[17]={0}; double Lmn_H[17][7]={{0}}; double Lmn_A[17][7][3]={{{0}}};    /* Latent, DR, new */
    double Lmp[17]={0}; double Lmp_H[17][7]={{0}}; double Lmp_A[17][7][3]={{{0}}};    /* Latent, DR, previous */
    double Nsn[17]={0}; double Nsn_H[17][7]={{0}}; double Nsn_A[17][7][3]={{{0}}};    /* Smear negative, DS, new */
    double Nsp[17]={0}; double Nsp_H[17][7]={{0}}; double Nsp_A[17][7][3]={{{0}}};    /* Smear negative, DS, previous */
    double Nmn[17]={0}; double Nmn_H[17][7]={{0}}; double Nmn_A[17][7][3]={{{0}}};    /* Smear negative, DR, new */
    double Nmp[17]={0}; double Nmp_H[17][7]={{0}}; double Nmp_A[17][7][3]={{{0}}};    /* Smear negative, DR, previous */
    double Isn[17]={0}; double Isn_H[17][7]={{0}}; double Isn_A[17][7][3]={{{0}}};    /* Smear positive, DS, new */
    double Isp[17]={0}; double Isp_H[17][7]={{0}}; double Isp_A[17][7][3]={{{0}}};    /* Smear positive, DS, previous */
    double Imn[17]={0}; double Imn_H[17][7]={{0}}; double Imn_A[17][7][3]={{{0}}};    /* Smear positive, DR, new */
    double Imp[17]={0}; double Imp_H[17][7]={{0}}; double Imp_A[17][7][3]={{{0}}};    /* Smear positive, DR, previous */
    double PTn[17]={0}; double PTn_H[17][7]={{0}}; double PTn_A[17][7][3]={{{0}}};     /* Post PT, new - also move people here if they are false positive for TB and receive Rx */
    double PTp[17]={0}; double PTp_H[17][7]={{0}}; double PTp_A[17][7][3]={{{0}}};     /* Post PT, previous - also move people here if they are false positive for TB and receive Rx */

    /* These are the rates of change (same names but prefixed with d) */
    double dS[17]={0};   double dS_H[17][7]={{0}};   double dS_A[17][7][3]={{{0}}};
    double dLsn[17]={0}; double dLsn_H[17][7]={{0}}; double dLsn_A[17][7][3]={{{0}}};
    double dLsp[17]={0}; double dLsp_H[17][7]={{0}}; double dLsp_A[17][7][3]={{{0}}};
    double dLmn[17]={0}; double dLmn_H[17][7]={{0}}; double dLmn_A[17][7][3]={{{0}}};
    double dLmp[17]={0}; double dLmp_H[17][7]={{0}}; double dLmp_A[17][7][3]={{{0}}};
    double dNsn[17]={0}; double dNsn_H[17][7]={{0}}; double dNsn_A[17][7][3]={{{0}}};
    double dNsp[17]={0}; double dNsp_H[17][7]={{0}}; double dNsp_A[17][7][3]={{{0}}};
    double dNmn[17]={0}; double dNmn_H[17][7]={{0}}; double dNmn_A[17][7][3]={{{0}}};
    double dNmp[17]={0}; double dNmp_H[17][7]={{0}}; double dNmp_A[17][7][3]={{{0}}};
    double dIsn[17]={0}; double dIsn_H[17][7]={{0}}; double dIsn_A[17][7][3]={{{0}}};
    double dIsp[17]={0}; double dIsp_H[17][7]={{0}}; double dIsp_A[17][7][3]={{{0}}};
    double dImn[17]={0}; double dImn_H[17][7]={{0}}; double dImn_A[17][7][3]={{{0}}};
    double dImp[17]={0}; double dImp_H[17][7]={{0}}; double dImp_A[17][7][3]={{{0}}};
    double dPTn[17]={0}; double dPTn_H[17][7]={{0}}; double dPTn_A[17][7][3]={{{0}}}; 
    double dPTp[17]={0}; double dPTp_H[17][7]={{0}}; double dPTp_A[17][7][3]={{{0}}}; 

    /* intergers to use as counters */ 
    int i,j,l,ij;
     
    /* Then need to assign the variables to the correct bit of y (we do the same for the rates of change after solving the equations) */
     
    /* HIV- */ 
    
    int n_age = 17;     /* Number of age groups */
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

    /* Create vectors of disease parameters (by age - now 5 year bins) to use in derivatives - includes BCG effect on risk of primary disease */
    double a_age[17];
    double sig_age[17];
    double v_age[17];
    double muN_age[17];
    double muI_age[17];
    double bcg = (BCG_cov*(1-BCG_eff)+(1-BCG_cov));  /* those on bcg (BCG_COV) have RR of (1-BCG_eff) */

    a_age[0] = a0*bcg;
    a_age[1] = a5*bcg;
    a_age[2] = a10*bcg;
      
    sig_age[0] = sig0;
    sig_age[1] = sig5;
    sig_age[2] = sig10;
      
    v_age[0] = v;
    v_age[1] = v;
    v_age[2] = v;
      
    muN_age[0] = mu_N0;
    muN_age[1] = mu_N;  /* note for mortality only 0-4 has different rate from 15+ */
    muN_age[2] = mu_N;
      
    muI_age[0] = mu_I0;
    muI_age[1] = mu_I;
    muI_age[2] = mu_I;
  
    for (i=3; i<n_age; i++){
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
    double a_age_H[17][7];
    double p_H[7];
    double v_age_H[17][7];
    double a_age_A[17][7][3];
    double p_A[7][3];
    double v_age_A[17][7][3];
    double muN_H_A[17][3];
    double muI_H_A[17][3];
    for (j=0; j<n_HIV; j++){
      p_H[j] = p*RR1p*pow(RR2p,-1*(500-mid_CD4[j])/100);
      for (i=0; i<n_age; i++){
        a_age_H[i][j] = fmin(a_age[i]*RR1a*pow(RR2a,(500-mid_CD4[j])/100),1); /* fmin ensures proportion developing active disease can't go above 1 - assume this cap is also applied in TIME */
        v_age_H[i][j] = v_age[i]*RR1v*pow(RR2v,(500-mid_CD4[j])/100);
        for (l=0; l<n_ART; l++){
          a_age_A[i][j][l] = fmax(a_age_H[i][j]*(1-ART_TB[l]),a_age[i]);   /* fmax or fmin ensures being on ART can't be better than being HIV- */
          v_age_A[i][j][l] = fmax(v_age_H[i][j]*(1-ART_TB[l]),v_age[i]);
          
          p_A[j][l] = fmin(1-(1-p_H[j])*(1-ART_TB[l]),p);                  /* protection term gets higher as ART is taken - now mimics TIME code */

          muN_H_A[i][l] = fmax(muN_H*(1-ART_mort[l]),muN_age[i]); /* make sure mortality can't go lower than HIV- */ 
          muI_H_A[i][l] = fmax(muI_H*(1-ART_mort[l]),muI_age[i]);
        
        }
      }
    }
  
    /* Set up parameters for HIV model - these are taken from AIM */

    double H_CD4[7][17] = {{0},{0},{0},{0},{0},{0},{0}}; /* Distribution of new HIV infections (age,CD4) - assume distribution of new child infections mirrors adults */
    for (i=0; i<5; i++) {
      H_CD4[0][i] = 0.643; H_CD4[1][i] = 0.357;
    }
    for (i=5; i<7; i++){
      H_CD4[0][i] = 0.607; H_CD4[1][i] = 0.393;
    }
    for (i=7; i<9; i++){
      H_CD4[0][i] = 0.585; H_CD4[1][i] = 0.415;
    }
    for (i=9; i<17; i++){
      H_CD4[0][i] = 0.552; H_CD4[1][i] = 0.448;
    }

    /* Have updated these values based on the durations in the AIM manual (rate = 1/duration) as they are different from rates in AIM editor in software */
    /* those are actually risks (i.e. 1-exp(-rate)) */
    double H_prog[8][17] = {{0},{0},{0},{0},{0},{0},{0},{0}};  /* Progression through CD4 categories (age, CD4) - has extra row to avoid progression in/out of first/last groups */
    for (i=0; i<3; i++){
      H_prog[1][i] = 0.298; H_prog[2][i] = 0.239; H_prog[3][i] = 0.183; H_prog[4][i] = 0.183; H_prog[5][i] = 0.130; H_prog[6][i] = 0.130;
    }
    for (i=3; i<5; i++){
      H_prog[1][i] = 0.117; H_prog[2][i] = 0.223; H_prog[3][i] = 0.294; H_prog[4][i] = 0.508; H_prog[5][i] = 0.214; H_prog[6][i] = 0.348;
    }
    for (i=5; i<7; i++){
      H_prog[1][i] = 0.147; H_prog[2][i] = 0.240; H_prog[3][i] = 0.452; H_prog[4][i] = 1.087; H_prog[5][i] = 0.637; H_prog[6][i] = 1.449;
    }
    for (i=7; i<9; i++){
      H_prog[1][i] = 0.183; H_prog[2][i] = 0.355; H_prog[3][i] = 0.581; H_prog[4][i] = 1.250; H_prog[5][i] = 0.676; H_prog[6][i] = 1.449;
    }
    for (i=9; i<17; i++){
      H_prog[1][i] = 0.213; H_prog[2][i] = 0.535; H_prog[3][i] = 0.855; H_prog[4][i] = 1.818; H_prog[5][i] = 0.952; H_prog[6][i] = 2.000;
    }

    double H_mort[7][17]; /* Mortality due to HIV (no ART) (age, CD4) */
    
    H_mort[0][0] = 0.312; H_mort[1][0] = 0.382; H_mort[2][0] = 0.466; H_mort[3][0] = 0.466; H_mort[4][0] = 0.569; H_mort[5][0] = 0.569; H_mort[6][0] = 0.569;
    
    for (i=1; i<3; i++){
      H_mort[0][i] = 0.039; H_mort[1][i] = 0.048; H_mort[2][i] = 0.058; H_mort[3][i] = 0.058; H_mort[4][i] = 0.071; H_mort[5][i] = 0.071; H_mort[6][i] = 0.071;
    }
    for (i=3; i<5; i++){
      H_mort[0][i] = 0.005; H_mort[1][i] = 0.011; H_mort[2][i] = 0.026; H_mort[3][i] = 0.061; H_mort[4][i] = 0.139; H_mort[5][i] = 0.321; H_mort[6][i] = 0.737;
    }    
    for (i=5; i<7; i++){
      H_mort[0][i] = 0.004; H_mort[1][i] = 0.010; H_mort[2][i] = 0.026; H_mort[3][i] = 0.069; H_mort[4][i] = 0.185; H_mort[5][i] = 0.499; H_mort[6][i] = 1.342;
    } 
    for (i=7; i<9; i++){
      H_mort[0][i] = 0.005; H_mort[1][i] = 0.013; H_mort[2][i] = 0.036; H_mort[3][i] = 0.096; H_mort[4][i] = 0.258; H_mort[5][i] = 0.691; H_mort[6][i] = 1.851;
    } 
    for (i=9; i<17; i++){
      H_mort[0][i] = 0.005; H_mort[1][i] = 0.013; H_mort[2][i] = 0.032; H_mort[3][i] = 0.080; H_mort[4][i] = 0.203; H_mort[5][i] = 0.513; H_mort[6][i] = 1.295;
    } 

    double A_mort[3][7][17]; /* On ART mortality (age,starting CD4, time on ART)  - this is an average of male and female values weigthed by sex of those on ART  */
    
    A_mort[0][0][0] = 0.062; A_mort[0][1][0] = 0.282;  A_mort[0][2][0] = 0.214; A_mort[0][3][0] = 0.214; A_mort[0][4][0] = 0.749;  A_mort[0][5][0] = 0.749;  A_mort[0][6][0] = 0.749; 
    A_mort[1][0][0] = 0.147; A_mort[1][1][0] = 0.210;  A_mort[1][2][0] = 0.199; A_mort[1][3][0] = 0.199; A_mort[1][4][0] = 0.376;  A_mort[1][5][0] = 0.376;  A_mort[1][6][0] = 0.376; 
    A_mort[2][0][0] = 0.062; A_mort[2][1][0] = 0.089;  A_mort[2][2][0] = 0.084; A_mort[2][3][0] = 0.084; A_mort[2][4][0] = 0.159;  A_mort[2][5][0] = 0.159;  A_mort[2][6][0] = 0.159; 
    
    A_mort[0][0][1] = 0.008; A_mort[0][1][1] = 0.035;  A_mort[0][2][1] = 0.027; A_mort[0][3][1] = 0.027; A_mort[0][4][1] = 0.094;  A_mort[0][5][1] = 0.094;  A_mort[0][6][1] = 0.094; 
    A_mort[1][0][1] = 0.018; A_mort[1][1][1] = 0.026;  A_mort[1][2][1] = 0.025; A_mort[1][3][1] = 0.025; A_mort[1][4][1] = 0.047;  A_mort[1][5][1] = 0.047;  A_mort[1][6][1] = 0.047; 
    A_mort[2][0][1] = 0.008; A_mort[2][1][1] = 0.011;  A_mort[2][2][1] = 0.011; A_mort[2][3][1] = 0.011; A_mort[2][4][1] = 0.020;  A_mort[2][5][1] = 0.020;  A_mort[2][6][1] = 0.020; 
    
    A_mort[0][0][2] = 0.007; A_mort[0][1][2] = 0.033;  A_mort[0][2][2] = 0.025; A_mort[0][3][2] = 0.025; A_mort[0][4][2] = 0.088;  A_mort[0][5][2] = 0.088;  A_mort[0][6][2] = 0.088; 
    A_mort[1][0][2] = 0.009; A_mort[1][1][2] = 0.012;  A_mort[1][2][2] = 0.012; A_mort[1][3][2] = 0.012; A_mort[1][4][2] = 0.022;  A_mort[1][5][2] = 0.022;  A_mort[1][6][2] = 0.022; 
    A_mort[2][0][2] = 0.004; A_mort[2][1][2] = 0.005;  A_mort[2][2][2] = 0.005; A_mort[2][3][2] = 0.005; A_mort[2][4][2] = 0.009;  A_mort[2][5][2] = 0.009;  A_mort[2][6][2] = 0.009; 
    
    for (i=3; i<5; i++){
      A_mort[0][0][i] = 0.0050; A_mort[0][1][i] = 0.0115;  A_mort[0][2][i] = 0.0264; A_mort[0][3][i] = 0.0557; A_mort[0][4][i] = 0.0953;  A_mort[0][5][i] = 0.1553;  A_mort[0][6][i] = 0.3418; 
      A_mort[1][0][i] = 0.0050; A_mort[1][1][i] = 0.0115;  A_mort[1][2][i] = 0.0204; A_mort[1][3][i] = 0.0216; A_mort[1][4][i] = 0.0272;  A_mort[1][5][i] = 0.0343;  A_mort[1][6][i] = 0.0496; 
      A_mort[2][0][i] = 0.0050; A_mort[2][1][i] = 0.0068;  A_mort[2][2][i] = 0.0073; A_mort[2][3][i] = 0.0079; A_mort[2][4][i] = 0.0105;  A_mort[2][5][i] = 0.0139;  A_mort[2][6][i] = 0.0211; 
    }
    for (i=5; i<7; i++){
      A_mort[0][0][i] = 0.0035; A_mort[0][1][i] = 0.0095;  A_mort[0][2][i] = 0.0256; A_mort[0][3][i] = 0.0477; A_mort[0][4][i] = 0.0792;  A_mort[0][5][i] = 0.1305;  A_mort[0][6][i] = 0.2897; 
      A_mort[1][0][i] = 0.0035; A_mort[1][1][i] = 0.0095;  A_mort[1][2][i] = 0.0242; A_mort[1][3][i] = 0.0278; A_mort[1][4][i] = 0.0351;  A_mort[1][5][i] = 0.0442;  A_mort[1][6][i] = 0.0640; 
      A_mort[2][0][i] = 0.0035; A_mort[2][1][i] = 0.0081;  A_mort[2][2][i] = 0.0093; A_mort[2][3][i] = 0.0100; A_mort[2][4][i] = 0.0134;  A_mort[2][5][i] = 0.0177;  A_mort[2][6][i] = 0.0271; 
    }
    for (i=7; i<9; i++){
      A_mort[0][0][i] = 0.0050; A_mort[0][1][i] = 0.0134;  A_mort[0][2][i] = 0.0359; A_mort[0][3][i] = 0.0495; A_mort[0][4][i] = 0.0833;  A_mort[0][5][i] = 0.1384;  A_mort[0][6][i] = 0.3094; 
      A_mort[1][0][i] = 0.0050; A_mort[1][1][i] = 0.0134;  A_mort[1][2][i] = 0.0267; A_mort[1][3][i] = 0.0283; A_mort[1][4][i] = 0.0362;  A_mort[1][5][i] = 0.0461;  A_mort[1][6][i] = 0.0675; 
      A_mort[2][0][i] = 0.0050; A_mort[2][1][i] = 0.0076;  A_mort[2][2][i] = 0.0083; A_mort[2][3][i] = 0.0091; A_mort[2][4][i] = 0.0128;  A_mort[2][5][i] = 0.0174;  A_mort[2][6][i] = 0.0275; 
    }
    for (i=9; i<17; i++){
      A_mort[0][0][i] = 0.0050; A_mort[0][1][i] = 0.0126;  A_mort[0][2][i] = 0.0319; A_mort[0][3][i] = 0.0489; A_mort[0][4][i] = 0.0872;  A_mort[0][5][i] = 0.1496;  A_mort[0][6][i] = 0.3431; 
      A_mort[1][0][i] = 0.0050; A_mort[1][1][i] = 0.0126;  A_mort[1][2][i] = 0.0306; A_mort[1][3][i] = 0.0366; A_mort[1][4][i] = 0.0481;  A_mort[1][5][i] = 0.0625;  A_mort[1][6][i] = 0.0935; 
      A_mort[2][0][i] = 0.0045; A_mort[2][1][i] = 0.0066;  A_mort[2][2][i] = 0.0076; A_mort[2][3][i] = 0.0087; A_mort[2][4][i] = 0.0141;  A_mort[2][5][i] = 0.0209;  A_mort[2][6][i] = 0.0355; 
    }
    
    double A_prog[4] = {0,2,2,0}; /* Progression through time on ART, 6 monthly time blocks - 0 ensure no progression into first catergory and no progression out of last category*/
    double A_start[3] = {1,0,0};  /* Used to make sure ART initiations are only added to the fist time on ART box */ 
    
    /* sum up various totals */

    /* Use sumsum function to add up HIV- */
    double Total_S = sumsum(S,0,16);                        /* Total susceptible */
    double Total_Ls = sumsum(Lsn,0,16)+sumsum(Lsp,0,16);    /* Total LTBI with drug susceptible (DS) strain */
    double Total_Lm = sumsum(Lmn,0,16)+sumsum(Lmp,0,16);    /* Total LTBI with drug resistant (DR) strain */
    double Total_Ns_N = sumsum(Nsn,0,16)+sumsum(Nsp,0,16);    /* Total DS smear negative TB */
    double Total_Nm_N = sumsum(Nmn,0,16)+sumsum(Nmp,0,16);    /* Total DR smear negative TB */
    double Total_Is_N = sumsum(Isn,0,16)+sumsum(Isp,0,16);    /* Total DS smear positive TB */
    double Total_Im_N = sumsum(Imn,0,16)+sumsum(Imp,0,16);    /* Total DR smear positive TB */
    double Total_PT = sumsum(PTn,0,16)+sumsum(PTp,0,16);      /* Post PT */
    
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
    double TB_deaths_neg[17];
    double TB_deaths_HIV[17][7];
    double TB_deaths_ART[17][7][3];
    double TB_deaths_HIV_age[17] = {0};
    double TB_deaths_ART_age[17] = {0};
    double TB_deaths[17];
    
    double HIV_deaths_HIV[17] = {0}; ;
    double HIV_deaths_ART[17] = {0};

    double up_H_mort[7][17];
    double up_A_mort[3][7][17];
    double m_b[17];
    double rate_dis_death[17];
    
    double tot_age[17] = {0};
    double tot_age_HIV[17][7];
    double tot_age_ART[17][7][3];
    
    /*double Tot_deaths = 0;*/
    double Tot_deaths_age[17];
    double ART_deaths_age[17] = {0};
    double Tot_deaths=0;
    double tot_age_neg[17] = {0};
    
 
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
          TB_deaths_ART[i][j][l] = (Nsn_A[i][j][l]+Nsp_A[i][j][l]+Nmn_A[i][j][l]+Nmp_A[i][j][l])*muN_H_A[i][l] + 
                                   (Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l])*muI_H_A[i][l]; 
          TB_deaths_ART_age[i] = TB_deaths_ART_age[i] + TB_deaths_ART[i][j][l];
          TB_deaths[i] = TB_deaths[i] + TB_deaths_ART[i][j][l];
          /* Calculate size of ART age group */
          tot_age_ART[i][j][l] = S_A[i][j][l]+Lsn_A[i][j][l]+Lsp_A[i][j][l]+Lmn_A[i][j][l]+Lmp_A[i][j][l]+Nsn_A[i][j][l]+Nsp_A[i][j][l]+
                                 Nmn_A[i][j][l]+Nmp_A[i][j][l]+Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l]+PTn_A[i][j][l]+PTp_A[i][j][l];
          /* Update size of age group */
          tot_age[i] = tot_age[i] + tot_age_ART[i][j][l];
          /* Adjust ART mortality probability to remove TB deaths (only if there is any ART yet (otherwise we get a divide by 0 error)) */
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
          ART_deaths_age[i] = ART_deaths_age[i] + TB_deaths_ART[i][j][l] + up_A_mort[l][j][i]*(S_A[i][j][l]+Lsn_A[i][j][l]+Lsp_A[i][j][l]+Lmn_A[i][j][l]+Lmp_A[i][j][l]+
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
          ART_deaths_age[i] = ART_deaths_age[i] + m_b[i]*(S_A[i][j][l]+Lsn_A[i][j][l]+Lsp_A[i][j][l]+Lmn_A[i][j][l]+Lmp_A[i][j][l]+Nsn_A[i][j][l]+Nsp_A[i][j][l]+
                                                  Nmn_A[i][j][l]+Nmp_A[i][j][l]+Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l]+PTn_A[i][j][l]+PTp_A[i][j][l]);
        }
      }
      
    } 
    double TB_deaths_tot = sumsum(TB_deaths,0,16);
    double ART_deaths = sumsum(ART_deaths_age,0,16);
    
    /* Sum up populations over CD4 categories, with and without ART and calculate rates of ART initiation by age */
    
    double ART_prop[17][7] = {{0}};     /* Proportion of CD4 category who should start ART by age */
    double temp1[17][7] = {{0}};
    double CD4_dist[17][7] = {{0}};     /* Not on ART by CD4 and age */
    double CD4_dist_ART[17][7] = {{0}}; /* On ART by CD4 and age*/
    double CD4_dist_all[7] = {0};       /* Not on ART by CD4 */
    double CD4_dist_ART_all[7] = {0};   /* On ART by CD4 */
    double CD4_deaths[17][7] = {{0}};   /* Deaths by CD4 (no ART) */
    double ART_new[17] = {0};           /* Number of new people to put on ART by age */
    double ART_el[17] = {0};            /* Number who are eligible but not on ART */
    double ART_el_deaths[17] = {0};     /* Number eligible who will die */
    double ART_on[17] = {0};            /* Number who should be on ART by age */
    double Tot_ART[17] = {0};           /* Number currently on ART by age */
    
    for (i=0; i<n_age; i++){

      for (j=0; j<n_HIV; j++){
        
        CD4_dist[i][j] = S_H[i][j]+Lsn_H[i][j]+Lsp_H[i][j]+Lmn_H[i][j]+Lmp_H[i][j]+
                         Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j]+ 
                         Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j]+
                         PTn_H[i][j]+PTp_H[i][j];
                      
        CD4_dist_all[j] = CD4_dist_all[j] + CD4_dist[i][j];
        
        CD4_deaths[i][j] = H_mort[j][i]*(S_H[i][j]+Lsn_H[i][j]+Lsp_H[i][j]+Lmn_H[i][j]+Lmp_H[i][j]+
                                         Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j]+ 
                                         Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j]+
                                         PTn_H[i][j]+PTp_H[i][j]);
                                                          
                                                          
        for (l=0; l<n_ART; l++){
          
          CD4_dist_ART[i][j] = CD4_dist_ART[i][j]+S_A[i][j][l]+Lsn_A[i][j][l]+Lsp_A[i][j][l]+Lmn_A[i][j][l]+Lmp_A[i][j][l]+
                                                  Nsn_A[i][j][l]+Nsp_A[i][j][l]+Nmn_A[i][j][l]+Nmp_A[i][j][l]+ 
                                                  Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l]+
                                                  PTn_A[i][j][l]+PTp_A[i][j][l];                
        
        } 
        
        CD4_dist_ART_all[j] = CD4_dist_ART_all[j] + CD4_dist_ART[i][j];
        
        Tot_ART[i] = Tot_ART[i] + CD4_dist_ART[i][j];  /* sum up number currently on ART */
         
        if (j>=Athresh) { /* If this CD4 is eligible for ART */
          ART_on[i] = ART_on[i] + (CD4_dist[i][j] + CD4_dist_ART[i][j])*forc[i+82]; /* number who should be on ART - HIV+ population times coverage (by CD4) */
        }

      }
      ART_new[i] = fmax(0,ART_on[i] - (Tot_ART[i] - ART_deaths_age[i]));   /* number who need to start is number who should be on minus those already on plus those on ART who will die in current time */ 
        
      /* Then work out where these should go by CD4 - based on proportion of eligible population in CD4 group and proportion of deaths which occur in CD4 group */
      for (j=Athresh; j<n_HIV; j++) {
        ART_el[i] = ART_el[i] + CD4_dist[i][j];
        ART_el_deaths[i] = ART_el_deaths[i] + CD4_deaths[i][j];
      }
      
      /* check that number to be put on isn't greater than number eligbile */
      ART_new[i] = fmin(ART_el[i],ART_new[i]);
      
      if (ART_el[i] > 0){
        for (j=Athresh; j<n_HIV; j++) {
          if (CD4_dist[i][j] > 0) {
            
            ART_prop[i][j] = (((CD4_dist[i][j]/ART_el[i])+(CD4_deaths[i][j]/ART_el_deaths[i]))/2)*(ART_new[i]/CD4_dist[i][j]); /* applies weighting and size of CD4 group to work out % of CD4 group that should move */
            /*ART_prop[i][j] = 0.2;*/
          }
        }
      }
    }

    /* Force of infection */
    double FS = beta*(Total_Ns_N*rel_inf + Total_Ns_H*rel_inf_H + Total_Is_N + Total_Is_H)/Total; 
    double FM = fit_cost*beta*(Total_Nm_N*rel_inf + Total_Nm_H*rel_inf_H + Total_Im_N + Total_Im_H)/Total; 
    
    /* Variables to store numbers of new cases */
    double TB_cases_age[17] = {0};
    double TB_cases_neg_age[17];
    double TB_cases_neg = 0;
    double TB_cases_pos_age[17][7];
    double TB_cases_pos = 0;
    double TB_cases_ART_age[17][7][3];
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
      dS[i] = - (FS + FM)*S[i] - forc[i+18]*S[i] - m_b[i]*S[i] + (S[i]/tot_age[i])*(forc[i+52]);
       
      /* Latent, ds, naive */
      dLsn[i] = - m_b[i]*Lsn[i] + S_to_Lsn + Lmn_to_Lsn + r*(Isn[i] + Nsn[i]) - Lsn_to_Lmn - Lsn_to_Nsn - Lsn_to_Isn  - Lsn_to_Nmn - Lsn_to_Imn -
                forc[i+18]*Lsn[i] + (Lsn[i]/tot_age[i])*(forc[i+52]) + PTn_to_Lsn -
                health*kneg*(1-sp_I_neg*sp_N_neg)*l_s*tneg_s*Lsn[i];                   
      
      /* Latent, ds, prev */
      dLsp[i] = - m_b[i]*Lsp[i] + Lmp_to_Lsp + r*(Isp[i] + Nsp[i]) - Lsp_to_Lmp - Lsp_to_Nsp - Lsp_to_Isp  - Lsp_to_Nmp - Lsp_to_Imp +     
                kneg*(l_s*(1-e)*tneg_s*((dstneg_p*sp_m_neg)+(1-dstneg_p)) + l_m*tneg_m*dstneg_p*(1-sp_m_neg))*(se_I_neg*Isp[i] + se_N_neg*rel_d*Nsp[i]) + 
                kneg*(l_s*(1-e)*tneg_s*((dstneg_n*sp_m_neg)+(1-dstneg_n)) + l_m*tneg_m*dstneg_n*(1-sp_m_neg))*(se_I_neg*Isn[i] + se_N_neg*rel_d*Nsn[i]) - /* Care */
                forc[i+18]*Lsp[i] + (Lsp[i]/tot_age[i])*(forc[i+52]) + PTp_to_Lsp -
                health*kneg*(1-sp_I_neg*sp_N_neg)*l_s*tneg_s*Lsp[i];         

      /* Latent, mdr, naive */ 
      dLmn[i] = - m_b[i]*Lmn[i] + S_to_Lmn + Lsn_to_Lmn + r*(Imn[i] + Nmn[i]) - Lmn_to_Lsn - Lmn_to_Nsn - Lmn_to_Isn - Lmn_to_Nmn - Lmn_to_Imn  -                
                forc[i+18]*Lmn[i] + (Lmn[i]/tot_age[i])*(forc[i+52]) + PTn_to_Lmn;              
      
      /* Latent, mdr, prev */
      dLmp[i] = - m_b[i]*Lmp[i] + Lsp_to_Lmp + r*(Imp[i] + Nmp[i]) - Lmp_to_Lsp - Lmp_to_Nsp - Lmp_to_Isp - Lmp_to_Nmp - Lmp_to_Imp +
                kneg*(dstneg_p*se_m_neg*l_m*tneg_m + l_s*tneg_s*eff_p*((1-dstneg_p)+dstneg_p*(1-se_m_neg)))*(se_I_neg*Imp[i] + se_N_neg*rel_d*Nmp[i]) +
                kneg*(dstneg_n*se_m_neg*l_m*tneg_m + l_s*tneg_s*eff_n*((1-dstneg_n)+dstneg_n*(1-se_m_neg)))*(se_I_neg*Imn[i] + se_N_neg*rel_d*Nmn[i]) - /* Care */
                forc[i+18]*Lmp[i] + (Lmp[i]/tot_age[i])*(forc[i+52]) + PTp_to_Lmp;  
                
      /* Smear neg, ds, new */
      dNsn[i] = - m_b[i]*Nsn[i] + S_to_Nsn + Lsn_to_Nsn + Lmn_to_Nsn - (theta + r + muN_age[i])*Nsn[i] -
                kneg*se_N_neg*rel_d*(l_s*(dstneg_n*sp_m_neg + (1-dstneg_n)) + dstneg_n*(1-sp_m_neg)*l_m)*Nsn[i] - /* care */ 
                forc[i+18]*Nsn[i] + (Nsn[i]/tot_age[i])*(forc[i+52]) + PTn_to_Nsn;  
      
      /* Smear neg, ds, prev */                             
      dNsp[i] = - m_b[i]*Nsp[i] + Lsp_to_Nsp + Lmp_to_Nsp - (theta + r + muN_age[i])*Nsp[i] -
                kneg*rel_d*se_N_neg*((1-e)*tneg_s*l_s*(sp_m_neg*dstneg_p + (1-dstneg_p)) + dstneg_p*(1-sp_m_neg)*l_m*tneg_m + l_s*e*(sp_m_neg*dstneg_p +(1-dstneg_p)))*Nsp[i] + 
                kneg*rel_d*se_N_neg*(l_s*(1-e)*(1-tneg_s)*(dstneg_n*sp_m_neg + (1-dstneg_n)) + dstneg_n*(1-sp_m_neg)*(1-tneg_m)*l_m)*Nsn[i] - /* Care */
                forc[i+18]*Nsp[i] + (Nsp[i]/tot_age[i])*(forc[i+52]) + PTp_to_Nsp;  

      /* Smear neg, mdr, new */
      dNmn[i] = - m_b[i]*Nmn[i] + S_to_Nmn + Lsn_to_Nmn + Lmn_to_Nmn - (theta + r + muN_age[i])*Nmn[i] -
                kneg*rel_d*se_N_neg*(l_m*dstneg_n*se_m_neg + l_s*((1-dstneg_n)+dstneg_n*(1-se_m_neg)))*Nmn[i] - /* Care */
                forc[i+18]*Nmn[i] + (Nmn[i]/tot_age[i])*(forc[i+52]) + PTn_to_Nmn;   
                   
      /* Smear neg, mdr, prev */
      dNmp[i] = - m_b[i]*Nmp[i] + Lsp_to_Nmp + Lmp_to_Nmp - (theta + r + muN_age[i])*Nmp[i] -
                kneg*rel_d*se_N_neg*(dstneg_p*se_m_neg*l_m*tneg_m + l_s*tneg_s*eff_p*((1-dstneg_p)+dstneg_p*(1-se_m_neg)))*Nmp[i] + 
                kneg*l_s*rel_d*e*se_N_neg*((dstneg_n*sp_m_neg+(1-dstneg_n))*Nsn[i]+(dstneg_p*sp_m_neg+(1-dstneg_p))*Nsp[i]) + 
                kneg*se_N_neg*rel_d*(dstneg_n*l_m*se_m_neg*(1-tneg_m) + l_s*(1-(tneg_s*eff_n))*((1-dstneg_n)+dstneg_n*(1-se_m_neg)))*Nmn[i] - /* Care */
                forc[i+18]*Nmp[i] + (Nmp[i]/tot_age[i])*(forc[i+52]) + PTp_to_Nmp;  

      /* Smear pos, ds, new */
      dIsn[i] = - m_b[i]*Isn[i] + S_to_Isn + Lsn_to_Isn + Lmn_to_Isn + theta*Nsn[i] - (r + muI_age[i])*Isn[i] - 
                kneg*se_I_neg*((dstneg_n*sp_m_neg + (1-dstneg_n))*l_s +(dstneg_n*(1-sp_m_neg)*l_m))*Isn[i] - /* Care */  
                forc[i+18]*Isn[i] + (Isn[i]/tot_age[i])*(forc[i+52]) + PTn_to_Isn;  
      
      /* Smear pos, ds, prev */
      dIsp[i] = - m_b[i]*Isp[i] + Lsp_to_Isp + Lmp_to_Isp + theta*Nsp[i] - (r + muI_age[i])*Isp[i] -   
                kneg*se_I_neg*(l_s*(1-e)*tneg_s*(sp_m_neg*dstneg_p + (1-dstneg_p)) + dstneg_p*(1-sp_m_neg)*l_m*tneg_m + l_s*e*(sp_m_neg*dstneg_p + (1-dstneg_p)))*Isp[i] + 
                kneg*se_I_neg*(l_s*(1-e)*(1-tneg_s)*(dstneg_n*sp_m_neg + (1-dstneg_n)) + (dstneg_n*(1-sp_m_neg)*l_m*(1-tneg_m)))*Isn[i] - /* Care */
                forc[i+18]*Isp[i] + (Isp[i]/tot_age[i])*(forc[i+52]) + PTp_to_Isp; 

      /* Smear pos, mdr, new */
      dImn[i] = - m_b[i]*Imn[i] + S_to_Imn + Lsn_to_Imn + Lmn_to_Imn + theta*Nmn[i] - (r + muI_age[i])*Imn[i] -
                kneg*se_I_neg*(l_m*dstneg_n*se_m_neg + (1-dstneg_n)*l_s + dstneg_n*(1-se_m_neg)*l_s)*Imn[i] - /* Care */ 
                forc[i+18]*Imn[i] + (Imn[i]/tot_age[i])*(forc[i+52]) + PTn_to_Imn;  
                
                
      /* Smear pos, mdr, prev */
      dImp[i] = - m_b[i]*Imp[i] + Lsp_to_Imp + Lmp_to_Imp + theta*Nmp[i] - (r + muI_age[i])*Imp[i] -
                kneg*se_I_neg*(l_m*dstneg_p*tneg_m*se_m_neg + l_s*tneg_s*eff_p*((1-dstneg_p)+dstneg_p*(1-se_m_neg)))*Imp[i] + 
                kneg*se_I_neg*(se_m_neg*l_m*dstneg_n*(1-tneg_m) + l_s*(1-(tneg_s*eff_n))*((1-dstneg_n)+dstneg_n*(1-se_m_neg)))*Imn[i] + 
                kneg*se_I_neg*l_s*e*((dstneg_n*sp_m_neg + (1-dstneg_n))*Isn[i]+(dstneg_p*sp_m_neg+ (1-dstneg_p))*Isp[i]) -
                forc[i+18]*Imp[i] + (Imp[i]/tot_age[i])*(forc[i+52]) + PTp_to_Imp;  
                
      /* Post PT, ds, new */          
      dPTn[i] = - m_b[i]*PTn[i] - PTn_to_Lsn - PTn_to_Nsn - PTn_to_Isn - PTn_to_Lmn - PTn_to_Nmn - PTn_to_Imn +
                health*kneg*(1-sp_I_neg*sp_N_neg)*l_s*tneg_s*Lsn[i] -
                forc[i+18]*PTn[i] + (PTn[i]/tot_age[i])*(forc[i+52]);
      
      /* Post PT, ds, prev */ 
      dPTp[i] = - m_b[i]*PTp[i] - PTp_to_Lsp - PTp_to_Nsp - PTp_to_Isp - PTp_to_Lmp - PTp_to_Nmp - PTp_to_Imp +
                health*kneg*(1-sp_I_neg*sp_N_neg)*l_s*tneg_s*Lsp[i] -
                forc[i+18]*PTp[i] + (PTp[i]/tot_age[i])*(forc[i+52]);
                        
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
          double PTnH_to_ImnH = FM*a_age_H[i][j]*(1-p_H[j])*sig_H*PTn_H[i][j];               /* Post PT to smear positive DR disease (no disease history) */
      
          double PTpH_to_LspH = FS*(1-a_age_H[i][j])*(1-p_H[j])*PTp_H[i][j];            /* Post PT to latent DS (prior Rx) */
          double PTpH_to_NspH = FS*a_age_H[i][j]*(1-p_H[j])*(1-sig_H)*PTp_H[i][j];      /* Post PT to smear negative DS disease (prior Rx) */
          double PTpH_to_IspH = FS*a_age_H[i][j]*(1-p_H[j])*sig_H*PTp_H[i][j];          /* Post PT to smear positive DS disease (prior Rx) */
          double PTpH_to_LmpH = FM*(1-a_age_H[i][j])*(1-p_H[j])*g*PTp_H[i][j];          /* Post PT to latent DR (prior Rx) */
          double PTpH_to_NmpH = FM*a_age_H[i][j]*(1-p_H[j])*(1-sig_H)*PTp_H[i][j];      /* Post PT to smear negative DR disease (prior Rx) */
          double PTpH_to_ImpH = FM*a_age_H[i][j]*(1-p_H[j])*sig_H*PTp_H[i][j];          /* Post PT to smear positive DR disease (prior Rx) */

          dS_H[i][j] = - m_b[i]*S_H[i][j] - 
                      (FS + FM)*S_H[i][j] +
                      forc[i+18]*H_CD4[j][i]*S[i] - H_prog[j+1][i]*S_H[i][j] + H_prog[j][i]*S_H[i][j-1] - 
                      up_H_mort[j][i]*S_H[i][j] - ART_prop[i][j]*S_H[i][j] + (S_H[i][j]/tot_age[i])*(forc[i+52]);

          dLsn_H[i][j] = - m_b[i]*Lsn_H[i][j] +
                        FS*((1-a_age_H[i][j])*S_H[i][j] + (1-a_age_H[i][j])*(1-p_H[j])*(1-g)*Lmn_H[i][j]) -
                        (v_age_H[i][j] + FS*a_age_H[i][j]*(1-p_H[j]) + FM*a_age_H[i][j]*(1-p_H[j]) + FM*(1-a_age_H[i][j])*(1-p_H[j])*g)*Lsn_H[i][j] +
                        r_H*(Isn_H[i][j] + Nsn_H[i][j]) +
                        forc[i+18]*H_CD4[j][i]*Lsn[i] - H_prog[j+1][i]*Lsn_H[i][j] + H_prog[j][i]*Lsn_H[i][j-1] - 
                        up_H_mort[j][i]*Lsn_H[i][j] - ART_prop[i][j]*Lsn_H[i][j] + (Lsn_H[i][j]/tot_age[i])*(forc[i+52])-
                        health*kpos*(1-sp_I_pos*sp_N_pos)*l_s*tpos_s*Lsn_H[i][j] + PTnH_to_LsnH;

          dLsp_H[i][j] = - m_b[i]*Lsp_H[i][j] +
                        FS*(1-a_age_H[i][j])*(1-p_H[j])*(1-g)*Lmp_H[i][j] -
                        (v_age_H[i][j] + FS*a_age_H[i][j]*(1-p_H[j]) + FM*a_age_H[i][j]*(1-p_H[j]) + FM*(1-a_age_H[i][j])*(1-p_H[j])*g)*Lsp_H[i][j] +
                        r_H*(Isp_H[i][j]+Nsp_H[i][j]) +
                        kpos*(l_s*(1-e)*tpos_s*((dstpos_p*sp_m_pos)+(1-dstpos_p)) + l_m*tpos_m*dstpos_p*(1-sp_m_pos))*(se_I_pos*Isp_H[i][j] + se_N_pos*rel_d*Nsp_H[i][j]) + 
                        kpos*(l_s*(1-e)*tpos_s*((dstpos_n*sp_m_pos)+(1-dstpos_n)) + l_m*tpos_m*dstpos_n*(1-sp_m_pos))*(se_I_pos*Isn_H[i][j] + se_N_pos*rel_d*Nsn_H[i][j]) +
                        forc[i+18]*H_CD4[j][i]*Lsp[i] - H_prog[j+1][i]*Lsp_H[i][j] + H_prog[j][i]*Lsp_H[i][j-1] - 
                        up_H_mort[j][i]*Lsp_H[i][j] - ART_prop[i][j]*Lsp_H[i][j] + (Lsp_H[i][j]/tot_age[i])*(forc[i+52])-
                        health*kpos*(1-sp_I_pos*sp_N_pos)*l_s*tpos_s*Lsp_H[i][j] + PTpH_to_LspH;

          dLmn_H[i][j] = - m_b[i]*Lmn_H[i][j] +
                        FM*((1-a_age_H[i][j])*S_H[i][j] + (1-a_age_H[i][j])*(1-p_H[j])*g*Lsn_H[i][j]) -
                        (v_age_H[i][j] + FM*a_age_H[i][j]*(1-p_H[j]) + FS*a_age_H[i][j]*(1-p_H[j]) + FS*(1-a_age_H[i][j])*(1-p_H[j])*(1-g))*Lmn_H[i][j] +
                        r_H*(Imn_H[i][j]+Nmn_H[i][j]) +
                        forc[i+18]*H_CD4[j][i]*Lmn[i] - H_prog[j+1][i]*Lmn_H[i][j] + H_prog[j][i]*Lmn_H[i][j-1] - 
                        up_H_mort[j][i]*Lmn_H[i][j] - ART_prop[i][j]*Lmn_H[i][j] + (Lmn_H[i][j]/tot_age[i])*(forc[i+52]) + PTnH_to_LmnH;
                     
          dLmp_H[i][j] = - m_b[i]*Lmp_H[i][j] + 
                        FM*(1-a_age_H[i][j])*(1-p_H[j])*g*Lsp_H[i][j] -
                        (v_age_H[i][j] + FM*a_age_H[i][j]*(1-p_H[j]) + FS*a_age_H[i][j]*(1-p_H[j]) + FS*(1-a_age_H[i][j])*(1-p_H[j])*(1-g))*Lmp_H[i][j] +
                        r_H*(Imp_H[i][j]+Nmp_H[i][j]) +
                        kpos*(dstpos_p*se_m_pos*l_m*tpos_m + l_s*tpos_s*eff_p*((1-dstpos_p)+dstpos_p*(1-se_m_pos)))*(se_I_pos*Imp_H[i][j] + se_N_pos*rel_d*Nmp_H[i][j]) +
                        kpos*(dstpos_n*se_m_pos*l_m*tpos_m + l_s*tpos_s*eff_n*((1-dstpos_n)+dstpos_n*(1-se_m_pos)))*(se_I_pos*Imn_H[i][j] + se_N_pos*rel_d*Nmn_H[i][j]) + 
                        forc[i+18]*H_CD4[j][i]*Lmp[i] - H_prog[j+1][i]*Lmp_H[i][j] + H_prog[j][i]*Lmp_H[i][j-1] - 
                        up_H_mort[j][i]*Lmp_H[i][j] - ART_prop[i][j]*Lmp_H[i][j] + (Lmp_H[i][j]/tot_age[i])*(forc[i+52]) + PTpH_to_LmpH;

          dNsn_H[i][j] = - m_b[i]*Nsn_H[i][j] +
                        (v_age_H[i][j]*(1-sig_H) + FS*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lsn_H[i][j] +
                        FS*a_age_H[i][j]*(1-sig_H)*(S_H[i][j] + (1-p_H[j])*Lmn_H[i][j]) -
                        (theta_H + r_H + muN_H)*Nsn_H[i][j] -
                        kpos*se_N_pos*rel_d*(l_s*(dstpos_n*sp_m_pos + (1-dstpos_n)) + dstpos_n*(1-sp_m_pos)*l_m)*Nsn_H[i][j] + 
                        forc[i+18]*H_CD4[j][i]*Nsn[i] - H_prog[j+1][i]*Nsn_H[i][j] + H_prog[j][i]*Nsn_H[i][j-1] - 
                        up_H_mort[j][i]*Nsn_H[i][j] - ART_prop[i][j]*Nsn_H[i][j] + (Nsn_H[i][j]/tot_age[i])*(forc[i+52]) + PTnH_to_NsnH;
                                                
          dNsp_H[i][j] = - m_b[i]*Nsp_H[i][j] + 
                        (v_age_H[i][j]*(1-sig_H) + FS*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lsp_H[i][j] +
                        FS*a_age_H[i][j]*(1-sig_H)*(1-p_H[j])*Lmp_H[i][j] -
                        (theta_H + r_H + muN_H)*Nsp_H[i][j] -
                        kpos*rel_d*se_N_pos*((1-e)*tpos_s*l_s*(sp_m_pos*dstpos_p + (1-dstpos_p)) + dstpos_p*(1-sp_m_pos)*l_m*tpos_m + l_s*e*(sp_m_pos*dstpos_p +(1-dstpos_p)))*Nsp_H[i][j] + 
                        kpos*rel_d*se_N_pos*(l_s*(1-e)*(1-tpos_s)*(dstpos_n*sp_m_pos + (1-dstpos_n)) + dstpos_n*(1-sp_m_pos)*(1-tpos_m)*l_m)*Nsn_H[i][j] +                       
                        forc[i+18]*H_CD4[j][i]*Nsp[i] - H_prog[j+1][i]*Nsp_H[i][j] + H_prog[j][i]*Nsp_H[i][j-1] - 
                        up_H_mort[j][i]*Nsp_H[i][j] - ART_prop[i][j]*Nsp_H[i][j] + (Nsp_H[i][j]/tot_age[i])*(forc[i+52]) + PTpH_to_NspH;

          dNmn_H[i][j] = - m_b[i]*Nmn_H[i][j] +
                        (v_age_H[i][j]*(1-sig_H) + FM*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lmn_H[i][j] +
                        FM*a_age_H[i][j]*(1-sig_H)*(S_H[i][j] + (1-p_H[j])*Lsn_H[i][j]) -
                        (theta_H + r_H + muN_H)*Nmn_H[i][j] -
                        kpos*rel_d*se_N_pos*(l_m*dstpos_n*se_m_pos + l_s*((1-dstpos_n)+dstpos_n*(1-se_m_pos)))*Nmn_H[i][j] +
                        forc[i+18]*H_CD4[j][i]*Nmn[i] - H_prog[j+1][i]*Nmn_H[i][j] + H_prog[j][i]*Nmn_H[i][j-1] - 
                        up_H_mort[j][i]*Nmn_H[i][j] - ART_prop[i][j]*Nmn_H[i][j] + (Nmn_H[i][j]/tot_age[i])*(forc[i+52]) + PTnH_to_NmnH;
                       
          dNmp_H[i][j] = - m_b[i]*Nmp_H[i][j] +
                        (v_age_H[i][j]*(1-sig_H) + FM*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lmp_H[i][j] +
                        FM*a_age_H[i][j]*(1-sig_H)*(1-p_H[j])*Lsp_H[i][j] -
                        (theta_H + r_H + muN_H)*Nmp_H[i][j] -
                        kpos*rel_d*se_N_pos*(dstpos_p*se_m_pos*l_m*tpos_m + l_s*tpos_s*eff_p*((1-dstpos_p)+dstpos_p*(1-se_m_pos)))*Nmp_H[i][j] + 
                        kpos*l_s*rel_d*e*se_N_pos*((dstpos_n*sp_m_pos+(1-dstpos_n))*Nsn_H[i][j]+(dstpos_p*sp_m_pos+(1-dstpos_p))*Nsp_H[i][j]) + 
                        kpos*se_N_pos*rel_d*(dstpos_n*l_m*se_m_pos*(1-tpos_m) + l_s*(1-(tpos_s*eff_n))*((1-dstpos_n)+dstpos_n*(1-se_m_pos)))*Nmn_H[i][j] +
                        forc[i+18]*H_CD4[j][i]*Nmp[i] - H_prog[j+1][i]*Nmp_H[i][j] + H_prog[j][i]*Nmp_H[i][j-1] - 
                        up_H_mort[j][i]*Nmp_H[i][j] - ART_prop[i][j]*Nmp_H[i][j] + (Nmp_H[i][j]/tot_age[i])*(forc[i+52]) + PTpH_to_NmpH;

          dIsn_H[i][j] = - m_b[i]*Isn_H[i][j] +
                        (v_age_H[i][j]*sig_H + FS*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lsn_H[i][j] +
                        FS*a_age_H[i][j]*sig_H*S_H[i][j] +
                        FS*a_age_H[i][j]*(1-p_H[j])*sig_H*Lmn_H[i][j] +
                        theta_H*Nsn_H[i][j] -
                        (r_H + muI_H)*Isn_H[i][j] -
                        kpos*se_I_pos*((dstpos_n*sp_m_pos + (1-dstpos_n))*l_s +(dstpos_n*(1-sp_m_pos)*l_m))*Isn_H[i][j] +
                        forc[i+18]*H_CD4[j][i]*Isn[i] - H_prog[j+1][i]*Isn_H[i][j] + H_prog[j][i]*Isn_H[i][j-1] - 
                        up_H_mort[j][i]*Isn_H[i][j] - ART_prop[i][j]*Isn_H[i][j] + (Isn_H[i][j]/tot_age[i])*(forc[i+52]) + PTnH_to_IsnH;

          dIsp_H[i][j] = - m_b[i]*Isp_H[i][j] +
                        (v_age_H[i][j]*sig_H + FS*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lsp_H[i][j] +
                        FS*a_age_H[i][j]*sig_H*(1-p_H[j])*Lmp_H[i][j] +
                        theta_H*Nsp_H[i][j] +
                        kpos*se_I_pos*(l_s*(1-e)*(1-tpos_s)*(dstpos_n*sp_m_pos + (1-dstpos_n)) + (dstpos_n*(1-sp_m_pos)*l_m*(1-tpos_m)))*Isn_H[i][j] - 
                        (r_H + muI_H)*Isp_H[i][j] -
                        kpos*se_I_pos*(l_s*(1-e)*tpos_s*(sp_m_pos*dstpos_p + (1-dstpos_p)) + dstpos_p*(1-sp_m_pos)*l_m*tpos_m + l_s*e*(sp_m_pos*dstpos_p + (1-dstpos_p)))*Isp_H[i][j] +
                        forc[i+18]*H_CD4[j][i]*Isp[i] - H_prog[j+1][i]*Isp_H[i][j] + H_prog[j][i]*Isp_H[i][j-1] - 
                        up_H_mort[j][i]*Isp_H[i][j] - ART_prop[i][j]*Isp_H[i][j] + (Isp_H[i][j]/tot_age[i])*(forc[i+52]) + PTpH_to_IspH;       
     
          dImn_H[i][j] = - m_b[i]*Imn_H[i][j] +
                        (v_age_H[i][j]*sig_H + FM*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lmn_H[i][j] +
                        FM*a_age_H[i][j]*sig_H*S_H[i][j] +
                        FM*a_age_H[i][j]*(1-p_H[j])*sig_H*Lsn_H[i][j] +
                        theta_H*Nmn_H[i][j] -
                        kpos*se_I_pos*(l_m*dstpos_n*se_m_pos + (1-dstpos_n)*l_s + dstpos_n*(1-se_m_pos)*l_s)*Imn_H[i][j] -
                        (r_H + muI_H)*Imn_H[i][j] +
                        forc[i+18]*H_CD4[j][i]*Imn[i] - H_prog[j+1][i]*Imn_H[i][j] + H_prog[j][i]*Imn_H[i][j-1] - 
                        up_H_mort[j][i]*Imn_H[i][j] - ART_prop[i][j]*Imn_H[i][j] + (Imn_H[i][j]/tot_age[i])*(forc[i+52]) + PTnH_to_ImnH;

          dImp_H[i][j] = - m_b[i]*Imp_H[i][j] +
                        (v_age_H[i][j]*sig_H + FM*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lmp_H[i][j] +
                        FM*a_age_H[i][j]*sig_H*(1-p_H[j])*Lsp_H[i][j] +
                        theta_H*Nmp_H[i][j] + 
                        kpos*se_I_pos*(se_m_pos*l_m*dstpos_n*(1-tpos_m) + l_s*(1-(tpos_s*eff_n))*((1-dstpos_n)+dstpos_n*(1-se_m_pos)))*Imn_H[i][j] + 
                        kpos*se_I_pos*l_s*e*((dstpos_n*sp_m_pos + (1-dstpos_n))*Isn_H[i][j]+(dstpos_p*sp_m_pos+ (1-dstpos_p))*Isp_H[i][j]) -
                        (r_H + muI_H)*Imp_H[i][j] -
                        kpos*se_I_pos*(l_m*dstpos_p*tpos_m*se_m_pos + l_s*tpos_s*eff_p*((1-dstpos_p)+dstpos_p*(1-se_m_pos)))*Imp_H[i][j]+
                        forc[i+18]*H_CD4[j][i]*Imp[i] - H_prog[j+1][i]*Imp_H[i][j] + H_prog[j][i]*Imp_H[i][j-1] - 
                        up_H_mort[j][i]*Imp_H[i][j] - ART_prop[i][j]*Imp_H[i][j] + (Imp_H[i][j]/tot_age[i])*(forc[i+52]) + PTpH_to_ImpH;       
                                                      
          dPTn_H[i][j] = -m_b[i]*PTn_H[i][j] - PTnH_to_LsnH - PTnH_to_NsnH - PTnH_to_IsnH - PTnH_to_LmnH - PTnH_to_NmnH - PTnH_to_ImnH +
                         health*kpos*(1-sp_I_pos*sp_N_pos)*l_s*tpos_s*Lsn_H[i][j] +
                         forc[i+18]*H_CD4[j][i]*PTn[i] - H_prog[j+1][i]*PTn_H[i][j] + H_prog[j][i]*PTn_H[i][j-1] - 
                         up_H_mort[j][i]*PTn_H[i][j] - ART_prop[i][j]*PTn_H[i][j] + (PTn_H[i][j]/tot_age[i])*(forc[i+52]);
          
          dPTp_H[i][j] = - m_b[i]*PTp_H[i][j] - PTpH_to_LspH - PTpH_to_NspH - PTpH_to_IspH - PTpH_to_LmpH - PTpH_to_NmpH - PTpH_to_ImpH +
                         health*kpos*(1-sp_I_pos*sp_N_pos)*l_s*tpos_s*Lsp_H[i][j] +
                         forc[i+18]*H_CD4[j][i]*PTp[i] - H_prog[j+1][i]*PTp_H[i][j] + H_prog[j][i]*PTp_H[i][j-1] - 
                         up_H_mort[j][i]*PTp_H[i][j] - ART_prop[i][j]*PTp_H[i][j] + (PTp_H[i][j]/tot_age[i])*(forc[i+52]);
                   
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
       
       /* HIV+ on ART: loop through time on ART, CD4 at initiation, age */
          for (l=0; l<n_ART; l++){
        
        
            double PTnA_to_LsnA = FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*PTn_A[i][j][l];            /* Post PT to latent DS (no disease history) */
            double PTnA_to_NsnA = FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H)*PTn_A[i][j][l];      /* Post PT to smear negative DS disease (no disease history) */
            double PTnA_to_IsnA = FS*a_age_A[i][j][l]*(1-p_A[j][l])*sig_H*PTn_A[i][j][l];          /* Post PT to smear positive DS disease (no disease history) */
            double PTnA_to_LmnA = FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g*PTn_A[i][j][l];          /* Post PT to latent DR (no disease history) */
            double PTnA_to_NmnA = FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H)*PTn_A[i][j][l];      /* Post PT to smear negative DR disease (no disease history) */
            double PTnA_to_ImnA = FM*a_age_A[i][j][l]*(1-p_A[j][l])*sig_H*PTn_A[i][j][l];          /* Post PT to smear positive DR disease (no disease history) */
      
            double PTpA_to_LspA = FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*PTp_A[i][j][l];            /* Post PT to latent DS (prior Rx) */
            double PTpA_to_NspA = FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H)*PTp_A[i][j][l];      /* Post PT to smear negative DS disease (prior Rx) */
            double PTpA_to_IspA = FS*a_age_A[i][j][l]*(1-p_A[j][l])*sig_H*PTp_A[i][j][l];          /* Post PT to smear positive DS disease (prior Rx) */
            double PTpA_to_LmpA = FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g*PTp_A[i][j][l];          /* Post PT to latent DR (prior Rx) */
            double PTpA_to_NmpA = FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H)*PTp_A[i][j][l];      /* Post PT to smear negative DR disease (prior Rx) */
            double PTpA_to_ImpA = FM*a_age_A[i][j][l]*(1-p_A[j][l])*sig_H*PTp_A[i][j][l];          /* Post PT to smear positive DR disease (prior Rx) */
        
            dS_A[i][j][l] = - m_b[i]*S_A[i][j][l] - (FS + FM)*S_A[i][j][l] +
                            ART_prop[i][j]*A_start[l]*S_H[i][j] + A_prog[l]*S_A[i][j][l-1] - A_prog[l+1]*S_A[i][j][l] - up_A_mort[l][j][i]*S_A[i][j][l] + 
                            (S_A[i][j][l]/tot_age[i])*(forc[i+52]);

            dLsn_A[i][j][l] = - m_b[i]*Lsn_A[i][j][l] +
                            FS*((1-a_age_A[i][j][l])*S_A[i][j][l] + (1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g)*Lmn_A[i][j][l]) -
                            (v_age_A[i][j][l] + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g)*Lsn_A[i][j][l] +
                            r_H*(Isn_A[i][j][l] + Nsn_A[i][j][l]) +
                            ART_prop[i][j]*A_start[l]*Lsn_H[i][j] + A_prog[l]*Lsn_A[i][j][l-1] - A_prog[l+1]*Lsn_A[i][j][l] - up_A_mort[l][j][i]*Lsn_A[i][j][l] +
                            (Lsn_A[i][j][l]/tot_age[i])*(forc[i+52]) -
                            health*kpos*(1-sp_I_pos*sp_N_pos)*l_s*tART_s*Lsn_A[i][j][l] + PTnA_to_LsnA;

            dLsp_A[i][j][l] = - m_b[i]*Lsp_A[i][j][l] +
                            FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g)*Lmp_A[i][j][l] -
                            (v_age_A[i][j][l] + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g)*Lsp_A[i][j][l] +
                            r_H*(Isp_A[i][j][l]+Nsp_A[i][j][l]) +
                            kpos*(l_s*(1-e)*tART_s*((dstpos_p*sp_m_pos)+(1-dstpos_p)) + l_m*tART_m*dstpos_p*(1-sp_m_pos))*(se_I_pos*Isp_A[i][j][l] + se_N_pos*rel_d*Nsp_A[i][j][l]) + 
                            kpos*(l_s*(1-e)*tART_s*((dstpos_n*sp_m_pos)+(1-dstpos_n)) + l_m*tART_m*dstpos_n*(1-sp_m_pos))*(se_I_pos*Isn_A[i][j][l] + se_N_pos*rel_d*Nsn_A[i][j][l]) +
                            ART_prop[i][j]*A_start[l]*Lsp_H[i][j] + A_prog[l]*Lsp_A[i][j][l-1] - A_prog[l+1]*Lsp_A[i][j][l] - up_A_mort[l][j][i]*Lsp_A[i][j][l] +
                            (Lsp_A[i][j][l]/tot_age[i])*(forc[i+52]) -
                            health*kpos*(1-sp_I_pos*sp_N_pos)*l_s*tART_s*Lsp_A[i][j][l] + PTpA_to_LspA;
                             
            dLmn_A[i][j][l] = - m_b[i]*Lmn_A[i][j][l] +
                            FM*((1-a_age_A[i][j][l])*S_A[i][j][l] + (1-a_age_A[i][j][l])*(1-p_A[j][l])*g*Lsn_A[i][j][l]) -
                            (v_age_A[i][j][l] + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g))*Lmn_A[i][j][l] +
                            r_H*(Imn_A[i][j][l]+Nmn_A[i][j][l]) +
                            ART_prop[i][j]*A_start[l]*Lmn_H[i][j] + A_prog[l]*Lmn_A[i][j][l-1] - A_prog[l+1]*Lmn_A[i][j][l] - up_A_mort[l][j][i]*Lmn_A[i][j][l] +
                            (Lmn_A[i][j][l]/tot_age[i])*(forc[i+52]) + PTnA_to_LmnA;
      
            dLmp_A[i][j][l] = - m_b[i]*Lmp_A[i][j][l] + 
                            FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g*Lsp_A[i][j][l] -
                            (v_age_A[i][j][l] + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g))*Lmp_A[i][j][l] +
                            r_H*(Imp_A[i][j][l]+Nmp_A[i][j][l]) +
                            kpos*(dstpos_p*se_m_pos*l_m*tART_m + l_s*tART_s*eff_p*((1-dstpos_p)+dstpos_p*(1-se_m_pos)))*(se_I_pos*Imp_A[i][j][l] + se_N_pos*rel_d*Nmp_A[i][j][l]) +
                            kpos*(dstpos_n*se_m_pos*l_m*tART_m + l_s*tART_s*eff_n*((1-dstpos_n)+dstpos_n*(1-se_m_pos)))*(se_I_pos*Imn_A[i][j][l] + se_N_pos*rel_d*Nmn_A[i][j][l]) + 
                            ART_prop[i][j]*A_start[l]*Lmp_H[i][j] + A_prog[l]*Lmp_A[i][j][l-1] - A_prog[l+1]*Lmp_A[i][j][l] - up_A_mort[l][j][i]*Lmp_A[i][j][l] +
                            (Lmp_A[i][j][l]/tot_age[i])*(forc[i+52]) + PTpA_to_LmpA;

            dNsn_A[i][j][l] = - m_b[i]*Nsn_A[i][j][l]+
                            (v_age_A[i][j][l]*(1-sig_H) + FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lsn_A[i][j][l] +
                            FS*a_age_A[i][j][l]*(1-sig_H)*(S_A[i][j][l] + (1-p_A[j][l])*Lmn_A[i][j][l]) -
                            (theta_H + r_H + muN_H_A[i][l])*Nsn_A[i][j][l] -
                            kpos*se_N_pos*rel_d*(l_s*(dstpos_n*sp_m_pos + (1-dstpos_n)) + dstpos_n*(1-sp_m_pos)*l_m)*Nsn_A[i][j][l] +
                            ART_prop[i][j]*A_start[l]*Nsn_H[i][j] + A_prog[l]*Nsn_A[i][j][l-1] - A_prog[l+1]*Nsn_A[i][j][l] - up_A_mort[l][j][i]*Nsn_A[i][j][l] +
                            (Nsn_A[i][j][l]/tot_age[i])*(forc[i+52]) + PTnA_to_NsnA;
                                    
            dNsp_A[i][j][l] = - m_b[i]*Nsp_A[i][j][l] + 
                            (v_age_A[i][j][l]*(1-sig_H) + FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lsp_A[i][j][l] +
                            FS*a_age_A[i][j][l]*(1-sig_H)*(1-p_A[j][l])*Lmp_A[i][j][l] -
                            (theta_H + r_H + muN_H_A[i][l])*Nsp_A[i][j][l] -
                            kpos*rel_d*se_N_pos*((1-e)*tART_s*l_s*(sp_m_pos*dstpos_p + (1-dstpos_p)) + dstpos_p*(1-sp_m_pos)*l_m*tART_m + l_s*e*(sp_m_pos*dstpos_p +(1-dstpos_p)))*Nsp_A[i][j][l] + 
                            kpos*rel_d*se_N_pos*(l_s*(1-e)*(1-tART_s)*(dstpos_n*sp_m_pos + (1-dstpos_n)) + dstpos_n*(1-sp_m_pos)*(1-tART_m)*l_m)*Nsn_A[i][j][l] +       
                            ART_prop[i][j]*A_start[l]*Nsp_H[i][j] + A_prog[l]*Nsp_A[i][j][l-1] - A_prog[l+1]*Nsp_A[i][j][l] - up_A_mort[l][j][i]*Nsp_A[i][j][l] +
                            (Nsp_A[i][j][l]/tot_age[i])*(forc[i+52]) + PTpA_to_NspA;
          
            dNmn_A[i][j][l] = - m_b[i]*Nmn_A[i][j][l] +
                            (v_age_A[i][j][l]*(1-sig_H) + FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lmn_A[i][j][l] +
                            FM*a_age_A[i][j][l]*(1-sig_H)*(S_A[i][j][l] + (1-p_A[j][l])*Lsn_A[i][j][l]) -
                            (theta_H + r_H + muN_H_A[i][l])*Nmn_A[i][j][l] -
                            kpos*rel_d*se_N_pos*(l_m*dstpos_n*se_m_pos + l_s*((1-dstpos_n)+dstpos_n*(1-se_m_pos)))*Nmn_A[i][j][l] +
                            ART_prop[i][j]*A_start[l]*Nmn_H[i][j] + A_prog[l]*Nmn_A[i][j][l-1] - A_prog[l+1]*Nmn_A[i][j][l] - up_A_mort[l][j][i]*Nmn_A[i][j][l] +
                            (Nmn_A[i][j][l]/tot_age[i])*(forc[i+52]) + PTnA_to_NmnA;
          
            dNmp_A[i][j][l] = - m_b[i]*Nmp_A[i][j][l] +
                            (v_age_A[i][j][l]*(1-sig_H) + FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lmp_A[i][j][l] +
                            FM*a_age_A[i][j][l]*(1-sig_H)*(1-p_A[j][l])*Lsp_A[i][j][l] -
                            (theta_H + r_H + muN_H_A[i][l])*Nmp_A[i][j][l] -
                            kpos*rel_d*se_N_pos*(dstpos_p*se_m_pos*l_m*tART_m + l_s*tART_s*eff_p*((1-dstpos_p)+dstpos_p*(1-se_m_pos)))*Nmp_A[i][j][l] + 
                            kpos*l_s*rel_d*e*se_N_pos*((dstpos_n*sp_m_pos+(1-dstpos_n))*Nsn_A[i][j][l]+(dstpos_p*sp_m_pos+(1-dstpos_p))*Nsp_A[i][j][l]) + 
                            kpos*se_N_pos*rel_d*(dstpos_n*l_m*se_m_pos*(1-tART_m) + l_s*(1-(tART_s*eff_n))*((1-dstpos_n)+dstpos_n*(1-se_m_pos)))*Nmn_A[i][j][l] +
                            ART_prop[i][j]*A_start[l]*Nmp_H[i][j] + A_prog[l]*Nmp_A[i][j][l-1] - A_prog[l+1]*Nmp_A[i][j][l] - up_A_mort[l][j][i]*Nmp_A[i][j][l] +
                            (Nmp_A[i][j][l]/tot_age[i])*(forc[i+52]) + PTpA_to_NmpA;

            dIsn_A[i][j][l] = - m_b[i]*Isn_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_H + FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lsn_A[i][j][l] +
                            FS*a_age_A[i][j][l]*sig_H*S_A[i][j][l] +
                            FS*a_age_A[i][j][l]*(1-p_A[j][l])*sig_H*Lmn_A[i][j][l] +
                            theta_H*Nsn_A[i][j][l] -
                            (r_H + muI_H_A[i][l])*Isn_A[i][j][l] -
                            kpos*se_I_pos*((dstpos_n*sp_m_pos + (1-dstpos_n))*l_s +(dstpos_n*(1-sp_m_pos)*l_m))*Isn_A[i][j][l] +
                            ART_prop[i][j]*A_start[l]*Isn_H[i][j] + A_prog[l]*Isn_A[i][j][l-1] - A_prog[l+1]*Isn_A[i][j][l] - up_A_mort[l][j][i]*Isn_A[i][j][l] +
                            (Isn_A[i][j][l]/tot_age[i])*(forc[i+52]) + PTnA_to_IsnA;

            dIsp_A[i][j][l] = - m_b[i]*Isp_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_H + FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lsp_A[i][j][l] +
                            FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l])*Lmp_A[i][j][l] +
                            theta_H*Nsp_A[i][j][l] +
                            kpos*se_I_pos*(l_s*(1-e)*(1-tART_s)*(dstpos_n*sp_m_pos + (1-dstpos_n)) + (dstpos_n*(1-sp_m_pos)*l_m*(1-tART_m)))*Isn_A[i][j][l] - 
                            (r_H + muI_H_A[i][l])*Isp_A[i][j][l] -
                            kpos*se_I_pos*(l_s*(1-e)*tART_s*(sp_m_pos*dstpos_p + (1-dstpos_p)) + dstpos_p*(1-sp_m_pos)*l_m*tART_m + l_s*e*(sp_m_pos*dstpos_p + (1-dstpos_p)))*Isp_A[i][j][l] +
                            ART_prop[i][j]*A_start[l]*Isp_H[i][j] + A_prog[l]*Isp_A[i][j][l-1] - A_prog[l+1]*Isp_A[i][j][l] - up_A_mort[l][j][i]*Isp_A[i][j][l] +
                            (Isp_A[i][j][l]/tot_age[i])*(forc[i+52]) + PTpA_to_IspA;
                       
            dImn_A[i][j][l] = - m_b[i]*Imn_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_H + FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lmn_A[i][j][l] +
                            FM*a_age_A[i][j][l]*sig_H*S_A[i][j][l] +
                            FM*a_age_A[i][j][l]*(1-p_A[j][l])*sig_H*Lsn_A[i][j][l] +
                            theta_H*Nmn_A[i][j][l] -
                            kpos*se_I_pos*(l_m*dstpos_n*se_m_pos + (1-dstpos_n)*l_s + dstpos_n*(1-se_m_pos)*l_s)*Imn_A[i][j][l] -
                            (r_H + muI_H_A[i][l])*Imn_A[i][j][l] +
                            ART_prop[i][j]*A_start[l]*Imn_H[i][j] + A_prog[l]*Imn_A[i][j][l-1] - A_prog[l+1]*Imn_A[i][j][l] - up_A_mort[l][j][i]*Imn_A[i][j][l] +
                            (Imn_A[i][j][l]/tot_age[i])*(forc[i+52]) + PTnA_to_ImnA;

            dImp_A[i][j][l] = - m_b[i]*Imp_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_H + FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lmp_A[i][j][l] +
                            FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l])*Lsp_A[i][j][l] +
                            theta_H*Nmp_A[i][j][l] +
                            kpos*se_I_pos*(se_m_pos*l_m*dstpos_n*(1-tART_m) + l_s*(1-(tART_s*eff_n))*((1-dstpos_n)+dstpos_n*(1-se_m_pos)))*Imn_A[i][j][l] + 
                            kpos*se_I_pos*l_s*e*((dstpos_n*sp_m_pos + (1-dstpos_n))*Isn_A[i][j][l]+(dstpos_p*sp_m_pos+ (1-dstpos_p))*Isp_A[i][j][l]) -
                            (r_H + muI_H_A[i][l])*Imp_A[i][j][l] -                     
                            kpos*se_I_pos*(l_m*dstpos_p*tART_m*se_m_pos + l_s*tART_s*eff_p*((1-dstpos_p)+dstpos_p*(1-se_m_pos)))*Imp_A[i][j][l] +
                            ART_prop[i][j]*A_start[l]*Imp_H[i][j] + A_prog[l]*Imp_A[i][j][l-1] - A_prog[l+1]*Imp_A[i][j][l] - up_A_mort[l][j][i]*Imp_A[i][j][l] +
                            (Imp_A[i][j][l]/tot_age[i])*(forc[i+52]) + PTpA_to_ImpA;             

            dPTn_A[i][j][l] = -m_b[i]*PTn_A[i][j][j] - PTnA_to_LsnA - PTnA_to_NsnA - PTnA_to_IsnA - PTnA_to_LmnA - PTnA_to_NmnA - PTnA_to_ImnA +
                              health*kpos*(1-sp_I_pos*sp_N_pos)*l_s*tART_s*Lsn_A[i][j][l] +
                              ART_prop[i][j]*A_start[l]*PTn_H[i][j] + A_prog[l]*PTn_A[i][j][l-1] - A_prog[l+1]*PTn_A[i][j][l] - up_A_mort[l][j][i]*PTn_A[i][j][l] +
                              (PTn_A[i][j][l]/tot_age[i])*(forc[i+52]);  

            dPTp_A[i][j][l] = - m_b[i]*PTp_A[i][j][l] - PTpA_to_LspA - PTpA_to_NspA - PTpA_to_IspA - PTpA_to_LmpA - PTpA_to_NmpA - PTpA_to_ImpA +
                              health*kpos*(1-sp_I_pos*sp_N_pos)*l_s*tART_s*Lsp_A[i][j][l] +   
                              ART_prop[i][j]*A_start[l]*PTp_H[i][j] + A_prog[l]*PTp_A[i][j][l-1] - A_prog[l+1]*PTp_A[i][j][l] - up_A_mort[l][j][i]*PTp_A[i][j][l] +
                              (PTp_A[i][j][l]/tot_age[i])*(forc[i+52]);  
        
            TB_cases_ART_age[i][j][l] = (v_age_A[i][j][l]*(1-sig_H) + FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lsn_A[i][j][l] + FS*a_age_A[i][j][l]*(1-sig_H)*(S_A[i][j][l] + (1-p_A[j][l])*(Lmn_A[i][j][l] + PTn_A[i][j][l])) +
                                      (v_age_A[i][j][l]*(1-sig_H) + FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lsp_A[i][j][l] + FS*a_age_A[i][j][l]*(1-sig_H)*(1-p_A[j][l])*(Lmp_A[i][j][l] + PTp_A[i][j][l]) +
                                      (v_age_A[i][j][l]*(1-sig_H) + FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lmn_A[i][j][l] + FM*a_age_A[i][j][l]*(1-sig_H)*(S_A[i][j][l] + (1-p_A[j][l])*(Lsn_A[i][j][l] + PTn_A[i][j][l])) +
                                      (v_age_A[i][j][l]*(1-sig_H) + FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lmp_A[i][j][l] + FM*a_age_A[i][j][l]*(1-sig_H)*(1-p_A[j][l])*(Lsp_A[i][j][l] + PTp_A[i][j][l]) +
                                      (v_age_A[i][j][l]*sig_H + FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lsn_A[i][j][l] + FS*a_age_A[i][j][l]*sig_H*(S_A[i][j][l] + (1-p_A[j][l])*(Lmn_A[i][j][l] + PTn_A[i][j][l])) +
                                      (v_age_A[i][j][l]*sig_H + FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lsp_A[i][j][l] + FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l])*(Lmp_A[i][j][l] + PTp_A[i][j][l]) +
                                      (v_age_A[i][j][l]*sig_H + FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lmn_A[i][j][l] + FM*a_age_A[i][j][l]*sig_H*(S_A[i][j][l] + (1-p_A[j][l])*(Lsn_A[i][j][l] + PTn_A[i][j][l])) +
                                      (v_age_A[i][j][l]*sig_H + FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lmp_A[i][j][l] + FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l])*(Lsp_A[i][j][l] + PTp_A[i][j][l]);

            TB_cases_ART = TB_cases_ART + TB_cases_ART_age[i][j][l];
        
            TB_cases_age[i] = TB_cases_age[i] + TB_cases_ART_age[i][j][l];
        
          }  /* end loop on ART */

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
      
          for (l=0; l<n_ART; l++){
            
            DS_correct = DS_correct + rel_d*kpos*se_N_pos*(dstpos_n*sp_m_pos + (1-dstpos_n))*l_s*Nsn_A[i][j][l] +
                                      rel_d*kpos*se_N_pos*(dstpos_p*sp_m_pos + (1-dstpos_p))*l_s*Nsp_A[i][j][l] +                                         
                                      kpos*se_I_pos*(dstpos_n*sp_m_pos + (1-dstpos_n))*l_s*Isn_A[i][j][l] +       
                                      kpos*se_I_pos*(dstpos_p*sp_m_pos + (1-dstpos_p))*l_s*Isp_A[i][j][l];
                 
                 
            DS_incorrect = DS_incorrect + rel_d*kpos*se_N_pos*(dstpos_n*(1-se_m_pos) + (1-dstpos_n))*l_s*Nmn_A[i][j][l] + 
                                          rel_d*kpos*se_N_pos*(dstpos_p*(1-se_m_pos) + (1-dstpos_p))*l_s*Nmp_A[i][j][l] +                                        
                                          kpos*se_I_pos*(dstpos_n*(1-se_m_pos) + (1-dstpos_n))*l_s*Imn_A[i][j][l] +       
                                          kpos*se_I_pos*(dstpos_p*(1-se_m_pos) + (1-dstpos_p))*l_s*Imp_A[i][j][l];  
                 
            MDR_correct = MDR_correct + rel_d*kpos*se_N_pos*l_m*(dstpos_n*Nmn_A[i][j][l] + dstpos_p*Nmp_A[i][j][l]) + 
                                        kpos*se_I_pos*l_m*(dstpos_n*Imn_A[i][j][l] + dstpos_p*Imp_A[i][j][l]);  
                                      
            MDR_incorrect = MDR_incorrect + rel_d*kpos*se_N_pos*(1-sp_m_pos)*l_m*(dstpos_n*Nsn_A[i][j][l] + dstpos_p*Nsp_A[i][j][l]) + 
                                            kpos*se_I_pos*(1-sp_m_pos)*l_m*(dstpos_n*Isn_A[i][j][l] + dstpos_p*Isp_A[i][j][l]);  
                 
            FP = FP + health*kpos*(1-sp_I_pos*sp_N_pos)*l_s*(S_A[i][j][l]+Lsn_A[i][j][l]+Lsp_A[i][j][l]+Lmn_A[i][j][l]+Lmp_A[i][j][l]+PTn_A[i][j][l]+PTp_A[i][j][l]);  
            
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
    yout[5] = Total_Ns_N + Total_Ns_H;
    yout[6] = Total_Nm_N + Total_Nm_H;
    yout[7] = Total_N;
    yout[8] = Total_Is_N + Total_Is_H;
    yout[9] = Total_Im_N + Total_Im_H;
    yout[10] = Total_I;
    yout[11] = Total_DS;
    yout[12] = Total_MDR;
    yout[13] = FS*100;
    yout[14] = FM*100;
    yout[15] = CD4_dist_all[0];
    yout[16] = CD4_dist_all[1];
    yout[17] = CD4_dist_all[2];
    yout[18] = CD4_dist_all[3];
    yout[19] = CD4_dist_all[4];
    yout[20] = CD4_dist_all[5];
    yout[21] = CD4_dist_all[6];
    yout[22] = CD4_dist_ART_all[0];
    yout[23] = CD4_dist_ART_all[1];
    yout[24] = CD4_dist_ART_all[2];
    yout[25] = CD4_dist_ART_all[3];
    yout[26] = CD4_dist_ART_all[4];
    yout[27] = CD4_dist_ART_all[5];
    yout[28] = CD4_dist_ART_all[6];
    yout[29] = TB_deaths_tot;
    yout[30] = TB_cases_neg;
    yout[31] = TB_cases_pos;
    yout[32] = TB_cases_ART;
    yout[33] = births;
    yout[34] = Tot_deaths;
    yout[35] = DS_correct;
    yout[36] = DS_incorrect;
    yout[37] = MDR_correct;
    yout[38] = MDR_incorrect;
    yout[39] = FP;
}



