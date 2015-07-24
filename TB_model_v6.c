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
static double forc[176];

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
#define k forc[173]       /* detection rate */
#define dst_n forc[174]   /* dst coverage in new cases */
#define dst_p forc[175]   /* dst coverage in previoulsy treated cases */

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
    int N=176;
    odeforcs(&N, forc);
}

/* ##### EVENTS ARE USED TO ADD BIRTHS AND SHIFT THE POPULATION BY AGE - EQUIVALENT TO THE METHOD OF SCHENZLE ###### */

/* ~~~~~~ NEED TO UPDATE THIS TO DEAL WITH ALL DISEASE STATES ~~~~~~~~~~ */
void event(int *n, double *t, double *y) 
{
  int i;
  
  /* Store current population in temp */
  double temp[30537];
  for (i=0; i<30537; i++) temp[i] = y[i];
  
  /* Shift every age group forward one - currently just doing for susceptibles */
  for (i=1; i<30537; i++) y[i] = temp[i-1];
  
  /* Set every 0 age group to zero */
  for (i=0; i<30537; i+=81) y[i] = 0;
  /* Then add births into group 0 - only susceptibles get born */ 
  y[0] = birth_rate*sumsum(temp,0,30536)/1000;

  /* Set >80 age groups to the previous age group plus those still surviving */
  for (i=80; i<30537; i+=81) y[i] = temp[i] + temp[i-1]; /* the final age group is just the previous age plus those still left in */

}

/* ###### DERIVATIVE FUNCTION ###### */

void derivsc(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if (ip[0] <2) error("nout should be at least 2");
    
    /* Expand state variables so can use more meaningful names than y and ydot (the variables and rates of change are vectors y and ydot here) */
    /* There are 81 age groups, 7 HIV positive categories and 3 times on ART */
    /* _H = HIV+; _A = HIV+ on ART */
    
    /* These are the variables */
    double S[81];   double S_H[81][7];   double S_A[81][7][3];      /* Susceptible */
    double Lsn[81]; double Lsn_H[81][7]; double Lsn_A[81][7][3];    /* Latent, DS, new */
    double Lsp[81]; double Lsp_H[81][7]; double Lsp_A[81][7][3];    /* Latent, DS, previous */
    double Lmn[81]; double Lmn_H[81][7]; double Lmn_A[81][7][3];    /* Latent, DR, new */
    double Lmp[81]; double Lmp_H[81][7]; double Lmp_A[81][7][3];    /* Latent, DR, previous */
    double Nsn[81]; double Nsn_H[81][7]; double Nsn_A[81][7][3];    /* Smear negative, DS, new */
    double Nsp[81]; double Nsp_H[81][7]; double Nsp_A[81][7][3];    /* Smear negative, DS, previous */
    double Nmn[81]; double Nmn_H[81][7]; double Nmn_A[81][7][3];    /* Smear negative, DR, new */
    double Nmp[81]; double Nmp_H[81][7]; double Nmp_A[81][7][3];    /* Smear negative, DR, previous */
    double Isn[81]; double Isn_H[81][7]; double Isn_A[81][7][3];    /* Smear positive, DS, new */
    double Isp[81]; double Isp_H[81][7]; double Isp_A[81][7][3];    /* Smear positive, DS, previous */
    double Imn[81]; double Imn_H[81][7]; double Imn_A[81][7][3];    /* Smear positive, DR, new */
    double Imp[81]; double Imp_H[81][7]; double Imp_A[81][7][3];    /* Smear positive, DR, previous */

    /* These are the rates of change (same names but prefixed with d) */
    double dS[81];   double dS_H[81][7];   double dS_A[81][7][3];
    double dLsn[81]; double dLsn_H[81][7]; double dLsn_A[81][7][3];
    double dLsp[81]; double dLsp_H[81][7]; double dLsp_A[81][7][3];
    double dLmn[81]; double dLmn_H[81][7]; double dLmn_A[81][7][3];
    double dLmp[81]; double dLmp_H[81][7]; double dLmp_A[81][7][3];
    double dNsn[81]; double dNsn_H[81][7]; double dNsn_A[81][7][3];
    double dNsp[81]; double dNsp_H[81][7]; double dNsp_A[81][7][3];
    double dNmn[81]; double dNmn_H[81][7]; double dNmn_A[81][7][3];
    double dNmp[81]; double dNmp_H[81][7]; double dNmp_A[81][7][3];
    double dIsn[81]; double dIsn_H[81][7]; double dIsn_A[81][7][3];
    double dIsp[81]; double dIsp_H[81][7]; double dIsp_A[81][7][3];
    double dImn[81]; double dImn_H[81][7]; double dImn_A[81][7][3];
    double dImp[81]; double dImp_H[81][7]; double dImp_A[81][7][3];

    /* intergers to use as counters */ 
    int i,j,l,ij;
     
    /* Then need to assign the variables to the correct bit of y (we do the same for the rates of change after solving the equations) */
    /* SEE IF CAN FIND A NEATER WAY TO DO THIS  - MAYBE JUST HAVE IT ALL IN A LOOP WITH A COUNTER TO KEEP TRACK OF WHERE WE ARE ??*/
     
    /* HIV- */ 
    
    int n_age = 81;     /* Number of age groups */
    int n_HIV = 7;      /* Number of HIV pos groups */
    int n_ART = 3;      /* Number of ART groups */
    int n_disease = 13; /* Number of disease states */
    
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
          ij = ij+1;
        }
      }
    }
    
    /* Adjust TB model parameters for age, HIV and ART */

    /* Create vectors of disease parameters (by age - now single year bins) to use in derivatives - includes BCG effect */
    double a_age[81];
    double sig_age[81];
    double v_age[81];
    double muN_age[81];
    double muI_age[81];
    double bcg = (BCG_cov*(1-BCG_eff)+(1-BCG_cov));  /* those on bcg (BCG_COV) have RR of (1-BCG_eff) */

    for (i=0; i<5; i++){
      a_age[i] = a0*bcg;
      sig_age[i] = sig0;
      v_age[i] = v*bcg;
      muN_age[i] = mu_N0;
      muI_age[i] = mu_I0;
    }
    for (i=5; i<10; i++){
      a_age[i] = a5*bcg;
      sig_age[i] = sig5;
      v_age[i] = v*bcg;
      muN_age[i] = mu_N;
      muI_age[i] = mu_I;
    }
    for (i=10; i<15; i++){
      a_age[i] = a10*bcg;
      sig_age[i] = sig10;
      v_age[i] = v*bcg;
      muN_age[i] = mu_N;
      muI_age[i] = mu_I;
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
    for (j=0; j<n_HIV; j++){
      p_H[j] = p*RR1p*pow(RR2p,-1*(500-mid_CD4[j])/100);
      for (i=0; i<n_age; i++){
        a_age_H[i][j] = a_age[i]*RR1a*pow(RR2a,(500-mid_CD4[j])/100);
        v_age_H[i][j] = v_age[i]*RR1v*pow(RR2v,(500-mid_CD4[j])/100);
        for (l=0; l<n_ART; l++){
          a_age_A[i][j][l] = fmax(a_age_H[i][j]*ART_TB[l],a_age[i]);   /* fmax or fmin ensures being on ART can't be better than being HIV- */
          v_age_A[i][j][l] = fmax(v_age_H[i][j]*ART_TB[l],v_age[i]);
          p_A[j][l] = fmin(p_H[j]/ART_TB[l],p);
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
      A_mort[1][0][i] = 0.009; A_mort[1][1][i] = 0.012;  A_mort[1][2][i] = 0.012; A_mort[1][3][i] = 0.022; A_mort[1][4][i] = 0.022;  A_mort[1][5][i] = 0.022;  A_mort[1][6][i] = 0.022; 
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
    double Total_Ns = sumsum(Nsn,0,80)+sumsum(Nsp,0,80);    /* Total DS smear negative TB */
    double Total_Nm = sumsum(Nmn,0,80)+sumsum(Nmp,0,80);    /* Total DR smear negative TB */
    double Total_Is = sumsum(Isn,0,80)+sumsum(Isp,0,80);    /* Total DS smear positive TB */
    double Total_Im = sumsum(Imn,0,80)+sumsum(Imp,0,80);    /* Total DR smear positive TB */
    
    /* Now loop through HIV and ART and add them in */
    for (j=0; j<n_HIV; j++){
      for (i=0; i<n_age; i++){
        Total_S = Total_S + S_H[i][j];
        Total_Ls = Total_Ls + Lsn_H[i][j]+Lsp_H[i][j];
        Total_Lm = Total_Lm + Lmn_H[i][j]+Lmp_H[i][j];
        Total_Ns = Total_Ns + Nsn_H[i][j]+Nsp_H[i][j];
        Total_Nm = Total_Nm + Nmn_H[i][j]+Nmp_H[i][j];
        Total_Is = Total_Is + Isn_H[i][j]+Isp_H[i][j];
        Total_Im = Total_Im + Imn_H[i][j]+Imp_H[i][j];
        for (l=0; l<n_ART; l++){
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
    double TB_deaths[81];
    double TB_ART_deaths = 0;
    double HIV_deaths[81];
    double ART_deaths[81];
    for (i=0; i<n_age; i++){
      HIV_deaths[i] = 0;
      ART_deaths[i] = 0;
    }
    
    for (i=0; i<n_age; i++) {
      TB_deaths[i] = (Nsn[i]+Nsp[i]+Nmn[i]+Nmp[i])*muN_age[i] + (Isn[i]+Isp[i]+Imn[i]+Imp[i])*muI_age[i];
      for(j=0; j<n_HIV; j++){
        TB_deaths[i] = TB_deaths[i] + (Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j])*muN_H + (Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j])*muI_H;
                                      
        HIV_deaths[i] = HIV_deaths[i]+H_mort[j][i]*(S_H[i][j]+Lsn_H[i][j]+Lsp_H[i][j]+Lmn_H[i][j]+Lmp_H[i][j]+
                        Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j]+Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j]);                              
                          
        for (l=0; l<n_ART; l++){
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
    double TB_deaths_tot = sumsum(TB_deaths,0,80);
    double ART_deaths_tot = sumsum(ART_deaths,0,80) + TB_ART_deaths; 
    
    /* Calculate total population by age */
    double tot_age[81];
    for (i=0; i<n_age; i++){
      tot_age[i] = S[i]+Lsn[i]+Lsp[i]+Lmn[i]+Lmp[i]+Nsn[i]+Nsp[i]+Nmn[i]+Nmp[i]+Isn[i]+Isp[i]+Imn[i]+Imp[i];
      for (j=0; j<n_HIV; j++){
        tot_age[i] = tot_age[i]+S_H[i][j]+
                     Lsn_H[i][j]+Lsp_H[i][j]+Lmn_H[i][j]+Lmp_H[i][j]+
                     Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j]+
                     Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j];
        for (l=0; l<n_ART; l++){
          tot_age[i] = tot_age[i]+S_A[i][j][l]+Lsn_A[i][j][l]+Lsp_A[i][j][l]+Lmn_A[i][j][l]+Lmp_A[i][j][l]+
                            Nsn_A[i][j][l]+Nsp_A[i][j][l]+Nmn_A[i][j][l]+Nmp_A[i][j][l]+ 
                            Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l];
        }
      }
    }
    /* Then calculate % of pop that die of disease by age - this is only calculated pre 2015 (determined by value of pop_ad) so that the same constant proportional correction is then applied in the future */
    double prop_dis_death[81];
    if(pop_ad > 0){
      for (i=0; i<n_age; i++){
        prop_dis_death[i] = (TB_deaths[i]+HIV_deaths[i]+ART_deaths[i])/tot_age[i];
      }
    }

    /* Then reduce background mortality by the disease mortality */
    double m_b[81];
    for (i=0; i<n_age; i++){
      m_b[i] = forc[i+1]-prop_dis_death[i];
    }
    
    /* Calculate total deaths */
    double Tot_deaths = 0;
    double Tot_deaths_age[81];
    for (i=0; i<n_age; i++){
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

    for (j=0; j<n_HIV; j++){
      for (i=0; i<n_age; i++){
        CD4_dist[j] = CD4_dist[j]+S_H[i][j]+Lsn_H[i][j]+Lsp_H[i][j]+Lmn_H[i][j]+Lmp_H[i][j]+
                      Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j]+ 
                      Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j];
        CD4_deaths[j] = CD4_deaths[j]+
                        H_mort[j][i]*(S_H[i][j]+Lsn_H[i][j]+Lsp_H[i][j]+Lmn_H[i][j]+Lmp_H[i][j]+
                        Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j]+ 
                        Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j]);
        for (l=0; l<n_ART; l++){
          CD4_dist_ART[j] = CD4_dist_ART[j]+S_A[i][j][l]+Lsn_A[i][j][l]+Lsp_A[i][j][l]+Lmn_A[i][j][l]+Lmp_A[i][j][l]+
                            Nsn_A[i][j][l]+Nsp_A[i][j][l]+Nmn_A[i][j][l]+Nmp_A[i][j][l]+ 
                            Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l];         
                      
        }            
      } 
      Tot_ART = Tot_ART + CD4_dist_ART[j];  /* sum up number currently on ART */
      ART_on = ART_on + (CD4_dist[j] + CD4_dist_ART[j])*forc[j+163]; /* number who should be on ART - HIV+ population times coverage (by CD4) */
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
      for (j=Athresh; j<n_ART; j++) {
        if (CD4_dist[j] > 0) {
          ART_prop[j] = (((CD4_dist[j]/ART_el)+(CD4_deaths[j]/ART_el_deaths))/2)*(ART_new/CD4_dist[j]); /* applies weighting and size of CD4 group to work out % of CD4 group that should move */
           
        }
      }
    }
    
    /* Force of infection */
    double FS = beta*(Total_Ns*rel_inf + Total_Is)/Total; 
    double FM = fit_cost*beta*(Total_Nm*rel_inf + Total_Im)/Total; 
    
    /* Variables to store numbers of new cases */
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
      
          
      /* Susceptible - NOTE BIRTHS ARE ADDED TO HERE IN THE EVENTS FUNCTION*/
      dS[i] = - (FS + FM)*S[i] - forc[i+82]*S[i] - m_b[i]*S[i];
       
      /* Latent, ds, naive */
      dLsn[i] = - m_b[i]*Lsn[i] + S_to_Lsn + Lmn_to_Lsn + r*(Isn[i] + Nsn[i]) - Lsn_to_Lmn - Lsn_to_Nsn - Lsn_to_Isn  - Lsn_to_Nmn - Lsn_to_Imn -
                forc[i+82]*Lsn[i];                   
      /* Latent, ds, prev */
      dLsp[i] = - m_b[i]*Lsp[i] + Lmp_to_Lsp + r*(Isp[i] + Nsp[i]) - Lsp_to_Lmp - Lsp_to_Nsp - Lsp_to_Isp  - Lsp_to_Nmp - Lsp_to_Imp + 
                k*l_s*(1-e)*tau_s*(Isn[i]+Isp[i]) + k*l_s*(1-e)*tau_s*d*(Nsn[i]+Nsp[i]) - /* Care */
                forc[i+82]*Lsp[i]; /* HIV */          
      /* Latent, mdr, naive */ 
      dLmn[i] = - m_b[i]*Lmn[i] + S_to_Lmn + Lsn_to_Lmn + r*(Imn[i] + Nmn[i]) - Lmn_to_Lsn - Lmn_to_Nsn - Lmn_to_Isn - Lmn_to_Nmn - Lmn_to_Imn  -                
                forc[i+82]*Lmn[i]; /* HIV */             
      /* Latent, mdr, prev */
      dLmp[i] = - m_b[i]*Lmp[i] + Lsp_to_Lmp + r*(Imp[i] + Nmp[i]) - Lmp_to_Lsp - Lmp_to_Nsp - Lmp_to_Isp - Lmp_to_Nmp - Lmp_to_Imp +
                k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(Imn[i]+d*Nmn[i]) + k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(Imp[i]+d*Nmp[i]) - /* Care */
                forc[i+82]*Lmp[i]; /* HIV */          
      /* Smear neg, ds, new */
      dNsn[i] = - m_b[i]*Nsn[i] + S_to_Nsn + Lsn_to_Nsn + Lmn_to_Nsn - (theta + r + muN_age[i])*Nsn[i] -
                k*l_s*d*Nsn[i] - /* care */
                forc[i+82]*Nsn[i]; /* HIV */
      /* Smear neg, ds, prev */                             
      dNsp[i] = - m_b[i]*Nsp[i] + Lsp_to_Nsp + Lmp_to_Nsp - (theta + r + muN_age[i])*Nsp[i] -
                (k*l_s*d*(1-e)*tau_s + k*l_s*d*e)*Nsp[i] + k*l_s*d*(1-e)*(1-tau_s)*Nsn[i] - /* Care */
                forc[i+82]*Nsp[i]; /* HIV */
      /* Smear neg, mdr, new */
      dNmn[i] = - m_b[i]*Nmn[i] + S_to_Nmn + Lsn_to_Nmn + Lmn_to_Nmn - (theta + r + muN_age[i])*Nmn[i] -
                k*(l_m*d*dst_n + l_s*d*(1-dst_n))*Nmn[i] - /* Care */
                forc[i+82]*Nmn[i]; /* HIV */        
      /* Smear neg, mdr, prev */
      dNmp[i] = - m_b[i]*Nmp[i] + Lsp_to_Nmp + Lmp_to_Nmp - (theta + r + muN_age[i])*Nmp[i] -
                k*(l_m*d*dst_p*tau_m + l_s*d*(1-dst_p)*tau_s*eff_p)*Nmp[i] + k*l_s*d*e*(Nsn[i]+Nsp[i]) + (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*Nmn[i] - /* Care */
                forc[i+82]*Nmp[i]; /* HIV */
      /* Smear pos, ds, new */
      dIsn[i] = - m_b[i]*Isn[i] + S_to_Isn + Lsn_to_Isn + Lmn_to_Isn + theta*Nsn[i] - (r + muI_age[i])*Isn[i] - 
                k*l_s*Isn[i] - /* Care */
                forc[i+82]*Isn[i]; /* HIV */
      /* Smear pos, ds, prev */
      dIsp[i] = - m_b[i]*Isp[i] + Lsp_to_Isp + Lmp_to_Isp + theta*Nsp[i] - (r + muI_age[i])*Isp[i] - 
                (k*l_s*(1-e)*tau_s + k*l_s*e)*Isp[i] + k*l_s*(1-e)*(1-tau_s)*Isn[i] - /* Care */
                forc[i+82]*Isp[i]; /* HIV */
      /* Smear pos, mdr, new */
      dImn[i] = - m_b[i]*Imn[i] + S_to_Imn + Lsn_to_Imn + Lmn_to_Imn + theta*Nmn[i] - (r + muI_age[i])*Imn[i] -
                (k*l_m*dst_n + k*l_s*(1-dst_n))*Imn[i] - /* Care */ 
                forc[i+82]*Imn[i]; /* HIV */
      /* Smear pos, mdr, prev */
      dImp[i] = - m_b[i]*Imp[i] + Lsp_to_Imp + Lmp_to_Imp + theta*Nmp[i] - (r + muI_age[i])*Imp[i] -
                k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*Imp[i] + (k*l_m*dst_n*(1-tau_m) + k*l_s*(1-dst_n)*(1-(tau_s*eff_n)))*Imn[i] + k*l_s*e*(Isn[i]+Isp[i]) -
                forc[i+82]*Imp[i];
                        
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
    
      for (j=0; j<n_HIV; j++){      /* CD4 */

        dS_H[i][j] = - m_b[i]*S_H[i][j] - 
                     (FS + FM)*S_H[i][j] +
                     forc[i+82]*H_CD4[j][i]*S[i] - H_prog[j+1][i]*S_H[i][j] + H_prog[j][i]*S_H[i][j-1] - 
                     H_mort[j][i]*S_H[i][j] - ART_prop[j]*S_H[i][j];

        dLsn_H[i][j] = - m_b[i]*Lsn_H[i][j] +
                       FS*((1-a_age_H[i][j])*S_H[i][j] + (1-a_age_H[i][j])*(1-p_H[j])*(1-g)*Lmn_H[i][j]) -
                       (v_age_H[i][j] + FS*a_age_H[i][j]*(1-p_H[j]) + FM*a_age_H[i][j]*(1-p_H[j]) + FM*(1-a_age_H[i][j])*(1-p_H[j])*g)*Lsn_H[i][j] +
                       r*(Isn_H[i][j] + Nsn_H[i][j]) +
                       forc[i+82]*H_CD4[j][i]*Lsn[i] - H_prog[j+1][i]*Lsn_H[i][j] + H_prog[j][i]*Lsn_H[i][j-1] - 
                       H_mort[j][i]*Lsn_H[i][j] - ART_prop[j]*Lsn_H[i][j];
        
        dLsp_H[i][j] = - m_b[i]*Lsp_H[i][j] +
                       FS*(1-a_age_H[i][j])*(1-p_H[j])*(1-g)*Lmp_H[i][j] -
                       (v_age_H[i][j] + FS*a_age_H[i][j]*(1-p_H[j]) + FM*a_age_H[i][j]*(1-p_H[j]) + FM*(1-a_age_H[i][j])*(1-p_H[j])*g)*Lsp_H[i][j] +
                       r_H*(Isp_H[i][j]+Nsp_H[i][j]) +
                       k*l_s*(1-e)*tau_s*(Isn_H[i][j]+Isp_H[i][j]) + 
                       k*l_s*(1-e)*tau_s*d*(Nsn_H[i][j]+Nsp_H[i][j]) +
                       forc[i+82]*H_CD4[j][i]*Lsp[i] - H_prog[j+1][i]*Lsp_H[i][j] + H_prog[j][i]*Lsp_H[i][j-1] - 
                       H_mort[j][i]*Lsp_H[i][j] - ART_prop[j]*Lsp_H[i][j];
                             
        dLmn_H[i][j] = - m_b[i]*Lmn_H[i][j] +
                       FM*((1-a_age_H[i][j])*S_H[i][j] + (1-a_age_H[i][j])*(1-p_H[j])*g*Lsn_H[i][j]) -
                       (v_age_H[i][j] + FM*a_age_H[i][j]*(1-p_H[j]) + FS*a_age_H[i][j]*(1-p_H[j]) + FS*(1-a_age_H[i][j])*(1-p_H[j])*(1-g))*Lmn_H[i][j] +
                       r_H*(Imn_H[i][j]+Nmn_H[i][j]) +
                       forc[i+82]*H_CD4[j][i]*Lmn[i] - H_prog[j+1][i]*Lmn_H[i][j] + H_prog[j][i]*Lmn_H[i][j-1] - 
                       H_mort[j][i]*Lmn_H[i][j] - ART_prop[j]*Lmn_H[i][j];
                     
        dLmp_H[i][j] = - m_b[i]*Lmp_H[i][j] + 
                       FM*(1-a_age_H[i][j])*(1-p_H[j])*g*Lsp_H[i][j] -
                       (v_age_H[i][j] + FM*a_age_H[i][j]*(1-p_H[j]) + FS*a_age_H[i][j]*(1-p_H[j]) + FS*(1-a_age_H[i][j])*(1-p_H[j])*(1-g))*Lmp_H[i][j] +
                       r_H*(Imp_H[i][j]+Nmp_H[i][j]) +
                       k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(Imn_H[i][j]+d*Nmn_H[i][j]) +
                       k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(Imp_H[i][j]+d*Nmp_H[i][j]) +
                       forc[i+82]*H_CD4[j][i]*Lmp[i] - H_prog[j+1][i]*Lmp_H[i][j] + H_prog[j][i]*Lmp_H[i][j-1] - 
                       H_mort[j][i]*Lmp_H[i][j] - ART_prop[j]*Lmp_H[i][j];
           
        dNsn_H[i][j] = - m_b[i]*Nsn_H[i][j] +
                       (v_age_H[i][j]*(1-sig_H) + FS*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lsn_H[i][j] +
                       FS*a_age_H[i][j]*(1-sig_H)*(S_H[i][j] + (1-p_H[j])*Lmn_H[i][j]) -
                       (theta + r_H + k*l_s*d + muN_H)*Nsn_H[i][j] +
                       forc[i+82]*H_CD4[j][i]*Nsn[i] - H_prog[j+1][i]*Nsn_H[i][j] + H_prog[j][i]*Nsn_H[i][j-1] - 
                       H_mort[j][i]*Nsn_H[i][j] - ART_prop[j]*Nsn_H[i][j];
                                                
        dNsp_H[i][j] = - m_b[i]*Nsp_H[i][j] + 
                       (v_age_H[i][j]*(1-sig_H) + FS*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lsp_H[i][j] +
                       FS*a_age_H[i][j]*(1-sig_H)*(1-p_H[j])*Lmp_H[i][j] -
                       (theta + r_H + k*l_s*d*(1-e)*tau_s + k*l_s*d*e + muN_H)*Nsp_H[i][j] +
                       k*l_s*d*(1-e)*(1-tau_s)*Nsn_H[i][j] +
                       forc[i+82]*H_CD4[j][i]*Nsp[i] - H_prog[j+1][i]*Nsp_H[i][j] + H_prog[j][i]*Nsp_H[i][j-1] - 
                       H_mort[j][i]*Nsp_H[i][j] - ART_prop[j]*Nsp_H[i][j];
                       
        dNmn_H[i][j] = - m_b[i]*Nmn_H[i][j] +
                       (v_age_H[i][j]*(1-sig_H) + FM*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lmn_H[i][j] +
                       FM*a_age_H[i][j]*(1-sig_H)*(S_H[i][j] + (1-p_H[j])*Lsn_H[i][j]) -
                       (theta + r_H + k*l_m*d*dst_n + k*l_s*d*(1-dst_n) + muN_H)*Nmn_H[i][j] +
                       forc[i+82]*H_CD4[j][i]*Nmn[i] - H_prog[j+1][i]*Nmn_H[i][j] + H_prog[j][i]*Nmn_H[i][j-1] - 
                       H_mort[j][i]*Nmn_H[i][j] - ART_prop[j]*Nmn_H[i][j];
                       
        dNmp_H[i][j] = - m_b[i]*Nmp_H[i][j] +
                       (v_age_H[i][j]*(1-sig_H) + FM*a_age_H[i][j]*(1-p_H[j])*(1-sig_H))*Lmp_H[i][j] +
                       FM*a_age_H[i][j]*(1-sig_H)*(1-p_H[j])*Lsp_H[i][j] -
                       (theta + r_H + k*l_m*d*dst_p*tau_m + k*l_s*d*(1-dst_p)*tau_s*eff_p + muN_H)*Nmp_H[i][j] +
                       k*l_s*d*e*(Nsn_H[i][j]+Nsp_H[i][j]) +
                       (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*Nmn_H[i][j] +
                       forc[i+82]*H_CD4[j][i]*Nmp[i] - H_prog[j+1][i]*Nmp_H[i][j] + H_prog[j][i]*Nmp_H[i][j-1] - 
                       H_mort[j][i]*Nmp_H[i][j] - ART_prop[j]*Nmp_H[i][j];
           
        dIsn_H[i][j] = - m_b[i]*Isn_H[i][j] +
                       (v_age_H[i][j]*sig_H + FS*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lsn_H[i][j] +
                       FS*a_age_H[i][j]*sig_H*S_H[i][j] +
                       FS*a_age_H[i][j]*(1-p_H[j])*sig_H*Lmn_H[i][j] +
                       theta*Nsn_H[i][j] -
                       (r_H + k*l_s + muI_H)*Isn_H[i][j] +
                       forc[i+82]*H_CD4[j][i]*Isn[i] - H_prog[j+1][i]*Isn_H[i][j] + H_prog[j][i]*Isn_H[i][j-1] - 
                       H_mort[j][i]*Isn_H[i][j] - ART_prop[j]*Isn_H[i][j];

        dIsp_H[i][j] = - m_b[i]*Isp_H[i][j] +
                       (v_age_H[i][j]*sig_H + FS*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lsp_H[i][j] +
                       FS*a_age_H[i][j]*sig_H*(1-p_H[j])*Lmp_H[i][j] +
                       theta*Nsp_H[i][j] +
                       k*l_s*(1-e)*(1-tau_s)*Isn_H[i][j] -
                       (r_H + k*l_s*(1-e)*tau_s + k*l_s*e + muI_H)*Isp_H[i][j] +
                       forc[i+82]*H_CD4[j][i]*Isp[i] - H_prog[j+1][i]*Isp_H[i][j] + H_prog[j][i]*Isp_H[i][j-1] - 
                       H_mort[j][i]*Isp_H[i][j] - ART_prop[j]*Isp_H[i][j];
                                   
        dImn_H[i][j] = - m_b[i]*Imn_H[i][j] +
                       (v_age_H[i][j]*sig_H + FM*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lmn_H[i][j] +
                       FM*a_age_H[i][j]*sig_H*S_H[i][j] +
                       FM*a_age_H[i][j]*(1-p_H[j])*sig_H*Lsn_H[i][j] +
                       theta*Nmn_H[i][j] -
                       (r_H + k*l_m*dst_n + k*l_s*(1-dst_n) + muI_H)*Imn_H[i][j] +
                       forc[i+82]*H_CD4[j][i]*Imn[i] - H_prog[j+1][i]*Imn_H[i][j] + H_prog[j][i]*Imn_H[i][j-1] - 
                       H_mort[j][i]*Imn_H[i][j] - ART_prop[j]*Imn_H[i][j];

        dImp_H[i][j] = - m_b[i]*Imp_H[i][j] +
                       (v_age_H[i][j]*sig_H + FM*a_age_H[i][j]*sig_H*(1-p_H[j]))*Lmp_H[i][j] +
                       FM*a_age_H[i][j]*sig_H*(1-p_H[j])*Lsp_H[i][j] +
                       theta*Nmp_H[i][j] +
                       (k*l_m*dst_n*(1-tau_m) + k*l_s*(1-dst_n)*(1-(tau_s*eff_n)))*Imn_H[i][j] +
                       k*l_s*e*(Isn_H[i][j]+Isp_H[i][j]) -
                       (r_H + k*l_m*dst_p*tau_m + k*l_s*(1-dst_p)*tau_s*eff_p + muI_H)*Imp_H[i][j] +
                       forc[i+82]*H_CD4[j][i]*Imp[i] - H_prog[j+1][i]*Imp_H[i][j] + H_prog[j][i]*Imp_H[i][j-1] - 
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
        for (l=0; l<n_ART; l++){
        
          dS_A[i][j][l] = - m_b[i]*S_A[i][j][l] - (FS + FM)*S_A[i][j][l] +
                          ART_prop[j]*A_start[l]*S_H[i][j] + A_prog[l]*S_A[i][j][l-1] - A_prog[l+1]*S_A[i][j][l] - A_mort[l][j][i]*S_A[i][j][l];

          dLsn_A[i][j][l] = - m_b[i]*Lsn_A[i][j][l] +
                            FS*((1-a_age_A[i][j][l])*S_A[i][j][l] + (1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g)*Lmn_A[i][j][l]) -
                            (v_age_A[i][j][l] + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g)*Lsn_A[i][j][l] +
                            r_H*(Isn_A[i][j][l] + Nsn_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lsn_H[i][j] + A_prog[l]*Lsn_A[i][j][l-1] - A_prog[l+1]*Lsn_A[i][j][l] - A_mort[l][j][i]*Lsn_A[i][j][l];
        
          dLsp_A[i][j][l] = - m_b[i]*Lsp_A[i][j][l] +
                            FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g)*Lmp_A[i][j][l] -
                            (v_age_A[i][j][l] + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g)*Lsp_A[i][j][l] +
                            r_H*(Isp_A[i][j][l]+Nsp_A[i][j][l]) +
                            k*l_s*(1-e)*tau_s*(Isn_A[i][j][l]+Isp_A[i][j][l]) + 
                            k*l_s*(1-e)*tau_s*d*(Nsn_A[i][j][l]+Nsp_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lsp_H[i][j] + A_prog[l]*Lsp_A[i][j][l-1] - A_prog[l+1]*Lsp_A[i][j][l] - A_mort[l][j][i]*Lsp_A[i][j][l];
                             
          dLmn_A[i][j][l] = - m_b[i]*Lmn_A[i][j][l] +
                            FM*((1-a_age_A[i][j][l])*S_A[i][j][l] + (1-a_age_A[i][j][l])*(1-p_A[j][l])*g*Lsn_A[i][j][l]) -
                            (v_age_A[i][j][l] + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g))*Lmn_A[i][j][l] +
                            r_H*(Imn_A[i][j][l]+Nmn_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lmn_H[i][j] + A_prog[l]*Lmn_A[i][j][l-1] - A_prog[l+1]*Lmn_A[i][j][l] - A_mort[l][j][i]*Lmn_A[i][j][l];
                     
          dLmp_A[i][j][l] = - m_b[i]*Lmp_A[i][j][l] + 
                            FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g*Lsp_A[i][j][l] -
                            (v_age_A[i][j][l] + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g))*Lmp_A[i][j][l] +
                            r_H*(Imp_A[i][j][l]+Nmp_A[i][j][l]) +
                            k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(Imn_A[i][j][l]+d*Nmn_A[i][j][l]) +
                            k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(Imp_A[i][j][l]+d*Nmp_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lmp_H[i][j] + A_prog[l]*Lmp_A[i][j][l-1] - A_prog[l+1]*Lmp_A[i][j][l] - A_mort[l][j][i]*Lmp_A[i][j][l];
           
          dNsn_A[i][j][l] = - m_b[i]*Nsn_A[i][j][l]+
                            (v_age_A[i][j][l]*(1-sig_H) + FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lsn_A[i][j][l] +
                            FS*a_age_A[i][j][l]*(1-sig_H)*(S_A[i][j][l] + (1-p_A[j][l])*Lmn_A[i][j][l]) -
                            (theta + r_H + k*l_s*d + muN_H*ART_mort[l])*Nsn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Nsn_H[i][j] + A_prog[l]*Nsn_A[i][j][l-1] - A_prog[l+1]*Nsn_A[i][j][l] - A_mort[l][j][i]*Nsn_A[i][j][l];
                                                
          dNsp_A[i][j][l] = - m_b[i]*Nsp_A[i][j][l] + 
                            (v_age_A[i][j][l]*(1-sig_H) + FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lsp_A[i][j][l] +
                            FS*a_age_A[i][j][l]*(1-sig_H)*(1-p_A[j][l])*Lmp_A[i][j][l] -
                            (theta + r_H + k*l_s*d*(1-e)*tau_s + k*l_s*d*e + muN_H*ART_mort[l])*Nsp_A[i][j][l] +
                            k*l_s*d*(1-e)*(1-tau_s)*Nsn_A[i][j][l] +      
                            ART_prop[j]*A_start[l]*Nsp_H[i][j] + A_prog[l]*Nsp_A[i][j][l-1] - A_prog[l+1]*Nsp_A[i][j][l] - A_mort[l][j][i]*Nsp_A[i][j][l];
                       
          dNmn_A[i][j][l] = - m_b[i]*Nmn_A[i][j][l] +
                            (v_age_A[i][j][l]*(1-sig_H) + FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lmn_A[i][j][l] +
                            FM*a_age_A[i][j][l]*(1-sig_H)*(S_A[i][j][l] + (1-p_A[j][l])*Lsn_A[i][j][l]) -
                            (theta + r_H + k*l_m*d*dst_n + k*l_s*d*(1-dst_n) + muN_H*ART_mort[l])*Nmn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Nmn_H[i][j] + A_prog[l]*Nmn_A[i][j][l-1] - A_prog[l+1]*Nmn_A[i][j][l] - A_mort[l][j][i]*Nmn_A[i][j][l];
                       
          dNmp_A[i][j][l] = - m_b[i]*Nmp_A[i][j][l] +
                            (v_age_A[i][j][l]*(1-sig_H) + FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_H))*Lmp_A[i][j][l] +
                            FM*a_age_A[i][j][l]*(1-sig_H)*(1-p_A[j][l])*Lsp_A[i][j][l] -
                            (theta + r_H + k*l_m*d*dst_p*tau_m + k*l_s*d*(1-dst_p)*tau_s*eff_p + muN_H*ART_mort[l])*Nmp_A[i][j][l] +
                            k*l_s*d*e*(Nsn_A[i][j][l]+Nsp_A[i][j][l]) +
                            (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*Nmn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Nmp_H[i][j] + A_prog[l]*Nmp_A[i][j][l-1] - A_prog[l+1]*Nmp_A[i][j][l] - A_mort[l][j][i]*Nmp_A[i][j][l];
           
          dIsn_A[i][j][l] = - m_b[i]*Isn_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_H + FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lsn_A[i][j][l] +
                            FS*a_age_A[i][j][l]*sig_H*S_A[i][j][l] +
                            FS*a_age_A[i][j][l]*(1-p_A[j][l])*sig_H*Lmn_A[i][j][l] +
                            theta*Nsn_A[i][j][l] -
                            (r_H + k*l_s + muI_H*ART_mort[l])*Isn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Isn_H[i][j] + A_prog[l]*Isn_A[i][j][l-1] - A_prog[l+1]*Isn_A[i][j][l] - A_mort[l][j][i]*Isn_A[i][j][l];

          dIsp_A[i][j][l] = - m_b[i]*Isp_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_H + FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lsp_A[i][j][l] +
                            FS*a_age_A[i][j][l]*sig_H*(1-p_A[j][l])*Lmp_A[i][j][l] +
                            theta*Nsp_A[i][j][l] +
                            k*l_s*(1-e)*(1-tau_s)*Isn_A[i][j][l] -
                            (r_H + k*l_s*(1-e)*tau_s + k*l_s*e + muI_H*ART_mort[l])*Isp_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Isp_H[i][j] + A_prog[l]*Isp_A[i][j][l-1] - A_prog[l+1]*Isp_A[i][j][l] - A_mort[l][j][i]*Isp_A[i][j][l];
                                   
          dImn_A[i][j][l] = - m_b[i]*Imn_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_H + FM*a_age_A[i][j][l]*sig_H*(1-p_A[j][l]))*Lmn_A[i][j][l] +
                            FM*a_age_A[i][j][l]*sig_H*S_A[i][j][l] +
                            FM*a_age_A[i][j][l]*(1-p_A[j][l])*sig_H*Lsn_A[i][j][l] +
                            theta*Nmn_A[i][j][l] -
                            (r_H + k*l_m*dst_n + k*l_s*(1-dst_n) + muI_H*ART_mort[l])*Imn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Imn_H[i][j] + A_prog[l]*Imn_A[i][j][l-1] - A_prog[l+1]*Imn_A[i][j][l] - A_mort[l][j][i]*Imn_A[i][j][l];

          dImp_A[i][j][l] = - m_b[i]*Imp_A[i][j][l] +
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
    for (i=0; i<n_age; i++) ydot[i] = dS[i];             
    for (i=n_age; i<n_age*2; i++) ydot[i] = dLsn[i-n_age];      
    for (i=n_age*2; i<n_age*3; i++) ydot[i] = dLsp[i-n_age*2];       /* Lsp: 34-50 */
    for (i=n_age*3; i<n_age*4; i++) ydot[i] = dLmn[i-n_age*3];       /* Lmn: 51-67 */
    for (i=n_age*4; i<n_age*5; i++) ydot[i] = dLmp[i-n_age*4];       /* Lmp: 68-84 */
    for (i=n_age*5; i<n_age*6; i++) ydot[i] = dNsn[i-n_age*5];      /* Nsn: 85-101 */
    for (i=n_age*6; i<n_age*7; i++) ydot[i] = dNsp[i-n_age*6];    /* Nsp: 102-118 */
    for (i=n_age*7; i<n_age*8; i++) ydot[i] = dNmn[i-n_age*7];    /* Nmn: 119-135 */
    for (i=n_age*8; i<n_age*9; i++) ydot[i] = dNmp[i-n_age*8];    /* Nmp: 136-152 */ 
    for (i=n_age*9; i<n_age*10; i++) ydot[i] = dIsn[i-n_age*9];    /* Isn: 153-169 */
    for (i=n_age*10; i<n_age*11; i++) ydot[i] = dIsp[i-n_age*10];    /* Isp: 170-186 */
    for (i=n_age*11; i<n_age*12; i++) ydot[i] = dImn[i-n_age*11];    /* Imn: 187-203 */
    for (i=n_age*12; i<n_age*13; i++) ydot[i] = dImp[i-n_age*12];    /* Imp: 204-220 */
    
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
    yout[44] = births;
    yout[45] = Tot_deaths;
    
}



