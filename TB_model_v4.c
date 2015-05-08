/* TB model in C code to call from R */

/* To add: BCG
           ART effects
           Turning off population scaling
           Incidence and mortality outputs 
           Add in time varying parameters for care and control
*/

/* compile within R with system("R CMD SHLIB TB_model_v3.c") */

#include <R.h>
#include <math.h>

static double parms[47];
static double forc[46];

/* A trick to keep up with the parameters and forcings */
#define age1 parms[0]
#define age2 parms[1]
#define beta parms[2]
#define a_a parms[3]
#define a0 parms[4]
#define a5 parms[5]
#define a10 parms[6]
#define p parms[7]
#define v parms[8]
#define sig_a parms[9]
#define sig0 parms[10]
#define sig5 parms[11]
#define sig10 parms[12]
#define rel_inf parms[13]
#define theta parms[14]
#define r parms[15]
#define mu_N parms[16]
#define mu_N0 parms[17]
#define mu_I parms[18]
#define mu_I0 parms[19]
#define fit_cost parms[20]
#define e parms[21]
#define g parms[22]
#define l_s parms[23]
#define l_m parms[24]
#define d parms[25]
#define tau_s parms[26]
#define tau_m parms[27]
#define eff_n parms[28]
#define eff_p parms[29]
#define dst_n parms[30]
#define dst_p parms[31]
#define muN_H parms[32]
#define muI_H parms[33]
#define RR1a parms[34]
#define RR2a parms[35]
#define RR1v parms[36]
#define RR2v parms[37]
#define RR1p parms[38]
#define RR2p parms[39]
#define ART_TB1 parms[40]
#define ART_TB2 parms[41]
#define ART_TB3 parms[42]
#define ART_mort1 parms[43]
#define ART_mort2 parms[44]
#define ART_mort3 parms[45]
#define BCG_eff parms[46]

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

#define Ahigh forc[35]
#define A500 forc[36]
#define A349 forc[37]
#define A249 forc[38]
#define A199 forc[39]
#define A99 forc[40]
#define A50 forc[41]
#define Athresh forc[42]

#define BCG_cov forc[43]

#define pop_ad forc[44]

#define k forc[45]

/* Function to sum array from element i_start to i_end */
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

/* initializers*/
void parmsc(void (* odeparms)(int *, double *))
{
    int N=47;
    odeparms(&N, parms);
}

void forcc(void (* odeforcs)(int *, double *))
{
    int N=46;
    odeforcs(&N, forc);
}


/* derivative function */
void derivsc(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if (ip[0] <2) error("nout should be at least 2");
    
    /* Expand state variables so can use more meaningful names than y and ydot  */
    double S[17];
    double Lsn[17];
    double Lsp[17];
    double Lmn[17];
    double Lmp[17];
    double Nsn[17];
    double Nsp[17];
    double Nmn[17];
    double Nmp[17];
    double Isn[17];
    double Isp[17];
    double Imn[17];
    double Imp[17];
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
     
    int i;  
    int j;
    int l;
    int ij;
     
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
    
    /* Create vectors of aging rates to use in derivatives */
    double age_out[17];       /* age_out is lower for final age group because it is wider */
    double age_in[17];        /* age_in is set to zero for first age group as new borns only enter the S class */
    for (i=0; i<17; i++) {        
      age_out[i] = age1;  
      age_in[i] = age1;
    }
    age_out[16] = age2;
    age_in[0] = 0.0;
           
    /* Create vectors of disease parameters (by age) to use in derivatives - includes BCG effect */
    double a_age[17];
    a_age[0] = a0*(BCG_cov*(1-BCG_eff) + (1-BCG_cov));
    a_age[1] = a5*(BCG_cov*(1-BCG_eff) + (1-BCG_cov));
    a_age[2] = a10*(BCG_cov*(1-BCG_eff) + (1-BCG_cov));
    double v_age[17];
    v_age[0] = v*(BCG_cov*(1-BCG_eff) + (1-BCG_cov));
    v_age[1] = v*(BCG_cov*(1-BCG_eff) + (1-BCG_cov));
    v_age[2] = v*(BCG_cov*(1-BCG_eff) + (1-BCG_cov));
    double sig_age[17];
    sig_age[0] = sig0;
    sig_age[1] = sig5;
    sig_age[2] = sig10;
    for (i=3; i<17; i++) {
      a_age[i] = a_a;
      sig_age[i] = sig_a;
      v_age[i] = v;
    }
 
    double muN_age[17];
    double muI_age[17];
    muN_age[0] = mu_N0;
    muI_age[0] = mu_I0;
    for (i=1; i<17; i++) {
      muN_age[i] = mu_N;
      muI_age[i] = mu_I;
    }

    /* Now adjust parameters to for HIV and ART

    /* and adjust TB parameters for HIV - mortality rates are passed straight in*/
    double mid_CD4[7] = {500,425,300,225,150,75,25}; /* mid points of CD4 categories */
    /* and create vectors of ART RRs */
    double ART_TB[3] = {ART_TB1,ART_TB2,ART_TB3};
    double ART_mort[3] = {ART_mort1,ART_mort2,ART_mort3};
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
          a_age_A[i][j][l] = a_age_H[i][j]*ART_TB[l];
          v_age_A[i][j][l] = v_age_H[i][j]*ART_TB[l];
          p_A[j][l] = p_H[j]/ART_TB[l];
        }
      }
    }
    
    /* Set up parameters for HIV model - these are taken from AIM */

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
      {0,0,0,0.2,0.2,0.213,0.213,0.299,0.299,0.414,0.414,0.414,0.414,0.414,0.414,0.414,0.414},
      {0,0,0,0.255,0.255,0.364,0.364,0.441,0.441,0.575,0.575,0.575,0.575,0.575,0.575,0.575,0.575},
      {0,0,0,0.398,0.398,0.663,0.663,0.713,0.713,0.838,0.838,0.838,0.838,0.838,0.838,0.838,0.838},
      {0,0,0,0.193,0.193,0.471,0.471,0.491,0.491,0.614,0.614,0.614,0.614,0.614,0.614,0.614,0.614},
      {0,0,0,0.294,0.294,0.765,0.765,0.765,0.765,0.865,0.865,0.865,0.865,0.865,0.865,0.865,0.865},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    };
    
    double H_mort[7][17] = { /* Mortality due to HIV (no ART) (age, CD4) */
      {0,0,0,0.005,0.005,0.004,0.004,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005},
      {0,0,0,0.011,0.011,0.01,0.01,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013},
      {0,0,0,0.026,0.026,0.026,0.026,0.036,0.036,0.032,0.032,0.032,0.032,0.032,0.032,0.032,0.032},
      {0,0,0,0.061,0.061,0.069,0.069,0.096,0.096,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08},
      {0,0,0,0.139,0.139,0.185,0.185,0.258,0.258,0.203,0.203,0.203,0.203,0.203,0.203,0.203,0.203},
      {0,0,0,0.321,0.321,0.499,0.499,0.691,0.691,0.513,0.513,0.513,0.513,0.513,0.513,0.513,0.513},
      {0,0,0,0.737,0.737,1.342,1.342,1.851,1.851,1.295,1.295,1.295,1.295,1.295,1.295,1.295,1.295}
    };
    
    double A_mort[3][7][17] = { /* On ART mortality (age,starting CD4, time on ART) */
      {{0,0,0,0.0050,0.0050,0.0035,0.0035,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050},
       {0,0,0,0.0115,0.0115,0.0095,0.0095,0.0134,0.0134,0.0126,0.0126,0.0126,0.0126,0.0126,0.0126,0.0126,0.0126},
       {0,0,0,0.0264,0.0264,0.0256,0.0256,0.0359,0.0359,0.0319,0.0319,0.0319,0.0319,0.0319,0.0319,0.0319,0.0319},
       {0,0,0,0.0607,0.0607,0.0563,0.0563,0.0587,0.0587,0.0594,0.0594,0.0594,0.0594,0.0594,0.0594,0.0594,0.0594},
       {0,0,0,0.1113,0.1113,0.0929,0.0929,0.0980,0.0980,0.1039,0.1039,0.1039,0.1039,0.1039,0.1039,0.1039,0.1039},
       {0,0,0,0.1810,0.1810,0.1525,0.1525,0.1619,0.1619,0.1762,0.1762,0.1762,0.1762,0.1762,0.1762,0.1762,0.1762},
       {0,0,0,0.3974,0.3974,0.3373,0.3373,0.3605,0.3605,0.4009,0.4009,0.4009,0.4009,0.4009,0.4009,0.4009,0.4009}},
      {{0,0,0,0.0050,0.0050,0.0035,0.0035,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050},
       {0,0,0,0.0115,0.0115,0.0095,0.0095,0.0134,0.0134,0.0126,0.0126,0.0126,0.0126,0.0126,0.0126,0.0126,0.0126},
       {0,0,0,0.0244,0.0244,0.0256,0.0256,0.0323,0.0323,0.0319,0.0319,0.0319,0.0319,0.0319,0.0319,0.0319,0.0319},
       {0,0,0,0.0258,0.0258,0.0332,0.0332,0.0341,0.0341,0.0451,0.0451,0.0451,0.0451,0.0451,0.0451,0.0451,0.0451},
       {0,0,0,0.0323,0.0323,0.0417,0.0417,0.0433,0.0433,0.0584,0.0584,0.0584,0.0584,0.0584,0.0584,0.0584,0.0584},
       {0,0,0,0.0405,0.0405,0.0523,0.0523,0.0548,0.0548,0.0751,0.0751,0.0751,0.0751,0.0751,0.0751,0.0751,0.0751},
       {0,0,0,0.0583,0.0583,0.0753,0.0753,0.0797,0.0797,0.1112,0.1112,0.1112,0.1112,0.1112,0.1112,0.1112,0.1112}},
      {{0,0,0,0.0050,0.0050,0.0035,0.0035,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050,0.0050},
       {0,0,0,0.0086,0.0086,0.0095,0.0095,0.0101,0.0101,0.0102,0.0102,0.0102,0.0102,0.0102,0.0102,0.0102,0.0102},
       {0,0,0,0.0092,0.0092,0.0117,0.0117,0.0109,0.0109,0.0114,0.0114,0.0114,0.0114,0.0114,0.0114,0.0114,0.0114},
       {0,0,0,0.0098,0.0098,0.0125,0.0125,0.0118,0.0118,0.0127,0.0127,0.0127,0.0127,0.0127,0.0127,0.0127,0.0127},
       {0,0,0,0.0129,0.0129,0.0165,0.0165,0.0161,0.0161,0.0190,0.0190,0.0190,0.0190,0.0190,0.0190,0.0190,0.0190},
       {0,0,0,0.0168,0.0168,0.0216,0.0126,0.0216,0.0216,0.0268,0.0268,0.0268,0.0268,0.0268,0.0268,0.0268,0.0268},
       {0,0,0,0.0251,0.0251,0.0324,0.0324,0.0333,0.0333,0.0438,0.0438,0.0438,0.0438,0.0438,0.0438,0.0438,0.0438}}
    };      
    
    double A_prog[4] = {0,2,2,0}; /* Progression through time on ART, 6 monthly time blocks */
    
    double A_start[3] = {1,0,0};  /* Used to make sure ART initiations are only added to the fist time on ART box */ 
    
    /* sum up various totals - uses function sum_array(array,i_start,i_end) and then for loops to add in HIV+ */ 
    
    double Total_S = sumsum(S,0,16);
    double Total_Ls = sumsum(Lsn,0,16)+sumsum(Lsp,0,16);
    double Total_Lm = sumsum(Lmn,0,16)+sumsum(Lmp,0,16);
    double Total_Ns = sumsum(Nsn,0,16)+sumsum(Nsp,0,16);
    double Total_Nm = sumsum(Nmn,0,16)+sumsum(Nmp,0,16);
    double Total_Is = sumsum(Isn,0,16)+sumsum(Isp,0,16);
    double Total_Im = sumsum(Imn,0,16)+sumsum(Imp,0,16);
    
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
    double Total_L = Total_Ls + Total_Lm;
    double Total_N = Total_Ns + Total_Nm;
    double Total_I = Total_Is + Total_Im;
    double Total_DS = Total_Ns + Total_Is;
    double Total_MDR = Total_Nm + Total_Im;
    double Total = Total_S+Total_L+Total_N+Total_I;
    
    /* Sum up over CD4 categories, with and without ART and calculate rates of ART initiation */
    double CD4_dist[7] = {0,0,0,0,0,0,0};
    double CD4_dist_ART[7] = {0,0,0,0,0,0,0};
    double Tot_ART = 0;                   /* Total on ART */ 
    double ART_need = 0;                  /* Number needing ART */
    double ART_new = 0;                   /* Number who need to start ART */
    double ART_el = 0;                    /* Number who are eligible but not on ART */
    double ART_prop[7] = {0,0,0,0,0,0,0};  /* Proportion of CD4 category who should start ART */
    for (j=0; j<7; j++){
      for (i=0; i<17; i++){
        CD4_dist[j] = CD4_dist[j]+S_H[i][j]+Lsn_H[i][j]+Lsp_H[i][j]+Lmn_H[i][j]+Lmp_H[i][j]+
                      Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j]+ 
                      Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j];
        for (l=0; l<3; l++){
          CD4_dist_ART[j] = CD4_dist_ART[j]+S_A[i][j][l]+Lsn_A[i][j][l]+Lsp_A[i][j][l]+Lmn_A[i][j][l]+Lmp_A[i][j][l]+
                            Nsn_A[i][j][l]+Nsp_A[i][j][l]+Nmn_A[i][j][l]+Nmp_A[i][j][l]+ 
                            Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l];         
                      
        }            
      } 
      Tot_ART = Tot_ART + CD4_dist_ART[j];
      ART_need = ART_need + (CD4_dist[j] + CD4_dist_ART[j])*forc[j+35];
    }
    ART_new = fmax(0,ART_need - Tot_ART);
    /* Then work out where these should go by CD4 */
    for (j=Athresh; j<7; j++) ART_el = ART_el + CD4_dist[j];
    if (ART_el > 0){
      for (j=Athresh; j<7; j++) {
        if (CD4_dist[j] > 0) {
          ART_prop[j] = (CD4_dist[j]/ART_el)*(ART_new/CD4_dist[j]);      
        }
      }
    }
    
    /* Force of infection */
    double FS = beta*(Total_Ns*rel_inf + Total_Is)/Total; 
    double FM = fit_cost*beta*(Total_Nm*rel_inf + Total_Im)/Total; 
    
    /* Calculate deaths due to TB by age group - these will be distributed across all states to avoid double counting deaths */
    double TB_deaths[17];
    for (i=0; i<17; i++) {
      TB_deaths[i] = (Nsn[i]+Nsp[i]+Nmn[i]+Nmp[i])*muN_age[i] + 
                     (Isn[i]+Isp[i]+Imn[i]+Imp[i])*muI_age[i];
      for(j=0; j<7; j++){
        TB_deaths[i] = TB_deaths[i] + (Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j])*muN_H + 
                                      (Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j])*muI_H;
        for (l=0; l<3; l++){
          TB_deaths[i] = TB_deaths[i] + (Nsn_A[i][j][l]+Nsp_A[i][j][l]+Nmn_A[i][j][l]+Nmp_A[i][j][l])*muN_H*ART_mort[l] + 
                                        (Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l])*muI_H*ART_mort[l];
        }                               
      }
    /*TB_deaths[i] = TB_deaths[i]*pop_ad;*/
    } 
    
    /* Calculate deaths due to HIV by age group */
    double HIV_deaths[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    for (i=0; i<17; i++){
      for (j=0; j<7; j++){
        HIV_deaths[i] = HIV_deaths[i]+H_mort[j][i]*(S_H[i][j]+
                        Lsn_H[i][j]+Lsp_H[i][j]+Lmn_H[i][j]+Lmp_H[i][j]+
                        Nsn_H[i][j]+Nsp_H[i][j]+Nmn_H[i][j]+Nmp_H[i][j]+
                        Isn_H[i][j]+Isp_H[i][j]+Imn_H[i][j]+Imp_H[i][j]);
      }
    }
    
   /* Calculate deaths due to HIV when on ART by age group */
    double ART_deaths[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    for (i=0; i<17; i++){
      for (j=0; j<7; j++){
        for (l=0; l<3; l++){
          ART_deaths[i] = ART_deaths[i]+A_mort[l][j][i]*(S_A[i][j][l]+Lsn_A[i][j][l]+Lsp_A[i][j][l]+Lmn_A[i][j][l]+Lmp_A[i][j][l]+
                            Nsn_A[i][j][l]+Nsp_A[i][j][l]+Nmn_A[i][j][l]+Nmp_A[i][j][l]+ 
                            Isn_A[i][j][l]+Isp_A[i][j][l]+Imn_A[i][j][l]+Imp_A[i][j][l]);
        }
      }
    }
     
    /* Calculate total population by age - this is used to distribute those dying of TB/HIV as mortality rates already include TB deaths */
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
 
    /* Derivatives */ 
 
    /* HIV-: loop through ages*/ 
    
    /* Susceptible - note this is done separately to deal with births */
    dS[0] = s_birth*birth_rate*Total/1000 - age1*S[0] - (FS + FM)*S[0] - forc[18]*S[0] + 
            (HIV_deaths[0] + ART_deaths[0] + TB_deaths[0])*S[0]/tot_age[0];
    for (i=1; i<17; i++) dS[i] = age_in[i]*forc[i+1]*S[i-1] - age_out[i]*S[i] - (FS + FM)*S[i] - forc[i+18]*S[i] +
                                 (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*S[i]/tot_age[i];
    /* Other TB states */
    for (i=0; i<17; i++){
      
      /* Latent, ds, naive */
      dLsn[i] = age_in[i]*forc[i+1]*Lsn[i-1] - age_out[i]*Lsn[i] +
                                   FS*((1-a_age[i])*S[i] + (1-a_age[i])*(1-p)*(1-g)*Lmn[i]) -
                                   (v_age[i] + FS*a_age[i]*(1-p) + FM*a_age[i]*(1-p) + FM*(1-a_age[i])*(1-p)*g)*Lsn[i] +
                                   r*(Isn[i] + Nsn[i]) - 
                                   forc[i+18]*Lsn[i] +
                                   (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Lsn[i]/tot_age[i];
      /* Latent, ds, prev */
      dLsp[i] = age_in[i]*forc[i+1]*Lsp[i-1] - age_out[i]*Lsp[i] +
                                   FS*(1-a_age[i])*(1-p)*(1-g)*Lmp[i] -
                                   (v_age[i] + FS*a_age[i]*(1-p) + FM*a_age[i]*(1-p) + FM*(1-a_age[i])*(1-p)*g)*Lsp[i] +
                                   r*(Isp[i]+Nsp[i]) +
                                   k*l_s*(1-e)*tau_s*(Isn[i]+Isp[i]) + 
                                   k*l_s*(1-e)*tau_s*d*(Nsn[i]+Nsp[i]) - 
                                   forc[i+18]*Lsp[i] +
                                   (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Lsp[i]/tot_age[i];
      /* Latent, mdr, naive */ 
      dLmn[i] = age_in[i]*forc[i+1]*Lmn[i-1] - age_out[i]*Lmn[i] +
                                   FM*((1-a_age[i])*S[i] + (1-a_age[i])*(1-p)*g*Lsn[i]) -
                                   (v_age[i] + FM*a_age[i]*(1-p) + FS*a_age[i]*(1-p) + FS*(1-a_age[i])*(1-p)*(1-g))*Lmn[i] +
                                   r*(Imn[i]+Nmn[i]) - 
                                   forc[i+18]*Lmn[i] +
                                   (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Lmn[i]/tot_age[i];
      /* Latent, mdr, prev */
      dLmp[i] = age_in[i]*forc[i+1]*Lmp[i-1] - age_out[i]*Lmp[i] + 
                                   FM*(1-a_age[i])*(1-p)*g*Lsp[i] -
                                   (v_age[i] + FM*a_age[i]*(1-p) + FS*a_age[i]*(1-p) + FS*(1-a_age[i])*(1-p)*(1-g))*Lmp[i] +
                                   r*(Imp[i]+Nmp[i]) +
                                   k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(Imn[i]+d*Nmn[i]) +
                                   k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(Imp[i]+d*Nmp[i]) - 
                                   forc[i+18]*Lmp[i] +
                                   (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Lmp[i]/tot_age[i];
      /* Smear neg, ds, new */
      dNsn[i] = age_in[i]*forc[i+1]*Nsn[i-1] - age_out[i]*Nsn[i] +
                                   (v_age[i]*(1-sig_age[i]) + FS*a_age[i]*(1-p)*(1-sig_age[i]))*Lsn[i] +
                                   FS*a_age[i]*(1-sig_age[i])*(S[i] + (1-p)*Lmn[i]) -
                                   (theta + r + k*l_s*d + muN_age[i])*Nsn[i] - 
                                   forc[i+18]*Nsn[i] +
                                   (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Nsn[i]/tot_age[i];
      /* Smear neg, ds, prev */                             
      dNsp[i] = age_in[i]*forc[i+1]*Nsp[i-1] - age_out[i]*Nsp[i] + 
                                   (v_age[i]*(1-sig_age[i]) + FS*a_age[i]*(1-p)*(1-sig_age[i]))*Lsp[i] +
                                   FS*a_age[i]*(1-sig_age[i])*(1-p)*Lmp[i] -
                                   (theta + r + k*l_s*d*(1-e)*tau_s + k*l_s*d*e + muN_age[i])*Nsp[i] +
                                   k*l_s*d*(1-e)*(1-tau_s)*Nsn[i] - 
                                   forc[i+18]*Nsp[i] +
                                   (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Nsp[i]/tot_age[i];
      /* Smear neg, mdr, new */
      dNmn[i] = age_in[i]*forc[i+1]*Nmn[i-1] - age_out[i]*Nmn[i] +
                                   (v_age[i]*(1-sig_age[i]) + FM*a_age[i]*(1-p)*(1-sig_age[i]))*Lmn[i] +
                                   FM*a_age[i]*(1-sig_age[i])*(S[i] + (1-p)*Lsn[i]) -
                                   (theta + r + k*l_m*d*dst_n + k*l_s*d*(1-dst_n) + muN_age[i])*Nmn[i] - 
                                   forc[i+18]*Nmn[i] +
                                   (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Nmn[i]/tot_age[i];
      /* Smear neg, mdr, prev */
      dNmp[i] = age_in[i]*forc[i+1]*Nmp[i-1] - age_out[i]*Nmp[i] +
                                   (v_age[i]*(1-sig_age[i]) + FM*a_age[i]*(1-p)*(1-sig_age[i]))*Lmp[i] +
                                   FM*a_age[i]*(1-sig_age[i])*(1-p)*Lsp[i] -
                                   (theta + r + k*l_m*d*dst_p*tau_m + k*l_s*d*(1-dst_p)*tau_s*eff_p + muN_age[i])*Nmp[i] +
                                   k*l_s*d*e*(Nsn[i]+Nsp[i]) +
                                   (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*Nmn[i] - 
                                   forc[i+18]*Nmp[i] +
                                   (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Nmp[i]/tot_age[i];
    
      /* Smear pos, ds, new */
      dIsn[i] = age_in[i]*forc[i+1]*Isn[i-1] - age_out[i]*Isn[i] +
                                   (v_age[i]*sig_age[i] + FS*a_age[i]*sig_age[i]*(1-p))*Lsn[i] +
                                   FS*a_age[i]*sig_age[i]*S[i] +
                                   FS*a_age[i]*(1-p)*sig_age[i]*Lmn[i] +
                                   theta*Nsn[i] -
                                   (r + k*l_s + muI_age[i])*Isn[i] - 
                                   forc[i+18]*Isn[i] +
                                   (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Isn[i]/tot_age[i];

      /* Smear pos, ds, prev */
      dIsp[i] = age_in[i]*forc[i+1]*Isp[i-1] - age_out[i]*Isp[i] +
                                   (v_age[i]*sig_age[i] + FS*a_age[i]*sig_age[i]*(1-p))*Lsp[i] +
                                   FS*a_age[i]*sig_age[i]*(1-p)*Lmp[i] +
                                   theta*Nsp[i] +
                                   k*l_s*(1-e)*(1-tau_s)*Isn[i] -
                                   (r + k*l_s*(1-e)*tau_s + k*l_s*e + muI_age[i])*Isp[i] - 
                                   forc[i+18]*Isp[i] +
                                   (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Isp[i]/tot_age[i];
                                   
      /* Smear pos, mdr, new */
      dImn[i] = age_in[i]*forc[i+1]*Imn[i-1] - age_out[i]*Imn[i] +
                                   (v_age[i]*sig_age[i] + FM*a_age[i]*sig_age[i]*(1-p))*Lmn[i] +
                                   FM*a_age[i]*sig_age[i]*S[i] +
                                   FM*a_age[i]*(1-p)*sig_age[i]*Lsn[i] +
                                   theta*Nmn[i] -
                                   (r + k*l_m*dst_n + k*l_s*(1-dst_n) + muI_age[i])*Imn[i] - 
                                   forc[i+18]*Imn[i] +
                                   (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Imn[i]/tot_age[i];
    
      /* Smear pos, mdr, prev */

      dImp[i] = age_in[i]*forc[i+1]*Imp[i-1] - age_out[i]*Imp[i] +
                                      (v_age[i]*sig_age[i] + FM*a_age[i]*sig_age[i]*(1-p))*Lmp[i] +
                                      FM*a_age[i]*sig_age[i]*(1-p)*Lsp[i] +
                                      theta*Nmp[i] +
                                      (k*l_m*dst_n*(1-tau_m) + k*l_s*(1-dst_n)*(1-(tau_s*eff_n)))*Imn[i] +
                                      k*l_s*e*(Isn[i]+Isp[i]) -
                                      (r + k*l_m*dst_p*tau_m + k*l_s*(1-dst_p)*tau_s*eff_p + muI_age[i])*Imp[i] - 
                                      forc[i+18]*Imp[i] +
                                      (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Imp[i]/tot_age[i];
    }

    /* HIV+: Loop through CD4 categories and age*/
  
    for (j=0; j<7; j++){      /* CD4 */
      for (i=0; i<17; i++){   /* age */
        dS_H[i][j] = age_in[i]*forc[i+1]*S_H[i-1][j] - age_out[i]*S_H[i][j] - (FS + FM)*S_H[i][j] +
                     forc[i+18]*H_CD4[j][i]*S[i] - H_prog[j+1][i]*S_H[i][j] + H_prog[j][i]*S_H[i][j-1] - 
                     H_mort[j][i]*S_H[i][j] - ART_prop[j]*S_H[i][j] + 
                     (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*S_H[i][j]/tot_age[i];

        dLsn_H[i][j] = age_in[i]*forc[i+1]*Lsn_H[i-1][j] - age_out[i]*Lsn_H[i][j] +
                       FS*((1-a_age_H[i][j])*S_H[i][j] + (1-a_age_H[i][j])*(1-p_H[j])*(1-g)*Lmn_H[i][j]) -
                       (v_age_H[i][j] + FS*a_age_H[i][j]*(1-p_H[j]) + FM*a_age_H[i][j]*(1-p_H[j]) + FM*(1-a_age_H[i][j])*(1-p_H[j])*g)*Lsn_H[i][j] +
                       r*(Isn_H[i][j] + Nsn_H[i][j]) +
                       forc[i+18]*H_CD4[j][i]*Lsn[i] - H_prog[j+1][i]*Lsn_H[i][j] + H_prog[j][i]*Lsn_H[i][j-1] - 
                       H_mort[j][i]*Lsn_H[i][j] - ART_prop[j]*Lsn_H[i][j] + 
                       (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Lsn_H[i][j]/tot_age[i];
        
        dLsp_H[i][j] = age_in[i]*forc[i+1]*Lsp_H[i-1][j] - age_out[i]*Lsp_H[i][j] +
                       FS*(1-a_age_H[i][j])*(1-p_H[j])*(1-g)*Lmp_H[i][j] -
                       (v_age_H[i][j] + FS*a_age_H[i][j]*(1-p_H[j]) + FM*a_age_H[i][j]*(1-p_H[j]) + FM*(1-a_age_H[i][j])*(1-p_H[j])*g)*Lsp_H[i][j] +
                       r*(Isp_H[i][j]+Nsp_H[i][j]) +
                       k*l_s*(1-e)*tau_s*(Isn_H[i][j]+Isp_H[i][j]) + 
                       k*l_s*(1-e)*tau_s*d*(Nsn_H[i][j]+Nsp_H[i][j]) +
                       forc[i+18]*H_CD4[j][i]*Lsp[i] - H_prog[j+1][i]*Lsp_H[i][j] + H_prog[j][i]*Lsp_H[i][j-1] - 
                       H_mort[j][i]*Lsp_H[i][j] - ART_prop[j]*Lsp_H[i][j] + 
                       (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Lsp_H[i][j]/tot_age[i];
                             
        dLmn_H[i][j] = age_in[i]*forc[i+1]*Lmn_H[i-1][j] - age_out[i]*Lmn_H[i][j] +
                       FM*((1-a_age_H[i][j])*S_H[i][j] + (1-a_age_H[i][j])*(1-p_H[j])*g*Lsn_H[i][j]) -
                       (v_age_H[i][j] + FM*a_age_H[i][j]*(1-p_H[j]) + FS*a_age_H[i][j]*(1-p_H[j]) + FS*(1-a_age_H[i][j])*(1-p_H[j])*(1-g))*Lmn_H[i][j] +
                       r*(Imn_H[i][j]+Nmn_H[i][j]) +
                       forc[i+18]*H_CD4[j][i]*Lmn[i] - H_prog[j+1][i]*Lmn_H[i][j] + H_prog[j][i]*Lmn_H[i][j-1] - 
                       H_mort[j][i]*Lmn_H[i][j] - ART_prop[j]*Lmn_H[i][j] + 
                       (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Lmn_H[i][j]/tot_age[i];
                     
        dLmp_H[i][j] = age_in[i]*forc[i+1]*Lmp_H[i-1][j] - age_out[i]*Lmp_H[i][j] + 
                       FM*(1-a_age_H[i][j])*(1-p_H[j])*g*Lsp_H[i][j] -
                       (v_age_H[i][j] + FM*a_age_H[i][j]*(1-p_H[j]) + FS*a_age_H[i][j]*(1-p_H[j]) + FS*(1-a_age_H[i][j])*(1-p_H[j])*(1-g))*Lmp_H[i][j] +
                       r*(Imp_H[i][j]+Nmp_H[i][j]) +
                       k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(Imn_H[i][j]+d*Nmn_H[i][j]) +
                       k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(Imp_H[i][j]+d*Nmp_H[i][j]) +
                       forc[i+18]*H_CD4[j][i]*Lmp[i] - H_prog[j+1][i]*Lmp_H[i][j] + H_prog[j][i]*Lmp_H[i][j-1] - 
                       H_mort[j][i]*Lmp_H[i][j] - ART_prop[j]*Lmp_H[i][j] + 
                       (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Lmp_H[i][j]/tot_age[i];
           
        dNsn_H[i][j] = age_in[i]*forc[i+1]*Nsn_H[i-1][j] - age_out[i]*Nsn_H[i][j] +
                       (v_age_H[i][j]*(1-sig_age[i]) + FS*a_age_H[i][j]*(1-p_H[j])*(1-sig_age[i]))*Lsn_H[i][j] +
                       FS*a_age_H[i][j]*(1-sig_age[i])*(S_H[i][j] + (1-p_H[j])*Lmn_H[i][j]) -
                       (theta + r + k*l_s*d + muN_H)*Nsn_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Nsn[i] - H_prog[j+1][i]*Nsn_H[i][j] + H_prog[j][i]*Nsn_H[i][j-1] - 
                       H_mort[j][i]*Nsn_H[i][j] - ART_prop[j]*Nsn_H[i][j] + 
                       (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Nsn_H[i][j]/tot_age[i];
                                                
        dNsp_H[i][j] = age_in[i]*forc[i+1]*Nsp_H[i-1][j] - age_out[i]*Nsp_H[i][j] + 
                       (v_age_H[i][j]*(1-sig_age[i]) + FS*a_age_H[i][j]*(1-p_H[j])*(1-sig_age[i]))*Lsp_H[i][j] +
                       FS*a_age_H[i][j]*(1-sig_age[i])*(1-p_H[j])*Lmp_H[i][j] -
                       (theta + r + k*l_s*d*(1-e)*tau_s + k*l_s*d*e + muN_H)*Nsp_H[i][j] +
                       k*l_s*d*(1-e)*(1-tau_s)*Nsn_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Nsp[i] - H_prog[j+1][i]*Nsp_H[i][j] + H_prog[j][i]*Nsp_H[i][j-1] - 
                       H_mort[j][i]*Nsp_H[i][j] - ART_prop[j]*Nsp_H[i][j] + 
                       (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Nsp_H[i][j]/tot_age[i];
                       
        dNmn_H[i][j] = age_in[i]*forc[i+1]*Nmn_H[i-1][j] - age_out[i]*Nmn_H[i][j] +
                       (v_age_H[i][j]*(1-sig_age[i]) + FM*a_age_H[i][j]*(1-p_H[j])*(1-sig_age[i]))*Lmn_H[i][j] +
                       FM*a_age_H[i][j]*(1-sig_age[i])*(S_H[i][j] + (1-p_H[j])*Lsn_H[i][j]) -
                       (theta + r + k*l_m*d*dst_n + k*l_s*d*(1-dst_n) + muN_H)*Nmn_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Nmn[i] - H_prog[j+1][i]*Nmn_H[i][j] + H_prog[j][i]*Nmn_H[i][j-1] - 
                       H_mort[j][i]*Nmn_H[i][j] - ART_prop[j]*Nmn_H[i][j] + 
                       (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Nmn_H[i][j]/tot_age[i];
                       
        dNmp_H[i][j] = age_in[i]*forc[i+1]*Nmp_H[i-1][j] - age_out[i]*Nmp_H[i][j] +
                       (v_age_H[i][j]*(1-sig_age[i]) + FM*a_age_H[i][j]*(1-p_H[j])*(1-sig_age[i]))*Lmp_H[i][j] +
                       FM*a_age_H[i][j]*(1-sig_age[i])*(1-p_H[j])*Lsp_H[i][j] -
                       (theta + r + k*l_m*d*dst_p*tau_m + k*l_s*d*(1-dst_p)*tau_s*eff_p + muN_H)*Nmp_H[i][j] +
                       k*l_s*d*e*(Nsn_H[i][j]+Nsp_H[i][j]) +
                       (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*Nmn_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Nmp[i] - H_prog[j+1][i]*Nmp_H[i][j] + H_prog[j][i]*Nmp_H[i][j-1] - 
                       H_mort[j][i]*Nmp_H[i][j] - ART_prop[j]*Nmp_H[i][j] + 
                       (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Nmp_H[i][j]/tot_age[i];
           
        dIsn_H[i][j] = age_in[i]*forc[i+1]*Isn_H[i-1][j] - age_out[i]*Isn_H[i][j] +
                       (v_age_H[i][j]*sig_age[i] + FS*a_age_H[i][j]*sig_age[i]*(1-p_H[j]))*Lsn_H[i][j] +
                       FS*a_age_H[i][j]*sig_age[i]*S_H[i][j] +
                       FS*a_age_H[i][j]*(1-p_H[j])*sig_age[i]*Lmn_H[i][j] +
                       theta*Nsn_H[i][j] -
                       (r + k*l_s + muI_H)*Isn_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Isn[i] - H_prog[j+1][i]*Isn_H[i][j] + H_prog[j][i]*Isn_H[i][j-1] - 
                       H_mort[j][i]*Isn_H[i][j] - ART_prop[j]*Isn_H[i][j] + 
                       (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Isn_H[i][j]/tot_age[i];

        dIsp_H[i][j] = age_in[i]*forc[i+1]*Isp_H[i-1][j] - age_out[i]*Isp_H[i][j] +
                       (v_age_H[i][j]*sig_age[i] + FS*a_age_H[i][j]*sig_age[i]*(1-p_H[j]))*Lsp_H[i][j] +
                       FS*a_age_H[i][j]*sig_age[i]*(1-p_H[j])*Lmp_H[i][j] +
                       theta*Nsp_H[i][j] +
                       k*l_s*(1-e)*(1-tau_s)*Isn_H[i][j] -
                       (r + k*l_s*(1-e)*tau_s + k*l_s*e + muI_H)*Isp_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Isp[i] - H_prog[j+1][i]*Isp_H[i][j] + H_prog[j][i]*Isp_H[i][j-1] - 
                       H_mort[j][i]*Isp_H[i][j] - ART_prop[j]*Isp_H[i][j] + 
                       (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Isp_H[i][j]/tot_age[i];
                                   
        dImn_H[i][j] = age_in[i]*forc[i+1]*Imn_H[i-1][j] - age_out[i]*Imn_H[i][j] +
                       (v_age_H[i][j]*sig_age[i] + FM*a_age_H[i][j]*sig_age[i]*(1-p_H[j]))*Lmn_H[i][j] +
                       FM*a_age_H[i][j]*sig_age[i]*S_H[i][j] +
                       FM*a_age_H[i][j]*(1-p_H[j])*sig_age[i]*Lsn_H[i][j] +
                       theta*Nmn_H[i][j] -
                       (r + k*l_m*dst_n + k*l_s*(1-dst_n) + muI_H)*Imn_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Imn[i] - H_prog[j+1][i]*Imn_H[i][j] + H_prog[j][i]*Imn_H[i][j-1] - 
                       H_mort[j][i]*Imn_H[i][j] - ART_prop[j]*Imn_H[i][j] + 
                       (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Imn_H[i][j]/tot_age[i];

        dImp_H[i][j] = age_in[i]*forc[i+1]*Imp_H[i-1][j] - age_out[i]*Imp_H[i][j] +
                       (v_age_H[i][j]*sig_age[i] + FM*a_age_H[i][j]*sig_age[i]*(1-p_H[j]))*Lmp_H[i][j] +
                       FM*a_age_H[i][j]*sig_age[i]*(1-p_H[j])*Lsp_H[i][j] +
                       theta*Nmp_H[i][j] +
                       (k*l_m*dst_n*(1-tau_m) + k*l_s*(1-dst_n)*(1-(tau_s*eff_n)))*Imn_H[i][j] +
                       k*l_s*e*(Isn_H[i][j]+Isp_H[i][j]) -
                       (r + k*l_m*dst_p*tau_m + k*l_s*(1-dst_p)*tau_s*eff_p + muI_H)*Imp_H[i][j] +
                       forc[i+18]*H_CD4[j][i]*Imp[i] - H_prog[j+1][i]*Imp_H[i][j] + H_prog[j][i]*Imp_H[i][j-1] - 
                       H_mort[j][i]*Imp_H[i][j] - ART_prop[j]*Imp_H[i][j] + 
                       (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Imp_H[i][j]/tot_age[i];
                     
      }
    }
    
    /* HIV+ on ART: loop through time on ART, CD4 at initiation, age - NEED TO ADD IN IMPACT OF ART ON TB PARAMETERS */
    for (l=0; l<3; l++){
      for (j=0; j<7; j++){
        for (i=0; i<17; i++){
          dS_A[i][j][l] = age_in[i]*forc[i+1]*S_A[i-1][j][l] - age_out[i]*S_A[i][j][l] - (FS + FM)*S_A[i][j][l] +
                          ART_prop[j]*A_start[l]*S_H[i][j] + A_prog[l]*S_A[i][j][l-1] - A_prog[l+1]*S_A[i][j][l] - A_mort[l][j][i]*S_A[i][j][l] +
                          (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*S_A[i][j][l]/tot_age[i];

          dLsn_A[i][j][l] = age_in[i]*forc[i+1]*Lsn_A[i-1][j][l] - age_out[i]*Lsn_A[i][j][l] +
                            FS*((1-a_age_A[i][j][l])*S_A[i][j][l] + (1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g)*Lmn_A[i][j][l]) -
                            (v_age_A[i][j][l] + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g)*Lsn_A[i][j][l] +
                            r*(Isn_A[i][j][l] + Nsn_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lsn_H[i][j] + A_prog[l]*Lsn_A[i][j][l-1] - A_prog[l+1]*Lsn_A[i][j][l] - A_mort[l][j][i]*Lsn_A[i][j][l] +
                            (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Lsn_A[i][j][l]/tot_age[i];
        
          dLsp_A[i][j][l] = age_in[i]*forc[i+1]*Lsp_A[i-1][j][l] - age_out[i]*Lsp_A[i][j][l] +
                            FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g)*Lmp_A[i][j][l] -
                            (v_age_A[i][j][l] + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g)*Lsp_A[i][j][l] +
                            r*(Isp_A[i][j][l]+Nsp_A[i][j][l]) +
                            k*l_s*(1-e)*tau_s*(Isn_A[i][j][l]+Isp_A[i][j][l]) + 
                            k*l_s*(1-e)*tau_s*d*(Nsn_A[i][j][l]+Nsp_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lsp_H[i][j] + A_prog[l]*Lsp_A[i][j][l-1] - A_prog[l+1]*Lsp_A[i][j][l] - A_mort[l][j][i]*Lsp_A[i][j][l] +
                            (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Lsp_A[i][j][j]/tot_age[i];
                             
          dLmn_A[i][j][l] = age_in[i]*forc[i+1]*Lmn_A[i-1][j][l] - age_out[i]*Lmn_A[i][j][l] +
                            FM*((1-a_age_A[i][j][l])*S_A[i][j][l] + (1-a_age_A[i][j][l])*(1-p_A[j][l])*g*Lsn_A[i][j][l]) -
                            (v_age_A[i][j][l] + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g))*Lmn_A[i][j][l] +
                            r*(Imn_A[i][j][l]+Nmn_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lmn_H[i][j] + A_prog[l]*Lmn_A[i][j][l-1] - A_prog[l+1]*Lmn_A[i][j][l] - A_mort[l][j][i]*Lmn_A[i][j][l] +
                            (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Lmn_A[i][j][l]/tot_age[i];
                     
          dLmp_A[i][j][l] = age_in[i]*forc[i+1]*Lmp_A[i-1][j][l] - age_out[i]*Lmp_A[i][j][l] + 
                            FM*(1-a_age_A[i][j][l])*(1-p_A[j][l])*g*Lsp_A[i][j][l] -
                            (v_age_A[i][j][l] + FM*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*a_age_A[i][j][l]*(1-p_A[j][l]) + FS*(1-a_age_A[i][j][l])*(1-p_A[j][l])*(1-g))*Lmp_A[i][j][l] +
                            r*(Imp_A[i][j][l]+Nmp_A[i][j][l]) +
                            k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(Imn_A[i][j][l]+d*Nmn_A[i][j][l]) +
                            k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(Imp_A[i][j][l]+d*Nmp_A[i][j][l]) +
                            ART_prop[j]*A_start[l]*Lmp_H[i][j] + A_prog[l]*Lmp_A[i][j][l-1] - A_prog[l+1]*Lmp_A[i][j][l] - A_mort[l][j][i]*Lmp_A[i][j][l] +
                            (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Lmp_A[i][j][l]/tot_age[i];
           
           /* GOT TO HERE REPLACING TB PARAMS WITH ART VALUES */
           
          dNsn_A[i][j][l] = age_in[i]*forc[i+1]*Nsn_A[i-1][j][l]- age_out[i]*Nsn_A[i][j][l] +
                            (v_age_A[i][j][l]*(1-sig_age[i]) + FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_age[i]))*Lsn_A[i][j][l] +
                            FS*a_age_A[i][j][l]*(1-sig_age[i])*(S_A[i][j][l] + (1-p_A[j][l])*Lmn_A[i][j][l]) -
                            (theta + r + k*l_s*d + muN_H*ART_mort[l])*Nsn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Nsn_H[i][j] + A_prog[l]*Nsn_A[i][j][l-1] - A_prog[l+1]*Nsn_A[i][j][l] - A_mort[l][j][i]*Nsn_A[i][j][l] +
                            (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Nsn_A[i][j][l]/tot_age[i];
                                                
          dNsp_A[i][j][l] = age_in[i]*forc[i+1]*Nsp_A[i-1][j][l] - age_out[i]*Nsp_A[i][j][l] + 
                            (v_age_A[i][j][l]*(1-sig_age[i]) + FS*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_age[i]))*Lsp_A[i][j][l] +
                            FS*a_age_A[i][j][l]*(1-sig_age[i])*(1-p_A[j][l])*Lmp_A[i][j][l] -
                            (theta + r + k*l_s*d*(1-e)*tau_s + k*l_s*d*e + muN_H*ART_mort[l])*Nsp_A[i][j][l] +
                            k*l_s*d*(1-e)*(1-tau_s)*Nsn_A[i][j][l] +      
                            ART_prop[j]*A_start[l]*Nsp_H[i][j] + A_prog[l]*Nsp_A[i][j][l-1] - A_prog[l+1]*Nsp_A[i][j][l] - A_mort[l][j][i]*Nsp_A[i][j][l] +       
                            (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Nsp_A[i][j][l]/tot_age[i];
                       
          dNmn_A[i][j][l] = age_in[i]*forc[i+1]*Nmn_A[i-1][j][l] - age_out[i]*Nmn_A[i][j][l] +
                            (v_age_A[i][j][l]*(1-sig_age[i]) + FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_age[i]))*Lmn_A[i][j][l] +
                            FM*a_age_A[i][j][l]*(1-sig_age[i])*(S_A[i][j][l] + (1-p_A[j][l])*Lsn_A[i][j][l]) -
                            (theta + r + k*l_m*d*dst_n + k*l_s*d*(1-dst_n) + muN_H*ART_mort[l])*Nmn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Nmn_H[i][j] + A_prog[l]*Nmn_A[i][j][l-1] - A_prog[l+1]*Nmn_A[i][j][l] - A_mort[l][j][i]*Nmn_A[i][j][l] + 
                            (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Nmn_A[i][j][l]/tot_age[i];
                       
          dNmp_A[i][j][l] = age_in[i]*forc[i+1]*Nmp_A[i-1][j][l] - age_out[i]*Nmp_A[i][j][l] +
                            (v_age_A[i][j][l]*(1-sig_age[i]) + FM*a_age_A[i][j][l]*(1-p_A[j][l])*(1-sig_age[i]))*Lmp_A[i][j][l] +
                            FM*a_age_A[i][j][l]*(1-sig_age[i])*(1-p_A[j][l])*Lsp_A[i][j][l] -
                            (theta + r + k*l_m*d*dst_p*tau_m + k*l_s*d*(1-dst_p)*tau_s*eff_p + muN_H*ART_mort[l])*Nmp_A[i][j][l] +
                            k*l_s*d*e*(Nsn_A[i][j][l]+Nsp_A[i][j][l]) +
                            (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*Nmn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Nmp_H[i][j] + A_prog[l]*Nmp_A[i][j][l-1] - A_prog[l+1]*Nmp_A[i][j][l] - A_mort[l][j][i]*Nmp_A[i][j][l] + 
                            (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Nmp_A[i][j][l]/tot_age[i];
           
          dIsn_A[i][j][l] = age_in[i]*forc[i+1]*Isn_A[i-1][j][l] - age_out[i]*Isn_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_age[i] + FS*a_age_A[i][j][l]*sig_age[i]*(1-p_A[j][l]))*Lsn_A[i][j][l] +
                            FS*a_age_A[i][j][l]*sig_age[i]*S_A[i][j][l] +
                            FS*a_age_A[i][j][l]*(1-p_A[j][l])*sig_age[i]*Lmn_A[i][j][l] +
                            theta*Nsn_A[i][j][l] -
                            (r + k*l_s + muI_H*ART_mort[l])*Isn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Isn_H[i][j] + A_prog[l]*Isn_A[i][j][l-1] - A_prog[l+1]*Isn_A[i][j][l] - A_mort[l][j][i]*Isn_A[i][j][l] +
                            (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Isn_A[i][j][l]/tot_age[i];

          dIsp_A[i][j][l] = age_in[i]*forc[i+1]*Isp_A[i-1][j][l] - age_out[i]*Isp_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_age[i] + FS*a_age_A[i][j][l]*sig_age[i]*(1-p_A[j][l]))*Lsp_A[i][j][l] +
                            FS*a_age_A[i][j][l]*sig_age[i]*(1-p_A[j][l])*Lmp_A[i][j][l] +
                            theta*Nsp_A[i][j][l] +
                            k*l_s*(1-e)*(1-tau_s)*Isn_A[i][j][l] -
                            (r + k*l_s*(1-e)*tau_s + k*l_s*e + muI_H*ART_mort[l])*Isp_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Isp_H[i][j] + A_prog[l]*Isp_A[i][j][l-1] - A_prog[l+1]*Isp_A[i][j][l] - A_mort[l][j][i]*Isp_A[i][j][l] +
                            (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Isp_A[i][j][l]/tot_age[i];
                                   
          dImn_A[i][j][l] = age_in[i]*forc[i+1]*Imn_A[i-1][j][l] - age_out[i]*Imn_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_age[i] + FM*a_age_A[i][j][l]*sig_age[i]*(1-p_A[j][l]))*Lmn_A[i][j][l] +
                            FM*a_age_A[i][j][l]*sig_age[i]*S_A[i][j][l] +
                            FM*a_age_A[i][j][l]*(1-p_A[j][l])*sig_age[i]*Lsn_A[i][j][l] +
                            theta*Nmn_A[i][j][l] -
                            (r + k*l_m*dst_n + k*l_s*(1-dst_n) + muI_H*ART_mort[l])*Imn_A[i][j][l] +
                            ART_prop[j]*A_start[l]*Imn_H[i][j] + A_prog[l]*Imn_A[i][j][l-1] - A_prog[l+1]*Imn_A[i][j][l] - A_mort[l][j][i]*Imn_A[i][j][l] + 
                            (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Imn_A[i][j][l]/tot_age[i];

          dImp_A[i][j][l] = age_in[i]*forc[i+1]*Imp_A[i-1][j][l] - age_out[i]*Imp_A[i][j][l] +
                            (v_age_A[i][j][l]*sig_age[i] + FM*a_age_A[i][j][l]*sig_age[i]*(1-p_A[j][l]))*Lmp_A[i][j][l] +
                            FM*a_age_A[i][j][l]*sig_age[i]*(1-p_A[j][l])*Lsp_A[i][j][l] +
                            theta*Nmp_A[i][j][l] +
                            (k*l_m*dst_n*(1-tau_m) + k*l_s*(1-dst_n)*(1-(tau_s*eff_n)))*Imn_A[i][j][l] +
                            k*l_s*e*(Isn_A[i][j][l]+Isp_A[i][j][l]) -
                            (r + k*l_m*dst_p*tau_m + k*l_s*(1-dst_p)*tau_s*eff_p + muI_H*ART_mort[l])*Imp_A[i][j][l] +                     
                            ART_prop[j]*A_start[l]*Imp_H[i][j] + A_prog[l]*Imp_A[i][j][l-1] - A_prog[l+1]*Imp_A[i][j][l] - A_mort[l][j][i]*Imp_A[i][j][l] + 
                            (HIV_deaths[i] + ART_deaths[i] + TB_deaths[i])*Imp_A[i][j][l]/tot_age[i];
               
        }
      }
    }

    /* Put function values into ydot */
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
    
    /* add variables to yout */
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
}




