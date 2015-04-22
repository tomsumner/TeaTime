/* Attempt to write TB model in C code to call from R */
/* Version 4 is first attempt to add in HIV - structure is based on AIM */

/* just add it for TB susceptible first and run with no TB to see if we can get the HIV dynamics correct 
have also dropped TB death adjustment for now, will need to add this back in together with adjustment for HIV mortality */

/* compile within R with system("R CMD SHLIB Demog_model.c") */

#include <R.h>

static double parms[33];
static double forc[35];

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
#define k parms[23]
#define l_s parms[24]
#define l_m parms[25]
#define d parms[26]
#define tau_s parms[27]
#define tau_m parms[28]
#define eff_n parms[29]
#define eff_p parms[30]
#define dst_n parms[31]
#define dst_p parms[32]

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
    int N=33;
    odeparms(&N, parms);
}

void forcc(void (* odeforcs)(int *, double *))
{
    int N=35;
    odeforcs(&N, forc);
}


/* derivative function */
void derivsc(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if (ip[0] <2) error("nout should be at least 2");
    
    /* Expand state variables so can use more meaningful names than y and ydot  */
    double S[17];
    double S_H[17][7];
   /* double S_A[357]; */
    
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
    
    double dS[17];
    double dS_H[17][7]; 
  /*  double dS_A[357];*/
    
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
     
    int i=0;  
    int j=0; 
     
    for (i=0; i<17; i++) S[i] = y[i];             /* Susceptible, HIV- */

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
    
    int ih = 221;
    for (j=0; j<7; j++){
      for (i=0; i<17; i++){  /* sus, HIV+ */
        S_H[i][j] = y[ih];
        ih = ih+1;
      }
    }
    
    /* Create vectors of aging rates to use in derivatives */
    double age_out[17];                         /* age_out is lower for final age group because it is wider */
    for (i=0; i<16; i++) age_out[i] = age1;            
    age_out[16] = age2;
    
    double age_in[17];                          /* age_in is set to zero for first age group as new borns only enter the S class */
    for (i=1; i<17; i++) age_in[i] = age1;            
    age_in[0] = 0.0;
 
    /* Create vectors of disease parameters (by age) to use in derivatives */
    double a_age[17];
    a_age[0] = a0;
    a_age[1] = a5;
    a_age[2] = a10;
    for (i=3; i<17; i++) a_age[i] = a_a;
 
    double sig_age[17];
    sig_age[0] = sig0;
    sig_age[1] = sig5;
    sig_age[2] = sig10;
    for (i=3; i<17; i++) sig_age[i] = sig_a;
 
    double muN_age[17];
    muN_age[0] = mu_N0;
    for (i=1; i<17; i++) muN_age[i] = mu_N;
 
    double muI_age[17];
    muI_age[0] = mu_I0;
    for (i=1; i<17; i++) muI_age[i] = mu_I;
 
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
    
  /*  double A_mort[357];      /* On ART mortality (age,starting CD4, time on ART) */
  /*  double A_prog_out[357];  /* Progression through time on ART */
  /*  double A_prog_in[357];   /* Progression through time on ART */
  /*  for (i=0; i<357; i++) {
      
      A_mort[i] = 0;
      A_prog_out[i] = 0;
      A_prog_in[i] = 0;
      
    }*/
 
   /* double Astart[3];   /* ensure ART starters only enter first ART category */
  /*  Astart[0]=1;
    Astart[2]=0;
    Astart[3]=0;*/
    
    
   /* for (i=340; i<697; i++) S_A[i-340] = y[i];    /* sus, ART */
       
    /* sum up various totals - uses function sum_array(array,i_start,i_end) */ 
    
    double Total_S = sumsum(S,0,16);
    double Total_SH = 0;
    for (j=0; j<7; j++){
      for (i=0; i<17; i++){
        Total_SH = Total_SH + S_H[i][j];
      }
    }
    
    double Total_Ls = sumsum(Lsn,0,16)+sumsum(Lsp,0,16);
    double Total_Lm = sumsum(Lmn,0,16)+sumsum(Lmp,0,16);
    double Total_L = Total_Ls + Total_Lm;
    
    double Total_Ns = sumsum(Nsn,0,16)+sumsum(Nsp,0,16);
    double Total_Nm = sumsum(Nmn,0,16)+sumsum(Nmp,0,16);
    double Total_N = Total_Ns + Total_Nm;
    double Total_Is = sumsum(Isn,0,16)+sumsum(Isp,0,16);
    double Total_Im = sumsum(Imn,0,16)+sumsum(Imp,0,16);
    double Total_I = Total_Is + Total_Im;
 
    double Total_DS = Total_Ns + Total_Is;
    double Total_MDR = Total_Nm + Total_Im;
 
    double Total = Total_S+Total_SH+Total_L+Total_N+Total_I;
    
    double CD4_dist[7] = {0,0,0,0,0,0,0};
    
    for (j=0; j<7; j++){
      for (i=0; i<17; i++){
        CD4_dist[j] = CD4_dist[j]+S_H[i][j]
      } 
    }
 
    /* Force of infection */
    double FS = beta*(Total_Ns*rel_inf + Total_Is)/Total; 
    double FM = fit_cost*beta*(Total_Nm*rel_inf + Total_Im)/Total; 
    
    /* Calculate deaths due to TB by age group - these will be distributed across all states to avoid double counting deaths */
    double TB_deaths[17];
    for (i=0; i<17; i++) TB_deaths[i] = (Nsn[i]+Nsp[i]+Nmn[i]+Nmp[i])*muN_age[i] + 
                                        (Isn[i]+Isp[i]+Imn[i]+Imp[i])*muI_age[i]; /* will need updating to include HIV_ TB deaths */
     
     
    /* Calculate deaths due to HIV by age group */
    double HIV_deaths[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    for (i=0; i<17; i++){
      for (j=0; j<7; j++){
        HIV_deaths[i] = HIV_deaths[i]+H_mort[j][i]*S_H[i][j];
      }
    }
     
    /* Calculate total population by age - this is used to distribute those dying of TB/HIV as mortality rates already include TB deaths */
    double tot_age[17];
    for (i=0; i<17; i++){
      tot_age[i] = S[i]+Lsn[i]+Lsp[i]+Lmn[i]+Lmp[i]+Nsn[i]+Nsp[i]+Nmn[i]+Nmp[i]+Isn[i]+Isp[i]+Imn[i]+Imp[i];
      for (j=0; j<7; j++){
        tot_age[i] = tot_age[i]+S_H[i][j];
      }
    }
 
    /* Derivatives */ 
 
    /* Susceptible */ 
    
    /* HIV- */
    dS[0] = s_birth*birth_rate*Total/1000 - age1*S[0] - (FS + FM)*S[0] - forc[18]*S[0] + 
            TB_deaths[0]*S[0]/tot_age[0] + HIV_deaths[0]*S[0]/tot_age[0];
    for (i=1; i<17; i++) dS[i] = age_in[i]*forc[i+1]*S[i-1] - age_out[i]*S[i] - (FS + FM)*S[i] - forc[i+18]*S[i] +
                                 + TB_deaths[i]*S[i]/tot_age[i] + HIV_deaths[i]*S[i]/tot_age[i];

    
   /* for (i=0; i<119; i++) dS_H[i] = 0;*/
   /* for (i=0; i<357; i++) dS_A[i] = 0;*/

    /* HIV+
    Loop through CD4 categories and through ages for each category
    Need to work out how to stop aging running into next CD4 catergory and vice versa - think just need to set the H_prog to zero in the appropriate places
    */
  
   for (j=0; j<7; j++){      /* CD4 */
      for (i=0; i<17; i++){   /* age */
        dS_H[i][j] = age_in[i]*forc[i+1]*S_H[i-1][j] - age_out[i]*S_H[i][j] - (FS + FM)*S_H[i][j] +
                   forc[i+18]*H_CD4[j][i]*S[i] - H_prog[j+1][i]*S_H[i][j] + H_prog[j][i]*S_H[i][j-1] -
                   H_mort[j][i]*S_H[i][j] + HIV_deaths[i]*S_H[i][j]/tot_age[i];
                   /*H_ART[j][i]*S_H[i][j];*/
        ih = ih+1;
      }
    }
    
    /* HIV+, ART */
    
  /*  ih = 0;
    for (jt=0; jt<3; jt++){     /* time on ART */
    /*  ia = 0;
      for (j=0; j<7; j++){      /* CD4 at ART start */
     /*   for (i=0; i<17; i++) {  /* age */

     /*     dS_HA[ih] = age_in[i]*forc[i+1]*S_HA[ih-1] - age_out[i]*S_HA[ih] - (FS + FM)*S_HA[ih] - 
                      A_mort[ih]*S_HA[ih] - A_prog_out[jt]*S_HA[ih] + A_prog_in[jt]*S_HA[ih-119] + H_ART[ia]*S_H[ia]*Astart[jt];
          ih = ih+1;
          ia = ia+1;
        }
      }
    }
    */

    /* Latent,s,n*/

    for (i=0; i<17; i++) dLsn[i] = age_in[i]*forc[i+1]*Lsn[i-1] - age_out[i]*Lsn[i] +
                                   FS*((1-a_age[i])*S[i] + (1-a_age[i])*(1-p)*(1-g)*Lmn[i]) -
                                   (v + FS*a_age[i]*(1-p) + FM*a_age[i]*(1-p) + FM*(1-a_age[i])*(1-p)*g)*Lsn[i] +
                                   r*(Isn[i] + Nsn[i]);
    
    /* Latent,s,p [ok] */
  
    for (i=0; i<17; i++) dLsp[i] = age_in[i]*forc[i+1]*Lsp[i-1] - age_out[i]*Lsp[i] +
                                   FS*(1-a_age[i])*(1-p)*(1-g)*Lmp[i] -
                                   (v + FS*a_age[i]*(1-p) + FM*a_age[i]*(1-p) + FM*(1-a_age[i])*(1-p)*g)*Lsp[i] +
                                   r*(Isp[i]+Nsp[i]) +
                                   k*l_s*(1-e)*tau_s*(Isn[i]+Isp[i]) + 
                                   k*l_s*(1-e)*tau_s*d*(Nsn[i]+Nsp[i]);
    
    /* Latent,m,n [ok] */
                
    for (i=0; i<17; i++) dLmn[i] = age_in[i]*forc[i+1]*Lmn[i-1] - age_out[i]*Lmn[i] +
                                   FM*((1-a_age[i])*S[i] + (1-a_age[i])*(1-p)*g*Lsn[i]) -
                                   (v + FM*a_age[i]*(1-p) + FS*a_age[i]*(1-p) + FS*(1-a_age[i])*(1-p)*(1-g))*Lmn[i] +
                                   r*(Imn[i]+Nmn[i]);  
    
    /* Latent,m,p [ok] */
                
    for (i=0; i<17; i++) dLmp[i] = age_in[i]*forc[i+1]*Lmp[i-1] - age_out[i]*Lmp[i] + 
                                   FM*(1-a_age[i])*(1-p)*g*Lsp[i] -
                                   (v + FM*a_age[i]*(1-p) + FS*a_age[i]*(1-p) + FS*(1-a_age[i])*(1-p)*(1-g))*Lmp[i] +
                                   r*(Imp[i]+Nmp[i]) +
                                   k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(Imn[i]+d*Nmn[i]) +
                                   k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(Imp[i]+d*Nmp[i]);

    /* Smear negative,s,n */
                
    for (i=0; i<17; i++) dNsn[i] = age_in[i]*forc[i+1]*Nsn[i-1] - age_out[i]*Nsn[i] +
                                   (v*(1-sig_age[i]) + FS*a_age[i]*(1-p)*(1-sig_age[i]))*Lsn[i] +
                                   FS*a_age[i]*(1-sig_age[i])*(S[i] + (1-p)*Lmn[i]) -
                                   (theta + r + k*l_s*d + muN_age[i])*Nsn[i];

    /* Smear negative,s,p [ok] */
                    
    for (i=0; i<17; i++) dNsp[i] = age_in[i]*forc[i+1]*Nsp[i-1] - age_out[i]*Nsp[i] + 
                                   (v*(1-sig_age[i]) + FS*a_age[i]*(1-p)*(1-sig_age[i]))*Lsp[i] +
                                   FS*a_age[i]*(1-sig_age[i])*(1-p)*Lmp[i] -
                                   (theta + r + k*l_s*d*(1-e)*tau_s + k*l_s*d*e + muN_age[i])*Nsp[i] +
                                   k*l_s*d*(1-e)*(1-tau_s)*Nsn[i];

    /* Smear negative,m,n [ok] */

    for (i=0; i<17; i++) dNmn[i] = age_in[i]*forc[i+1]*Nmn[i-1] - age_out[i]*Nmn[i] +
                                   (v*(1-sig_age[i]) + FM*a_age[i]*(1-p)*(1-sig_age[i]))*Lmn[i] +
                                   FM*a_age[i]*(1-sig_age[i])*(S[i] + (1-p)*Lsn[i]) -
                                   (theta + r + k*l_m*d*dst_n + k*l_s*d*(1-dst_n) + muN_age[i])*Nmn[i];
    
    /* Smear negative,m,p [ok] */

    for (i=0; i<17; i++) dNmp[i] = age_in[i]*forc[i+1]*Nmp[i-1] - age_out[i]*Nmp[i] +
                                   (v*(1-sig_age[i]) + FM*a_age[i]*(1-p)*(1-sig_age[i]))*Lmp[i] +
                                   FM*a_age[i]*(1-sig_age[i])*(1-p)*Lsp[i] -
                                   (theta + r + k*l_m*d*dst_p*tau_m + k*l_s*d*(1-dst_p)*tau_s*eff_p + muN_age[i])*Nmp[i] +
                                   k*l_s*d*e*(Nsn[i]+Nsp[i]) +
                                   (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*Nmn[i]; 

    /* Smear positive,s,n [ok] */

    for (i=0; i<17; i++) dIsn[i] = age_in[i]*forc[i+1]*Isn[i-1] - age_out[i]*Isn[i] +
                                   (v*sig_age[i] + FS*a_age[i]*sig_age[i]*(1-p))*Lsn[i] +
                                   FS*a_age[i]*sig_age[i]*S[i] +
                                   FS*a_age[i]*(1-p)*sig_age[i]*Lmn[i] +
                                   theta*Nsn[i] -
                                   (r + k*l_s + muI_age[i])*Isn[i];

    /* Smear positive,s,p [ok] */

    for (i=0; i<17; i++) dIsp[i] = age_in[i]*forc[i+1]*Isp[i-1] - age_out[i]*Isp[i] +
                                   (v*sig_age[i] + FS*a_age[i]*sig_age[i]*(1-p))*Lsp[i] +
                                   FS*a_age[i]*sig_age[i]*(1-p)*Lmp[i] +
                                   theta*Nsp[i] +
                                   k*l_s*(1-e)*(1-tau_s)*Isn[i] -
                                   (r + k*l_s*(1-e)*tau_s + k*l_s*e + muI_age[i])*Isp[i];
                                   
    /* Smear positive,m,n [ok] */

    for (i=0; i<17; i++) dImn[i] = age_in[i]*forc[i+1]*Imn[i-1] - age_out[i]*Imn[i] +
                                   (v*sig_age[i] + FM*a_age[i]*sig_age[i]*(1-p))*Lmn[i] +
                                   FM*a_age[i]*sig_age[i]*S[i] +
                                   FM*a_age[i]*(1-p)*sig_age[i]*Lsn[i] +
                                   theta*Nmn[i] -
                                   (r + k*l_m*dst_n + k*l_s*(1-dst_n) + muI_age[i])*Imn[i];
    
    /* Smear positive,m,n [OK] */

    for (i=0; i<17; i++) dImp[i] = age_in[i]*forc[i+1]*Imp[i-1] - age_out[i]*Imp[i] +
                                      (v*sig_age[i] + FM*a_age[i]*sig_age[i]*(1-p))*Lmp[i] +
                                      FM*a_age[i]*sig_age[i]*(1-p)*Lsp[i] +
                                      theta*Nmp[i] +
                                      (k*l_m*dst_n*(1-tau_m) + k*l_s*(1-dst_n)*(1-(tau_s*eff_n)))*Imn[i] +
                                      k*l_s*e*(Isn[i]+Isp[i]) -
                                      (r + k*l_m*dst_p*tau_m + k*l_s*(1-dst_p)*tau_s*eff_p + muI_age[i])*Imp[i];
                
    /* Put function values into ydot */

    for (i=0; i<17; i++) ydot[i] = dS[i];         /* S: 0-16 */
     
    for (i=17; i<34; i++) ydot[i] = dLsn[i-17];      /* Lsn: 17-33 */
    for (i=34; i<51; i++) ydot[i] = dLsp[i-34];      /* Lsp: 34-50 */
    for (i=51; i<68; i++) ydot[i] = dLmn[i-51];      /* Lmn: 51-67 */
    for (i=68; i<85; i++) ydot[i] = dLmp[i-68];      /* Lmp: 68-84 */
     
    for (i=85; i<102; i++) ydot[i] = dNsn[i-85];     /* Nsn: 85-101 */
    for (i=102; i<119; i++) ydot[i] = dNsp[i-102];    /* Nsp: 102-118 */
    for (i=119; i<136; i++) ydot[i] = dNmn[i-119];    /* Nmn: 119-135 */
    for (i=136; i<153; i++) ydot[i] = dNmp[i-136];    /* Nmp: 136-152 */
       
    for (i=153; i<170; i++) ydot[i] = dIsn[i-153];    /* Isn: 153-169 */
    for (i=170; i<187; i++) ydot[i] = dIsp[i-170];    /* Isp: 170-186 */
    for (i=187; i<204; i++) ydot[i] = dImn[i-187];    /* Imn: 187-203 */
    for (i=204; i<221; i++) ydot[i] = dImp[i-204];    /* Imp: 204-220 */

    ih = 221;
    for (j=0; j<7; j++){
      for (i=0; i<17; i++){
        ydot[ih] = dS_H[i][j];
        ih=ih+1;
      }
    }

  /*  for (i=340; i<697; i++) ydot[i] = dS_A[i-340];*/
  
    /* add variables to yout */

    yout[0] = Total;
    yout[1] = Total_S;
    yout[2] = Total_SH;
    yout[3] = Total_Ls;
    yout[4] = Total_Lm;
    yout[5] = Total_L;
    yout[6] = Total_Ns;
    yout[7] = Total_Nm;
    yout[8] = Total_N;
    yout[9] = Total_Is;
    yout[10] = Total_Im;
    yout[11] = Total_I;
    yout[12] = Total_DS;
    yout[13] = Total_MDR;
    yout[14] = FS;
    yout[15] = FM;
    yout[16] = CD4_dist[0];
    yout[17] = CD4_dist[1];
    yout[18] = CD4_dist[2];
    yout[19] = CD4_dist[3];
    yout[20] = CD4_dist[4];
    yout[21] = CD4_dist[5];
    yout[22] = CD4_dist[6];
    
}

                                                                    
  





