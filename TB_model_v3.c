/* Attempt to write TB model in C code to call from R */
/* Version 4 is first attempt to add in HIV - structure is based on AIM */

/* just add it for TB susceptible first and run with no TB to see if we can get the HIV dynamics correct
/* have also dropped TB death adjustment for now, will need to add this back in together with adjustment for HIV mortality

/* compile within R with system("R CMD SHLIB Demog_model.c") */

#include <R.h>

static double parms[33];
static double forc[18];

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
    int N=18;
    odeforcs(&N, forc);
}


/* derivative function */
void derivsc(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if (ip[0] <2) error("nout should be at least 2");
    
    /* Expand state variables so can use more meaningful names than y and ydot  */
    double S_S[17];
    double S_H[119];
    double S_HA[357];
    
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
    
    double dS_S[17];
    double dS_H[119];
    double dS_HA[357];
    
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
     
    for (i=0; i<17; i++) S_S[i] = y[i];          
    for (i=17; i<136; i++) S_H[i-17] = y[i];      
    for (i=136; i<493; i++) S_HA[i-136] = y[i]; 
     
    for (i=493; i<510; i++) Lsn[i-493] = y[i];      
    for (i=510; i<527; i++) Lsp[i-510] = y[i];      
    for (i=527; i<544; i++) Lmn[i-527] = y[i];      
    for (i=544; i<561; i++) Lmp[i-544] = y[i];     
    
    for (i=561; i<578; i++) Nsn[i-561] = y[i];      
    for (i=578; i<595; i++) Nsp[i-578] = y[i];    
    for (i=595; i<612; i++) Nmn[i-595] = y[i];    
    for (i=612; i<629; i++) Nmp[i-612] = y[i];    
       
    for (i=629; i<646; i++) Isn[i-629] = y[i];    
    for (i=646; i<663; i++) Isp[i-646] = y[i];    
    for (i=663; i<680; i++) Imn[i-663] = y[i];    
    for (i=680; i<697; i++) Imp[i-680] = y[i];    
       
    /* sum up various totals - uses function sum_array(array,i_start,i_end) */ 
    
    double Total_S = sumsum(S_S,0,16) + sumsum(S_H,0,118) + sumsum(S_HA,0,356);
    
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
 
    double Total = Total_S+Total_L+Total_N+Total_I;
 
    /* Force of infection */
    double FS = beta*(Total_Ns*rel_inf + Total_Is)/Total; 
    double FM = fit_cost*beta*(Total_Nm*rel_inf + Total_Im)/Total; 
    
    /* Create vectors of aging rates to use in derivatives */
    double age_out[17];                         /* age_out is lower for final age group because it is wider */
    for (i=0; i<16; i++) age_out[i] = age1;            
    age_out[16] = age2;
    
    double age_in[17];                          /* age_in is set to zero for first age group as new borns only enter the S class */
    for (i=1; i<17; i++) age_in[i] = age1;            
    age_in[0] = 0;
 
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
 
    /* Derivatives */ 
 
    /* Susceptible */ 
    
    /* HIV- */
    dS_S[0] = s_birth*birth_rate*Total/1000 - age_out[0]*S_S[0] - (FS + FM)*S_S[0] - 
              H[0]*S_S[0];
    for (i=1; i<17; i++) dS_S[i] = age_in[i]*forc[i+1]*S_S[i-1] - age_out[i]*S_S[i] - (FS + FM)*S_S[i] - 
                                   H[i]*S_S[i];

    /* HIV+
    Loop through CD4 categories and through ages for each category
    Need to work out how to stop aging running into next CD4 catergory and vice versa - think just need to set the H_prog to zero in the appropriate places
    */
    int ih = 0;
    for (j=0; j<7; j++){
      for (i=0; i<17; i++) {

        dS_H[ih] = age_in[i]*forc[i+1]*S_H[ih-1] - age_out[i]*S_H[ih] - (FS + FM)*S_H[ih] + 
                  H[i]*H_CD4[ih]*S_S[i] - H_mort[ih]*S_H[ih] - H_prog[ih]*S_H[ih] + H_prog[j-1]*S_H[ih-1] - H_ART[ih]*S_H[ih];
        ih = ih++;
      
      }
    }
    
    /* HIV+, ART */
    
    int ih = 0;
    for (j=0; j<7; j++){
      for (i=0; i<17; i++) {

        dS_H[i] = age_in[i]*forc[i+1]*S_H[ih-1] - age_out[i]*S_H[ih] - (FS + FM)*S_H[ih] + 
                  H[i]*H_CD4[j]*S_S[ih] - H_mort[ih]*S_H[ih] - H_prog[ih]*S_H[ih] + H_prog[ih-1]*S_H[ih-1] - H_ART[ih]*S_H[ih];
        ih = ih++;
      
      }
    }
    
    

    /* Latent,s,n*/

    for (i=0; i<17; i++) dLsn[i] = age_in[i]*forc[i+1]*Lsn[i-1] - age_out[i]*Lsn[i] +
                                   FS*((1-a_age[i])*S[i] + (1-a_age[i])*(1-p)*(1-g)*Lmn[i]) -
                                   (v + FS*a_age[i]*(1-p) + FM*a_age[i]*(1-p) + FM*(1-a_age[i])*(1-p)*g)*Lsn[i] +
                                   r*(Isn[i] + Nsn[i]) +
                                   TB_deaths[i]*Lsn[i]/tot_age[i];
    
    /* Latent,s,p [ok] */
  
    for (i=0; i<17; i++) dLsp[i] = age_in[i]*forc[i+1]*Lsp[i-1] - age_out[i]*Lsp[i] +
                                   FS*(1-a_age[i])*(1-p)*(1-g)*Lmp[i] -
                                   (v + FS*a_age[i]*(1-p) + FM*a_age[i]*(1-p) + FM*(1-a_age[i])*(1-p)*g)*Lsp[i] +
                                   r*(Isp[i]+Nsp[i]) +
                                   k*l_s*(1-e)*tau_s*(Isn[i]+Isp[i]) + 
                                   k*l_s*(1-e)*tau_s*d*(Nsn[i]+Nsp[i]) +
                                   TB_deaths[i]*Lsp[i]/tot_age[i];
    
    /* Latent,m,n [ok] */
                
    for (i=0; i<17; i++) dLmn[i] = age_in[i]*forc[i+1]*Lmn[i-1] - age_out[i]*Lmn[i] +
                                   FM*((1-a_age[i])*S[i] + (1-a_age[i])*(1-p)*g*Lsn[i]) -
                                   (v + FM*a_age[i]*(1-p) + FS*a_age[i]*(1-p) + FS*(1-a_age[i])*(1-p)*(1-g))*Lmn[i] +
                                   r*(Imn[i]+Nmn[i]) +
                                   TB_deaths[i]*Lmn[i]/tot_age[i];  
    
    /* Latent,m,p [ok] */
                
    for (i=0; i<17; i++) dLmp[i] = age_in[i]*forc[i+1]*Lmp[i-1] - age_out[i]*Lmp[i] + 
                                   FM*(1-a_age[i])*(1-p)*g*Lsp[i] -
                                   (v + FM*a_age[i]*(1-p) + FS*a_age[i]*(1-p) + FS*(1-a_age[i])*(1-p)*(1-g))*Lmp[i] +
                                   r*(Imp[i]+Nmp[i]) +
                                   k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(Imn[i]+d*Nmn[i]) +
                                   k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(Imp[i]+d*Nmp[i]) +
                                   TB_deaths[i]*Lmp[i]/tot_age[i];

    /* Smear negative,s,n */
                
    for (i=0; i<17; i++) dNsn[i] = age_in[i]*forc[i+1]*Nsn[i-1] - age_out[i]*Nsn[i] +
                                   (v*(1-sig_age[i]) + FS*a_age[i]*(1-p)*(1-sig_age[i]))*Lsn[i] +
                                   FS*a_age[i]*(1-sig_age[i])*(S[i] + (1-p)*Lmn[i]) -
                                   (theta + r + k*l_s*d + muN_age[i])*Nsn[i] +
                                   TB_deaths[i]*Nsn[i]/tot_age[i];

    /* Smear negative,s,p [ok] */
                    
    for (i=0; i<17; i++) dNsp[i] = age_in[i]*forc[i+1]*Nsp[i-1] - age_out[i]*Nsp[i] + 
                                   (v*(1-sig_age[i]) + FS*a_age[i]*(1-p)*(1-sig_age[i]))*Lsp[i] +
                                   FS*a_age[i]*(1-sig_age[i])*(1-p)*Lmp[i] -
                                   (theta + r + k*l_s*d*(1-e)*tau_s + k*l_s*d*e + muN_age[i])*Nsp[i] +
                                   k*l_s*d*(1-e)*(1-tau_s)*Nsn[i] +
                                   TB_deaths[i]*Nsp[i]/tot_age[i];

    /* Smear negative,m,n [ok] */

    for (i=0; i<17; i++) dNmn[i] = age_in[i]*forc[i+1]*Nmn[i-1] - age_out[i]*Nmn[i] +
                                   (v*(1-sig_age[i]) + FM*a_age[i]*(1-p)*(1-sig_age[i]))*Lmn[i] +
                                   FM*a_age[i]*(1-sig_age[i])*(S[i] + (1-p)*Lsn[i]) -
                                   (theta + r + k*l_m*d*dst_n + k*l_s*d*(1-dst_n) + muN_age[i])*Nmn[i] +
                                   TB_deaths[i]*Nmn[i]/tot_age[i];
    
    /* Smear negative,m,p [ok] */

    for (i=0; i<17; i++) dNmp[i] = age_in[i]*forc[i+1]*Nmp[i-1] - age_out[i]*Nmp[i] +
                                   (v*(1-sig_age[i]) + FM*a_age[i]*(1-p)*(1-sig_age[i]))*Lmp[i] +
                                   FM*a_age[i]*(1-sig_age[i])*(1-p)*Lsp[i] -
                                   (theta + r + k*l_m*d*dst_p*tau_m + k*l_s*d*(1-dst_p)*tau_s*eff_p + muN_age[i])*Nmp[i] +
                                   k*l_s*d*e*(Nsn[i]+Nsp[i]) +
                                   (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*Nmn[i] +
                                   TB_deaths[i]*Nmp[i]/tot_age[i]; 

    /* Smear positive,s,n [ok] */

    for (i=0; i<17; i++) dIsn[i] = age_in[i]*forc[i+1]*Isn[i-1] - age_out[i]*Isn[i] +
                                   (v*sig_age[i] + FS*a_age[i]*sig_age[i]*(1-p))*Lsn[i] +
                                   FS*a_age[i]*sig_age[i]*S[i] +
                                   FS*a_age[i]*(1-p)*sig_age[i]*Lmn[i] +
                                   theta*Nsn[i] -
                                   (r + k*l_s + muI_age[i])*Isn[i] +
                                   TB_deaths[i]*Isn[i]/tot_age[i];

    /* Smear positive,s,p [ok] */

    for (i=0; i<17; i++) dIsp[i] = age_in[i]*forc[i+1]*Isp[i-1] - age_out[i]*Isp[i] +
                                   (v*sig_age[i] + FS*a_age[i]*sig_age[i]*(1-p))*Lsp[i] +
                                   FS*a_age[i]*sig_age[i]*(1-p)*Lmp[i] +
                                   theta*Nsp[i] +
                                   k*l_s*(1-e)*(1-tau_s)*Isn[i] -
                                   (r + k*l_s*(1-e)*tau_s + k*l_s*e + muI_age[i])*Isp[i] +
                                   TB_deaths[i]*Isp[i]/tot_age[i];
                                   
    /* Smear positive,m,n [ok] */

    for (i=0; i<17; i++) dImn[i] = age_in[i]*forc[i+1]*Imn[i-1] - age_out[i]*Imn[i] +
                                   (v*sig_age[i] + FM*a_age[i]*sig_age[i]*(1-p))*Lmn[i] +
                                   FM*a_age[i]*sig_age[i]*S[i] +
                                   FM*a_age[i]*(1-p)*sig_age[i]*Lsn[i] +
                                   theta*Nmn[i] -
                                   (r + k*l_m*dst_n + k*l_s*(1-dst_n) + muI_age[i])*Imn[i] +
                                   TB_deaths[i]*Imn[i]/tot_age[i];
    
    /* Smear positive,m,n [OK] */

    for (i=0; i<17; i++) dImp[i] = age_in[i]*forc[i+1]*Imp[i-1] - age_out[i]*Imp[i] +
                                      (v*sig_age[i] + FM*a_age[i]*sig_age[i]*(1-p))*Lmp[i] +
                                      FM*a_age[i]*sig_age[i]*(1-p)*Lsp[i] +
                                      theta*Nmp[i] +
                                      (k*l_m*dst_n*(1-tau_m) + k*l_s*(1-dst_n)*(1-(tau_s*eff_n)))*Imn[i] +
                                      k*l_s*e*(Isn[i]+Isp[i]) -
                                      (r + k*l_m*dst_p*tau_m + k*l_s*(1-dst_p)*tau_s*eff_p + muI_age[i])*Imp[i] +
                                      TB_deaths[i]*Imp[i]/tot_age[i];
                
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
    
}

                                                                    
  




