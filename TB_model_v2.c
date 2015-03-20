/* Attempt to write TB model in C code to call from R */

/* compile within R with system("R CMD SHLIB Demog_model.c") */

#include <R.h>

static double parms[25];
static double forc[18];

/* A trick to keep up with the parameters and forcings */
#define age1 parms[0]
#define age2 parms[1]
#define beta parms[2]
#define a parms[3]
#define p parms[4]
#define v parms[5]
#define sig parms[6]
#define rel_inf parms[7]
#define theta parms[8]
#define r parms[9]
#define mu_N parms[10]
#define mu_I parms[11]
#define fit_cost parms[12]
#define e parms[13]
#define g parms[14]
#define k parms[15]
#define l_s parms[16]
#define l_m parms[17]
#define d parms[18]
#define tau_s parms[19]
#define tau_m parms[20]
#define eff_n parms[21]
#define eff_p parms[22]
#define dst_n parms[23]
#define dst_p parms[24]

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
    int N=25;
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
     
    int i=0;  
     
    for (i=0; i<17; i++) S[i] = y[i];            /* S: 0-16 */
     
    for (i=17; i<34; i++) Lsn[i-17] = y[i];      /* Lsn: 17-33 */
    for (i=34; i<51; i++) Lsp[i-34] = y[i];      /* Lsp: 34-50 */
    for (i=51; i<68; i++) Lmn[i-51] = y[i];      /* Lmn: 51-67 */
    for (i=68; i<85; i++) Lmp[i-68] = y[i];      /* Lmp: 68-84 */
    
    for (i=85; i<102; i++) Nsn[i-85] = y[i];      /* Nsn: 85-101 */
    for (i=102; i<119; i++) Nsp[i-102] = y[i];    /* Nsp: 102-118 */
    for (i=119; i<136; i++) Nmn[i-119] = y[i];    /* Nmn: 119-135 */
    for (i=136; i<153; i++) Nmp[i-136] = y[i];    /* Nmp: 136-152 */
       
    for (i=153; i<170; i++) Isn[i-153] = y[i];    /* Isn: 153-169 */
    for (i=170; i<187; i++) Isp[i-170] = y[i];    /* Isp: 170-186 */
    for (i=187; i<204; i++) Imn[i-187] = y[i];    /* Imn: 187-203 */
    for (i=204; i<221; i++) Imp[i-204] = y[i];    /* Imp: 204-220 */
       
    /* sum up various totals - uses function sum_array(array,i_start,i_end) */ 
    
    double Total_S = sumsum(S,0,16);
    
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
 
    /* Derivatives */ 
 
    /* Susceptible */ 
    
    dS[0] = s_birth*birth_rate*Total/1000 - age1*S[0] - (FS + FM)*S[0];
    for (i=1; i<16; i++) dS[i] = age1*forc[i+1]*S[i-1] - age1*S[i] - (FS + FM)*S[i];  
    dS[16] = age1*s80*S[15] - age2*S[16] - (FS + FM)*S[16];

    /* Latent,s,n*/
  
    dLsn[0] =  -age1*Lsn[0] 
               + FS*((1-a)*S[0] + (1-a)*(1-p)*(1-g)*Lmn[0])
               - (v + FS*a*(1-p) + FM*a*(1-p) + FM*(1-a)*(1-p)*g)*Lsn[0]
               + r*(Isn[0] + Nsn[0]);
    
    for (i=1; i<16; i++) dLsn[i] = age1*forc[i+1]*Lsn[i-1] - age1*Lsn[i]
                                   + FS*((1-a)*S[i] + (1-a)*(1-p)*(1-g)*Lmn[i])
                                   - (v + FS*a*(1-p) + FM*a*(1-p) + FM*(1-a)*(1-p)*g)*Lsn[i]
                                   + r*(Isn[i] + Nsn[i]); 
    
    dLsn[16] = age1*s80*Lsn[15] - age2*Lsn[16] 
               + FS*((1-a)*S[16] + (1-a)*(1-p)*(1-g)*Lmn[16])
               - (v + FS*a*(1-p) + FM*a*(1-p) + FM*(1-a)*(1-p)*g)*Lsn[16]
               + r*(Isn[16] + Nsn[16]);
    
    /* Latent,s,p [ok] */

    dLsp[0] =  - age1*Lsp[0] 
               + FS*(1-a)*(1-p)*(1-g)*Lmp[0]
               - (v + FS*a*(1-p) + FM*a*(1-p) + FM*(1-a)*(1-p)*g)*Lsp[0]
               + r*(Isp[0]+Nsp[0])
               + k*l_s*(1-e)*tau_s*(Isn[0]+Isp[0]) 
               + k*l_s*(1-e)*tau_s*d*(Nsn[0]+Nsp[0]);
  
    for (i=1; i<16; i++) dLsp[i] = age1*forc[i+1]*Lsp[i-1] - age1*Lsp[i] 
                                   + FS*(1-a)*(1-p)*(1-g)*Lmp[i]
                                   - (v + FS*a*(1-p) + FM*a*(1-p) + FM*(1-a)*(1-p)*g)*Lsp[i]
                                   + r*(Isp[i]+Nsp[i])
                                   + k*l_s*(1-e)*tau_s*(Isn[i]+Isp[i]) 
                                   + k*l_s*(1-e)*tau_s*d*(Nsn[i]+Nsp[i]);
    
    dLsp[16] = age1*s80*Lsp[15] - age2*Lsp[16] 
               + FS*(1-a)*(1-p)*(1-g)*Lmp[16]
               - (v + FS*a*(1-p) + FM*a*(1-p) + FM*(1-a)*(1-p)*g)*Lsp[16]
               + r*(Isp[16]+Nsp[16])
               + k*l_s*(1-e)*tau_s*(Isn[16]+Isp[16]) 
               + k*l_s*(1-e)*tau_s*d*(Nsn[16]+Nsp[16]);
    
    /* Latent,m,n [ok] */
  
    dLmn[0] =  - age1*Lmn[0] 
               + FM*((1-a)*S[0] + (1-a)*(1-p)*g*Lsn[0])
               - (v + FM*a*(1-p) + FS*a*(1-p) + FS*(1-a)*(1-p)*(1-g))*Lmn[0]
               + r*(Imn[0]+Nmn[0]);
                
    for (i=1; i<16; i++) dLmn[i] = age1*forc[i+1]*Lmn[i-1] - age1*Lmn[i] 
                                   + FM*((1-a)*S[i] + (1-a)*(1-p)*g*Lsn[i])
                                   - (v + FM*a*(1-p) + FS*a*(1-p) + FS*(1-a)*(1-p)*(1-g))*Lmn[i]
                                   + r*(Imn[i]+Nmn[i]);  
    
    dLmn[16] = age1*s80*Lmn[15] - age2*Lmn[16] 
               + FM*((1-a)*S[16] + (1-a)*(1-p)*g*Lsn[16])
               - (v + FM*a*(1-p) + FS*a*(1-p) + FS*(1-a)*(1-p)*(1-g))*Lmn[16]
               + r*(Imn[16]+Nmn[16]);

    /* Latent,m,p [ok] */
    
    dLmp[0] = - age1*Lmp[0] 
              + FM*(1-a)*(1-p)*g*Lsp[0]
              - (v + FM*a*(1-p) + FS*a*(1-p) + FS*(1-a)*(1-p)*(1-g))*Lmp[0]
              + r*(Imp[0]+Nmp[0])
              + k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(Imn[0]+d*Nmn[0])
              + k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(Imp[0]+d*Nmp[0]);
                
    for (i=1; i<16; i++) dLmp[i] = age1*forc[i+1]*Lmp[i-1] - age1*Lmp[i] 
                                   + FM*(1-a)*(1-p)*g*Lsp[i]
                                   - (v + FM*a*(1-p) + FS*a*(1-p) + FS*(1-a)*(1-p)*(1-g))*Lmp[i]
                                   + r*(Imp[i]+Nmp[i])
                                   + k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(Imn[i]+d*Nmn[i])
                                   + k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(Imp[i]+d*Nmp[i]);
    
    dLmp[16] = age1*s80*Lmp[15] - age2*Lmp[16] 
               + FM*(1-a)*(1-p)*g*Lsp[16]
               - (v + FM*a*(1-p) + FS*a*(1-p) + FS*(1-a)*(1-p)*(1-g))*Lmp[16]
               + r*(Imp[16]+Nmp[16])
               + k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(Imn[16]+d*Nmn[16])
               + k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(Imp[16]+d*Nmp[16]);


    /* Smear negative,s,n */

    dNsn[0] = - age1*Nsn[0] 
              + (v*(1-sig) + FS*a*(1-p)*(1-sig))*Lsn[0]
              + FS*a*(1-sig)*(S[0] + (1-p)*Lmn[0])
              - (theta + r + k*l_s*d + mu_N)*Nsn[0];
                
    for (i=1; i<16; i++) dNsn[i] = age1*forc[i+1]*Nsn[i-1] - age1*Nsn[i] 
                                   + (v*(1-sig) + FS*a*(1-p)*(1-sig))*Lsn[i]
                                   + FS*a*(1-sig)*(S[i] + (1-p)*Lmn[i])
                                   - (theta + r + k*l_s*d + mu_N)*Nsn[i];
    
    dNsn[16] = age1*s80*Nsn[15] - age2*Nsn[16] 
               + (v*(1-sig) + FS*a*(1-p)*(1-sig))*Lsn[16]
               + FS*a*(1-sig)*(S[16] + (1-p)*Lmn[16])
               - (theta + r + k*l_s*d + mu_N)*Nsn[16];

    /* Smear negative,s,p [ok] */

    dNsp[0] = - age1*Nsp[0] 
              + (v*(1-sig) + FS*a*(1-p)*(1-sig))*Lsp[0]
              + FS*a*(1-sig)*(1-p)*Lmp[0]
              - (theta + r + k*l_s*d*(1-e)*tau_s + k*l_s*d*e + mu_N)*Nsp[0]
              + k*l_s*d*(1-e)*(1-tau_s)*Nsn[0];
    
                
    for (i=1; i<16; i++) dNsp[i] = age1*forc[i+1]*Nsp[i-1] - age1*Nsp[i] 
                                   + (v*(1-sig) + FS*a*(1-p)*(1-sig))*Lsp[i]
                                   + FS*a*(1-sig)*(1-p)*Lmp[i]
                                   - (theta + r + k*l_s*d*(1-e)*tau_s + k*l_s*d*e + mu_N)*Nsp[i]
                                   + k*l_s*d*(1-e)*(1-tau_s)*Nsn[i];
    
    dNsp[16] = age1*s80*Nsp[15] - age2*Nsp[16] 
               + (v*(1-sig) + FS*a*(1-p)*(1-sig))*Lsp[16]
               + FS*a*(1-sig)*(1-p)*Lmp[16]
               - (theta + r + k*l_s*d*(1-e)*tau_s + k*l_s*d*e + mu_N)*Nsp[16]
               + k*l_s*d*(1-e)*(1-tau_s)*Nsn[16];

    /* Smear negative,m,n [ok] */

    dNmn[0] = - age1*Nmn[0] 
              + (v*(1-sig) + FM*a*(1-p)*(1-sig))*Lmn[0]
              + FM*a*(1-sig)*(S[0] + (1-p)*Lsn[0])
              - (theta + r + k*l_m*d*dst_n + k*l_s*d*(1-dst_n) + mu_N)*Nmn[0];

    for (i=1; i<16; i++) dNmn[i] = age1*forc[i+1]*Nmn[i-1] - age1*Nmn[i] 
                                   + (v*(1-sig) + FM*a*(1-p)*(1-sig))*Lmn[i]
                                   + FM*a*(1-sig)*(S[i] + (1-p)*Lsn[i])
                                   - (theta + r + k*l_m*d*dst_n + k*l_s*d*(1-dst_n) + mu_N)*Nmn[i];
    
    dNmn[16] = age1*s80*Nmn[15] - age2*Nmn[16] 
               + (v*(1-sig) + FM*a*(1-p)*(1-sig))*Lmn[16]
               + FM*a*(1-sig)*(S[16] + (1-p)*Lsn[16])
               - (theta + r + k*l_m*d*dst_n + k*l_s*d*(1-dst_n) + mu_N)*Nmn[16];
    
    /* Smear negative,m,p [ok] */
    
    dNmp[0] = - age1*Nmp[0]
              + (v*(1-sig) + FM*a*(1-p)*(1-sig))*Lmp[0]
              + FM*a*(1-sig)*(1-p)*Lsp[0]
              - (theta + r + k*l_m*d*dst_p*tau_m + k*l_s*d*(1-dst_p)*tau_s*eff_p + mu_N)*Nmp[0]
              + k*l_s*d*e*(Nsn[0]+Nsp[0])
              + (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*Nmn[0]; 

    for (i=1; i<16; i++) dNmp[i] = age1*forc[i+1]*Nmp[i-1] - age1*Nmp[i] 
                                   + (v*(1-sig) + FM*a*(1-p)*(1-sig))*Lmp[i]
                                   + FM*a*(1-sig)*(1-p)*Lsp[i]
                                   - (theta + r + k*l_m*d*dst_p*tau_m + k*l_s*d*(1-dst_p)*tau_s*eff_p + mu_N)*Nmp[i]
                                   + k*l_s*d*e*(Nsn[i]+Nsp[i])
                                   + (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*Nmn[i]; 

    dNmp[16] = age1*s80*Nmp[15] - age2*Nmp[16] 
               + (v*(1-sig) + FM*a*(1-p)*(1-sig))*Lmp[16]
               + FM*a*(1-sig)*(1-p)*Lsp[16]
               - (theta + r + k*l_m*d*dst_p*tau_m + k*l_s*d*(1-dst_p)*tau_s*eff_p + mu_N)*Nmp[16]
               + k*l_s*d*e*(Nsn[16]+Nsp[16])
               + (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*Nmn[16]; 

    /* Smear positive,s,n [ok] */
    
    dIsn[0] = - age1*Isn[0]
              + (v*sig + FS*a*sig*(1-p))*Lsn[0]
              + FS*a*sig*S[0]
              + FS*a*(1-p)*sig*Lmn[0]
              + theta*Nsn[0]
              - (r + k*l_s +mu_I)*Isn[0];

    for (i=1; i<16; i++) dIsn[i] = age1*forc[i+1]*Isn[i-1] - age1*Isn[i] 
                                   + (v*sig + FS*a*sig*(1-p))*Lsn[i]
                                   + FS*a*sig*S[i]
                                   + FS*a*(1-p)*sig*Lmn[i]
                                   + theta*Nsn[i]
                                   - (r + k*l_s +mu_I)*Isn[i];

    dIsn[16] = age1*s80*Isn[15] - age2*Isn[16] 
               + (v*sig + FS*a*sig*(1-p))*Lsn[16]
               + FS*a*sig*S[16]
               + FS*a*(1-p)*sig*Lmn[16]
               + theta*Nsn[16]
               - (r + k*l_s +mu_I)*Isn[16];

    /* Smear positive,s,p [ok] */

    dIsp[0] = - age1*Isp[0]
              + (v*sig + FS*a*sig*(1-p))*Lsp[0]
              + FS*a*sig*(1-p)*Lmp[0]
              + theta*Nsp[0]
              + k*l_s*(1-e)*(1-tau_s)*Isn[0]
              - (r + k*l_s*(1-e)*tau_s + k*l_s*e + mu_I)*Isp[0];

    for (i=1; i<16; i++) dIsp[i] = age1*forc[i+1]*Isp[i-1] - age1*Isp[i] 
                                   + (v*sig + FS*a*sig*(1-p))*Lsp[i]
                                   + FS*a*sig*(1-p)*Lmp[i]
                                   + theta*Nsp[i]
                                   + k*l_s*(1-e)*(1-tau_s)*Isn[i]
                                   - (r + k*l_s*(1-e)*tau_s + k*l_s*e + mu_I)*Isp[i];

    dIsp[16] = age1*s80*dIsp[15] - age2*Isp[16] 
               + (v*sig + FS*a*sig*(1-p))*Lsp[16]
               + FS*a*sig*(1-p)*Lmp[16]
               + theta*Nsp[16]
               + k*l_s*(1-e)*(1-tau_s)*Isn[16]
               - (r + k*l_s*(1-e)*tau_s + k*l_s*e + mu_I)*Isp[16];

    /* Smear positive,m,n [ok] */
    
    dImn[0] = - age1*Imn[0]
              + (v*sig + FM*a*sig*(1-p))*Lmn[0]
              + FM*a*sig*S[0]
              + FM*a*(1-p)*sig*Lsn[0]
              + theta*Nmn[0]
              - (r + k*l_m*dst_n + k*l_s*(1-dst_n) + mu_I)*Imn[0];

    for (i=1; i<16; i++) dImn[i] = age1*forc[i+1]*Imn[i-1] - age1*Imn[i] 
                                   + (v*sig + FM*a*sig*(1-p))*Lmn[i]
                                   + FM*a*sig*S[i]
                                   + FM*a*(1-p)*sig*Lsn[i]
                                   + theta*Nmn[i]
                                   - (r + k*l_m*dst_n + k*l_s*(1-dst_n) + mu_I)*Imn[i];

    dImn[16] = age1*s80*Imn[15] - age2*Imn[16] 
                + (v*sig + FM*a*sig*(1-p))*Lmn[16]
                + FM*a*sig*S[16]
                + FM*a*(1-p)*sig*Lsn[16]
                + theta*Nmn[16]
                - (r + k*l_m*dst_n + k*l_s*(1-dst_n) + mu_I)*Imn[16];
    
    /* Smear positive,m,n [OK] */
       
    dImp[0] = - age1*Imp[0]
              + (v*sig + FM*a*sig*(1-p))*Lmp[0]
              + FM*a*sig*(1-p)*Lsp[0]
              + theta*Nmp[0]
              + (k*l_m*dst_n*(1-tau_m) + k*l_s*(1-dst_n)*(1-(tau_s*eff_n)))*Imn[0]
              + k*l_s*e*(Isn[0]+Isp[0])
              - (r + k*l_m*dst_p*tau_m + k*l_s*(1-dst_p)*tau_s*eff_p + mu_I)*Imp[0];

    for (i=1; i<16; i++) dImp[i] = age1*forc[i+1]*Imp[i-1] - age1*Imp[i] 
                                      + (v*sig + FM*a*sig*(1-p))*Lmp[i]
                                      + FM*a*sig*(1-p)*Lsp[i]
                                      + theta*Nmp[i]
                                      + (k*l_m*dst_n*(1-tau_m) + k*l_s*(1-dst_n)*(1-(tau_s*eff_n)))*Imn[i]
                                      + k*l_s*e*(Isn[i]+Isp[i])
                                      - (r + k*l_m*dst_p*tau_m + k*l_s*(1-dst_p)*tau_s*eff_p + mu_I)*Imp[i];

    dImp[16] = age1*s80*Imp[15] - age2*Imp[16] 
                + (v*sig + FM*a*sig*(1-p))*Lmp[16]
                + FM*a*sig*(1-p)*Lsp[16]
                + theta*Nmp[16]
                + (k*l_m*dst_n*(1-tau_m) + k*l_s*(1-dst_n)*(1-(tau_s*eff_n)))*Imn[16]
                + k*l_s*e*(Isn[16]+Isp[16])
                - (r + k*l_m*dst_p*tau_m + k*l_s*(1-dst_p)*tau_s*eff_p + mu_I)*Imp[16];
                
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

                                                                    
  





