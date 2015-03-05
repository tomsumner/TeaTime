/* Attempt to write TB model in C code to call from R

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
    
    /* S:     0-16
       L_s_n: 17-33
       L_s_p: 34-50
       L_m_n: 51-67
       L_m_p: 68-84
       
       N_s_n: 85-101
       N_s_p: 102-118
       N_m_n: 119-135
       N_m_p: 136-152
       
       I_s_n: 153-169
       I_s_p: 170-186
       I_m_n: 187-203
       I_m_p: 204-220*/
       
    /* sum up various totals - uses function sum_array(array,i_start,i_end) */ 
    double Total_S = sumsum(y,0,16);
    double Total_L = sumsum(y,17,33);
    
    double Total_Ns = sumsum(y,85,118);
    double Total_Nm = sumsum(y,119,152);
    double Total_N = Total_Ns + Total_Nm;
    double Total_Is = sumsum(y,153,186);
    double Total_Im = sumsum(y,187,220);
    double Total_I = Total_Is + Total_Im;
 
    double Total_DS = Total_Ns + Total_Is;
    double Total_MDR = Total_Nm + Total_Im;
 
    double Total = sumsum(y,0,220);
 
    /* Force of infection */
    double FS = beta*(Total_Ns*rel_inf + Total_Is)/Total; 
    double FM = fit_cost*beta*(Total_Nm*rel_inf + Total_Im)/Total;
 
    /*Susceptible: 0-16 */ 
    int i=0; 
 
    ydot[0] = s_birth*birth_rate*Total/1000 - age1*y[0] - (FS + FM)*y[0];
    for (i=1; i<16; i++) ydot[i] = age1*forc[i+1]*y[i-1] - age1*y[i] - (FS + FM)*y[i];  
    ydot[16] = age1*s80*y[15] - age2*y[16] - (FS + FM)*y[16];

    /*Latent,s,n: 17-33*/
    
    ydot[17] =  -age1*y[17] 
                + FS*((1-a)*y[0] + (1-a)*(1-p)*(1-g)*y[51])
                - (v + FS*a*(1-p) + FM*a*(1-p) + FM*(1-a)*(1-p)*g)*y[17]
                + r*(y[153] + y[85]);
    
    for (i=18; i<33; i++) ydot[i] = age1*forc[i-16]*y[i-1] - age1*y[i]
                                    + FS*((1-a)*y[i-17] + (1-a)*(1-p)*(1-g)*y[i])
                                    - (v + FS*a*(1-p) + FM*a*(1-p) + FM*(1-a)*(1-p)*g)*y[i]
                                    + r*(y[i+136] + y[i+68]); 
    
    ydot[33] = age1*s80*y[32] - age2*y[33] 
               + FS*((1-a)*y[16] + (1-a)*(1-p)*(1-g)*y[67])
               - (v + FS*a*(1-p) + FM*a*(1-p) + FM*(1-a)*(1-p)*g)*y[33]
               + r*(y[169] + y[101]);
    
    /*Latent,s,p: 34-50*/
    
    ydot[34] =  - age1*y[34] 
                + FS*(1-a)*(1-p)*(1-g)*y[68]
                - (v + FS*a*(1-p) + FM*a*(1-p) + FM*(1-a)*(1-p)*g)*y[34]
                + r*(y[170]+y[102])
                + k*l_s*(1-e)*tau_s*(y[170]+y[153]) 
                + k*l_s*(1-e)*tau_s*d*(y[102]+y[85]);
  
    for (i=35; i<50; i++) ydot[i] = age1*forc[i-33]*y[i-1] - age1*y[i] 
                                    + FS*(1-a)*(1-p)*(1-g)*y[i+34]
                                    - (v + FS*a*(1-p) + FM*a*(1-p) + FM*(1-a)*(1-p)*g)*y[i]
                                    + r*(y[i+136]+y[i+68])
                                    + k*l_s*(1-e)*tau_s*(y[i+136]+y[i+119]) 
                                    + k*l_s*(1-e)*tau_s*d*(y[i+68]+y[i+51]);
    
    ydot[50] = age1*s80*y[49] - age2*y[50] 
               + FS*(1-a)*(1-p)*(1-g)*y[84]
               - (v + FS*a*(1-p) + FM*a*(1-p) + FM*(1-a)*(1-p)*g)*y[50]
               + r*(y[186]+y[118])
               + k*l_s*(1-e)*tau_s*(y[186]+y[169]) 
               + k*l_s*(1-e)*tau_s*d*(y[118]+y[101]);
    
    /*Latent,m,n: 51-67*/
  
    ydot[51] =  - age1*y[51] 
                + FM*((1-a)*y[0] + (1-a)*(1-p)*g*y[17])
                - (v + FM*a*(1-p) + FS*a*(1-p) + FS*(1-a)*(1-p)*(1-g))*y[51]
                + r*(y[187]+y[119]);
                
    for (i=52; i<67; i++) ydot[i] = age1*forc[i-50]*y[i-1] - age1*y[i] 
                                    + FM*((1-a)*y[i-51] + (1-a)*(1-p)*g*y[i-34])
                                    - (v + FM*a*(1-p) + FS*a*(1-p) + FS*(1-a)*(1-p)*(1-g))*y[i]
                                    + r*(y[i+136]+y[i+68]);  
    
    ydot[67] = age1*s80*y[66] - age2*y[67] 
               + FM*((1-a)*y[16] + (1-a)*(1-p)*g*y[33])
               - (v + FM*a*(1-p) + FS*a*(1-p) + FS*(1-a)*(1-p)*(1-g))*y[67]
               + r*(y[203]+y[135]);

    
    /*Latent,m,p: 68-84*/

    
    ydot[68] = - age1*y[68] 
               + FM*(1-a)*(1-p)*g*y[34]
               - (v + FM*a*(1-p) + FS*a*(1-p) + FS*(1-a)*(1-p)*(1-g))*y[68]
               + r*(y[204]+y[136])
               + k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(y[187]+d*y[119])
               + k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(y[204]+d*y[136]);
                
    for (i=69; i<84; i++) ydot[i] = age1*forc[i-67]*y[i-1] - age1*y[i] 
                                    + FM*(1-a)*(1-p)*g*y[i-34]
                                    - (v + FM*a*(1-p) + FS*a*(1-p) + FS*(1-a)*(1-p)*(1-g))*y[i]
                                    + r*(y[i+136]+y[i+68])
                                    + k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(y[i+119]+d*y[i+51])
                                    + k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(y[i+136]+d*y[i+136]);
    
    ydot[84] = age1*s80*y[83] - age2*y[84] 
               + FM*(1-a)*(1-p)*g*y[50]
               - (v + FM*a*(1-p) + FS*a*(1-p) + FS*(1-a)*(1-p)*(1-g))*y[84]
               + r*(y[220]+y[152])
               + k*(l_m*dst_n*tau_m + l_s*(1-dst_n)*tau_s*eff_n)*(y[203]+d*y[135])
               + k*(l_m*dst_p*tau_m + l_s*(1-dst_p)*tau_s*eff_p)*(y[220]+d*y[152]);


    /*Smear negative,s,n: 85-101*/

    ydot[85] = - age1*y[85] 
               + (v*(1-sig) + FS*a*(1-p)*(1-sig))*y[17]
               + FS*a*(1-sig)*(y[0] + (1-p)*y[51])
               - (theta + r + k*l_s*d + mu_N)*y[85];
                
    for (i=86; i<101; i++) ydot[i] = age1*forc[i-84]*y[i-1] - age1*y[i] 
                                     + (v*(1-sig) + FS*a*(1-p)*(1-sig))*y[i-68]
                                     + FS*a*(1-sig)*(y[i-85] + (1-p)*y[i-34])
                                     - (theta + r + k*l_s*d + mu_N)*y[i];
    
    ydot[101] = age1*s80*y[100] - age2*y[101] 
                + (v*(1-sig) + FS*a*(1-p)*(1-sig))*y[33]
                + FS*a*(1-sig)*(y[16] + (1-p)*y[67])
                - (theta + r + k*l_s*d + mu_N)*y[101];


    /*Smear negative,s,p: 102-118*/

    ydot[102] = - age1*y[102] 
                + (v*(1-sig) + FS*a*(1-p)*(1-sig))*y[34]
                + FS*a*(1-sig)*(1-p)*y[68]
                - (theta + r + k*l_s*d*(1-e)*tau_s + k*l_s*d*e + mu_N)*y[102]
                + k*l_s*d*(1-e)*(1-tau_s)*y[85];
    
                
    for (i=103; i<118; i++) ydot[i] = age1*forc[i-101]*y[i-1] - age1*y[i] 
                                      + (v*(1-sig) + FS*a*(1-p)*(1-sig))*y[i-68]
                                      + FS*a*(1-sig)*(1-p)*y[i-34]
                                      - (theta + r + k*l_s*d*(1-e)*tau_s + k*l_s*d*e + mu_N)*y[i]
                                      + k*l_s*d*(1-e)*(1-tau_s)*y[i-17];
    
    ydot[118] = age1*s80*y[117] - age2*y[118] 
                + (v*(1-sig) + FS*a*(1-p)*(1-sig))*y[50]
                + FS*a*(1-sig)*(1-p)*y[84]
                - (theta + r + k*l_s*d*(1-e)*tau_s + k*l_s*d*e + mu_N)*y[118]
                + k*l_s*d*(1-e)*(1-tau_s)*y[101];

    /*Smear negative,m,n: 119-135*/

    ydot[119] = - age1*y[119] 
                + (v*(1-sig) + FM*a*(1-p)*(1-sig))*y[51]
                + FM*a*(1-sig)*(y[0] + (1-p)*y[17])
                - (theta + r + k*l_m*d*dst_n + k*l_s*d*(1-dst_n) + mu_N)*y[119];

    for (i=120; i<135; i++) ydot[i] = age1*forc[i-118]*y[i-1] - age1*y[i] 
                                      + (v*(1-sig) + FM*a*(1-p)*(1-sig))*y[i-68]
                                      + FM*a*(1-sig)*(y[0] + (1-p)*y[i-102])
                                      - (theta + r + k*l_m*d*dst_n + k*l_s*d*(1-dst_n) + mu_N)*y[i];
    
    ydot[135] = age1*s80*y[134] - age2*y[135] 
                + (v*(1-sig) + FM*a*(1-p)*(1-sig))*y[67]
                + FM*a*(1-sig)*(y[0] + (1-p)*y[33])
                - (theta + r + k*l_m*d*dst_n + k*l_s*d*(1-dst_n) + mu_N)*y[135];
    
    /*Smear negative,m,p: 136-152*/
    
    ydot[136] = - age1*y[136]
                + (v*(1-sig) + FM*a*(1-p)*(1-sig))*y[68]
                + FM*a*(1-sig)*(1-p)*y[34]
                - (theta + r + k*l_m*d*dst_p*tau_m + k*l_s*d*(1-dst_p)*tau_s*eff_p + mu_N)*y[136]
                + k*l_s*d*e*(y[85]+y[102])
                + (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*y[119]; 

    for (i=137; i<152; i++) ydot[i] = age1*forc[i-135]*y[i-1] - age1*y[i] 
                                      + (v*(1-sig) + FM*a*(1-p)*(1-sig))*y[i-68]
                                      + FM*a*(1-sig)*(1-p)*y[i-102]
                                      - (theta + r + k*l_m*d*dst_p*tau_m + k*l_s*d*(1-dst_p)*tau_s*eff_p + mu_N)*y[i]
                                      + k*l_s*d*e*(y[i-51]+y[i-34])
                                      + (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*y[i-17]; 

    ydot[152] = age1*s80*y[151] - age2*y[152] 
                + (v*(1-sig) + FM*a*(1-p)*(1-sig))*y[84]
                + FM*a*(1-sig)*(1-p)*y[50]
                - (theta + r + k*l_m*d*dst_p*tau_m + k*l_s*d*(1-dst_p)*tau_s*eff_p + mu_N)*y[152]
                + k*l_s*d*e*(y[101]+y[118])
                + (k*l_m*d*dst_n*(1-tau_m) + k*l_s*d*(1-dst_n)*(1-(tau_s*eff_n)))*y[135]; 

    /*Smear positive,s,n: 153-169*/
    
    ydot[153] = - age1*y[153]
                + (v*sig + FS*a*sig*(1-p))*y[17]
                + FS*a*sig*y[0]
                + FS*a*(1-p)*sig*y[51]
                + theta*y[85]
                - (r + k*l_s +mu_I)*y[153];

    for (i=154; i<169; i++) ydot[i] = age1*forc[i-152]*y[i-1] - age1*y[i] 
                                    + (v*sig + FS*a*sig*(1-p))*y[i-136]
                                    + FS*a*sig*y[i-153]
                                    + FS*a*(1-p)*sig*y[i-102]
                                    + theta*y[i-68]
                                    - (r + k*l_s +mu_I)*y[i];

    ydot[169] = age1*s80*y[168] - age2*y[169] 
                + (v*sig + FS*a*sig*(1-p))*y[33]
                + FS*a*sig*y[16]
                + FS*a*(1-p)*sig*y[67]
                + theta*y[101]
                - (r + k*l_s +mu_I)*y[169];

    /*Smear positive,s,p: 170-186*/

    ydot[170] = - age1*y[170]
                + (v*sig + FS*a*sig*(1-p))*y[34]
                + FS*a*sig*(1-p)*y[68]
                + theta*y[102]
                + k*l_s*(1-e)*(1-tau_s)*y[153]
                - (r + k*l_s*(1-e)*tau_s + k*l_s*e + mu_I)*y[170];

    for (i=171; i<186; i++) ydot[i] = age1*forc[i-169]*y[i-1] - age1*y[i] 
                                      + (v*sig + FS*a*sig*(1-p))*y[i-136]
                                      + FS*a*sig*(1-p)*y[i-102]
                                      + theta*y[i-68]
                                      + k*l_s*(1-e)*(1-tau_s)*y[i-17]
                                      - (r + k*l_s*(1-e)*tau_s + k*l_s*e + mu_I)*y[i];

    ydot[186] = age1*s80*y[185] - age2*y[186] 
                + (v*sig + FS*a*sig*(1-p))*y[50]
                + FS*a*sig*(1-p)*y[84]
                + theta*y[118]
                + k*l_s*(1-e)*(1-tau_s)*y[169]
                - (r + k*l_s*(1-e)*tau_s + k*l_s*e + mu_I)*y[186];

    /*Smear positive,m,n: 187-203*/
    
    ydot[187] = - age1*y[187]
                + (v*sig + FM*a*sig*(1-p))*y[51]
                + FM*a*sig*y[0]
                + FM*a*(1-p)*sig*y[17]
                + theta*y[119]
                - (r + k*l_m*dst_n + k*l_s*(1-dst_n) + mu_I)*y[187];

    for (i=188; i<203; i++) ydot[i] = age1*forc[i-186]*y[i-1] - age1*y[i] 
                                      + (v*sig + FM*a*sig*(1-p))*y[i-136]
                                      + FM*a*sig*y[i-187]
                                      + FM*a*(1-p)*sig*y[i-170]
                                      + theta*y[i-68]
                                      - (r + k*l_m*dst_n + k*l_s*(1-dst_n) + mu_I)*y[i];

    ydot[203] = age1*s80*y[202] - age2*y[203] 
                + (v*sig + FM*a*sig*(1-p))*y[67]
                + FM*a*sig*y[16]
                + FM*a*(1-p)*sig*y[33]
                + theta*y[135]
                - (r + k*l_m*dst_n + k*l_s*(1-dst_n) + mu_I)*y[203];
    
    /*Smear positive,m,n: 204-220*/
       
    ydot[204] = - age1*y[204]
                + (v*sig + FM*a*sig*(1-p))*y[68]
                + FM*a*sig*(1-p)*y[34]
                + theta*y[136]
                + (k*l_m*dst_n*(1-tau_m) + k*l_s*(1-dst_n)*(1-(tau_s*eff_n)))*y[187]
                + k*l_s*e*(y[153]+y[170])
                - (r + k*l_m*dst_p*tau_m + k*l_s*(1-dst_p)*tau_s*eff_p + mu_I)*y[204];


    for (i=205; i<220; i++) ydot[i] = age1*forc[i-203]*y[i-1] - age1*y[i] 
                                      + (v*sig + FM*a*sig*(1-p))*y[i-136]
                                      + FM*a*sig*(1-p)*y[i-170]
                                      + theta*y[i-68]
                                      + (k*l_m*dst_n*(1-tau_m) + k*l_s*(1-dst_n)*(1-(tau_s*eff_n)))*y[i-17]
                                      + k*l_s*e*(y[i-51]+y[i-34])
                                      - (r + k*l_m*dst_p*tau_m + k*l_s*(1-dst_p)*tau_s*eff_p + mu_I)*y[i];

    ydot[220] = age1*s80*y[219] - age2*y[220] 
                + (v*sig + FM*a*sig*(1-p))*y[84]
                + FM*a*sig*(1-p)*y[50]
                + theta*y[152]
                + (k*l_m*dst_n*(1-tau_m) + k*l_s*(1-dst_n)*(1-(tau_s*eff_n)))*y[203]
                + k*l_s*e*(y[169]+y[186])
                - (r + k*l_m*dst_p*tau_m + k*l_s*(1-dst_p)*tau_s*eff_p + mu_I)*y[220];

    yout[0] = Total;
    yout[1] = Total_S;
    yout[2] = Total_L;
    yout[3] = Total_Ns;
    yout[4] = Total_Nm;
    yout[5] = Total_N;
    yout[6] = Total_Is;
    yout[7] = Total_Im;
    yout[8] = Total_I;
    yout[9] = Total_DS;
    yout[10] = Total_MDR;
    yout[11] = FS;
    yout[12] = FM;
    
}



                                                                    
  





