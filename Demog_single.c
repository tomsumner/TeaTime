/* Demog model - testing with single year age bins */

/* C libraries needed */
#include <R.h>
#include <math.h>

/* You need to define number of parameters and forcing functions passed to the model here */
/* These must match number in intializer functions below */
static double parms[2];
static double forc[82];

/* ###### A TRICK TO KEEP UP WITH THE PARAMETERS AND FORCINGS ###### */

/* !!!!!!! NOTE IN C INDEXING STARTS AT 0 AND NOT 1 (i.e. parms[0] is the first value in parms) !!!!!! */

/* Parameters are constant over time */
#define age_1 parms[0]       /* rate of entry/exit to/from age group = 1/width */
#define age_2 parms[1]

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
    int N=2;
    odeparms(&N, parms);
}

/* ###### FUNCTION TO INITIALIZE FORCINGS PASSED FROM R - if the number of parameters is changed you must update N here ###### */
void forcc(void (* odeforcs)(int *, double *))
{
    int N=82;
    odeforcs(&N, forc);
}

/* ##### EVENTS ARE USED TO ADD BIRTHS AND SHIFT THE POPULATION BY AGE - EQUIVALENT TO THE METHOD OF SCHENZLE ###### */
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

void derivsc(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if (ip[0] <2) error("nout should be at least 2");
    
    /* Expand state variables so can use more meaningful names than y and ydot (the variables and rates of change are vectors y and ydot here) */
    /* There are 17 age groups, 7 HIV positive categories and 3 times on ART */
    /* _H = HIV+; _A = HIV+ on ART */
    
    /* These are the variables */
    double S[81];       
   
    /* These are the rates of change (same names but prefixed with d) */
    double dS[81];
   
    /* intergers to use as counters */ 
    int i;  

    /* map y to S */
    for (i=0; i<81; i++) S[i] = y[i];            

    /* Use sumsum function to add up population */
    double Total = sumsum(S,0,80);                

    /* Assign mortality rates */
    double m_b[81];
    for (i=0; i<81; i++){
      m_b[i] = forc[i+1];
    }
    
    /* Derivatives */ 

    for (i=0; i<81; i++) dS[i] = - m_b[i]*S[i];
    
    /* Calculate total births */
    double births;
    births = birth_rate*Total/1000;
    
    /* Calculate total deaths */
    double Tot_deaths = 0;
    for (i=0; i<81; i++){
      Tot_deaths =  Tot_deaths + m_b[i]*S[i];
    }
    
    /* map ds to ydot */

    for (i=0; i<81; i++) ydot[i] = dS[i];             

    /* Finally assign the things we want to use in R (in addition to the state variables) to yout */
    /* The ode call in R needs to define the number of these and give names */
    yout[0] = Total;
    yout[1] = births;
    yout[2] = Tot_deaths;
    
}



