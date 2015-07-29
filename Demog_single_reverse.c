/* Demog model - testing with single year age bins */

/* C libraries needed */
#include <R.h>
#include <math.h>

/* You need to define number of parameters and forcing functions passed to the model here */
/* These must match number in intializer functions below */
static double parms[2];
static double forc[18];

/* ###### A TRICK TO KEEP UP WITH THE PARAMETERS AND FORCINGS ###### */

/* !!!!!!! NOTE IN C INDEXING STARTS AT 0 AND NOT 1 (i.e. parms[0] is the first value in parms) !!!!!! */

/* Parameters are constant over time */
#define d parms[0]       /* rate of entry/exit to/from age group = 1/width */
#define d2 parms[1]

/* Forcings are time dependant functions */
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
    int N=18;
    odeforcs(&N, forc);
}

/* ##### EVENTS ARE USED TO ADD BIRTHS AND SHIFT THE POPULATION BY AGE - EQUIVALENT TO THE METHOD OF SCHENZLE ###### */
void event(int *n, double *t, double *y) 
{
  int i;
  double temp[17];
  for (i=0; i<17; i++) temp[i] = y[i];
  y[0] = (4.0/5.0)*y[0]+birth_rate*sumsum(temp,0,16)/1000;
  for (i=1; i<16; i++) y[i] = (4.0/5.0)*temp[i] + (1.0/5.0)*temp[i-1];
  y[16] = y[16] + (1/5)*temp[15];
}


/* ###### DERIVATIVE FUNCTION ###### */

void derivsc(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if (ip[0] <2) error("nout should be at least 2");
    
    /* Expand state variables so can use more meaningful names than y and ydot (the variables and rates of change are vectors y and ydot here) */
    /* There are 17 age groups, 7 HIV positive categories and 3 times on ART */
    /* _H = HIV+; _A = HIV+ on ART */
    
    /* These are the variables */
    double S[17];       
   
    /* These are the rates of change (same names but prefixed with d) */
    double dS[17];
   
    /* intergers to use as counters */ 
    int i;  

    /* map y to S */
    for (i=0; i<17; i++) S[i] = y[i];            

    /* Use sumsum function to add up population */
    double Total = sumsum(S,0,16);                

    /* Assign mortality rates */
    double m_b[17];
    for (i=0; i<17; i++){
      m_b[i] = forc[i+1];
    }
    
    /* Derivatives */ 

    for (i=0; i<17; i++) dS[i] = - m_b[i]*S[i];
    
    /* Calculate total births */
    double births;
    births = birth_rate*Total/1000;
    
    /* Calculate total deaths */
    double Tot_deaths = 0;
    for (i=0; i<17; i++){
      Tot_deaths =  Tot_deaths + m_b[i]*S[i];
    }
    
    /* map ds to ydot */

    for (i=0; i<17; i++) ydot[i] = dS[i];             

    /* Finally assign the things we want to use in R (in addition to the state variables) to yout */
    /* The ode call in R needs to define the number of these and give names */
    yout[0] = Total;
    yout[1] = births;
    yout[2] = Tot_deaths;
    
}
