/* Attempt to write TB model in C code to call from R

/* compile within R with system("R CMD SHLIB Demog_model.c") */

#include <R.h>

static double parms[2];
static double forc[18];

/* A trick to keep up with the parameters and forcings */
/* progression through age groups */
#define d parms[0]
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

/* initializers*/
void parmsc(void (* odeparms)(int *, double *))
{
    int N=2;
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

    int i=0; 
    double Total = sumsum(y,0,16);

    /* Calculate the number of births */
    double births=birth_rate*Total/1000;

    double d_age[17] = {d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d2};

    ydot[0] = births - d_age[0]*y[0] - forc[i+1]*y[0];

    for (i=1; i<17; i++) ydot[i] = d_age[i-1]*y[i-1] - d_age[i]*y[i] - forc[i+1]*y[i];

    /* Calculate total deaths */
    double Tot_deaths = 0;
    for (i=0; i<17; i++){
      Tot_deaths =  Tot_deaths + forc[i+1]*y[i];
    }
    
    yout[0] = Total;
    yout[1] = births;
    yout[2] = Tot_deaths;

}



