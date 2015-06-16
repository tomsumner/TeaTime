/* Attempt to write TB model in C code to call from R

/* compile within R with system("R CMD SHLIB Demog_model.c") */

#include <R.h>

static double parms[2];
static double forc[42];

/* A trick to keep up with the parameters and forcings */
/* progression through age groups */
#define d parms[0]
#define d2 parms[1]

#define b20 forc[0]
#define b25 forc[1]
#define b30 forc[2]
#define b35 forc[3]
#define b40 forc[4]
#define b45 forc[5]
#define b50 forc[6]

#define s5f forc[7]
#define s10f forc[8]
#define s15f forc[9]
#define s20f forc[10]
#define s25f forc[11]
#define s30f forc[12]
#define s35f forc[13]
#define s40f forc[14]
#define s45f forc[15]
#define s50f forc[16]
#define s55f forc[17]
#define s60f forc[18]
#define s65f forc[19]
#define s70f forc[20]
#define s75f forc[21]
#define s80f forc[22]
#define s100f forc[23]

#define s5m forc[24]
#define s10m forc[25]
#define s15m forc[26]
#define s20m forc[27]
#define s25m forc[28]
#define s30m forc[29]
#define s35m forc[30]
#define s40m forc[31]
#define s45m forc[32]
#define s50m forc[33]
#define s55m forc[34]
#define s60m forc[35]
#define s65m forc[36]
#define s70m forc[37]
#define s75m forc[38]
#define s80m forc[39]
#define s100m forc[40]

#define TFR forc[41]

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
    int N=42;
    odeforcs(&N, forc);
}

/* derivative function */
void derivsc(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if (ip[0] <2) error("nout should be at least 2");

    int i=0; 
    double Total[3];

    /* sum up totals by sex and overall */
    Total[0] = sumsum(y,0,16);
    Total[1] = sumsum(y,17,33);
    Total[2] = Total[0] + Total[1];

    /* Calculate the number of births by gender */
    double births=0;
    for(i=3; i<10; i++){
      births = births + y[i]*forc[i-3]*TFR/(100*5);
    }
    double m_to_f = 1.03;
    double births_m = births*1.03/2.03;
    double births_f = births/2.03;

    /* Then calculate the proportion of people who should die per year by age */
    double m_b[34];
    double d_age[34] = {d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d2,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d2};
    for (i=0; i<34; i++){
      m_b[i] = ((1-forc[i+7])/(1/d_age[i]));
    }
 
    /* Female */
 
    ydot[0] = births_f - d_age[0]*y[0] - m_b[0]*y[0];

    for (i=1; i<17; i++) ydot[i] = d_age[i-1]*y[i-1] - d_age[i]*y[i] - m_b[i]*y[i];

    /* Male */

    ydot[17] = births_m - d_age[17]*y[17] - m_b[17]*y[17];

    for (i=18; i<34; i++) ydot[i] = d_age[i-1]*y[i-1] - d_age[i]*y[i] - m_b[i]*y[i];

    yout[0] = Total[0];
    yout[1] = Total[1];
    yout[2] = Total[2];

}



