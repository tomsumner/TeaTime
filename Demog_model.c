/* Attempt to write TB model in C code to call from R

/* compile within R with system("R CMD SHLIB Demog_model.c") */

#include <R.h>

static double parms[2];
static double forc[18];

/* A trick to keep up with the parameters and forcings */
#define d parms[0]
#define d2 parms[1]

#define birth_rate forc[0]
#define s_birth forc[1]
#define s5 forc[2]
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
/*void derivsc(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if (ip[0] <2) error("nout should be at least 2");
    
    double Total = y[0]+y[1]+y[2]+y[3]+y[4]+y[5]+y[6]+y[7]+y[8]+y[9]+y[10]+y[11]+y[12]+y[13]+y[14]+y[15]+y[16];
 
    ydot[0] = s_birth*birth_rate*Total/1000 - d*y[0];
    ydot[1] = d*s5*y[0] - d*y[1];
    ydot[2] = d*s10*y[1] - d*y[2];
    ydot[3] = d*s15*y[2] - d*y[3];
    ydot[4] = d*s20*y[3] - d*y[4];
    ydot[5] = d*s25*y[4] - d*y[5];
    ydot[6] = d*s30*y[5] - d*y[6];
    ydot[7] = d*s35*y[6] - d*y[7];
    ydot[8] = d*s40*y[7] - d*y[8];
    ydot[9] = d*s45*y[8] - d*y[9];
    ydot[10] = d*s50*y[9] - d*y[10];
    ydot[11] = d*s55*y[10] - d*y[11];
    ydot[12] = d*s60*y[11] - d*y[12];
    ydot[13] = d*s65*y[12] - d*y[13];
    ydot[14] = d*s70*y[13] - d*y[14];
    ydot[15] = d*s75*y[14] - d*y[15];
    ydot[16] = d*s80*y[15] - d2*y[16];

    yout[0] = Total;

}*/

/* derivative function */
void derivsc(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if (ip[0] <2) error("nout should be at least 2");
    
    double Total = y[0]+y[1]+y[2]+y[3]+y[4]+y[5]+y[6]+y[7]+y[8]+y[9]+y[10]+y[11]+y[12]+y[13]+y[14]+y[15]+y[16];
 
    ydot[0] = s_birth*birth_rate*Total/1000 - d*y[0];
    ydot[1] = d*s5*y[0] - d*y[1];
    ydot[2] = d*s10*y[1] - d*y[2];
    ydot[3] = d*s15*y[2] - d*y[3];
    ydot[4] = d*s20*y[3] - d*y[4];
    ydot[5] = d*s25*y[4] - d*y[5];
    ydot[6] = d*s30*y[5] - d*y[6];
    ydot[7] = d*s35*y[6] - d*y[7];
    ydot[8] = d*s40*y[7] - d*y[8];
    ydot[9] = d*s45*y[8] - d*y[9];
    ydot[10] = d*s50*y[9] - d*y[10];
    ydot[11] = d*s55*y[10] - d*y[11];
    ydot[12] = d*s60*y[11] - d*y[12];
    ydot[13] = d*s65*y[12] - d*y[13];
    ydot[14] = d*s70*y[13] - d*y[14];
    ydot[15] = d*s75*y[14] - d*y[15];
    ydot[16] = d*s80*y[15] - d2*y[16];

    yout[0] = Total;

}



