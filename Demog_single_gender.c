/* Demog model - testing with single year age bins and including gender */

/* !!!!!!! NOTE IN C INDEXING STARTS AT 0 AND NOT 1 (i.e. parms[0] is the first value in parms) !!!!!! */

/* C libraries needed */
#include <R.h>
#include <math.h>

/* You need to define number of parameters and forcing functions passed to the model here */
/* These must match number in intializer functions below */
static double parms[2];
static double forc[171];

/* header file including all the define statements - these are used to keep up with the parameters and forcings */
/* Parameters are constant over time  - NOTE THESE AREN'T CURRENTLY USED IN THIS MODEL */
#define age_1 parms[0]       /* rate of entry/exit to/from age group = 1/width */
#define age_2 parms[1]

/* Forcings are time dependant functions */

/* Age specific birth rates (by 5 year age group) */

#define b20 forc[0]
#define b25 forc[1]
#define b30 forc[2]
#define b35 forc[3]
#define b40 forc[4]
#define b45 forc[5]
#define b50 forc[6]

/* Total fertility rate */

#define TFR forc[7]

/* Ratio of male to female births */

#define m_to_f forc[8]

/* Age specific mortality (single year age groups, by gender) */

#define sf1 forc[9]         
#define sf2 forc[10] 
#define sf3 forc[11] 
#define sf4 forc[12] 
#define sf5 forc[13]
#define sf6 forc[14]
#define sf7 forc[15]
#define sf8 forc[16]
#define sf9 forc[17]
#define sf10 forc[18]
#define sf11 forc[19]
#define sf12 forc[20]
#define sf13 forc[21]
#define sf14 forc[22]
#define sf15 forc[23]
#define sf16 forc[24]
#define sf17 forc[25]
#define sf18 forc[26]       
#define sf19 forc[27] 
#define sf20 forc[28] 
#define sf21 forc[29] 
#define sf22 forc[30]
#define sf23 forc[31]
#define sf24 forc[32]
#define sf25 forc[33]
#define sf26 forc[34]
#define sf27 forc[35]
#define sf28 forc[36]
#define sf29 forc[37]
#define sf30 forc[38]
#define sf31 forc[39]
#define sf32 forc[40]
#define sf33 forc[41]
#define sf34 forc[42]
#define sf35 forc[43]
#define sf36 forc[44]
#define sf37 forc[45]
#define sf38 forc[46]         
#define sf39 forc[47] 
#define sf40 forc[48] 
#define sf41 forc[49] 
#define sf42 forc[50]
#define sf43 forc[51]
#define sf44 forc[52]
#define sf45 forc[53]
#define sf46 forc[54]
#define sf47 forc[55]
#define sf48 forc[56]
#define sf49 forc[57]
#define sf50 forc[58]
#define sf51 forc[59]
#define sf52 forc[60]
#define sf53 forc[61]
#define sf54 forc[62]
#define sf55 forc[63]
#define sf56 forc[64]
#define sf57 forc[65]
#define sf58 forc[66]        
#define sf59 forc[67] 
#define sf60 forc[68] 
#define sf61 forc[69] 
#define sf62 forc[70]
#define sf63 forc[71]
#define sf64 forc[72]
#define sf65 forc[73]
#define sf66 forc[74]
#define sf67 forc[75]
#define sf68 forc[76]
#define sf69 forc[77]
#define sf70 forc[78]
#define sf71 forc[79] 
#define sf72 forc[80]
#define sf73 forc[81]
#define sf74 forc[82]
#define sf75 forc[83]
#define sf76 forc[84]
#define sf77 forc[85]
#define sf78 forc[86]
#define sf79 forc[87]
#define sf80 forc[88]
#define sf81 forc[89]

#define sm1 forc[90]
#define sm2 forc[91]
#define sm3 forc[92] 
#define sm4 forc[93] 
#define sm5 forc[94]
#define sm6 forc[95]
#define sm7 forc[96]
#define sm8 forc[97]
#define sm9 forc[98]
#define sm10 forc[99]
#define sm11 forc[100]
#define sm12 forc[101]
#define sm13 forc[102]
#define sm14 forc[103]
#define sm15 forc[104]
#define sm16 forc[105]
#define sm17 forc[106]
#define sm18 forc[107]       
#define sm19 forc[108] 
#define sm20 forc[109] 
#define sm21 forc[110] 
#define sm22 forc[111]
#define sm23 forc[112]
#define sm24 forc[113]
#define sm25 forc[114]
#define sm26 forc[115]
#define sm27 forc[116]
#define sm28 forc[117]
#define sm29 forc[118]
#define sm30 forc[119]
#define sm31 forc[120]
#define sm32 forc[121]
#define sm33 forc[122]
#define sm34 forc[123]
#define sm35 forc[124]
#define sm36 forc[125]
#define sm37 forc[126]
#define sm38 forc[127]         
#define sm39 forc[128] 
#define sm40 forc[129] 
#define sm41 forc[130] 
#define sm42 forc[131]
#define sm43 forc[132]
#define sm44 forc[133]
#define sm45 forc[134]
#define sm46 forc[135]
#define sm47 forc[136]
#define sm48 forc[137]
#define sm49 forc[138]
#define sm50 forc[139]
#define sm51 forc[140]
#define sm52 forc[141]
#define sm53 forc[142]
#define sm54 forc[143]
#define sm55 forc[144]
#define sm56 forc[145]
#define sm57 forc[146]
#define sm58 forc[147]        
#define sm59 forc[148] 
#define sm60 forc[149] 
#define sm61 forc[150] 
#define sm62 forc[151]
#define sm63 forc[152]
#define sm64 forc[153]
#define sm65 forc[154]
#define sm66 forc[155]
#define sm67 forc[156]
#define sm68 forc[157]
#define sm69 forc[158]
#define sm70 forc[159]
#define sm71 forc[160] 
#define sm72 forc[161]
#define sm73 forc[162]
#define sm74 forc[163]
#define sm75 forc[164]
#define sm76 forc[165]
#define sm77 forc[166]
#define sm78 forc[167]
#define sm79 forc[168]
#define sm80 forc[169]
#define sm81 forc[170]


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
    int N=171;
    odeforcs(&N, forc);
}

/* ##### EVENTS ARE USED TO ADD BIRTHS AND SHIFT THE POPULATION BY AGE - EQUIVALENT TO THE METHOD OF SCHENZLE ###### */
void event(int *n, double *t, double *y) 
{
  int i;
  int j;
  int k = 15;
  
  /* Calculate the number of births by gender */
  double births=0;
  
  for(i=0; i<7; i++){
    for (j=1; j<6; j++){
      births = births + y[k]*forc[i]*TFR/(100*5);
      k = k+1;
    }  
  }
  double births_m = births/(1+(1/m_to_f));
  double births_f = births_m/m_to_f;
  
  double temp[162];
  for (i=0; i<162; i++) temp[i] = y[i];
  
  y[0] = births_f;
  for (i=1; i<80; i++) y[i] = temp[i-1];
  y[80] = y[80] + temp[79];
  
  y[81] = births_m;
  for (i=82; i<161; i++) y[i] = temp[i-1];
  y[161] = y[161] + temp[160];
  
}

/* ###### DERIVATIVE FUNCTION ###### */

void derivsc(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if (ip[0] <2) error("nout should be at least 2");
    
    /* Expand state variables so can use more meaningful names than y and ydot (the variables and rates of change are vectors y and ydot here) */
    
    /* These are the variables */
    double Sf[81];
    double Sm[81];
   
    /* These are the rates of change (same names but prefixed with d) */
    double dSf[81];
    double dSm[81];
   
    /* intergers to use as counters */ 
    int i;  

    /* map y to S */
    for (i=0; i<81; i++) {
      Sf[i] = y[i];            
      Sm[i] = y[i+81];    
    }

    /* Use sumsum function to add up population */
    double Total_f = sumsum(Sf,0,80);                
    double Total_m = sumsum(Sm,0,80);
    double Total = Total_f + Total_m;

    /* Assign mortality rates */
    double m_f[81];
    double m_m[81];
    for (i=0; i<81; i++){
      m_f[i] = forc[i+9];
      m_m[i] = forc[i+90];
    }
    
    /* Derivatives */ 

    for (i=0; i<81; i++) { 
      dSf[i] = - m_f[i]*Sf[i];
      dSm[i] = - m_m[i]*Sm[i];
    }
    
    /* Calculate total births */
    
    /* Calculate total deaths */
    double Tot_deaths_f = 0;
    double Tot_deaths_m = 0;
    double Tot_deaths = 0;
    for (i=0; i<81; i++){
      Tot_deaths_f =  Tot_deaths_f + m_f[i]*Sf[i];
      Tot_deaths_m =  Tot_deaths_m + m_m[i]*Sm[i];
    }
    Tot_deaths = Tot_deaths_f + Tot_deaths_m;
    
    /* map ds to ydot */

    for (i=0; i<81; i++) {
      ydot[i] = dSf[i]; 
      ydot[i+81] = dSm[i];
    }

    /* Finally assign the things we want to use in R (in addition to the state variables) to yout */
    /* The ode call in R needs to define the number of these and give names */
    yout[0] = Total;
    yout[1] = Total;
    yout[2] = Tot_deaths;
    
}



