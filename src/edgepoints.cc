#include <cmath>
#include <cassert>
#include <Rmath.h>
#include <Rdefines.h>
#include "statist.h"
using namespace std;

//#define MATHLIB_STANDALONE
//#include <Rmath.h>

#define QUAD(x) ((x)*(x))
#define PI 3.141592654
#define i2v(i) ((i*100000))
#define d2v(d) (((int)(d*100000)))

#define BFAK 1.6 // Faktor zum Beschleunigen der Schrittweite
#define DFAK 0.7 // Faktor zum D?mpfen der Schrittweite (Armijo-Schrittweite)
#define AFAK 0.5 // Faktor f?r die Abbruchbedingun im Armijo-Verfahren

// Median der ersten n Werte des Arrays data
// ACHTUNG: Das Array wird teilsortiert!!!
double median (double *const data, const int n)
{
  qsort(data,n,0,n-1,2);
  return n%2==1?data[(n-1)/2]:(data[(n-2)/2]+data[n/2])/2;
}

// MAD der ersten n Werte des Arrays data
// ACHTUNG: Das Array wird teilsortiert!!!
double mad (double *const data, const int n)
{
  double med, result;
  double *diff;
  int i;
  
  med=median(data,n);
  
  diff= (double *) R_alloc((n+1), sizeof(double));
  //UL: was: malloc()
  for (i=0; i<n; i++) 
    diff[i]=fabs(data[i]-med);
  
  result = median(diff,n)/0.6744898;
  //UL: was: free (diff);
  return result;
}

// Unteres Quartil der ersten n Werte des Arrays data
// ACHTUNG: Das Array wird teilsortiert!!!
double quartil(double *const data, const int n) {
  qsort(data,n,0,n-1,4);
  return data[(int)floor(((double)n)/4)];
}

// Q-Sch?tzer (=unteres Quartil der paarweisen Diff.)
double q_estimate(double *const data, const int n) {
  double *diff;
  double result;
  int i, j, l;
  
  diff= (double *) R_alloc((n*(n-1)+1), sizeof(double));
  //UL: was: malloc()
  l=0;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      if (i!=j){
        diff[l]=fabs(data[i]-data[j]);
        l++;
      }
      
      result=quartil(diff,l)/0.44;
      //UL: was: free(diff);
      return result;
}

//--------------------------------------------------------------------
// RDMKE f?r ?quidistante St?tzstellen und vollst?ndige Beobachtungen
//--------------------------------------------------------------------

// Drehung
double delta(const double x){
  return (x<0?-1:1);
}
double r_1(const double theta, const double x, const double y){
  double arctan;
  
  if (x==0)
    arctan = y >= 0 ? PI/2 : -PI/2;
  else
    arctan = atan(y/x);
  return sqrt(QUAD(x)+QUAD(y))*delta(x)*cos(arctan-theta);
}
double r_2(const double theta, const double x, const double y){
  double arctan;
  
  if (x==0)
    arctan = y >= 0 ? PI/2 : -PI/2;
  else
    arctan = atan(y/x);
  return sqrt(QUAD(x)+QUAD(y))*delta(x)*sin(arctan-theta);
}

// Scorefunktionen

double SIG, SIG2, SIG22, SIG4, SPS;

void setgaussconst(double sigma){
  SIG=sigma;
  SIG2=sigma*sigma;
  SIG22=2*SIG2;
  SIG4=SIG2*SIG2;
  SPS=2.506628*SIG;
}

double gauss (double x, int ableitung){
  double f, x2;
  
  assert (0 <= ableitung && ableitung <= 3);
  
  x2 = x*x;
  f = -exp(-x2/SIG22)/SPS;
  
  if (ableitung == 0)
    return f;
  if (ableitung == 1)
    return -f*x/SIG2;
  if (ableitung == 2)
    return f*(x2/SIG2-1)/SIG2;
  
  return f*(3*x-x*x2/SIG2)/SIG4;
}

void sethuber(double grenze){
  SIG=grenze;
  SIG2=1/(2*grenze);
}

double huber (double x, int ableitung){
  
  if (fabs(x)<SIG) {
    if (ableitung == 0)
      return SIG2*x*x;
    if (ableitung == 1)
      return SIG2*2*x;
    if (ableitung == 2)
      return SIG2*2;
    return 0;
  }
  if (ableitung == 0)
    return fabs(x)/2;
  if (ableitung == 1)
    return x>0?1:-1;
  return 0;
}

double quad (double x, int ableitung){
  
  if (ableitung == 0)
    return x*x;
  if (ableitung == 1)
    return 2*x;
  if (ableitung == 2)
    return 2;
  return 0;
}

// Kernfunktionen

double rechteckkern (const double x, const double y){
  return (d2v(fabs(x))<=i2v(1)?1:0)*(d2v(fabs(y))<=i2v(1)?1:0);
}


double dreieckkern (const double x, const double y){
  return (d2v(fabs(x))<i2v(1)?1-x:0)*(d2v(fabs(y))<i2v(1)?1-y:0);
}

double gdreieckkern (const double x, const double y){
  double y2;
  
  y2 = fabs(fabs(y)-0.5)*2;
  return (d2v(fabs(x))<i2v(1)?1-x:0)*(d2v(y2)<i2v(1)?1-y2:0);
}

double gauss1dim (const double x){
  if (d2v(fabs(x))<=i2v(1))
    return exp(-QUAD(x/0.5)/2)/(0.5*sqrt(2*PI));
  else
    return 0;
}    

double gausskern (const double x, const double y){
  return gauss1dim(x)*gauss1dim(y);
}

/////////////////////////////////////////////////////////////////////

extern "C" {
  
  void c_edgepoints (
      double *const z,           // Beobachtungen
      int    *const nrow,        // Kantenlaengen
      int    *const ncol, 
      int    *const kernel,      // 0: Rechteckkern
      // 1: Dreieckkern
      // 2: Dreieckkern mit Mittelpunkten +-0.5
      // 3: Gausskern
      // 4: beliebige Gewichtsmatrix ('kernmat')
      double *const h1n,         // Bandbreiten
      double *const h2n,
      int    *const estimator,   // Sch?tzung in den Fenstern:
      // 0: Kernsch?tzer mit 'kernel' als Kern
      // 1: Median ('kernel' mu? 0 sein)
      // 2: M-Kernsch?tzung Mean als Startwert
      // 3: M-Kernsch?tzung Median als Startwert
      // 4: beliebige Sch?tzung mit 'estfkt'
      //    als Sch?tzfunktion
      // 5: multipler Test auf gleiche 
      //    Erwartungswerte (Mean/emp.Var.)
      // 6: multipler Test auf gleiche 
      //    Erwartungswerte (Median/Q-Sch?tzer)
      int    *const score,       // Scorefunktion der M-Kernsch?tzung:
      // 0: negative Gau?dichte mit Var. 'sigma'
      // 1: Huber mit Grenze 'sigma'
      // 2: 'scorefkt'
      // 3: negative Quadratfkt -> normale Kern-
      // sch?tzung (sinnvoll nur zu Testzwecken)
      
      double *const sigma,       // Parameter der Scorefunktion
      double *const kernmat,     // Gewichtsmatrix
      double *const max_schritt, // Maximale Schrittl?nge in der Minimierung
      // der M-Kernsch?tzer. Sinnvollerwert ist
      // z.B. "Sigma der Scorefkt"*
      //      "Maximalwert der Beobachtungen"
      // estfkt
      int    *const asteps,      // Anzahl der Startwerte f?r den Winkel
      // R?ckgabewerte
      double *const angle,       // Winkel des maximalen Sprungs
      double *const value)       // und zugehoehriger Wert der scorefkt. bzw.
  // p-Wert bei 'estimator' = 5,6
  {
    // Nur zul?ssige Kerne
    assert (((*kernel)>=0) && ((*kernel)<=4));
    // Nur zul?ssige Sch?tzer
    assert (((*estimator)>=0) && ((*estimator)<=6) && ((*estimator) != 4));
    // Zul?ssige Winkelzahl
    assert ((*asteps) > 0);
    // Nur Rechteckkern bei Median
    assert (((*estimator)!=1) || ((*kernel)==0));
    // Zul?ssige Scorefunktionen
    assert (((*estimator)!=2 && (*estimator)!=3) ||
      ((*score)==0 || (*score)==1 || (*score)==9));
    
    // Gr??e der Umgebung
    int    env;                
    env = (int)ceil(sqrt(QUAD((*h1n))+QUAD((*h2n)))*
      ((*ncol)>(*nrow)?(*ncol):(*nrow)));
    
    // Kopie der Daten mit vergr??ertem Rand
    const int data_nrow = (*nrow) + 2 * env;
    const int data_ncol = (*ncol) + 2 * env;
    double** data = new double*[data_nrow];
    for (int i = 0; i < data_nrow; i++)
      data[i] = new double[data_ncol];
    // was: double* data = new double[(*nrow)+2*env][(*ncol)+2*env];
    
    // 3D-Array f?r die Gewichte an jedem Punkt bei jedem Winkel
    //const int gewichte_asteps = (*asteps);
    const int gewichte_dim = 2 * env + 1;
    double*** gewichte = new double**[(*asteps)];
    for (int i = 0; i < (*asteps); ++i) {
      gewichte[i] = new double*[gewichte_dim];
      
      for (int j = 0; j < gewichte_dim; ++j)
        gewichte[i][j] = new double[gewichte_dim];
    }
    // was: double* gewichte = new double[(*asteps)][2*env+1][2*env+1];
    
    // 3D-Array f?r die y-Koordinaten. Diese werden ben?tigt, um zu 
    // entscheiden, zu welchem Kern der Punkt geh?rt
    double*** yrot = new double**[(*asteps)];
    for (int i = 0; i < (*asteps); ++i) {
      yrot[i] = new double*[gewichte_dim];
      
      for (int j = 0; j < gewichte_dim; ++j)
        yrot[i][j] = new double[gewichte_dim];
    }
    // was: double* yrot = new double[(*asteps)][2*env+1][2*env+1];	
    
    // Beobachtungen in den beiden Fenstern, deren gewichtete Summe und 
    // Z?hler f?r deren Gr??e
    const int fenster = QUAD(2 * env + 1);
    double* fenster_l = new double[fenster];
    double* fenster_r = new double[fenster];
    double s_l, s_r;
    int    n_l, n_r;
    double sgew_l, sgew_r;
    
    // Variablen f?r die M-Kernsch?tzungen
    double y_l=0, y_r=0;
    double y0l, y1l, y2l, y0r, y1r, y2r;
    double richtung_l, letzte_richtung_l, 
    richtung_r, letzte_richtung_r;
    double diff_l, diff_r;
    double (*scorefkt) (double, int);
    double (*kernelfkt) (double, double);
    
    // Variablen f?r die multiplen Tests
    double mean_l, mean_r, var;
    double q_l, q_r, m_l, m_r;
    
    // andere Variablen
    int    i, j, k, m, n;
    double theta;
    double x, y, gew;
    double diff;
    double fak_l, fak_r;
    double hilf;
    
    // vergr??erten Rand berechnen
    for(i=0; i<(*nrow)+2*env; i++)
      for(j=0; j<(*ncol)+2*env; j++)
        data[i][j]=
          z[(i<=env ? 0 : (i>=(*nrow)+env ? (*nrow)-1 : i-env))
        +((*nrow)*
          ((j<=env ? 0 : (j>=(*ncol)+env ? (*ncol)-1 : j-env))))];
    
    // Scorefunktion setzen
    if (((*estimator)==2 || (*estimator)==3)){
      if ((*score)==0){
        setgaussconst(*sigma);
        scorefkt=gauss;
      }
      else if ((*score)==1){
        sethuber(*sigma);
        scorefkt=huber;
      }
      else if ((*score)==9){
        scorefkt=quad;
      }
      else
        scorefkt=gauss; // Compiler beruhigen
    }
    else
      scorefkt=gauss; // Compiler beruhigen
    
    // Kernfunktion setzen
    if ((*kernel)==0)
      kernelfkt=rechteckkern;
    else if ((*kernel)==1)
      kernelfkt=dreieckkern;
    else if ((*kernel)==2)
      kernelfkt=gdreieckkern;
    else if ((*kernel)==3)
      kernelfkt=gausskern;
    else
      kernelfkt=gausskern; // Compiler beruhigen
    
    // Koordinaten und Gewichte berechnen
    for (i=-env; i<=env; i++){           
      for(j=-env; j<=env; j++){        
        for (k=0; k<(*asteps); k++){
          theta=-PI/2+(k*PI/(*asteps));
          
          x=r_1(theta,(double)i/(*ncol),(double)j/(*nrow))/(*h1n);
          y=r_2(theta,(double)i/(*ncol),(double)j/(*nrow))/(*h2n);
          
          // Der Kern hat einen Tr?ger von [-1,1]x[-1,1], die ein-
          // seitigen Kerne m?ssen aber einen Tr?ger von 
          // [-0.5,0.5]x[-1,0] bzw. [-0.5,0.5]x[0,1] haben.
          // Darum mu? die x-Koordinate skaliert werden.
          x=2*x;
          
          if ((*kernel)!=4)
            gewichte[k][i+env][j+env]=kernelfkt(x,y);
          else
            gewichte[k][i+env][j+env]=
              kernmat[k*QUAD(2*env+1)+(i+env)*(2*env+1)+(j+env)];
          
          // y-Koordinate merken.
          yrot[k][i+env][j+env]=y;
        }
      }
    }
    
    // Eigentliche Berechnung
    for (i=0; i<(*nrow); i++){                  // F?r jede Beobachtung
      Rprintf("row: %i/%i\n",i,*nrow);
      for (j=0; j<(*ncol); j++){
        for (k=0; k<(*asteps); k++){        // und jeden Winkel
          //---------------------------------------------------
          // Startwerte bzw. tats?chliche Sch?tzungen berechnen
          //---------------------------------------------------
          
          // Gewichtete Beobachtungen im Fenster bestimmen:
          n_l = 0;
          n_r = 0;
          s_l = 0;
          s_r = 0;
          sgew_l = 0;
          sgew_r = 0;
          for (m=-env; m<=env; m++){
            for(n=-env; n<=env; n++){
              y=  yrot[k][m+env][n+env];
              gew=gewichte[k][m+env][n+env];
              
              if (gew != 0){
                // Zu welcher Seite geh?rt der Punkt? Wenn der
                // Punkt zu beiden geh?rt (y=0),
                // wird er gar nicht summiert.
                if (d2v(y) < i2v(0)){
                  fenster_l[n_l] = 
                    data[i+m+env][j+n+env];
                  s_l = s_l + gew*fenster_l[n_l];
                  n_l = n_l + 1;
                  sgew_l = sgew_l + gew;
                }
                else if (i2v(0) < d2v(y)){
                  fenster_r[n_r] 
                  = data[i+m+env][j+n+env];
                  s_r = s_r + gew*fenster_r[n_r];
                  n_r = n_r + 1;
                  sgew_r = sgew_r + gew;
                }
              }
            }
          }
          // M-Kernsch?tzungen berechnen
          if ((*estimator==2 || (*estimator)==3)) { 
            // Startwerte:
            if ((*estimator)==2){ // Mean als Startwert
              y_l = mean(fenster_l,n_l);
              y_r = mean(fenster_r,n_r);
            }
            else {                // Median als Startwert
              y_l = median(fenster_l,n_l);
              y_r = median(fenster_r,n_r);
            }
            richtung_l = richtung_r = 0;
            diff_l     = diff_r     = 0;
            fak_l      = fak_r      = 0.7;
            // konkav     = konvex     = 0;
            do {
              y0l = y1l = y2l = y0r = y1r = y2r = 0;
              // Werte im linken und rechten Fenster bestimmen
              for (m=-env; m<=env; m++){
                for(n=-env; n<=env; n++){
                  y=  yrot[k][m+env][n+env];
                  gew=gewichte[k][m+env][n+env];
                  if (gew != 0){
                    // Zu welcher Seite geh?rt der Punkt?
                    if (d2v(y) < i2v(0)){
                      y0l += gew * scorefkt(
                        y_l-data[i+m+env][j+n+env],0);
                      y1l += gew * scorefkt(
                        y_l-data[i+m+env][j+n+env],1);
                      y2l += gew * scorefkt(
                        y_l-data[i+m+env][j+n+env],2);
                    }
                    if (i2v(0) < d2v(y)){
                      y0r += gew * scorefkt(
                        y_r-data[i+m+env][j+n+env],0);
                      y1r += gew * scorefkt(
                        y_r-data[i+m+env][j+n+env],1);
                      y2r += gew * scorefkt(
                        y_r-data[i+m+env][j+n+env],2);
                    }
                  }
                }
              }
              // Schrittl?nge im linken Fenster bestimmen
              letzte_richtung_l = richtung_l;
              if (y2l <= 0){     // konkaver Fall
                richtung_l = (y1l > 0 ? -1 : 1);
                diff_l = BFAK*(*max_schritt) * richtung_l;
                fak_l  = 1/BFAK;
              }
              else {             // konvexer Fall: Newtow
                diff_l     = -y1l/y2l;
                fak_l = DFAK;
              }
              
              if (fak_l < 1){
                fak_l  = fak_l/DFAK;
                do {
                  fak_l = fak_l * DFAK;
                  hilf = 0;
                  for (m=-env; m<=env; m++){
                    for(n=-env; n<=env; n++){
                      y=  yrot[k][m+env][n+env];
                      gew=gewichte[k][m+env][n+env];
                      if ((gew != 0) &&
                          (d2v(y) < i2v(0))){
                        hilf += gew * scorefkt(
                          (y_l+fak_l*diff_l)
                        -data[i+m+env][j+n+env]
                        ,0);
                      }
                    }
                  }
                } while ((hilf - 
                  (y0l + AFAK*fak_l*y1l*diff_l))
                           > 0.00001);
              }
              
              // Schrittl?nge im rechten Fenster bestimmen
              letzte_richtung_r = richtung_r;
              if (y2r <= 0){     // konkaver Fall
                richtung_r = (y1r > 0 ? -1 : 1);
                diff_r = (*max_schritt) * richtung_r;
                fak_r  = 1/BFAK;
              }
              else {             // konvexer Fall: Newtow
                diff_r     = -y1r/y2r;
                fak_r = DFAK;
              }
              
              if (fak_r < 1){
                fak_r  = fak_r/DFAK;
                do {
                  fak_r = fak_r * DFAK;
                  hilf = 0;
                  for (m=-env; m<=env; m++){
                    for(n=-env; n<=env; n++){
                      y=  yrot[k][m+env][n+env];
                      gew=gewichte[k][m+env][n+env];
                      if ((gew != 0) &&
                          (i2v(0) < d2v(y))){
                        hilf += gew * scorefkt(
                          (y_r+fak_r*diff_r)
                        -data[i+m+env][j+n+env]
                        ,0);
                      }
                    }
                  }
                } while ((hilf - 
                  (y0r + AFAK*fak_r*y1r*diff_r))
                           > 0.00001);
              }
              
              y_l = y_l + fak_l*diff_l;
              y_r = y_r + fak_r*diff_r;
            } while (fabs(fak_l*diff_l) > 0.0001 || 
              fabs(fak_r*diff_r) > 0.0001);
          }
          
          // Differenz der Sch?tzungen in den Fenstern bzw.
          // p-Werte berechnen
          if ((*estimator)==0)
            diff=fabs((s_l/sgew_l)-(s_r/sgew_r));
          else if ((*estimator)==1)
            diff=fabs(median(fenster_l,n_l)-median(fenster_r,n_r));
          else if ((*estimator)==2 || (*estimator)==3)
            diff = fabs(y_l-y_r);
          else if ((*estimator)==5) { // Test (Mean/emp.Var)
            // Empirische Varianz berechnen
            mean_l = s_l/n_l;
            mean_r = s_r/n_r;
            if (d2v(fabs(mean_l-mean_r))>0) {
              var    = 0;
              for (m=0; m<n_l; m++)
                var = var + QUAD(fenster_l[m]-mean_l);
              for (m=0; m<n_r; m++)
                var = var + QUAD(fenster_r[m]-mean_r);
              var = var/(n_l+n_r-2);
              if (var==0)
                diff = 0;
              else 
                diff = -(1-pt(fabs(((mean_l-mean_r)/sqrt(var))*
                  sqrt((double)(n_l*n_r)/
                         (n_l+n_r))),
                         n_l+n_r-2,1,0))*(*asteps)*2;
            }
            else
              diff = -1;
          }
          else if ((*estimator)==6) { // Test (Median/Q-Sch?tzer)
            m_l=median(fenster_l,n_l);
            m_r=median(fenster_r,n_r);
            if (d2v(fabs(m_l-m_r))>0) {
              q_l=mad(fenster_l,n_l);
              q_r=mad(fenster_r,n_r);
              if (q_l+q_r==0)
                diff = 0;
              else 
                diff = -(1-pnorm(sqrt(((double)(n_l+n_r))/2)*
                  fabs(m_l-m_r)/
                                   sqrt(PI/2*(QUAD(q_l)+
                                     QUAD(q_r))),
                                     0,1,1,0))*(*asteps)*2;
            }
            else
              diff = -1;
          }
          else
            diff=0; // Nur, um den Compiler zu beruhigen
          
          if (k==0 || diff>value[i+(*nrow)*j]){
            value[i+(*nrow)*j]=diff;
            angle[i+(*nrow)*j]=-PI/2+(k*PI/(*asteps));
          }
        }
      }
    }
    // delete data
    for (int i = 0; i < data_nrow; i++)
      delete [] data[i];
    delete [] data;
    return;
    // delete gewichte
    for (int i = 0; i < (*asteps); ++i) {
      for (int j = 0; j < gewichte_dim; ++j)
        delete [] gewichte[i][j];
      
      delete [] gewichte[i];
    }
    delete [] gewichte;
    // delete yrot
    for (int i = 0; i < (*asteps); ++i) {
      for (int j = 0; j < gewichte_dim; ++j)
        delete [] yrot[i][j];
      
      delete [] yrot[i];
    }
    delete [] yrot;
  } // c_edgepoints
  
} // extern
