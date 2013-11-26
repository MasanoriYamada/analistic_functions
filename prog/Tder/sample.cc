#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include "VM.h"

void usage(const char* msg=NULL,...)
{
  if(msg!=NULL){
    va_list args;
    va_start(args, msg);
    vfprintf(stderr, msg, args);
    va_end(args);

    fprintf(stderr, "\n");
  }
  fprintf(stderr,
	  "usage: sample (E[MeV]) (a[fm]) (t0[fm]) (Ndeg)\n"
	  "output:\n"
	  "       (1) analytic.txt: analytic derivatives\n"
	  "       (2) numeric.txt:  numeric  derivatives\n"
	  );
}

int main(int argc,char* argv[])
{
  using namespace VM;

  argc--; argv++;

  if(argc != 4) { usage(); exit(1); }

  //  double hbarc = 197.3;
    double hbarc = 1.0;

  double E    = atof(argv[0]);
  double a    = atof(argv[1])/hbarc;

  double t0   = atof(argv[2])/hbarc;
  int    Ndeg = atoi(argv[3]);

  vandermonde van(Ndeg);
  
  double f[Ndeg];
  
  for(int i = 0; i < Ndeg; i++){
    double t = t0 + double(i)*a;
    f[i] = exp(E * t);
  }

  double df[Ndeg][Ndeg];

  for(int nth = 0; nth < Ndeg; nth++){
    fmadd(&df[nth][0], van.Mder[nth], f, true);
  }

  //---------------------------------------------------------------------
  // numeric
  //---------------------------------------------------------------------
  {
    FILE *fp;
    if((fp = fopen("numeric.txt","w"))==NULL){
      fprintf(stderr, "cannot fopen 'numeric.txt'\n");
      exit(1);
    }
    fprintf(fp, "# E=%4.1g MeV, a=%1.4g fm, t0=%1.4g fm, Ndeg=%d\n",
	    E, a*hbarc, t0*hbarc, Ndeg);
    fprintf(fp, "#%-10s       ", "t [fm]");
    fprintf(fp, "%-10s    ",   "f(t)   [fm^{0}]");
    fprintf(fp, "%-10s    ",   "f'(t)  [fm^{-1}]");
    fprintf(fp, "%-10s    ",   "f''(t) [fm^{-2}]");
    fprintf(fp, "....  \n");
    for(int i = 0; i < Ndeg; i++){
      double t = t0 + double(i)*a;
      fprintf(fp,"%1.8e    ", t*hbarc);
      for(int nth = 0; nth < Ndeg; nth++){
	fprintf(fp,"%1.12e  ", df[nth][i]/pow(a*hbarc, nth));
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
  }

  //---------------------------------------------------------------------
  // analytic
  //---------------------------------------------------------------------
  {
    FILE *fp;
    if((fp = fopen("analytic.txt","w"))==NULL){
      fprintf(stderr, "cannot fopen 'analytic.txt'\n");
      exit(1);
    }
    fprintf(fp, "# E=%4.1g MeV, a=%1.4g fm, t0=%1.4g fm, Ndeg=%d\n",
	    E, a*197.3, t0*197.3, Ndeg);
    fprintf(fp, "#%-10s       ", "t [fm]");
    fprintf(fp, "%-10s    ",   "f(t)   [fm^{0}]");
    fprintf(fp, "%-10s    ",   "f'(t)  [fm^{-1}]");
    fprintf(fp, "%-10s    ",   "f''(t) [fm^{-2}]");
    fprintf(fp, "....  \n");
    double t1 = 0.5*t0;
    double t2 = 1.5*(t0 + Ndeg*a);
    double dt = (t2 - t1)/100.0;
    for(double t = t1; t <= t2; t+= dt){
      fprintf(fp, "%1.8e    ", t*hbarc);
      for(int nth = 0; nth < Ndeg; nth++){
	fprintf(fp, "%1.8e  ", exp(E*t)*pow(E/hbarc,nth));
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
  
  return 0;
}
