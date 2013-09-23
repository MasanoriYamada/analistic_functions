#ifndef FIT_H_
#define FIT_H_

#include <complex>
#include <string>

typedef std::complex<double> COMPLEX; 


static const double ascale = 0.1209;
static const double hbarc = 197.327;
static const double ToMev = hbarc/ascale;

//set in out info

#define binCor(iT,iconf) binCor[iT + datasize*iconf]



void fit(double* datain,double* datain_sigma, double& b1, double& b2,double& chisq);


extern "C" {
  void mrqmin_(double x[],double y[],double sig[], int& ndata,double a[],int ia[],int& ma,
	       double covar[], double alpha[],int& nca, double &chisq,
	       void (*func)(double& x,double a[],double& yfit,double dyda[],int& ma),double& alambda);	       
  void single_exp(double& x,double a[],double& yfit,double dyda[],int& ma);
}



#endif
