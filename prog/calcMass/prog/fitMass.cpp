//--------------------------------------------------------------------------
/**
 * @Filex main_fit_pot3g.cpp
 * @brief fiting potential 3 gauss
 * @ingroup YAMADA
 * @author  M.YAMADA * @date    Sat Jun 13 22:09:45 2013
 */
//--------------------------------------------------------------------------

#include "../include/io.h"
//#include "../../../include/analys.h"
#include "../include/analys.h"
#include "../include/jack.h"
#include "../include/fit.h"
#include <sstream>

using namespace std;


std::string physInfo = "RC16x32_B1830Kud013760Ks013710C1761";


bool inBinary = true;
bool outBinary = false; 
int datasize = TnodeSites;

main(int argc,char** argv){

  dir_path=argv[1];
  cout <<"Directory path  ::"<<dir_path<<endl;
  in_dir_path = dir_path;
  out_dir_path = dir_path;
  root_mkdir(out_dir_path.c_str());
  in_dir_path = in_dir_path + "/tcor2pt/bintcor2pt/";
  out_dir_path = out_dir_path + "/tcor2pt/mass";
  root_mkdir(out_dir_path.c_str());

  cout <<in_dir_path<<endl;

  IODATA inData;
  inData.setReadBinaryMode(inBinary);
  inData.setWriteBinaryMode(outBinary);
  inData.setConfSize(binnumber);


  JACK jackCor(Confsize,binsize,datasize);
  COMPLEX* binCor = new COMPLEX[datasize*binnumber]();
  for (int iconf=0; iconf< binnumber; iconf++) {
    double * in = new double[datasize]();
    for (int iT = T_in; iT<T_fi +1 ;iT++){
      stringstream ss;
      ss<< in_dir_path <<setw(2)<<setfill('0')<< iT ;
      std::string inPath = ss.str();
      cout<<ss.str()<<endl;
      inData.callData(&(in[iT]),inPath,"tcor2pt",physInfo,iconf,iT);
      binCor(iT,iconf)  = in[iT];
      jackCor.setBinData(in,iconf);
    }
    delete [] in;
  }

  double* err= new double[datasize]();
  err = jackCor.calcErr();
  JACK jackMass(Confsize,binsize,1);
  cout << "@Start fit"<<endl;
  for (int iconf=0; iconf< binnumber; iconf++) {
    double b1,b2,chisq;
    double* binDataIn = new double[datasize];
    for(int iT = 0;iT<datasize;iT++){
      binDataIn[iT] = binCor(iT,iconf).real();
    }
    fit(binDataIn,err,b1,b2,chisq);
    cout<<b1<<" "<<b2<<" "<<chisq<<endl;
    double mass[1]={b2};
    jackMass.setBinData(mass,iconf);
  }//iconf
  double aveErr[2];
  aveErr[0] = * jackMass.calcAve();
  aveErr[1] = * jackMass.calcErr();
  inData.outData(aveErr,out_dir_path,"Mass",physInfo,0,0,2);
    
  cout <<"@End all jobs"<<endl; 
}


  //---------------------fit function---------------------------------------//
  void single_exp(double& x,double a[],double& yfit,double dyda[],int& ma)
  {
    
    double b1=a[0];
    double b2=a[1];
    yfit = b1*exp(-b2*x);
    
    dyda[0] = exp(-b2*x);
    dyda[1] = -b1*x*exp(-b2*x);
  }


void fit(double dataIn[],double dataErr[], double& b1, double& b2,double& chisq)
{
  double x[T_fi-T_in+1];
  double y[T_fi-T_in+1];
  double dy[T_fi-T_in+1];
  for(int id = 0; id < T_fi-T_in+1 ; id ++)
    {	x[id] = ascale/hbarc*(double)(T_in + id);
      y[id] = dataIn[id+T_in];
      dy[id] = dataErr[id+T_in];
    }
    
  int ndata = T_fi-T_in + 1;
  static double a[2];
  static int    is_initial = 1;
  int   ia[]={1,1}, ma =2;
  double covar[4],alpha[4];
  int    nca = 2;
  double alambda;
	
  if(is_initial==1){
    is_initial=0;
    a[0] =3.9577e-07;
    a[1] =2000.25182 ;
  }
  // The initial call
  alambda = -1.0;
  mrqmin_(x,y,dy, ndata,a,ia, ma,	covar,alpha, nca,chisq, single_exp, alambda);
  // iterations
  for(int iter=0; iter<65536; iter++){
    mrqmin_(x,y,dy, ndata,a,ia, ma,	covar,alpha, nca,chisq, single_exp, alambda);
    if (alambda > 1.0e64) break;
  }
  if (alambda <= 1.0e64) {
    cerr << "convergence is not achieved\n";
    exit(1);
  }
  // The final call
  alambda = 0.0;
  mrqmin_(x,y,dy, ndata,a,ia, ma,	covar,alpha, nca,chisq, single_exp, alambda);
  b1 = a[0];
  b2 = a[1];
    
}
  //////////////////////////////END OF FITTING /////////////////////////////////

