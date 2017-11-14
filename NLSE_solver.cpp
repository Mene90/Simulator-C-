#include <cmath>
#include <complex>
#include <iostream>
#include <valarray>
#include "fft.h"

using namespace std;

class NLSE_solver{

  typedef std::complex  <double>   Complex;
  typedef std::valarray <Complex>  CArray;
  typedef std::valarray <double>   DArray;

  int nstep;
  const double PI = 3.141592653589793238460;
  const std::complex<double> i(0, 1);

public:
  NLSE_solver(int nstep)
              : nstep{nstep} {}

  void solve (CArray& ux, CArray& beta, DArray& geffz, double dz){
      ssfm(ux,beta,geffz,dz);
  }

  int getNsteps(){return nstep;}

private:

  void ssfm (CArray& ux, CArray& beta, DArray& geffz, double dz){

    CArray Fb  = exp(-i*beta*dz);
    CArray Fhb = exp(-i*beta*dz/2);

    ln_step(Fhb,ux);
    nl_step(geffz[0],ux)

    for (k=1,k<nstep,k++){
      ln_step(Fb,ux);
      nl_step(geffz[k],ux)
    }

    ln_step(Fhb,ux);
  }

  void ln_step(CArray& Fb, CArray& ux){
    ifft(fft(ux)*Fb);
  }

  void nl_step(Complex geffz, CArray& ux){
    power = pow(abs(ux),2);
       ux = ux*exp(-i*geffz*power);
  }

}
