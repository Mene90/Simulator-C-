#include <complex>
#include <iostream>
#include <valarray>
#include "fft.h"

using namespace std;

class NLSE_solver{

  int nstep;
  const double PI = 3.141592653589793238460;

  typedef std::complex<double> Complex;
  typedef std::valarray<Complex> CArray;

public:
  NLSE_solver(int nstep)
              : nstep{nstep} {}
  void ssfm (CArray& x){
    ln_step(betaz,x);
    nl_step(xi,x)
    for (i=1,i<nstep,i++){
      ln_step(betaz,x);
      nl_step(xi,x)
    }
    ln_step(betaz,x);
    nl_step(xi,x);
  }

  void ln_step(CArray& expbetaz, CArray& x){
    fft(x);
    x = x*betaz;
    ifft(x);
  }

  void nl_step(Complex xi, CArray& x){
    power = pow(abs(x),2);
    
  }

}
