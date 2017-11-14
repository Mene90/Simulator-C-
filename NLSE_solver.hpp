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
  // const double PI = 3.141592653589793238460;
  // const std::complex<double> i(int re = 0, int  im = 1);

public:
  NLSE_solver(int nstep);
  void solve (CArray& ux, CArray& beta, DArray& geffz, double dz);
  int getNsteps();
private:
  void ssfm (CArray& ux, CArray& beta, DArray& geffz, double dz);
  void ln_step(CArray& Fb, CArray& ux);
  void nl_step(Complex geffz, CArray& ux);
  
};
