#include <cmath>
#include <complex>
#include <iostream>
#include <valarray>


//Constants.hpp
#ifndef CONST_H
#define CONST_H
#include "Constants.hpp"
#endif
//fft.hpp
#ifndef FFT_H
#define FFT_H
#include "fft.hpp"
#endif

using namespace std;

class NLSE_solver{

  typedef std::complex  <double>   Complex;
  typedef std::valarray <Complex>  CArray;
  typedef std::valarray <double>   DArray;

  int nstep;

public:
  NLSE_solver(int nstep);
  void solve (CArray& ux, CArray& beta, DArray& geffz, double dz);
  int getNsteps();
private:
  void ssfm (CArray& ux, CArray& beta, DArray& geffz, double dz);
  void ln_step(CArray& Fb, CArray& ux);
  void nl_step(Complex geffz, CArray& ux);

};
