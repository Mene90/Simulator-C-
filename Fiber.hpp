#include <cmath>
#include <iostream>
#include <complex>
#include <valarray>

#include "NLSE_solver.hpp"

class Fiber {

private:
  typedef std::complex  <double>   Complex;
  typedef std::valarray <Complex>  CArray;
  typedef std::valarray <double>   DArray;

  double length, alphadB, dispersion, slope, nlindex, aeff;
  double alphalin;

public:
Fiber(double length, double alphadB, double dispersion, double slope, double nlindex, double aeff);
void propagation (CArray& ux, double gmp, CArray& beta, NLSE_solver num_method);

};
