#include <cmath>
#include <iostream>
#include <complex>
#include <valarray>

//Constants.hpp
#ifndef CONST_H
#define CONST_H
#include "Constants.hpp"
#endif

//NLSE_solver.hpp
#ifndef NLSE_SOLVER_H
#define NLSE_SOLVER_H
#include "NLSE_solver.hpp"
#endif

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
double getGamma(double lambda);
double getBeta2(double lambda);
};
