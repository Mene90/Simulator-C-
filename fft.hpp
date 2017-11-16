#include <complex>
#include <iostream>
#include <valarray>

//Constants.hpp
#ifndef CONST_H
#define CONST_H
#include "Constants.hpp"
#endif


typedef std::complex  <double>   Complex;
typedef std::valarray <Complex>  CArray;
typedef std::valarray <double>   DArray;

void fft(CArray &x);
void ifft(CArray& x);
