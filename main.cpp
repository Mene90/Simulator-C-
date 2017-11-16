#include <cmath>
#include <random>
#include <complex>
#include <iostream>
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

//Fiber.hpp
#ifndef FIBER_H
#define FIBER_H
#include "Fiber.hpp"
#endif

using namespace std;

typedef std::complex  <double>   Complex;
typedef std::valarray <Complex>  CArray;
typedef std::valarray <double>   DArray;


CArray gaussGen (int N, double sgn){
  Complex nc(sgn,sgn);
  CArray x(nc,N);

  x = x.apply([](Complex n)  -> Complex  {
    double sgn = n.real();                              // per ora perche valarray.apply non permette capture nel lambda function
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> nd(0,sgn);
    return Complex(n.real()*nd(generator),n.imag()*nd(generator));
  });

  return x;
}

DArray getNormFreqs(int Nsymb, int Nt){

    double stepf = 1/Nsymb;
    double FN_array[Nsymb*Nt];

    std::generate_n (FN_array,Nt/2-stepf,[Nt,stepf](){
      double current  = -Nt/2;
      double tmp      = current;
      current        += stepf;
      return current++;
    });

    DArray FN(FN_array,Nsymb*Nt);

    return FN;
}


int main(int argc, char const *argv[]) {

  Fiber       fb(1e5,0.2,17,0,2.5e-20,80);
  NLSE_solver pro_sim(100);

  int N  = pow(2,4);
  int Nt = 1;


  DArray FN = getNormFreqs(N,Nt);

  CArray ux = gaussGen(N,0.7071);

  // for(auto n : ux) {
  //      std::cout << n << '\n';
  //  }
  //  std::cout << '\n';

  return 0;
}
