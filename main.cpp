#include <cmath>
#include <random>
#include <complex>
#include <iostream>
#include <valarray>

#include "Fiber.hpp"

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


int main(int argc, char const *argv[]) {

  Fiber       fb(1e5,0.2,17,0,2.5e-20,80);
  NLSE_solver pro_sim(100);

  CArray ux = gaussGen(20,0.7071);

  // for(auto n : ux) {
  //      std::cout << n << '\n';
  //  }
  //  std::cout << '\n';

  return 0;
}
