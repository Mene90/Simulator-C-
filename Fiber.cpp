#include <cmath>



class Fiber {

double length, alphadB, dispersion, slope, nlindex, aeff;
double b2, b3, gm, alphalin;

typedef std::complex  <double>   Complex;
typedef std::valarray <Complex>  CArray;
typedef std::valarray <double>   DArray;

public:
Fiber(double length, double alphadB, double dispersion, double slope, double nlindex, double aeff)
      : length{length}
      , alphadB{alphadB}
      , dispersion{dispersion}
      , slope{slope}
      , nlindex{nlindex}
      , aeff{aeff}
    {alphalin = (log(10)*1e-4)*alphadB; }
}

void propagation (double gmp, CArray& beta, CArray& Ux, NLSE_solver num_method){

      dz      = length/num_method::nstep;
      geff01  = alphalin*dz > 1e-6 ? gmp*(1-exp(-alphalin*dz))/alphalin : gmp*dz;

      int current = 0;
      int UniqueNumeber () {return dz*current++;};
      int zarray[num_method::nstep];
      std::generate_n (zarray,num_method::nstep-1,UniqueNumeber);

      DArray z(zarray,num_method::nstep);
      geffz = geff01*exp(z*alphalin);

      num_method::solve(beta,geff,dz);

}
