#include "Fiber.hpp"

Fiber::Fiber(double length, double alphadB, double dispersion, double slope, double nlindex, double aeff)
      : length{length}
      , alphadB{alphadB}
      , dispersion{dispersion}
      , slope{slope}
      , nlindex{nlindex}
      , aeff{aeff}
    {alphalin = (log(10)*1e-4)*alphadB; }


void Fiber::propagation (CArray& ux, double gmp, CArray& beta, NLSE_solver num_method){
      int    nsteps  = num_method.getNsteps();
      double dz      = length/nsteps;
      double geff01  = alphalin*dz > 1e-6 ? gmp*(1-exp(-alphalin*dz))/alphalin : gmp*dz;
      double zarray[nsteps];

      std::generate_n (zarray,nsteps-1, [dz]()
        {
          int current = 0;
          return dz*current++;
        });

      DArray z(zarray,nsteps);

      DArray geffz = geff01*exp(z*alphalin);

      num_method.solve(ux,beta,geffz,dz);

}
