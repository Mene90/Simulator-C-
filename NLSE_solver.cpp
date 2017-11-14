#include "NLSE_solver.hpp"

  NLSE_solver::NLSE_solver(int nstep)
              : nstep{nstep} {}

  void NLSE_solver::solve (CArray& ux, CArray& beta, DArray& geffz, double dz){
      ssfm(ux,beta,geffz,dz);
  }

  int NLSE_solver::getNsteps(){return nstep;}

  void NLSE_solver::ssfm (CArray& ux, CArray& beta, DArray& geffz, double dz){
    Complex i(0,1);
    CArray Fb  = exp(i*dz*beta);
    CArray Fhb = exp(i*dz*beta/2);

    ln_step(Fhb,ux);
    nl_step(geffz[0],ux);

    for (int k=1;k<nstep;k++){
      ln_step(Fb,ux);
      nl_step(geffz[k],ux);
    }

    ln_step(Fhb,ux);
  }

  void NLSE_solver::ln_step(CArray& Fb, CArray& ux){
    fft(ux);
    ux = ux * Fb;
    ifft(ux);
  }

  void NLSE_solver::nl_step(Complex geffz, CArray& ux){
    Complex i(0,1);
    ux = ux*exp(-1.0*i*geffz*pow(abs(ux),2));
  }
