#include <complex>

using namespace std;

class NLSE_solver{
  int nstep;
public:
  NLSE_solver(int nstep)
              : nstep{nstep} {}
  void ssfm (std::complex<double>& x){}
}
