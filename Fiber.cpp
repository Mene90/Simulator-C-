#include <cmath>

class Fiber {
double length, alphadB, dispersion, slope, nlindex, aeff;
double b2, b3, gm, alphalin;
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
