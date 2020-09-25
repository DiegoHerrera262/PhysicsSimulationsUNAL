#include "Normal_Modes.hpp"

int main(void){

  double lmin = 0.1;
  double lmax = 15.0;
  double dl = (lmax-lmin)/500.0;

  for(double l = lmin; l <= lmax; l+=dl)
    std::cout << l << "\t" << R_at_One(l) << std::endl;

  return 0;
}
