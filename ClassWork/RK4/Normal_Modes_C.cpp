#include "Normal_Modes.hpp"

double f(double x){
  return x*x - 2;
}

int main(void){
  std::cout.precision(10);
  double lambda = RegulaFalsi(R_at_One,14.0,16.0);
  double rinf = 0.005;
  double rsup = 1.0;
  Print_Drum_Bessel(lambda,rinf,rsup);
  return 0;
}
