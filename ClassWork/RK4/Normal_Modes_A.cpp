#include "Normal_Modes.hpp"

int main(void){

  double lambda = 1.0;
  double rinf = 0.005;
  double rsup = 10.0;
  Print_Drum_Bessel(lambda,rinf,rsup);

  return 0;
}
