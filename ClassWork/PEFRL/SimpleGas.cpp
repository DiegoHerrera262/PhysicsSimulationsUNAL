#include "DiscreteElements.hpp"
#include <cmath>
#include <random>

const double k = 1.0e-4;
const double lx = 10;
const double ly = 10;

int main(){

  // Simulation a 2D Gas of NumBalls balls
  int NumBalls = 2;
  int DOFs = 2;
  System Gas2D(NumBalls,DOFs);

  // Parameters for gas elements
  double m = 1.0; double r = 1.0; double kT = 10.0;
  double theta = 0.0; double v = std::sqrt(2.0*kT/m);
  std::vector<double> x_init; x_init.assign(DOFs,0.0);
  std::vector<double> v_init; v_init.assign(DOFs,0.0);

  // Prepare Random number generator for velocity angle
  // and initial coordinates
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> ang(0.0,2.0*M_PI);
  std::uniform_real_distribution<> x(r,lx-r);
  std::uniform_real_distribution<> y(r,ly-r);

  // Set up parameters of balls
  for(int i = 0; i < NumBalls; i++){
    x_init[0] = (i+1)*3*r; x_init[1] = (i+1)*3*r;
    theta = 14*M_PI/180.0;
    v_init[0] = v*std::cos(theta); v_init[1] = v*std::sin(theta);
    Gas2D.ElementList[i].set_mass(m);
    Gas2D.ElementList[i].set_radius(r);
    Gas2D.ElementList[i].set_Coordinates(x_init);
    Gas2D.ElementList[i].set_Velocities((-1.0)*v_init);
  }

  // Perform Simulation
  Gas2D.Evolve(0.0,12*lx/v,0.001);

  // Return memory
  Gas2D.KillSystem();

  return 0;
}
