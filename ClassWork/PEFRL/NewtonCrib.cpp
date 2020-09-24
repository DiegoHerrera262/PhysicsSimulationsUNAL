#include "DiscreteElements.hpp"
#include <cmath>

const double g = 1.0;
const double l = 2.0;
const double separation = 1.0;

int main(void){

  // Simulation a Newton Crib of NumBalls balls
  int NumBalls = 3;
  int DOFs = 2;
  System NewtonCrib(NumBalls,DOFs);

  // Condition for elliptical orbits r0 v0^2 < 2*100
  // Condition for circular orbits r0 v0^2 = 100
  // Period for circular orbits: T = 2*pi*r0/v0

  // Definition of parameters for 2 planet system
  double m0 = 16.0; double theta = 40;
  double x0 = -l*std::sin(theta/180 * M_PI);
  double y0 = -l*std::cos(theta/180.0 * M_PI);

  // Initial conditions
  std::vector<double> x_init;
  x_init.assign(DOFs,0.0);
  std::vector<double> v_init;
  v_init.assign(DOFs,0.0);

  // Set up parameters for leftmost ball
  x_init[0] = x0; x_init[1] = y0;
  NewtonCrib.ElementList[0].set_mass(1.0);
  NewtonCrib.ElementList[0].set_radius(0.5 * separation);
  NewtonCrib.ElementList[0].set_Coordinates(x_init);
  NewtonCrib.ElementList[0].set_Velocities(v_init);

  // Set up parameters for rest of balls
  for(int i = 1; i < NumBalls; i++){
    x_init[0] = i * separation; x_init[1] = - l;
    NewtonCrib.ElementList[i].set_mass(1.0);
    NewtonCrib.ElementList[i].set_radius(0.5 * separation);
    NewtonCrib.ElementList[i].set_Coordinates(x_init);
    NewtonCrib.ElementList[i].set_Velocities(v_init);
  }

  // Set up parameters for large planet
  /*
  x_init[0] = -x0; x_init[1] = y0;
  NewtonCrib.ElementList[1].set_mass(1.0);
  NewtonCrib.ElementList[0].set_Coordinates(x_init);
  NewtonCrib.ElementList[0].set_Velocities(v_init);
  */
  // Perform Simulation
  NewtonCrib.Evolve(0.0,10.0*M_PI*std::sqrt(2.0),0.001);

  // Return memory
  NewtonCrib.KillSystem();

  return 0;
}
