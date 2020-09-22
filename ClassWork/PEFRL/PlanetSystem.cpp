#include "DiscreteElements.hpp"
#include <cmath>

int main(void){

  // Simulation of Only 1 Planet orbiting a Force center
  int DOFs = 2;
  System PlanetSystem(2,DOFs);

  // Condition for elliptical orbits r0 v0^2 < 2*100
  // Condition for circular orbits r0 v0^2 = 100
  // Period for circular orbits: T = 2*pi*r0/v0

  // Definition of parameters for 2 planet system
  double m0 = 16.0; double m1 = 2*4.0;
  double r0 = 1.0; double r1 = -m0/m1 * r0;
  double v0 = std::sqrt(m1*m1 * std::pow((m1 + m0),-2.0) * (m1/r0));
  double v1 = -std::sqrt(m0*m0 * std::pow((m1 + m0),-2.0) * (m1/r0));

  // Initial conditions
  std::vector<double> x_init;
  x_init.assign(DOFs,0.0);
  std::vector<double> v_init;
  v_init.assign(DOFs,0.0);

  // Set up parameters for large planet
  x_init[0] = r0;
  v_init[1] = v0;
  PlanetSystem.ElementList[0].set_mass(m0);
  PlanetSystem.ElementList[0].set_Coordinates(x_init);
  PlanetSystem.ElementList[0].set_Velocities(v_init);

  // Set up parameters for small planet
  x_init[0] = r1;
  v_init[1] = v1;
  PlanetSystem.ElementList[1].set_mass(m1);
  PlanetSystem.ElementList[1].set_Coordinates(x_init);
  PlanetSystem.ElementList[1].set_Velocities(v_init);

  // Perform Simulation
  PlanetSystem.Evolve(0.0,2.0*M_PI*r0/v0,0.0001);

  // Return memory
  PlanetSystem.KillSystem();

  return 0;
}
