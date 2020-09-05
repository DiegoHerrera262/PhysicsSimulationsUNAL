#include "DiscreteElements.hpp"
#include <cmath>

int main(void){

  /*
  Element *syst;
  syst = new Element[2];

  syst[0].set_mass(3.0);
  syst[1].set_mass(4.0);

  std::cout << syst[0].get_mass() << std::endl;
  std::cout << syst[1].get_mass() << std::endl;
  */

  System PlanetSystem(2,2);

  // Parameters of the sun
  PlanetSystem.ElementList[0].set_mass(1.0);
  // Parameters of orbiting planet
  PlanetSystem.ElementList[1].set_mass(1.0);

  // Set up initial conditions
  std::vector<double> r_init;
  r_init.assign(2,0.0);
  std::vector<double> v_init;
  v_init.assign(2,0.0);

  double x0 = 5.0;
  double y0 = 0.0;
  r_init[0] = x0; r_init[1] = y0;

  double Vx0 = 0.0;
  double Vy0 = 0.0;
  v_init[0] = Vx0; v_init[1] = Vy0;

  // Set up initial position
  PlanetSystem.ElementList[1].set_Coordinates(r_init);
  PlanetSystem.ElementList[1].set_Velocities(v_init);
  PlanetSystem.ElementList[0].set_Coordinates(-1.0*r_init);
  PlanetSystem.ElementList[0].set_Velocities(v_init);

  // Perform Simulation
  PlanetSystem.Evolve(0.0,1.0,0.0001);

  // Return memory
  PlanetSystem.KillSystem();

  return 0;
}
