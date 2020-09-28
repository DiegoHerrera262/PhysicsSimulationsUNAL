#include "DiscreteElements.hpp"
#include <random>

const double lx = 60;
const double ly = 60;

int main(){

  // Simulation a 2D Gas of NumBalls balls
  int Nx = 9;
  int Ny = 9;
  int NumBalls = Nx * Ny;
  int DOFs = 2;
  System Gas2D(NumBalls,DOFs);
  double dx =  dx=lx/(Nx+1), dy=ly/(Ny+1);
  // Parameters for gas elements
  double m = 1.0; double r = 1.0; double kT = 10.0;
  double theta = 0.0; double v = std::sqrt(2.0*kT/m);
  vector3D x_init; x_init.cargue(0.0,0.0,0.0);
  vector3D v_init; v_init.cargue(0.0,0.0,0.0);
  // Prepare Random number generator for velocity angle
  // and initial coordinates
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> ang(0.0,2.0*M_PI);

  // Set up parameters of balls
  for(int i = 0; i < Nx; i++){
    for(int j = 0; j < Ny; j++){
      x_init[0] = (i+1)*dx; x_init[1] = (j+1)*dx;
      theta = ang(gen);
      v_init[0] = v*std::cos(theta); v_init[1] = v*std::sin(theta);
      Gas2D.ElementList[Nx*j+i].set_mass(m);
      Gas2D.ElementList[Nx*j+i].set_radius(r);
      Gas2D.ElementList[Nx*j+i].set_Coordinates(x_init);
      Gas2D.ElementList[Nx*j+i].set_Velocities(v_init);
    }
  }

  // Perform Simulation
  Gas2D.Evolve(0.0,10*lx/v,0.001);

  // Return memory
  Gas2D.KillSystem();

  return 0;
}
