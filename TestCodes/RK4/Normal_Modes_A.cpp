#include "RungeKutta4.hpp"

int main(void){

  double t_begin = 0.005;                 // Params of integrator
  double t_end = 10.0;
  double dt = (t_end-t_begin)/10000.0;
  int dim = 2;

  std::vector<double> y;                // Vector of initial conditions
  y.assign(dim,0.0);
  y[0] = 1.0; y[1] = 0.0;

  std::vector<double> Parameters;       // Parameters for SIR Model
  Parameters.assign(1,0.0);
  Parameters[0] = 1.0;

  Function Drum_Model;
  Drum_Model.set_Params(Parameters);

  RK4_Solver MyDynamicSystem(t_begin,t_end,dt,dim);
  MyDynamicSystem.Integrate(Drum_Model,y);
  MyDynamicSystem.PrintData();

  return 0;
}
