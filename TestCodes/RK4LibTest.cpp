#include "RungeKutta4.hpp"
/*******************************************************************************
                            FUNCTION TO INTEGRATE
********************************************************************************
This function characterizes the dinamical system. It recives as input an
scalar that coresponds to evolution parameter, and the current state of
the system (characterized by a vector). It returns de rate of chage of The
state of the system.
*******************************************************************************/

std::vector<double> EvolFunction(double t, std::vector<double> curr_state){
  std::vector<double> force;
  force.assign(curr_state.size(),0.0);
  // Definition of force equations
  force[1] = -(2.0*M_PI/5.0)*(2.0*M_PI/5.0)*curr_state[0];
  force[0] = curr_state[1];
  // Return Force
  return force;
}

/******************************************************************************/

int main(void){

  double t_begin = 0.0;                 // Params of integrator
  double t_end = 5.0;
  double dt = (t_end-t_begin)/10000.0;
  int dim = 2;

  std::vector<double> y;                // Vector of initial conditions
  y.assign(dim,0.0);
  y[0] = 1.0; y[1] = 0.0;

  RK4_Solver MyDynamicSystem(t_begin,t_end,dt,dim);
  MyDynamicSystem.Integrate(EvolFunction,y);
  MyDynamicSystem.PrintData();

  return 0;
}
