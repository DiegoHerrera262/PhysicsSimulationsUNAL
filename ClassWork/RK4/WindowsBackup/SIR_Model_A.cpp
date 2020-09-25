#include "RungeKutta4.hpp"
/*******************************************************************************
                            FUNCTION TO INTEGRATE
********************************************************************************
This function characterizes the dinamical system. It recives as input an
scalar that coresponds to evolution parameter, and the current state of
the system (characterized by a vector). It returns de rate of chage of The
state of the system.
*******************************************************************************/

std::vector<double> SIR_Eqns(double t, std::vector<double> curr_state){
  std::vector<double> force;
  force.assign(curr_state.size(),0.0);
  // Definition of force equations
  force[0] = -0.35*curr_state[0]*curr_state[1];
  force[1] = 0.35*curr_state[0]*curr_state[1]-0.08*curr_state[1];
  force[2] = 0.08*curr_state[1];
  // Return Force
  return force;
}

/******************************************************************************/

int main(void){

  double t_begin = 0.0;                 // Params of integrator
  double t_end = 80.0;
  double dt = (t_end-t_begin)/10000.0;
  int dim = 3;

  std::vector<double> y;                // Vector of initial conditions
  y.assign(dim,0.0);
  y[0] = 0.999; y[1] = 0.001; y[2] = 0.0;

  std::vector<double> Parameters;       // Parameters for SIR Model
  Parameters.assign(2,0.0);
  Parameters[0] = 0.35; Parameters[1] = 0.08;

  Function SIR_Model;
  SIR_Model.set_Params(Parameters);

  RK4_Solver MyDynamicSystem(t_begin,t_end,dt,dim);
  MyDynamicSystem.Integrate(SIR_Model,y);
  MyDynamicSystem.PrintData();

  return 0;
}
