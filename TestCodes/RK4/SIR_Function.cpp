#include "RungeKutta4.hpp"

/*******************************************************************************
                              FUNCTION OBJECT
********************************************************************************
This object is used for definig functions with user defined parameters
*******************************************************************************/
std::vector<double> Function::Value(double t, std::vector<double> curr_state){
  std::vector<double> force;
  force.assign(curr_state.size(),0.0);
  // Definition of force equations
  force[0] = -Params[0]*curr_state[0]*curr_state[1];
  force[1] = Params[0]*curr_state[0]*curr_state[1]-Params[1]*curr_state[1];
  force[2] = Params[1]*curr_state[1];
  // Return Force
  return force;
}
