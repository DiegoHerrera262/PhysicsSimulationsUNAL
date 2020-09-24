#include "DiscreteElements.hpp"
#include <cmath>

/******************************************************************************/
/*           THIS IS THE ONLY PART THAT SHOULD BE ACCESSIBLE TO USER          */
/******************************************************************************
                Interaction Force of the System (VERY IMPORTANT)
*******************************************************************************/
std::vector<double> System::InteractionForce(int i, int j){
  // Model Gravitational Attraction with elastic repulsion
  std::vector<double> Rij;              // Relative vector of two elements
  Rij = ElementList[i].Coordinates + (-1.0) * ElementList[j].Coordinates;
  double factor = -ElementList[i].mass * ElementList[j].mass;
  factor *= pow(norm_squared(Rij),-1.5);
  return factor * Rij;
}
/******************************************************************************
                Constraint Force of the System (VERY IMPORTANT)
*******************************************************************************/
std::vector<double> System::ConstraintForce(int i){
  // Assuming gravitational interaction on free space
  std::vector<double> Force;
  Force.assign(Elem_DOF,0.0);
  return Force;
}
