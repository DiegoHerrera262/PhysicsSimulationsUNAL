#include "DiscreteElements.hpp"
#include <cmath>
#include <algorithm>

const double Kwalls = 5e3;
const double Kballs = 1e8;
const double lx = 100;
const double ly = 100;
/******************************************************************************/
/*           THIS IS THE ONLY PART THAT SHOULD BE ACCESSIBLE TO USER          */
/******************************************************************************
                Interaction Force of the System (VERY IMPORTANT)
*******************************************************************************/
std::vector<double> System::InteractionForce(int i, int j){
  // Model Gravitational Attraction with elastic repulsion
  std::vector<double> Rij;              // Relative vector of two elements
  Rij = ElementList[i].Coordinates + (-1.0) * ElementList[j].Coordinates;
  double d = ElementList[i].radius + ElementList[j].radius;
  d = d - std::sqrt(norm_squared(Rij));
  d = std::max(0.0,d);
  double factor = pow(norm_squared(Rij),-0.5) * Kballs * pow(d,1.5);
  return factor * Rij;
}
/******************************************************************************
                Constraint Force of the System (VERY IMPORTANT)
*******************************************************************************/
// Computation of total force over element i
std::vector<double> System::ConstraintForce(int i){
  std::vector<double> Force;
  Force.assign(Elem_DOF,0.0);
  // Compute distances to walls
  double leftwallDist = ElementList[i].Coordinates[0];
  double rightwallDist = lx - ElementList[i].Coordinates[0];
  double bottomwallDist = ElementList[i].Coordinates[1];
  double topwallDist = ly - ElementList[i].Coordinates[1];
  // Compare if needed to correct
  leftwallDist = std::min(0.0,leftwallDist);
  rightwallDist = std::min(0.0,rightwallDist);
  bottomwallDist = std::min(0.0,bottomwallDist);
  topwallDist = std::min(0.0,topwallDist);
  // Compute force on x
  Force[0] = -Kwalls * (leftwallDist - rightwallDist);
  Force[1] = Kwalls * (topwallDist - bottomwallDist);
  return Force;
}
