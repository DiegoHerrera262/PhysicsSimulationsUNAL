#include "DiscreteElements.hpp"
#include <algorithm>

const double Kwalls = 5e3;
const double Kballs = 1e8;
const double lx = 60;
const double ly = 60;
/******************************************************************************/
/*           THIS IS THE ONLY PART THAT SHOULD BE ACCESSIBLE TO USER          */
/******************************************************************************
                Interaction Force of the System (VERY IMPORTANT)
*******************************************************************************/
void System::InteractionForce(int i, int j){
  // Model Gravitational Attraction with elastic repulsion
  // Relative vector of two elements
  vector3D Rij = ElementList[i].Coordinates;  // Relative vector of two elements
  Rij -= ElementList[j].Coordinates;
  double d = ElementList[i].radius + ElementList[j].radius;
  double n = norma(Rij);
  d = d - n;
  if(d > 0){
    double factor = pow(n,-1.0) * Kballs * pow(d,1.5);
    ElementList[i].Force += factor * Rij;
    ElementList[j].Force -= ElementList[i].Force;
  }
}
/******************************************************************************
                Constraint Force of the System (VERY IMPORTANT)
*******************************************************************************/
// Computation of total force over element i
void System::ConstraintForce(Element &ball){
  // Compute distances to walls
  double leftwallDist = ball.Coordinates[0];
  double rightwallDist = lx - ball.Coordinates[0];
  double bottomwallDist = ball.Coordinates[1];
  double topwallDist = ly - ball.Coordinates[1];
  // Compare if needed to correct
  leftwallDist = std::min(0.0,leftwallDist);
  rightwallDist = std::min(0.0,rightwallDist);
  bottomwallDist = std::min(0.0,bottomwallDist);
  topwallDist = std::min(0.0,topwallDist);
  // Compute force on x
  ball.Force[0] += -Kwalls * (leftwallDist - rightwallDist);
  ball.Force[1] += Kwalls * (topwallDist - bottomwallDist);
}
