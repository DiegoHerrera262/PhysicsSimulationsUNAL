#include "DiscreteElements.hpp"
#include <cmath>
#include <algorithm>

const double g = 1.0;
const double l = 2.0;
const double separation = 1.0;
const double Kstring = 1e5;
const double Kballs = 1e8;

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
std::vector<double> PendullumForce(\
  std::vector<double> &pinCoord, Element &Ball){
  // Read coordinate of
  std::vector<double> Rpin = Ball.get_Coordinates();
  if(pinCoord.size() != Rpin.size())
    throw std::length_error("pinCoord is not of right size");
  Rpin += (-1.0) * pinCoord;             // Compute rad vec to pinpoint
  // Compute constraint force using linear approximation
  double d = std::sqrt(norm_squared(Rpin));
  d = -Kstring * (1.0 - l/d);
  Rpin = d * Rpin;
  // Distance to pinpoint squared
  std::vector<double> mg;               // Compute gravitational force
  mg.assign(Rpin.size(),0.0);
  mg[1] = -Ball.get_mass() * g;

  /*
  mg[1] = -Ball.get_mass() * g;
  if(Rc2 - l*l >= 0.0){
    Force = std::pow(Rc2,-0.5) * Force;   // Normalize rad vec
    double centrifuc = norm_squared(Ball.get_Velocities());
    // Compute centrifugal component
    centrifuc *= Ball.get_mass() * std::pow(Rc2,-0.5);
    Force = (-dot(Force,mg) + centrifuc) * Force;
    return Force + mg;
  }
  else
    return mg;              // STRING IS NOT TENSE !!!
  */
  return mg + Rpin;
}
// Computation of total force over element i
void System::SetTotalForce(int i){
  std::vector<double> pinCoord;
  pinCoord.assign(Elem_DOF,0.0);
  pinCoord[0] = i * separation;
  std::vector<double> Force;
  Force = PendullumForce(pinCoord,ElementList[i]);
  for(int j = 0; j<NumElems; j++)
    if(j != i)
      Force += InteractionForce(i,j);

  ElementList[i].Force = Force;
}
