#include "DiscreteElements.hpp"
#include <cmath>

// Definition of parameters for 2 planet system
const double m0 = 16.0; double m1 = 2*4.0;
const double r0 = 1.0; double r1 = -m0/m1 * r0;
const double v0 = std::sqrt(m1*m1 * std::pow((m1 + m0),-2.0) * (m1/r0));
const double v1 = -std::sqrt(m0*m0 * std::pow((m1 + m0),-2.0) * (m1/r0));

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
  double n = norma2(Rij);
  double factor = -pow(n,-1.5) * ElementList[i].mass * ElementList[j].mass;
  ElementList[i].Force += factor * Rij;
  ElementList[j].Force -= ElementList[i].Force;
}
/******************************************************************************
                Constraint Force of the System (VERY IMPORTANT)
*******************************************************************************/
// Computation of total force over element i
void System::ConstraintForce(Element &ball){

}
