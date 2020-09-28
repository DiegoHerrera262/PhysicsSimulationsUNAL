#include "DiscreteElements.hpp"
/*******************************************************************************
                              VECTOR CLASS
********************************************************************************
This class was written by a programming god that really knew what he or she was
doing, not me. It was provided by Jose Daniel Munoz, teacher of the class.
*******************************************************************************/
void vector3D::cargue(double x0, double y0, double z0){
  v[0]=x0; v[1]=y0; v[2]=z0;
}
void vector3D::show(void){
  cout << "(" <<v[0]<< "," <<v[1]<< "," <<v[2]<< ")" << endl;
}
vector3D vector3D::operator=(vector3D v2){
  for(int i=0;i<3;i++)
    v[i] = v2.v[i];
  return *this;
}
vector3D vector3D::operator+(vector3D v2){
  vector3D total;
  for(int i=0;i<3;i++)
    total.v[i] = v[i] + v2.v[i];
  return total;
}
vector3D vector3D::operator+=(vector3D v2){
  *this = *this + v2;
  return *this;
}
vector3D vector3D::operator*(double a){
  vector3D total;
  for(int i=0;i<3;i++)
    total.v[i] = a*v[i];
  return total;
}
vector3D vector3D::operator*=(double a){
  *this = (*this)*a;
  return *this;
}
vector3D vector3D::operator/(double a){
  double inver = 1.0/a;
  vector3D total;
  for(int i=0;i<3;i++)
    total.v[i] = inver*v[i];
  return total;
}
vector3D vector3D::operator-(vector3D v2){
  return *this + v2*(-1);
}
vector3D vector3D::operator-=(vector3D v2){
  *this = *this - v2;
  return *this;
}
double vector3D::operator*(vector3D v2){
  double p=0;
  for(int i=0;i<3;i++)
    p += v[i]*v2.v[i];
  return p;
}
vector3D vector3D::operator^(vector3D v2){
  vector3D c;
  c.v[0] = v[1]*v2.v[2]-v[2]*v2.v[1];
  c.v[1] = v[2]*v2.v[0]-v[0]*v2.v[2];
  c.v[2] = v[0]*v2.v[1]-v[1]*v2.v[0];
  return c;
}
vector3D operator*(double a,vector3D v1){
  vector3D total;
  total = v1*a;
  return total;
}
double norma2(vector3D v1){
  double n=0;
  for(int i=0;i<3;i++)
    n += v1.v[i]*v1.v[i];
  return n;
}
double norma(vector3D v1){
  return sqrt(norma2(v1));
}
/*******************************************************************************
                              ELEMENT CLASS
********************************************************************************
Used as a discrete element for a particular many-body system. It contains the
properties of the class that are determinant for the interaction between
elements of the system
*******************************************************************************/
void Element::update_Coordinates(double dt, double coeff){
  Coordinates += Velocities * (coeff * dt);
}
void Element::update_Velocities(double dt, double coeff){
  Velocities += Force * (coeff * dt/mass);
}
/*******************************************************************************
                                SYSTEM CLASS
********************************************************************************
Used to model a system of elements described by the same set of properties.
This class should contain pair interactions and external interactions that
model the accelaration of the elements. It should contain an integration step
that computes the time evolution of the system, and a way to record and print
the simulation data.
*******************************************************************************/
// More efficient computation of pairwise interaction
void System::SetInternalForce(int i){
  for(int j = i+1; j < NumElems; j++){
    InteractionForce(i,j);
  }
}
// Updates interaction forces pairwise and applies constraint forces if needed
void System::SetTotalForce(int i){
  // Update pairwise interaction force (This is fixed for all simlations)
  SetInternalForce(i);
  // Include constraints if necessary here
  ConstraintForce(ElementList[i]);
}
// Print Current state of system
void System::PrintCurrState(){
  std::vector<double> data;
  for(int i = 0; i < NumElems; i++){
    for(int j = 0; j < Elem_DOF; j++)
      cout << ElementList[i].Coordinates[j] << "\t";
  }
  cout << endl;
}
// Perform Step Evolution of System Using Velocity Verlet
void System::update_all_Coordinates(double dt, double param){
  // Update position (FOR ALL PARTICLES)
  for(int i = 0; i<NumElems; i++)
    ElementList[i].update_Coordinates(dt,param);
}
void System::update_all_Forces(){
  // Update forces
  for(int i = 0; i<NumElems; i++)
    SetTotalForce(i);
}
void System::update_all_Velocities(double dt, double param){
  update_all_Forces();
  // Update velocities with coeff1 (FOR ALL PARTICLES)
  for(int i = 0; i<NumElems; i++){
    ElementList[i].update_Velocities(dt,param);
    // Reset Forces
    ElementList[i].Force.cargue(0.0,0.0,0.0);
  }
}
void System::StepEvolution(double dt){

  update_all_Coordinates(dt,IntegParams[0]);
  update_all_Velocities(dt,IntegParams[3]);
  update_all_Coordinates(dt,IntegParams[2]);
  update_all_Velocities(dt,IntegParams[1]);
  update_all_Coordinates(dt,IntegParams[4]);
  update_all_Velocities(dt,IntegParams[1]);
  update_all_Coordinates(dt,IntegParams[2]);
  update_all_Velocities(dt,IntegParams[3]);
  update_all_Coordinates(dt,IntegParams[0]);

}
// Perform Time evolution
void System::Evolve(double t_begin, double t_end, double dt){
  for(double t = t_begin; t <= t_end; t += dt){
    std::cout << t << '\t';
    PrintCurrState();
    StepEvolution(dt);
  }
}
// Return Memory
void System::KillSystem(){
  delete [] ElementList;
}
