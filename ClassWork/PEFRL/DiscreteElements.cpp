#include "DiscreteElements.hpp"
#include <cmath>

/*******************************************************************************
                            TIME SERIES CLASS
********************************************************************************
Used for storing integration data. It has a beginning time, an ending time,
an std::vector<double> for storing data and dimension. It must have methods for:
- Inserting a data vector (of size dimension) at specific time.
- Retrive data vector stored at specific time.
- Print chosen DOFs as function of time.
- Print all data series.
Data is stored in a linear array, and datapoints are accessed modulo DOFs.
*******************************************************************************/
// Inserting Element in specified position of TimeSeries (TESTED)
void TimeSeries::InsertDataPoint(double t, std::vector<double> d){
  int posData = int((NumPoints-1) * (t-t_begin)/(t_end-t_begin));
  if (posData <= NumPoints-1){
    if (d.size() == NumDOFs){
      double idx = posData * (NumDOFs + 1);
      DataVals[idx] = t;
      for(int i = 0; i<d.size(); i++)
        DataVals[idx+i+1] = d[i];
      }
      else
        throw std::length_error("Data inconsistent with DOFs");
  }
  else
    throw std::length_error("Point Outside Time Range");
}
// Retrieve Datapoint at specified position of TimeSeries (TESTED)
std::vector<double> TimeSeries::GetDataPoint(double t){
  int posData = int((NumPoints-1) * (t-t_begin)/(t_end-t_begin));
  if (posData <= NumPoints-1){
    std::vector<double> data;
    data.assign(NumDOFs+1,0.0);
    double idx = posData * (NumDOFs + 1);
    for(int i = 0; i<= NumDOFs; i++)
      data[i] = DataVals[i+idx];
    return data;
  }
  else
    throw std::length_error("Point Outside Time Range");
}
// Print whole TimeSeries (TESTED)
void TimeSeries::PrintDataSeries(){
  int idx = 0;
  for(int i=0; i < NumPoints; i++){
    idx = i*(NumDOFs + 1);
    if(DataVals[idx] != 0){
      std::cout << "Pt" << i << "\t";
      for(int j = 0; j <= NumDOFs; j++)
        std::cout << DataVals[idx + j] << "\t";
      std::cout << std::endl;
    }
  }
}
// Retrieve certain DOFs
std::vector<double> TimeSeries::GetDOFs(double t, std::vector<int> dofs){
  int posData = int((NumPoints-1) * (t-t_begin)/(t_end-t_begin));
  bool dofs_in_range = true;
  for(int i = 0; i < dofs.size(); i++)
    dofs_in_range = dofs_in_range && (dofs[i] <= NumDOFs);
  bool input_is_consistent = dofs_in_range && posData <= NumPoints-1;
  if (input_is_consistent){
    std::vector<double> data;
    data.assign(dofs.size(),0.0);
    double idx = posData * (NumDOFs + 1);
    for(int i = 0; i < data.size(); i++)
      data[i] = DataVals[dofs[i]+idx];
    return data;
  }
  else
    throw std::length_error("Point Outside Time Range");
}
/*******************************************************************************
                              ELEMENT CLASS
********************************************************************************
Used as a discrete element for a particular many-body system. It contains the
properties of the class that are determinant for the interaction between
elements of the system
*******************************************************************************/
void Element::update_Coordinates(double dt, double coeff){
  Coordinates += (coeff * dt) * Velocities;
}
void Element::update_Velocities(double dt, double coeff){
  Velocities += (coeff * dt/mass) * Force;
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
/*           THIS IS THE ONLY PART THAT SHOULD BE ACCESSIBLE TO USER          */
/******************************************************************************
                Interaction Force of the System (VERY IMPORTANT)
*******************************************************************************/
/*
std::vector<double> System::InteractionForce(int i, int j){
  // Model Gravitational Attraction with elastic repulsion
  std::vector<double> Rij;              // Relative vector of two elements
  Rij = ElementList[i].Coordinates + (-1.0) * ElementList[j].Coordinates;
  double factor = -ElementList[i].mass * ElementList[j].mass;
  factor *= pow(norm_squared(Rij),-1.5);
  return factor * Rij;
}
*/
// Computation of total force over element i
/*
void System::SetTotalForce(int i){

  std::vector<double> Force;
  Force.assign(Elem_DOF,0.0);
  for(int j = 0; j<NumElems; j++)
    if(j != i)
      Force += InteractionForce(i,j);

  ElementList[i].Force = Force;
}
*/
// Print Current state of system
void System::PrintCurrState(){
  std::vector<double> data;
  for(int i = 0; i < NumElems; i++){
    data = ElementList[i].get_Coordinates();
    for(int j = 0; j < Elem_DOF; j++)
      std::cout << data[j] << "\t";
  }
  std::cout << std::endl;
}
// Perform Step Evolution of System Using Velocity Verlet
void System::update_all_Coordinates(double dt, double param){
  // Update position (FOR ALL PARTICLES)
  for(int i = 0; i<NumElems; i++)
    ElementList[i].update_Coordinates(dt,param);
  // Update forces
  for(int i = 0; i<NumElems; i++)
    SetTotalForce(i);
}
void System::update_all_Velocities(double dt, double param){
  // Update velocities with coeff1 (FOR ALL PARTICLES)
  for(int i = 0; i<NumElems; i++)
    ElementList[i].update_Velocities(dt,param);
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
