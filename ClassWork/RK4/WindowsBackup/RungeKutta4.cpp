#include "RungeKutta4.hpp"
/*******************************************************************************
       OVERLOAD OF ADDITION & SCALAR MULTIPLICATION FOR SIMPLICITY
********************************************************************************
Overload of vector addition to simplify sintax
Template is used to define a function for all types of vectors
The arguments are passed by reference to avoid overuse of memory
throw is used to produce an exception in case of inconsistency

Overloaded operators:
- += for simplified addition
- + for addition of 2 vectors
- * for scalar multiplication
- * for direct multiplication
*******************************************************************************/
// Overload of += operator for vectors
template<typename T>
std::vector<T> operator +=(std::vector<T> &lhs, const std::vector<T> &rhs) {
    // Check if the vectors to add are of the same size
    if (lhs.size() != rhs.size())
      throw std::length_error("vectors must be same size to add");
    // Perform simple addition
    for(int i = 0; i<lhs.size(); i++)
      lhs[i] += rhs[i];
    return lhs;
}
// Overload of * operator for vectors and scalars
template<typename T>
std::vector<T> operator *(const T &scalar, std::vector<T> v){
  for(int i = 0; i<v.size(); i++)
    v[i] *= scalar;
  return v;
}
// Overload of * operator for direct product of vectors
template<typename T>
std::vector<T> operator *(const std::vector<T> &w, std::vector<T> v){
  // Check if the vectors to add are of the same size
  if (w.size() != v.size())
    throw std::length_error("vectors must be same size to add");
  // Perform direct vector multiplication
  for(int i = 0; i<v.size(); i++)
    v[i] *= w[i];
  return v;
}
// Overload of addition for vectors
template<typename T>
std::vector<T> operator +(std::vector<T> lhs, const std::vector<T> &rhs) {
    if (lhs.size() != rhs.size())
        throw std::length_error("vectors must be same size to add");
    return lhs += rhs;
}
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

/*******************************************************************************
                      RUNGE - KUTTA 4TH ORDER INTEGRATOR
********************************************************************************
Model for solver of ODEs using Runge Kutta 4 algorithm. It has as parameters:
- Beginning time of simulation
- Ending time of simulation
- Time step of simulation
- Number of DOFs
- TimeSeries to store integration data
*******************************************************************************/
// Integration and Storage
void RK4_Solver::Integrate(Function &DynEqns, std::vector<double> y0){

  // Parameters of RK4 Integration
  double t0 = t_begin;            // Initial time of simulation
  double tf = t_end;              // Final time of simulation
  double dt = time_step;          // Time step of integration

  std::vector<double> y;          // Vector of initial conditions
  y = y0;

  // Auxiliar variables for integration
  std::vector<double> k1;         // Stores f(t,y)
  k1.assign(DegOfFreed,0.0);
  std::vector<double> k2;         // Stores 1st correction
  k2.assign(DegOfFreed,0.0);
  std::vector<double> k3;         // Stores 2nd correction
  k3.assign(DegOfFreed,0.0);
  std::vector<double> k4;         // Stores 3rd correction
  k4.assign(DegOfFreed,0.0);

  for(double t = t0; t <= tf; t += dt){
    // Store DataPoint in TimeSeries
    SimulationData.InsertDataPoint(t,y);

    // Update auxiliar kvecs
    k1 = DynEqns.Value(t,y);
    k2 = DynEqns.Value(t+0.5*dt,y+0.5*dt*k1);
    k3 = DynEqns.Value(t+0.5*dt,y+0.5*dt*k2);
    k4 = DynEqns.Value(t+dt,y+dt*k3);

    // Update State of the system
    y = y + (1.0/6.0 * dt) * (k1 + 2.0 * (k2 + k3) + k4);
  }
}
// Wrapper for returning data drom solver
void RK4_Solver::PrintData(){
  SimulationData.PrintDataSeries();
}
