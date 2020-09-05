#include <iostream>
#include <vector>
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
class TimeSeries{

  private:
    double t_begin;                 // Beginning time of series
    double t_end;                   // Ending time of series
    int NumPoints;                  // Number of points in series
    int NumDOFs;                    // Number of Degrees of Freedom
    std::vector<double> DataVals;   // Values of series data

  public:
    // Constructor of the TimeSeries Class (TESTED)
    TimeSeries(double tBeg = 0.0, double tEnd = 1.0, int DOFs = 1, int NumP = 1000){
      t_begin = tBeg;
      t_end = tEnd;
      NumPoints = NumP;
      NumDOFs = DOFs;
      DataVals.assign(NumPoints*(NumDOFs+1),\
                  0.0);           // Need 1 extra slot for time !!!
    }
    // F**king Setters
    void set_t_begin(double tB){
      t_begin = tB;
    }
    void set_t_end(double tB){
      t_end = tB;
    }
    void set_NumPoints(int Num){
      NumPoints = Num;
    }
    void set_NumDOFs(int Num){
      NumDOFs = Num;
    }
    void set_DataVals(){
      DataVals.assign(NumPoints*(1+NumDOFs),0.0);
    }
    // F**king Getters
    double get_t_begin(){
      return t_begin;
    }
    double get_t_end(){
      return t_end;
    }
    int get_NumPoints(){
      return NumPoints;
    }
    int get_NumDOFs(){
      return NumDOFs;
    }
    std::vector<double> get_DataVals(){
      return DataVals;
    }
/******************************************************************************/
    // Inserting Element in specified position of TimeSeries (TESTED)
    void InsertDataPoint(double t, std::vector<double> d);
    // Retrieve Datapoint at specified position of TimeSeries (TESTED)
    std::vector<double> GetDataPoint(double t);
    // Print whole TimeSeries (TESTED)
    void PrintDataSeries();
    // Retrieve certain DOFs
    std::vector<double> GetDOFs(double t, std::vector<int> dofs);
};
/*******************************************************************************
                              ELEMENT CLASS
********************************************************************************
Used as a discrete element for a particular many-body system. It contains the
properties of the class that are determinant for the interaction between
elements of the system
*******************************************************************************/
class Element{

  private:
    // This properties change according to physical mmodelling
    double mass;
    double radius;
    int NumDOFs;
    // These following properties of class are mandatory
    std::vector<double> Coordinates;
    std::vector<double> Velocities;

  public:
    Element(double m = 1.0, double rad = 1.0, int DOFs = 3){
      mass = m;
      radius = rad;
      NumDOFs = DOFs;
      // Coordinates.assign(DOFs,0.0);
      // Velocities.assign(DOFs,0.0);
    }
    // F**king setters
    void set_mass(double m){
      mass = m;
    }
    void set_radius(double r){
      radius = r;
    }
    void set_NumDOFs(int dofs){
      NumDOFs = dofs;
    }
    void set_Coordinates(std::vector<double> c){
      if (c.size() == NumDOFs)
        Coordinates = c;
    }
    void set_Velocities(std::vector<double> v){
      if (v.size() == NumDOFs)
        Velocities = v;
    }
    // F**king getters
    double get_mass(){
      return mass;
    }
    double get_radius(){
      return radius;
    }
    int get_NumDOFs(){
      return NumDOFs;
    }
    std::vector<double> get_Coordinates(){
      return Coordinates;
    }
    std::vector<double> get_Velocities(){
      return Velocities;
    }
};
/*******************************************************************************
                                SYSTEM CLASS
********************************************************************************
Used to model a system of elements described by the same set of properties.
This class should contain pair interactions and external interactions that
model the accelaration of the elements. It should contain an integration step
that computes the time evolution of the system, and a way to record and print
the simulation data.
*******************************************************************************/
class System{

  private:
    int NumElems;
    int Elem_DOF;

  public:
    Element *ElementList;
    System(int numel = 1, int ElemDOFs = 3){
      // Initialize system parameters
      NumElems = numel;
      Elem_DOF = ElemDOFs;
      // Create array of elements
      ElementList = new Element[NumElems];
      // Initialize each elelement parameters
      std::vector<double> aux;
      aux.assign(Elem_DOF,0.0);
      for(int ii = 0; ii<NumElems; ii++){
        ElementList[ii].set_NumDOFs(Elem_DOF);
        ElementList[ii].set_Coordinates(aux);
        ElementList[ii].set_Velocities(aux);
      }
    }
    // Element interaction force
    std::vector<double> InteractionForce(int i, int j);
    // Net Force experienced by an element
    std::vector<double> TotalForce(int i);
    // Print current state of system
    void PrintCurrState();
    // Time Step evolution
    void StepEvolution(double dt);
    // Time Evolution for log register
    void Evolve(double t_begin, double t_end, double dt);
    // Return Memory
    void KillSystem();

};
