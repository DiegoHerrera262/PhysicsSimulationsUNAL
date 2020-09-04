#include <iostream>
#include <vector>
#include <cmath>
#include <exception>
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
};

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

class RK4_Solver{

  private:
    double t_begin;
    double t_end;
    double time_step;
    double DegOfFreed;
    TimeSeries SimulationData;

  public:
    // Constructor of the RK4_Solver
    RK4_Solver(double tBeg = 0.0, double tEnd = 1.0,\
      double dt = 0.1, int DOFs = 1){
            t_begin = tBeg;
            t_end = tEnd;
            time_step = dt;
            int NumPts = int((tEnd - tBeg)/dt + 1);
            // Initialize SimulationData (like constructor)
            SimulationData.set_t_begin(tBeg);
            SimulationData.set_t_end(tEnd);
            SimulationData.set_NumDOFs(DOFs);
            SimulationData.set_NumPoints(NumPts);
            SimulationData.set_DataVals();

    }
    // F**king Getters
    double get_t_begin(){
      return t_begin;
    }
    double get_t_end(){
      return t_end;
    }
    double get_time_step(){
      return time_step;
    }
    double get_DegOfFreed(){
      return DegOfFreed;
    }
    TimeSeries get_SimulationData(){
      return SimulationData;
    }
    // F**king Setters
    void set_t_begin(double tB){
      t_begin = tB;
    }
    void set_t_end(double tE){
      t_end = tE;
    }
    void set_time_step(double dt){
      time_step = dt;
    }
    void set_DegOfFreed(double dofs){
      DegOfFreed = dofs;
    }
/******************************************************************************/
    // Integration and Storage
    void Integrate(std::vector<double> (*func)(double, std::vector<double>),std::vector<double> y0);
    // Wrapper for returning data drom solver
    void PrintData();

};
