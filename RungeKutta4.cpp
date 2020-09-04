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

template<typename T>
std::vector<T> operator *(const T &scalar, std::vector<T> v){
  for(int i = 0; i<v.size(); i++)
    v[i] *= scalar;
  return v;
}

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
    TimeSeries(double tBeg = 0.0, double tEnd = 1.0, int DOFs = 1, \
            int NumP = 1000){
      t_begin = tBeg;
      t_end = tEnd;
      NumPoints = NumP;
      NumDOFs = DOFs;
      DataVals.assign(NumPoints*(NumDOFs+1),\
                  0.0);           // Need 1 extra slot for time !!!
    }
    // Inserting Element in specified position of TimeSeries (TESTED)
    void InsertDataPoint(double t, std::vector<double> d){
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
    std::vector<double> GetDataPoint(double t){
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
    void PrintDataSeries(){
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

};

/******************************************************************************/

/*******************************************************************************
                      RUNGE - KUTTA 4TH ORDER INTEGRATOR
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
    // Integration and Storage
    void Integrate(\
      std::vector<double> (*func)(double, std::vector<double>),\
      std::vector<double> y0){

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
          k1 = func(t,y);
          k2 = func(t+0.5*dt,y+0.5*dt*k1);
          k3 = func(t+0.5*dt,y+0.5*dt*k2);
          k4 = func(t+dt,y+dt*k3);

          // Update State of the system
          y = y + (1.0/6.0 * dt) * (k1 + 2.0 * (k2 + k3) + k4);
        }
    }
    // Wrapper for returning data drom solver
    void PrintData(){
      SimulationData.PrintDataSeries();
    }
};

/******************************************************************************/

void RK4_print(int dim, \
  std::vector<double> (*func)(double, std::vector<double>),\
  std::vector<double> params, std::vector<double> y0){

  // Parameters of RK4 Integration
  double t0 = params[0];                // Initial time of simulation
  double tf = params[1];                // Final time of simulation
  double dt = params[2];                // Time step of integration

  std::vector<double> y;                // Vector of initial conditions
  y = y0;

  // Auxiliar variables for integration
  std::vector<double> k1;               // Stores f(t,y)
  k1.assign(dim,0.0);
  std::vector<double> k2;               // Stores 1st correction
  k2.assign(dim,0.0);
  std::vector<double> k3;               // Stores 2nd correction
  k3.assign(dim,0.0);
  std::vector<double> k4;               // Stores 3rd correction
  k4.assign(dim,0.0);

  for(double t = t0; t <= tf; t += dt){
    // Print current time & state of system
    std::cout << t << "\t";
    for(double &el : y)
      std::cout << el << "\t";
    std::cout << "\n";

    // Update auxiliar kvecs
    k1 = func(t,y);
    k2 = func(t+0.5*dt,y+0.5*dt*k1);
    k3 = func(t+0.5*dt,y+0.5*dt*k2);
    k4 = func(t+dt,y+dt*k3);

    // Update State of the system
    y = y + (1.0/6.0 * dt) * (k1 + 2.0 * (k2 + k3) + k4);
  }
}

/******************************************************************************/
/*******************************************************************************
                            FUNCTION TO INTEGRATE
********************************************************************************
This function characterizes the dinamical system. It recives as input an
scalar that coresponds to evolution parameter, and the current state of
the system (characterized by a vector). It returns de rate of chage of The
state of the system.
*******************************************************************************/

std::vector<double> EvolFunction(double t, std::vector<double> curr_state){
  std::vector<double> force;
  force.assign(curr_state.size(),0.0);
  // Definition of force equations
  force[1] = -(2.0*M_PI/5.0)*(2.0*M_PI/5.0)*curr_state[0];
  force[0] = curr_state[1];
  // Return Force
  return force;
}

/******************************************************************************/

int main(void){


  double t_begin = 0.0;                 // Params of integrator
  double t_end = 5.0;
  double dt = (t_end-t_begin)/10000.0;
  int dim = 2;

  /*
  int NumPts = (int) (t_end - t_begin)/dt + 1;

  TimeSeries Test(t_begin,t_end,dim,NumPts);
  Test.PrintDataSeries();
  std::cout << "\n" << std::endl;

  std::vector<double> y; y.assign(dim,0.9);
  Test.InsertDataPoint(t_end,y);

  Test.PrintDataSeries();

  std::vector<double> point = Test.GetDataPoint(t_end);

  for(int i = 0; i<point.size(); i++)
    std::cout << point[i] << "\t";
  std::cout << "\n" << std::endl;
  */


  std::vector<double> y;                // Vector of initial conditions
  y.assign(dim,0.0);
  y[0] = 0.0; y[1] = 1.0;

  RK4_Solver MyDynamicSystem(t_begin,t_end,dt,dim);
  MyDynamicSystem.Integrate(EvolFunction,y);
  MyDynamicSystem.PrintData();

  /*
  int dim = 2;                          // Dimension of Integration problem

  // Parameters of RK4 Integration
  std::vector<double> params;
  params.assign(3,0.0);
  params[0] = 0.0;                      // Initial time of simulation
  params[1] = 5.0;                      // Final time of simulation
  params[2] = 0.01;                     // Time step of integration

  std::vector<double> y;                // Vector of initial conditions
  y.assign(dim,0.0);
  y[0] = 1.0; y[1] = 0.0;

  RK4_print(dim,EvolFunction,params,y);
  */

  return 0;
}
