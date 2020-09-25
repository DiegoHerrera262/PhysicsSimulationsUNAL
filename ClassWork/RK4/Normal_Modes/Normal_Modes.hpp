#include "../RungeKutta4.hpp"

/*******************************************************************************
                            FOR PRINTING BESSEL
********************************************************************************
Recieves eigenvalue as parameter
*******************************************************************************/
void Print_Drum_Bessel(double lambda, double r_inf, double r_sup){
  double t_begin = r_inf;                 // Params of integrator
  double t_end = r_sup;
  double dt = (t_end-t_begin)/10000.0;
  int dim = 2;

  std::vector<double> y;                // Vector of initial conditions
  y.assign(dim,0.0);
  y[0] = 1.0; y[1] = 0.0;               // Free end at center

  std::vector<double> Parameters;       // Parameters for SIR Model
  Parameters.assign(1,0.0);
  Parameters[0] = lambda;

  Function Drum_Model;
  Drum_Model.set_Params(Parameters);

  RK4_Solver MyDynamicSystem(t_begin,t_end,dt,dim);
  MyDynamicSystem.Integrate(Drum_Model,y);
  MyDynamicSystem.PrintData();
}

/*******************************************************************************
                      FOR IMPOSING BOUNDARY CONDITIONS
********************************************************************************
Boundary fixed at r = 1.0, Boundary free at r = 0.
*******************************************************************************/
double R_at_One(double lambda){
  double t_begin = 0.005;                 // Params of integrator
  double t_end = 1.0;
  double dt = (t_end-t_begin)/10000.0;
  int dim = 2;

  std::vector<double> y;                // Vector of initial conditions
  y.assign(dim,0.0);
  y[0] = 1.0; y[1] = 0.0;               // Free end at center

  std::vector<double> Parameters;       // Parameters for SIR Model
  Parameters.assign(1,0.0);
  Parameters[0] = lambda;

  Function Drum_Model;
  Drum_Model.set_Params(Parameters);

  RK4_Solver Bessel_Equation(t_begin,t_end,dt,dim);
  Bessel_Equation.Integrate(Drum_Model,y);

  std::vector<double> StateAtEnd;
  StateAtEnd = Bessel_Equation.SimulationData.GetDataPoint(0.9999999);

  return StateAtEnd[1];
}
/*******************************************************************************
                              FOR ROOT FINDING
********************************************************************************
Implemented Regula Falsi
*******************************************************************************/
double RegulaFalsi(double (*f)(double), double x_inf, double x_sup){
  double xi = x_inf; double xs = x_sup;
  double fi = f(xi); double fs = f(xs);
  double xr = 1.0; double fr = 1.0;
  do {
    if(fi*fs < 0.0){
      xr = (fs*xi - fi*xs)/(fs-fi);
      fr = f(xr);
      if(fr*fs > 0.0){
        xs = xr; fs = fr;
      }
      else if(fr*fi > 0.0){
        xi = xr; fi = fr;
      }
      else{
        std::cout << "Problem computing real numbers" << std::endl;
        break;
      }
    }
    else{
      std::cout << "Chose bad interval" << std::endl;
      break;
    }
  } while(std::abs(fr) > 1e-8);
  return xr;
}
