#include <iostream>
#include <cmath>
#include <vector>
/*******************************************************************************
                              VECTOR CLASS
********************************************************************************
This class was written by a programming god that really knew what he or she was
doing, not me. It was provided by Jose Daniel Munoz, teacher of the class.
*******************************************************************************/
using namespace std;
//---------------------- class vector3D --------------------
class vector3D{
private:
  double v[3];
 public:
  void   cargue(double x0, double y0, double z0);
  void   show(void);
  // Funciones de salida de componentes
  double x(void){return v[0];};
  double y(void){return v[1];};
  double z(void){return v[2];};
  //Lectura de Elementos
  double & operator[](int i){return v[i];};

  // Operaciones vectoriales
  vector3D    operator= (vector3D v2);
  vector3D    operator+ (vector3D v2);
  vector3D    operator+=(vector3D v2);
  vector3D    operator- (vector3D v2);
  vector3D    operator-=(vector3D v2);
  // Producto por escalar
  vector3D    operator* (double a);
  vector3D    operator*=(double a);
  friend  vector3D    operator* (double a,vector3D v1);
  // Division por escalar
  vector3D    operator/ (double a);
  // Producto cruz
  vector3D    operator^ (vector3D v2);
  // Producto punto
  double operator* (vector3D v2);
  // Norma
  friend  double norma2(vector3D v1);
  friend  double norma(vector3D v1);
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
    vector3D Coordinates;
    vector3D Velocities;
    vector3D Force;
    // Allows access to private properties & methods
    friend class System;

  public:
    Element(double m = 1.0, double rad = 1.0, int DOFs = 3){
      mass = m;
      radius = rad;
      NumDOFs = DOFs;
      Coordinates.cargue(0.0,0.0,0.0);
      Velocities.cargue(0.0,0.0,0.0);
      Force.cargue(0.0,0.0,0.0);
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
    void set_Coordinates(vector3D c){
      Coordinates = c;
    }
    void set_Velocities(vector3D v){
      Velocities = v;
    }
    void set_Force(vector3D f){
      Force = f;
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
    vector3D get_Coordinates(){
      return Coordinates;
    }
    vector3D get_Velocities(){
      return Velocities;
    }
    vector3D get_Force(){
      return Force;
    }
    void update_Coordinates(double dt, double coeff);
    void update_Velocities(double dt, double coeff);
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
    std::vector<double> IntegParams;

  public:
    Element *ElementList;
    System(int numel = 1, int ElemDOFs = 3){
      // Initialize system parameters
      NumElems = numel;
      Elem_DOF = ElemDOFs;
      // Create array of elements
      ElementList = new Element[NumElems];
      // Initialize each elelement parameters
      vector3D aux;
      aux.cargue(0.0,0.0,0.0);
      for(int ii = 0; ii<NumElems; ii++){
        ElementList[ii].set_NumDOFs(Elem_DOF);
        ElementList[ii].set_Coordinates(aux);
        ElementList[ii].set_Velocities(aux);
        ElementList[ii].set_Force(aux);
      }
      // Initailize integrator parameters
      IntegParams.assign(5,0.0);
      IntegParams[0] = 0.1786178958448091e00;                   // zeta
      IntegParams[1] = -0.2123418310626054e0;                   // lambda
      IntegParams[2] = -0.6626458266981849e-1;                  // chi
      IntegParams[3] = (1-2*IntegParams[1])/2.0;                // Coeff1
      IntegParams[4] = 1-2*(IntegParams[2] + IntegParams[0]);   // Coeff2
    }
    // Element interaction force
    void InteractionForce(int i, int j);
    void ConstraintForce(Element &ball);
    // Net Force experienced by an element
    void SetInternalForce(int i);
    void SetTotalForce(int i);
    // Print current state of system
    void PrintCurrState();
    // Time Step evolution
    void update_all_Coordinates(double dt, double param);
    void update_all_Forces();
    void update_all_Velocities(double dt, double param);
    void StepEvolution(double dt);
    // Time Evolution for log register
    void Evolve(double t_begin, double t_end, double dt);
    // Return Memory
    void KillSystem();

};
