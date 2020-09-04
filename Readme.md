# Welcome to PhysicsSimulationsUNAL

This repo is set up to be a tool for those physics students aiming at learning how to model Physical Systems using techniques such as:

* ODE integrators for dynamical Systems
* Discrete Elements
* Statistical methods (i.e. Monte Carlo)
* Cellular automata
* Lattice-Boltzman simulation

In addition to the C++ codes used in the course, examples an exercises are included to illustrate applications of the techniques studied. Also, a documentation of the use and scope of the routines implemented is included.

## General Objective

The general objective is to wrap up the techniques learnt in the course in handy C++ libraries that allow the user to perform complex physics simulations with little code. The routines should keep **OOP** in mind.

## Specific Objectives

This section is in constant work depending on the development of the course.

### Integrators and Discrete Elements

1. **Create a C++ library with routines that allow user to simulate interesting dynamical systems that are commonly modelled using simple ODEs**: This library should support different types of solver, such as RK4. The focus should be on the use of the solver as an object to solve a particular system characterised by a set of ODEs.

2. **Create a C++ library that allows simulation of a collection of arbitrary objects**: This library should be able to model a particular system as an object, composed of discrete elements that also should be modelled as an object. Some way of using the integrator library in the system class should be devised.
