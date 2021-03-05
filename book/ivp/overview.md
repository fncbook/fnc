# Initial-value problems for ODEs

```{epigraph}
Without precise calculations we could fly right through a star or bounce too close to a supernova and that'd end your trip real quick, wouldn't it?

--Han Solo, *Star Wars: A New Hope*
```

Quantities that change continuously in time or space are often modeled by differential equations. When everything depends on just one independent variable, we call the model an **ordinary differential equation** (ODE).  Differential equations need supplemental conditions to define both the modeling situation and the theoretical solutions uniquely. The **initial-value problem** (IVP), in which all of the conditions are given at a single value of the independent variable, is the simplest situation. Often the independent variable in this case represents time.

Methods for IVPs usually start from the known initial value and iterate or "march" forward from there. There is a large number of them, owing in part to differences in accuracy, stability, and convenience. The most broadly important methods fall into one of two camps: **Runge--Kutta** and **linear multistep** formulas. Each type introduces its own complications, and we will consider them separately.

**Important terms**

```{glossary}
Euler's method
  Prototype of all IVP solution methods, obtained by assuming constant derivatives for the solution over short time intervals.

initial-value problem
  An ordinary differential equation (possibly vector-valued) together with an initial condition.

implicit
  Describes a formula that defines a new solution value as the solution of an implicit equation.

local truncation error
  Discretization error made in one time step of a solution method.

generating polynomials
  A pair of polynomials whose coefficients match those of a multistep method for IVPs.

global error
  Error made by a method over the entire time interval of the solution.

multistep
  Formula using information over more than a single time step to advance the solution.

order of accuracy
  Leading exponent of discretization size in the local truncation error, also equal to the leading exponent in the global truncation error.

Runge--Kutta
  Describing a multistage method that evaluates the derivative of the solution more than once to advance a single step.

step size
  Increment in time between successive solution values in a numerical IVP solver.

stiff
  Describes an IVP in which stability is a greater restriction than accuracy for many solution methods.

```

**Important Julia commands and keywords**

```{glossary}
`DifferentialEquations`
  Package for solving ordinary differential equations in Julia.

`ODEProblem`
  Defines an initial-value problem from the time derivative function, initial value, and time domain.

`solve`
  Solves an initial-value problem (in `DifferentialEquations` package)
```
