# Roots of nonlinear equations

```{epigraph}
He says "I found her," and keeps repeating, "She's here."

-- C3PO, *Star Wars: A New Hope*
```

In this chapter we extend from linear algebra to deal with *nonlinear* algebraic problems. This kind of problem arises when there is a parameter or variable that can be changed in order to satisfy a constraint or achieve some goal. We start with scalar functions of a single variable, then generalize to $n$ variables and $n$ nonlinear equations. Finally, we generalize the problem of linear least squares to situations with more nonlinear constraints to satisfy than there are variables. In every case the strategy used is one of the cornerstones of numerical computing: *replace a problem you can't solve with an approximate one that you can.* In the context of nonlinear algebraic problems, the particular tactic is to set up and solve a sequence of linear problems of the types covered in the two previous chapters.

## Important terms

```{glossary}
fixed point iteration
  Repeated application of a function in hopes of converging to a fixed point.

fixed point problem
  Finding a value of a given function where the input and output values are the same; equivalent to rootfinding.

Gauss--Newton method
  Generalization of Newton's method for nonlinear least squares.

Jacobian matrix
  Matrix of first partial derivatives that defines the linearization of a vector-valued function.

linear convergence
  Sequence for which the difference between values and the limit asymptotically decreases by a constant factor at each term, making a straight line on a log--linear graph.

Newton's method
  Rootfinding iteration that uses the tangent line (linearization) of the given function in order to define the next root approximation.

nonlinear least squares problem
  Minimization of the 2-norm of the residual of a function that depends nonlinearly on the free parameters.

Quasi-Newton methods
  Rootfinding methods that overcome the issues of Jacobian computation and nonglobal convergence in Newton's method.

quadratic convergence
  Sequence for which the difference between values and the limit asymptotically decreases by a constant times the difference squared of the preceding difference.

rootfinding problem
  Finding a value of a given function that makes the function zero.

secant method
  Variant of Newton's method that uses a secant line rather than a tangent line to define a root estimate.

simple root
  Root at which the derivative of the function is nonzero.

superlinear convergence
  Sequence for which the convergence is asymptotically faster than any linear rate.
```

## Important Julia terms

```{glossary}
`nlsolve`
  Function in the `NLsolve` package for finding roots of scalar and vector-valued functions.
```
