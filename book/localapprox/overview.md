# Piecewise interpolation

```{epigraph}
You must feel the Force around you. Here, between you...me...the tree...the rock...everywhere!

--Yoda, *The Empire Strikes Back*
```

In many scientific problems the solution is a function. Accordingly, our next task is to represent functions numerically. This task is more difficult and complicated than the one we faced in representing real numbers. With numbers it's intuitively clear how one real value can stand for a small interval around it. But designating representatives for sets of functions is less straightforward—in fact, it's one of the core topics in computing. The process of converting functions into numerical representations of finite length is known as {term}`discretization`.

Once we have selected a method of discretization, we can define numerical analogs of our two favorite operations on functions, differentiation and integration. These are linear operations, so the most natural numerical analogs are linear operations too. As we will see in many of the chapters following this one, a lot of numerical computing boils down to converting calculus to algebra, with discretization as the link between them.

**Important terms**

```{glossary}
cardinal functions
  Functions that are one at one interpolation node and zero at all the others, forming a useful basis for interpolation.

cubic spline
  Piecwise cubic function with two globally continuous derivatives, often used for interpolation.

discretization
  Replacement of functions by approximations of finite length.

extrapolation
  Use of multiple discretization values to cancel out terms in an error expansion.

finite difference
  Linear combination of function values that approximates the value of a derivative of the function at a point.

hat functions
  Piecewise linear cardinal functions.

interpolation
  Construction of a function that passes through a given set of data points.

Newton--Cotes formula
  Linear combination of function values that approximates the definite integral of the function.

Nodes
  Values of the independent variable where an interpolant's values are prescribed.

order of accuracy
  Leading power of the truncation error as a function of discretization size.

trapezoid formula
  Newton-Cotes integration formula of order 2.

truncation error
  Difference between an exact value and an approximation, typically one that truncates an infinite series.


```

**Important Julia terms**

```{glossary}
`CubicSplineInterpolation`
  Piecewise cubic spline interpolant for given data points (in `Interpolations`).

`fit`
  Polynomial interpolant for given data points (in `Polynomials`).

`LinearInterpolation`
  Piecewise linear interpolant for given data points (in `Interpolations`).
```
