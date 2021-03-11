# Global function approximation

```{index} Solo, Han
```
```{index} The Empire Strikes Back
```
::::{epigraph}
Not entirely stable? I'm glad you're here to tell us these things.
-- Han Solo, *The Empire Strikes Back*
::::

In {ref}`../localfuncapprox/overview` we considered a few ways to map data values to functions via interpolation. The methods we deemed successful were piecewise low-degree polynomials. In this chapter we deal with approximations that are globally defined over the entire interval, not piecewise. 

The conditioning of a global polynomial is unacceptable for high degree interpolants of equally spaced data. We'll remedy that issue by changing how the interpolation nodes are distributed. With that change, polynomial interpolation becomes extremely accurate and fast. Then we will look beyond interpolation and beyond polynomials a bit, and consider the application of these global methods to numerical integration.

**Important terms in this chapter**

::::{glossary}
barycentric formula
  Computationally useful expression for the interpolating polynomial as a ratio of rational terms.

inner product
  Extension of the vector scalar product to a pair of functions.

Lagrange formula
  Theoretically useful expression for the interpolating polynomial.

Lagrange polynomial
  Cardinal function for polynomial interpolation on a given node set.

orthogonal polynomials
  Family of polynomials whose distinct members have an integral inner product equal to zero, as with Legendre and Chebyshev polynomials.

quasimatrix
  Collection of functions (such as orthogonal polynomials) that have algebraic parallels to columns of a matrix.

Runge phenomenon
  Manifestation of the instability of polynomial interpolation at equally spaced nodes as degree increases.

spectral convergence
  Exponentially rapid decrease in error as the number of interpolation nodes increases, as observed in Chebyshev polynomial and trigonometric interpolation.

trigonometric interpolation
  Interpolation of a periodic function by a linear combination of real or complex trigonometric functions.
::::
