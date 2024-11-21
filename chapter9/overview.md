# Global function approximation

```{index} Han Solo, The Empire Strikes Back
```
::::{epigraph}
Not entirely stable? I'm glad you're here to tell us these things.

— Han Solo, *The Empire Strikes Back*
::::

In [Chapter 5](../localapprox/overview) we considered a few ways to map data values to functions via interpolation. The methods we deemed successful were piecewise low-degree polynomials. In this chapter we deal with approximations that are globally defined over the entire interval, not piecewise. 

The conditioning of a global polynomial is unacceptable for high degree interpolants of equally spaced data. We'll remedy that issue by changing how the interpolation nodes are distributed. With that change, polynomial interpolation becomes extremely accurate and fast. Then we will look beyond interpolation and beyond polynomials a bit, and consider the application of these global methods to numerical integration.

