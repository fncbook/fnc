---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---

We create two vectors for data about the population of China. The first has the years of census data and the other has the population, in millions of people.

```{code-cell}
year = [1982, 2000, 2010, 2015]; 
pop = [1008.18, 1262.64, 1337.82, 1374.62];
```

:::{index} ! Julia; .-, ! Julia; .+
:::

:::{index} Julia; broadcasting
:::

It's convenient to measure time in years since 1980. We use `.-` to subtract a scalar from every element of a vector. We will also use a floating-point value in the subtraction, so the result is also in double precision.
```{tip}
A dotted operator such as `.-` or `.*` acts elementwise, broadcasting scalar values to match up with elements of an array.
```

```{code-cell}
t = year .- 1980.0
y = pop;
```

:::{index} ! Julia; comprehension
:::

Now we have four data points $(t_1,y_1),\dots,(t_4,y_4)$, so $n=4$ and we seek an interpolating cubic polynomial. We construct the associated Vandermonde matrix:
```{tip}
An expression inside square brackets and ending with a `for` statement is called a **comprehension**. It's often an easy and readable way to construct vectors and matrices. 
```

```{code-cell}
V = [ t[i]^j for i in 1:4, j in 0:3 ]
```

:::{index} ! Julia; \\
:::

To solve for the vector of polynomial coefficients, we use a backslash to solve the linear system:
```{tip}
A **backslash** `\` is used to solve a linear system of equations.
```

```{code-cell}
c = V \ y
```

The algorithms used by the backslash operator are the main topic of this chapter. As a check on the solution, we can compute the *residual*.

```{code-cell} julia
y - V * c
```

Using floating-point arithmetic, it is not realistic to expect exact equality of quantities; a relative difference comparable to $\macheps$ is all we can look for.

By our definitions, the elements of `c` are coefficients in ascending-degree order for the interpolating polynomial. We can use the polynomial to estimate the population of China in 2005:
```{tip}
The `Polynomials` package has functions to make working with polynomials easy and efficient.
```

```{code-cell}
using Polynomials
p = Polynomial(c)    # construct a polynomial
p(2005-1980)         # include the 1980 time shift
```

The official population value for 2005 was 1303.72, so our result is rather good. 

:::{index} ! Julia; scatter
:::

We can visualize the interpolation process. First, we plot the data as points.
```{tip}
The `scatter` function creates a scatter plot of points; you can specify a line connecting the points as well.
```
We want to superimpose a plot of the polynomial. We do that by evaluating it at a vector of points in the interval. The dot after the name of the polynomial is a universal way to apply a function to every element of an array, a technique known as **broadcasting**.
:::{card}
The `range` function constructs evenly spaced values given the endpoints and either the number of values, or the step size between them.

Adding a dot to the end of a function name causes it to be broadcast, i.e., applied to every element of a vector or matrix.
:::
::::

```{code-cell}
# Choose 500 times in the interval [0,35].
tt = range(0, 35, 500)

# Evaluate the polynomial at all the vector components.
yy = p.(tt)

foreach(println, yy[1:4])
```

:::{index} ! Julia; \!
:::

Now we use `plot!` to add to the current plot, rather than replacing it.
```{tip}
The `plot` function plots lines connecting the given $x$ and $y$ values; you can also specify markers at the points.
By convention, functions whose names end with the bang `!` change the value or state of something, in addition to possibly returning output.
```

```{code-cell}
plot!(tt, yy, label="interpolant")
```
