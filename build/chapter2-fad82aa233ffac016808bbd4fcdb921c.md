---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---
```{code-cell}
:tags: [remove-cell]
import Pkg; Pkg.activate("/Users/driscoll/Documents/GitHub/fnc")
using FundamentalsNumericalComputation
FNC.init_format()
```

(demo-fl-julia)=
``````{dropdown} Absolute and relative accuracy
:open: false
Recall the grade-school approximation to the number $\pi$.

```{code-cell}
@show p = 22/7;
```
::::{grid} 1 1 2 2
Not all the digits displayed for `p` are the same as those of $\pi$. 
:::{card}
:columns: 7

The value of `pi` is predefined and equivalent to `π`, which is entered by typing `\pi` followed immediately by the <kbd>Tab</kbd> key.
:::
::::

```{code-cell}
@show float(π);
```

```{index} ! Julia; string interpolation
```

::::{grid} 1 1 2 2

The absolute and relative accuracies of the approximation are as follows.
:::{card}
:columns: 7

A dollar sign `$` in a string substitutes (or *interpolates*) the named variable or expression into the string.

:::
::::

```{code-cell}
acc = abs(p-π)
println("absolute accuracy = $acc")
println("relative accuracy = $(acc/π)")
```

::::{grid} 1 1 2 2
Here we calculate the number of accurate digits in `p`.
:::{card}
:columns: 7

The `log` function is for the natural log. For other common bases, use `log10` or `log2`.

:::
::::

```{code-cell}
println("Number of accurate digits = $(-log10(acc/π))")
```
This last value could be rounded down by using `floor`.

``````



(function-ho-julia)=
`````{dropdown} **Horner's algorithm for evaluating a polynomial**
:open: true
```{code-block} julia
:lineno-start: 1
"""
    horner(c,x)

Evaluate a polynomial whose coefficients are given in ascending
order in `c`, at the point `x`, using Horner's rule.
"""
function horner(c,x)
    n = length(c)
    y = c[n]
    for k in n-1:-1:1
        y = x*y + c[k]
    end
    return y
end
```
`````

<!-- SECTION 1 -->
(demo-interp-vander-julia)=
``````{dropdown} **Vandermonde matrix and polynomial interpolation**
:open: false


We create two vectors for data about the population of China. The first has the years of census data and the other has the population, in millions of people.

```{code-cell}
year = [1982,2000,2010,2015]; 
pop = [1008.18, 1262.64, 1337.82, 1374.62];
```

:::{index} ! Julia; .-, ! Julia; .+
:::

:::{index} Julia; broadcasting
:::

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item}
:columns: 7


It's convenient to measure time in years since 1980. We use `.-` to subtract a scalar from every element of a vector. We will also use a floating-point value in the subtraction, so the result is also in double precision.


:::
:::{grid-item-card}
:columns: 5
 

A dotted operator such as `.-` or `.*` acts elementwise, broadcasting scalar values to match up with elements of an array.

:::
::::

```{code-cell}
t = year .- 1980.0
y = pop;
```

:::{index} ! Julia; comprehension
:::

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item}
:columns: 7


Now we have four data points $(t_1,y_1),\dots,(t_4,y_4)$, so $n=4$ and we seek an interpolating cubic polynomial. We construct the associated Vandermonde matrix:


:::
:::{grid-item-card}
:columns: 5
 


An expression inside square brackets and ending with a `for` statement is called a **comprehension**. It's often an easy and readable way to construct vectors and matrices. 

:::
::::

```{code-cell}
V = [ t[i]^j for i=1:4, j=0:3 ]
```

:::{index} ! Julia; \\
:::

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item}
:columns: 7


To solve for the vector of polynomial coefficients, we use a backslash to solve the linear system:

:::
:::{grid-item-card}
:columns: 5


A **backslash** `\` is used to solve a linear system of equations.

:::
::::

```{code-cell}
c = V \ y
```

The algorithms used by the backslash operator are the main topic of this chapter. As a check on the solution, we can compute the *residual*.

```{code-cell} julia
y - V*c
```

Using floating-point arithmetic, it is not realistic to expect exact equality of quantities; a relative difference comparable to $\macheps$ is all we can look for.

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item}
:columns: 7


By our definitions, the elements of `c` are coefficients in ascending-degree order for the interpolating polynomial. We can use the polynomial to estimate the population of China in 2005:


:::
:::{grid-item-card}
:columns: 5


The `Polynomials` package has functions to make working with polynomials easy and efficient.

:::
::::

```{code-cell}
p = Polynomial(c)    # construct a polynomial
p(2005-1980)         # include the 1980 time shift
```

The official population value for 2005 was 1303.72, so our result is rather good. 

:::{index} ! Julia; scatter
:::

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item}
:columns: 5


We can visualize the interpolation process. First, we plot the data as points.


:::
:::{grid-item-card}
:columns: 7
 

The `scatter` function creates a scatter plot of points; you can specify a line connecting the points as well.

:::
::::

```{code-cell}
scatter(t,y, label="actual", legend=:topleft,
    xlabel="years since 1980", ylabel="population (millions)", 
    title="Population of China")
```

:::{index} Julia; range
:::

```{index} ! Julia; broadcasting
```

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item}
:columns: 5


We want to superimpose a plot of the polynomial. We do that by evaluating it at a vector of points in the interval. The dot after the name of the polynomial is a universal way to apply a function to every element of an array, a technique known as **broadcasting**.


:::
:::{grid-item-card}
:columns: 7


The `range` function constructs evenly spaced values given the endpoints and either the number of values, or the step size between them.

Adding a dot to the end of a function name causes it to be broadcast, i.e., applied to every element of a vector or matrix.

:::
::::

```{code-cell}
# Choose 500 times in the interval [0,35].
tt = range(0,35,length=500)   
# Evaluate the polynomial at all the vector components.
yy = p.(tt)
foreach(println,yy[1:4])
```

:::{index} ! Julia; \!
:::

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item}
:columns: 5


Now we use `plot!` to add to the current plot, rather than replacing it.


:::
:::{grid-item-card}
:columns: 7


The `plot` function plots lines connecting the given $x$ and $y$ values; you can also specify markers at the points.

By convention, functions whose names end with the bang `!` change the value or state of something, in addition to possibly returning output.

:::
::::

```{code-cell}
plot!(tt,yy,label="interpolant")
```
``````



<!-- SECTION 2 -->
<!-- SECTION 3 -->
<!-- SECTION 4 -->
<!-- SECTION 5 -->
<!-- SECTION 6 -->
<!-- SECTION 7 -->
<!-- SECTION 8 -->
<!-- SECTION 9 -->
<!-- SECTION 10 -->
