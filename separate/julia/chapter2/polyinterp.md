---
numbering:
  enumerator: 2.1.%s
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.12
---
```{code-cell}
:tags: [remove-cell]
import Pkg; Pkg.activate("/Users/driscoll/Documents/GitHub/fnc")

using FNCFunctions

using Plots
default(
    titlefont=(11,"Helvetica"),
    guidefont=(11,"Helvetica"),
    linewidth = 2,
    markersize = 3,
    msa = 0,
    size=(500,320),
    label="",
    html_output_format = "svg"
)

using PrettyTables, LaTeXStrings, Printf
using LinearAlgebra

@ptconf backend = Val(:html) tf = tf_html_simple
```

(section-linsys-polyinterp)=

# Polynomial interpolation

```{index} interpolation
```

The United States carries out a census of its population every 10 years. Suppose we want to know the population at times in between the census years, or to estimate future populations. One technique is to find a polynomial that passes through all the data points.[^census]

::::{prf:definition} Polynomial interpolation
:label: definition-polyinterp
Given $n$ points $(t_1,y_1),\ldots,(t_n,y_n)$, where the $t_i$ are all distinct, the **polynomial interpolation** problem is to find a polynomial $p$ of degree less than $n$ such that $p(t_i)=y_i$ for all $i$.
::::

As posed in {numref}`Definition {number} <definition-polyinterp>`, the polynomial interpolation problem has a unique solution. Once the interpolating polynomial is found, it can be evaluated anywhere to estimate or predict values.

[^census]: We're quite certain that the U.S. Census Bureau uses more sophisticated modeling techniques than the one we present here!

## Interpolation as a linear system

Given data $(t_i,y_i)$ for $i=1,\ldots,n$, we seek a polynomial

```{math}
:label: vanderinterp
p(t) = c_1 + c_{2} t + c_3t^2 +  \cdots + c_{n} t^{n-1},
```

such that $y_i=p(t_i)$ for all $i$. These conditions are used to determine the coefficients $c_1\ldots,c_n$:

```{math}
\begin{split}
 c_1 + c_2 t_1 + \cdots + c_{n-1}t_1^{n-2} + c_nt_1^{n-1} &= y_1, \\
 c_1 + c_2 t_2 + \cdots + c_{n-1}t_2^{n-2} + c_nt_2^{n-1} &= y_2, \\
 c_1 + c_2 t_3 + \cdots + c_{n-1}t_3^{n-2} + c_nt_3^{n-1} &= y_3, \\
 \vdots \qquad & \\
 c_1 + c_2 t_n + \cdots + c_{n-1}t_n^{n-2} + c_nt_n^{n-1} &= y_n.
 \end{split}
```

These equations form a linear system for the coefficients $c_i$:

```{math}
:label: vandersystem
  \begin{bmatrix}
    1 & t_1 & \cdots & t_1^{n-2} & t_1^{n-1} \\
    1 & t_2 & \cdots & t_2^{n-2} & t_2^{n-1} \\
    1 & t_3 & \cdots & t_3^{n-2} & t_3^{n-1} \\
    \vdots & \vdots &  & \vdots & \vdots \\
    1 & t_n & \cdots & t_n^{n-2} & t_n^{n-1} \\
  \end{bmatrix}
  \begin{bmatrix}
    c_1  \\
    c_2  \\
    c_3 \\
    \vdots \\
    c_n
  \end{bmatrix}
  =
  \begin{bmatrix}
    y_1  \\
    y_2  \\
    y_3 \\
    \vdots \\
    y_n
  \end{bmatrix},
```

or more simply, $\mathbf{V} \mathbf{c} = \mathbf{y}$. The matrix $\mathbf{V}$ is of a special type.

```{index} ! Vandermonde matrix
```

::::{prf:definition} Vandermonde matrix
:label: definition-vandermonde
Given distinct values $t_1,\ldots,t_n$, a {term}`Vandermonde matrix` for these values is the $n\times n$ matrix appearing in {eq}`vandersystem`.
::::

Polynomial interpolation can therefore be posed as a linear system of equations with a Vandermonde matrix.

:::{index} ! Julia; \\
:::

:::{index} ! MATLAB; \\
:::

:::{index} ! Python; solve
:::

``````{important}
The backslash operator `\` is used to solve linear systems in Julia. It is mathematically equivalent to multiplication by the inverse of the matrix on its left.
``````

::::{prf:example} Linear system for polynomial interpolation
:label: demo-interp-vander


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
:class: dropdown
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
:class: dropdown
An expression inside square brackets and ending with a `for` statement is called a **comprehension**. It's often an easy and readable way to construct vectors and matrices. 
```

```{code-cell}
V = [ t[i]^j for i in 1:4, j in 0:3 ]
```

:::{index} Julia; \\
:::

To solve for the vector of polynomial coefficients, we use a backslash to solve the linear system:
```{tip}
:class: dropdown
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
:class: dropdown
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
:class: dropdown
The `scatter` function creates a scatter plot of points; you can specify a line connecting the points as well.
```

```{code-cell}
using Plots
scatter(t, y;
    label="actual",  legend=:topleft,
    xlabel="years since 1980",  ylabel="population (millions)", 
    title="Population of China")
```

:::{index} Julia; range
:::

```{index} ! Julia; broadcasting
```

::::{grid} 1 1 2 2
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
:class: dropdown
The `plot` function plots lines connecting the given $x$ and $y$ values; you can also specify markers at the points.
By convention, functions whose names end with the bang `!` change the value or state of something, in addition to possibly returning output.
```

```{code-cell}
plot!(tt, yy, label="interpolant")
```

::::

## Exercises

``````{exercise}
:label: problem-polyinterp-conditions

Suppose you want to interpolate the points $(-1,0)$, $(0,1)$, $(2,0)$, $(3,1)$, and $(4,2)$ by a polynomial of as low a degree as possible.

**(a)** ✍ Write out a linear system of equations for the coefficients of the interpolating polynomial.

**(b)** ⌨ Solve the system in (b) numerically.
``````

``````{exercise}
:label: problem-polyinterp-hermite

**(a)** ✍ Say you want to find a cubic polynomial $p$ such that $p(-1) =-2$, $p'(-1) =1$, $p(1) = 0$, and $p'(1) =-1$. (This is known as a *Hermite interpolant.*) Write out a linear system of equations for the coefficients of $p$. 

**(b)** ⌨ Numerically solve the linear system in part (a) and make a plot of $p$ over $-1 \le x \le 1$. 
``````

``````{exercise}
:label: problem-polyinterp-census

⌨ Here are population figures (in millions) for three countries over a 30-year period (from *United Nations World Population Prospects*, 2019).

|          | 1990  |  2000  | 2010 | 2020 |
|:------:|:----------:|:---------:|:---------:|:---------:|
| United States | 252.120  | 281.711 | 309.011 | 331.003 |
| India | 873.278 | 1,056.576 | 1,234.281 | 1,380.004 |
| Poland | 37.960 | 38.557 | 38.330 | 37.847 |

**(a)** Use cubic polynomial interpolation to estimate the population of the USA in 2005.

**(b)** Use cubic polynomial interpolation to estimate when the population of Poland peaked during this time period.

**(c)** Use cubic polynomial interpolation to make a plot of the Indian population over this period. Your plot should be well labeled and show a smooth curve as well as the original data points.
``````

``````{exercise}
:label: problem-polyinterp-delaware

⌨ Here are the official population figures for the state of Delaware, USA, every ten years from 1790 to 1900: 59096, 64273, 72674, 72749, 76748, 78085, 91532, 112216, 125015, 146608, 168493, 184735. For this problem, use 

```{math}
:numbered: false
t = \frac{\text{year} - 1860}{10}
```

as the independent (time) variable.

**(a)** Using only the data from years 1860 to 1900, plot the interpolating polynomial over the same range of years. Add the original data points to your plot as well.

**(b)** You might assume that adding more data will make the interpolation better. But this is not always the case. Use all the data above to create an interpolating polynomial of degree 11, and then plot that polynomial over the range 1860 to 1900. In what way is this fit clearly inferior to the previous one? (This phenomenon is studied in [Chapter 9](../globalapprox/overview).)
``````
