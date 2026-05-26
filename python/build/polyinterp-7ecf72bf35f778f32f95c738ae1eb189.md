---
numbering:
  enumerator: 2.1.%s
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---
```{code-cell}
:tags: [remove-cell]
from numpy import *
from scipy import linalg
from scipy.linalg import norm
from matplotlib.pyplot import *
from prettytable import PrettyTable
from timeit import default_timer as timer
import sys
sys.path.append('fncbook/')
import fncbook as FNC

# This (optional) block is for improving the display of plots.
# from IPython.display import set_matplotlib_formats
# set_matplotlib_formats("svg","pdf")
# %config InlineBackend.figure_format = 'svg'
rcParams["figure.figsize"] = [7, 4]
rcParams["lines.linewidth"] = 2
rcParams["lines.markersize"] = 4
rcParams['animation.html'] = "jshtml"  # or try "html5"
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
The function `numpy.linalg.solve` is used to solve linear systems in Python. It is mathematically equivalent to multiplication by the inverse of the matrix in its first argument.
``````

::::{prf:example} Linear system for polynomial interpolation
:label: demo-interp-vander

We create two vectors for data about the population of China. The first has the years of census data and the other has the population, in millions of people.

```{code-cell} 
year = arange(1980, 2020, 10)   # from 1980 to 2020 by 10
pop = array([984.736, 1148.364, 1263.638, 1330.141])
```

It's convenient to measure time in years since 1980. 

```{code-cell} 
t = year - 1980
y = pop
```

Now we have four data points $(t_1,y_1),\dots,(t_4,y_4)$, so $n=4$ and we seek an interpolating cubic polynomial. We construct the associated Vandermonde matrix: 

```{code-cell} 
V = vander(t)
print(V)
```

:::{index} Python; solve
:::

To solve a linear system $\mathbf{V} \mathbf{c} = \mathbf{y}$ for the vector of polynomial coefficients, we use `solve`, located in `scipy.linalg`:

```{code-cell}
from scipy import linalg
c = linalg.solve(V, y)
print(c)
```
The algorithms used by `solve` are the main topic of this chapter. As a check on the solution, we can compute the *residual* $\mathbf{y} - \mathbf{V} \mathbf{c}$, which should be small (near machine precision).

```{tip}
:class: dropdown
Matrix multiplication in NumPy is done with `@` or `matmul`.
```

```{code-cell} 
print(y - V @ c)
```

By our definitions, the coefficients in `c` are given in descending order of power in $t$. We can use the resulting polynomial to estimate the population of China in 2005:

```{code-cell} 
p = poly1d(c)          # construct a polynomial
print(p(2005 - 1980))     # apply the 1980 time shift
```

The official figure was 1303.72, so our result is rather good.

We can visualize the interpolation process. First, we plot the data as points. Then we add a plot of the interpolant, taking care to shift the $t$ variable back to actual years.

```{code-cell} 
scatter(year, y, color="k", label="data");
tt = linspace(0, 30, 300)   # 300 times from 1980 to 2010
plot(1980 + tt, p(tt), label="interpolant");
xlabel("year");
ylabel("population (millions)");
title("Population of China");
legend();
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
