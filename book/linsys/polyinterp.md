# Polynomial interpolation

```{index} interpolation
```

The United States carries out a census of its population every 10 years. Suppose we want to know the population at times in-between the census years, or to estimate future populations.[^census] One technique is to use {term}`interpolation`. Interpolation is the process of constructing a mathematical function of a continuous variable that passes through a given set of points.  Once an interpolating function, or *interpolant*, is found, it can be evaluated to estimate or predict values. By definition, the interpolant "predicts" the correct values at all of the known data points.

[^census]: We're quite certain that the U.S. Census Bureau uses more sophisticated modeling techniques than the one we present here!

## Interpolation as a linear system

Polynomials are one of the most popular and useful function types for interpolation. Suppose data $(t_i,y_i)$ are given for $i=1,\ldots,n$. If we assume that $n$ values should be enough to uniquely determine $n$ unknown coefficients, then we can use a polynomial of degree $n-1$ to interpolate the data. So we assume interpolation by $y=f(t)$, where

```{math}
:label: vanderinterp
f(t) = c_1 + c_{2} t + c_3t^2 +  \cdots + c_{n} t^{n-1}.
```

We determine the coefficients $c_1\ldots,c_n$ by setting $y_i=f(t_i)$ for all the values of $i$. Writing out these equations gives

```{math}
\begin{split}
 c_1 + c_2 t_1 + \cdots + c_{n-1}t_1^{n-2} + c_nt_1^{n-1} &= y_1 \\
 c_1 + c_2 t_2 + \cdots + c_{n-1}t_2^{n-2} + c_nt_2^{n-1} &= y_2 \\
 c_1 + c_2 t_3 + \cdots + c_{n-1}t_3^{n-2} + c_nt_3^{n-1} &= y_3 \\
 \vdots \qquad & \\
 c_1 + c_2 t_n + \cdots + c_{n-1}t_n^{n-2} + c_nt_n^{n-1} &= y_n 
 \end{split}
```

These equations are not linear in the $t_i$. However, they do form a linear system for the coefficients $c_i$:

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

```{index} Vandermonde matrix
```

```{margin}
Polynomial interpolation can be posed as a linear system of equations with a Vandermonde matrix.
```

or more simply, $\mathbf{V} \mathbf{c} = \mathbf{y}$. The matrix $\mathbf{V}$ is of a
special type known as a {term}`Vandermonde matrix`. Polynomial interpolation can be posed as a linear system of equations with a Vandermonde matrix.

```{proof:example} Julia demo
:class: demo
{doc}`demos/interp-vander`
```

The construction of the Vandermonde matrix in the example {doc}`demos/interp-vander` was done in the line

```julia
V = [ t[i]^j for i in 1:4, j in 0:3 ]
```

```{index} comprehension
```

This syntax uses a Julia {term}`comprehension`. For the most part comprehensions act like loops, but they are more compact and do not require an advance allocation of space for the result.

## Exercises

1. Suppose you want to interpolate the points $(-1,0)$, $(0,1)$, $(2,0)$, $(3,1)$, and $(4,2)$ by a polynomial of as low a degree as possible.
  
   **(a)** ✍ What is the maximum necessary degree of this polynomial?

   **(b)** ✍ Write out a linear system of equations for the coefficients of the interpolating polynomial.

   **(c)** ⌨ Use Julia to solve the system in (b) numerically.
  
    ````{only} solutions
    **(a)** The degree is one less than the number of points, so here it is 4.

    **(b)--(c)** Here are the matrix and vector of the linear system.
    ``` julia
    t = [-1,0,2,3,4]
    A = [ t^j for t in t, j in 0:4]
    y = [0,1,0,1,2]
    A\y
    ```
    ````

2. ✍ Say you want to find a cubic polynomial $p$ such that $p(0)=-2$, $p'(0)=1$, $p(1)=2$, and $p'(1)=-1$. (This is known as a *Hermite interpolant.*) Write out a linear system of equations for the coefficients of $p$. You do not need to solve the resulting system.

3. ⌨ Here are population figures for three countries over the same 30-year period as in {doc}`demos/interp-vander`.

    | Year | United States | China    | Germany |
    |:------:|:----------:|:---------:|:---------:|
    | 1980 | 227.225      | 984.736   | 78.298 |
    | 1990 | 249.623      | 1,148.364 | 79.380 |
    | 2000 | 282.172      | 1,263.638 | 82.184 |
    | 2010 | 308.282      | 1,330.141 | 81.644 |

    **(a)** Use cubic polynomial interpolation to estimate the population of China in 1992.

    **(b)** Use cubic polynomial interpolation to estimate the population of the USA in 1984.
  
    **(c)** Use cubic polynomial interpolation to make a plot of the German population from 1985 to 2000. Your plot should show a smooth curve and be well annotated.
  
4. ⌨ The shifting of years in {doc}`demos/interp-vander`, so that $t=0$ means 1980, was more significant than it might seem.

    **(a)** Using {doc}`demos/interp-vander` as a model, find the ratio $a_4/a_1$ using the original code.

    **(b)** Repeat the computation without shifting by 1980; that is, the entries of the time vector are the actual years. When you solve for the coefficients of the cubic polynomial, you will probably get a cryptic warning. Find $a_4/a_1$ again.

    **(c)** Continuing part (b), again estimate the population of China in 2005. How different, in relative terms, is the new answer?

    The phenomena observed in these experiments will be explained later in this chapter.
