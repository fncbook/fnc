# The interpolation problem

Formally, we now want to solve the following.

```{index} interpolation
```

```{proof:definition} Interpolation problem
Given $n+1$ distinct points $(t_0,y_0)$, $(t_1,y_1),\ldots,(t_n,y_n)$, with $t_0<t_1<\ldots <t_n$, find a function $p(x)$, called the **interpolant**, such that $p(t_k)=y_k$ for $k=0,\dots,n$.
```

The values $t_0,\ldots,t_n$ are called the {term}`nodes` of the interpolant. In this chapter, we use $t_k$ for the nodes and $x$ to denote the continuous independent variable. Take note that the nodes are numbered from zero to $n$. This is convenient for many of our mathematical statements, but less so in a language such as Julia in which all vectors are starting with one. Our approach is that *indices in a computer code have the same meaning as those identically named in the mathematical formulas*, and therefore must be increased by one whenever used in an indexing context.

## Polynomials

```{index} interpolation; by polynomials
```

````{sidebar} Demo
:class: demo
{doc}`demos/interp-global`
````

Polynomials are the obvious first candidate to serve as interpolating functions. They are easy to work with, and in \secref{linsysinterp} we saw that a linear system of equations can be used to determine the coefficients of a polynomial that passes through every member of a set of given points in the plane. However, it's not hard to find examples for which polynomial interpolation leads to unusable results.

```{margin}
Interpolation by a polynomial at equally spaced nodes is ill-conditioned as the degree of the polynomial grows.
```

In Chapter 9 we explore the large oscillations in the last figure of {doc}`demos/interp-global`; it turns out that one must abandon either equally spaced nodes or $n\to\infty$ for polynomials. In the rest of this chapter we will keep $n$ fairly small and let the nodes be unrestricted.

## Piecewise polynomials

```{index} interpolation; by piecewise polnyomials
```

In order to keep polynomial degrees small while interpolating large data sets, we will choose interpolants from the **piecewise polynomials**. Specifically, the interpolant $p$ must be a polynomial on each subinterval $[t_{k-1},t_k]$ for $k=1,\ldots,n$.

````{proof:example}
  Some examples of piecewise polynomials for the nodes  $t_0=-2$, $t_1=0$, $t_2=1$, and $t_3=4$ are $p_1(x)=x+1$, $p_2(x)=\operatorname{sign}(x)$, $p_3(x)=|x-1|^{3}$, and $p_4(x)=(\max\{0,x\})^{4}$. Note that $p_{1}$, $p_{2}$, and $p_4$ would also be piecewise polynomial on the node set $\{t_0,t_1,t_3\}$, but $p_3$ would not.
````

````{sidebar} Demo
:class: demo
{doc}`demos/interp-pwise`
````

Usually we designate in advance a maximum degree $d$ for each polynomial piece of $p(x)$. An important property of the piecewise polynomials of degree $d$ is that they form a vector space: that is, any linear combination of piecewise polynomials of degree $d$ is another piecewise polynomial of degree $d$. If $p$ and $q$ share the same node set, then the combination is piecewise polynomial on that node set.

We will consider piecewise linear interpolation in more detail in {doc}`pwlin`, and we look at this type of piecewise cubic interpolation in {doc}`splines`.

## Conditioning of interpolation

```{index} condition number; of interpolation
```

In the interpolation problem we are given the values $(t_k,y_k)$ for $k=0,\ldots,n$. Let us consider the nodes $t_k$ of the problem to be fixed, and let $a=t_0$, $b=t_n$. Then the data for the interpolation problem consists of a vector $\mathbf{y}$, and the result of the problem is a function on $[a,b]$.

Let $\mathcal{I}$ be a prescription for producing the interpolant from a data vector.  That is, $\mathcal{I}(\mathbf{y})=p$, where $p(t_k)=y_k$ for all $k$. The interpolation methods we will consider are all *linear*, in the sense that

```{math}
:label: interp-linearity
\mathcal{I}(\alpha\mathbf{y} + \beta\mathbf{z}) = \alpha \mathcal{I}(\mathbf{y}) + \beta \mathcal{I}(\mathbf{z})
```

for all vectors $\mathbf{y},\mathbf{z}$ and scalars $\alpha,\beta$.

Linearity greatly simplifies the analysis of the conditioning of interpolation. If the data are changed from $\mathbf{y}$ to $\mathbf{y}+ \Delta \mathbf{y}$, then

```{math}
  :label: lininterp-perturb
  \Delta \mathcal{I} = \mathcal{I}(\mathbf{y} + \Delta \mathbf{y}) - \mathcal{I}(\mathbf{y}) =  \mathcal{I}(\Delta \mathbf{y})
 = \sum_{k=0}^{n} (\Delta y)_k\mathcal{I}(\mathbf{e}_k),
```

where as always $\mathbf{e}_k$ is a column of an identity matrix. We use $\|\Delta \mathbf{y}\|_\infty$ to measure the size of the perturbation to the data, and for $\Delta \mathcal{I}$, which is a function, we use the functional infinity-norm or max-norm defined by

```{math}
:label: norm-function-inf
\| f\|_{\infty} = \max_{x \in [a,b]} |f(x)|.
```

The absolute condition number $\kappa(\mathbf{y})$ relates $\|\Delta \mathcal{I} \|_\infty$ to $\|\Delta \mathbf{y}\|_\infty$.

(theorem-interp-conditioning)=

````{proof:theorem} Interpolation conditioning
Suppose that $\mathcal{I}$ is a linear interpolation method. Then the absolute condition number of $\mathcal{I}$ satisfies
  
```{math}
:label: interp-conditioning
\max_{0\le k \le n} \bigl\| \mathcal{I}(\mathbf{e}_k) \bigr\|_\infty \le \kappa(\mathbf{y}) \le  \sum_{k=0}^n  \bigl\| \mathcal{I}(\mathbf{e}_k) \bigr\|_\infty,
```

if vectors and functions are measured in the infinity norm.
````

````{proof:proof}
Because of {eq}`lininterp-perturb`, we have

```{math}
\frac{\bigl\| \Delta \mathcal{I} \bigr\|_{\infty}}{\| \Delta\mathbf{y}\|_{\infty}} =
\left\|\, \sum_{k=0}^{n} \frac{\Delta y_k}{\|\Delta \mathbf{y}\|_{\infty}} \mathcal{I}(\mathbf{e}_k)\,  \right\|_{\infty}.
```

The absolute condition number maximizes this quantity over all $\Delta \mathbf{y}$. Suppose $j$ is such that $\|\mathcal{I}(\mathbf{e}_j)\|$ is maximal. Then let $\Delta \mathbf{y}=\mathbf{e}_j$ and the first inequality in {eq}`interp-conditioning` follows. The other inequality follows from the triangle inequality and the definition of $\|\Delta \mathbf{y}\|_{\infty}$.
````

```{index} cardinal functions
```

```{margin}
The condition number of a linear interpolation method is essentially that of its cardinal functions.
```

````{sidebar} Demo
:class: demo
{doc}`demos/interp-cond`
````

The [interpolation conditioning theorem](theorem-interp-conditioning) says that assessing the condition number of interpolation for any data can be simplified to measuring the effect of interpolation with each node "switched on" one at a time. The results of such interpolations are known as {term}`cardinal functions`.

## Exercises

1. ⌨ Create data by entering

    ``` julia
    t = -2:4;  y = tanh.(t);
    ```

    **(a)** Use `fit` to construct and plot the polynomial interpolant of the data.

    **(b)** Use `CubicSplineInterpolation` to construct and plot a piecewise cubic interpolant of the data.

    ````{only} solutions
    ````

2. ⌨ The following table gives the life expectancy in the U.S. by year of birth.

      | 1980 | 1985 | 1990 | 1995 | 2000 | 2005 | 2010 |
      |:---:|:----:|:-----:|:----:|:----:|:----:|:----:|
      | 73.7 | 74.7 | 75.4 | 75.8 | 77.0 | 77.8 | 78.7 |

    **(a)** Defining "year since 1980" as the independent variable, use `fit` to construct and plot the polynomial interpolant of the data.

    **(b)** Use `CubicSplineInterpolation` to construct and plot a piecewise cubic interpolant of the data.

    **(c)** Use both methods to estimate the life expectancy for a person born in 2007. Which value is more believable?
  
    ````{only} solutions
    t = (1980:5:2010)'
    y = [ 73.7 74.7 75.4 75.8 77.0 77.8 78.7]';

    clf
    t = t-1980;
    plot(t,y,'o')
    c = polyfit(t,y,length(t)-1);
    p = @(x) polyval(c,x);
    hold on
    fplot(p,[0 30]);
    x = linspace(0,30,300);
    plot(x,interp1(t,y,x,'cubic'))

    %% (b)
    polyval(p,27)
    interp1(t,y,27)
    ````

3. ⌨ The following two point sets define the top and bottom of a flying saucer shape.
    Top:

    ``` julia
    x = [ 0 , 0.51 , 0.96 , 1.06 , 1.29 , 1.55 , 1.73 , 2.13 , 2.61 ]
    y = [ 0 , 0.16 , 0.16 , 0.43 , 0.62 , 0.48 , 0.19 , 0.18 , 0    ]
    ```

    ``` julia
    x = [ 0 ,  0.58 ,  1.04 ,  1.25 ,  1.56 ,  1.76 ,  2.19 , 2.61 ]
    y = [ 0 , -0.16 , -0.15 , -0.30 , -0.29 , -0.12 , -0.12 , 0    ]
    ```

    Use `CubicSplineInterpolation` to make a picture of the flying saucer. 

    ````{only} solutions
    x1 = [ 0.00  0.51  0.96  1.06  1.29  1.55  1.73  2.13 2.61 ];
    y1 = [ 0  0.16  0.16  0.43  0.62  0.48  0.19  0.18 0 ]; 
    x = linspace(0,2.61,200);
    clf
    plot(x,interp1(x1,y1,x,'cubic'))
    hold on

    x2 = [ 0.00  0.58  1.04  1.25  1.56  1.76  2.19  2.61];
    y2 = [ 0.00  -0.16  -0.15 -0.30  -0.29  -0.12  -0.12  0];
    plot(x,interp1(x2,y2,x,'cubic'))
    ````

    (problem-quadratic-interpolant)=
4. ✍ Define

    ```{math}
    q(x) = a\frac{x(x-1)}{2} - b (x-1)(x+1) + c \frac{x(x+1)}{2}.
    ```

    **(a)** Show that $q$ is a polynomial interpolant of the points $(-1,a)$, $(0,b)$, $(1,c)$.

    **(b)** Use a change of variable to find a quadratic polynomial interpolant for the points $(x_0-h,a)$, $(x_0,b)$, $(x_0+h,c)$.
  
    ````{only} solutions
    ````

5. ✍ Use the formula of [this previous exercise](problem-quadratic-interpolant) and the [interpolation conditioning theorem](theorem-interp-conditioning)` to derive bounds on the condition number of quadratic polynomial interpolation at the nodes $-1$, $0$, $1$.

    ````{only} solutions
    ````
