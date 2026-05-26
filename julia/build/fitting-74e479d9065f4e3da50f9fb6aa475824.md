---
numbering:
  enumerator: 3.1.%s
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

(section-leastsq-fitting)=

# Fitting functions to data

```{index} interpolation; by polynomials
```

In {numref}`section-linsys-polyinterp` we saw how a polynomial can be used to interpolate data—that is, derive a continuous function that evaluates to give a set of prescribed values. But interpolation may not be appropriate in many applications.

::::{prf:example} Interpolating temperature data
:label: demo-fitting-tempinterp

Here are 5-year averages of the worldwide temperature anomaly as compared to the 1951–1980 average (source: NASA).

```{code-cell}
year = 1955:5:2000
temp = [ -0.0480, -0.0180, -0.0360, -0.0120, -0.0040,
       0.1180, 0.2100, 0.3320, 0.3340, 0.4560 ]
    
scatter(year, temp, label="data",
    xlabel="year", ylabel="anomaly (degrees C)", 
    legend=:bottomright)
```

A polynomial interpolant can be used to fit the data. Here we build one using a Vandermonde matrix. First, though, we express time as decades since 1950, as it improves the condition number of the matrix.

```{code-cell}
t = @. (year - 1950) / 10
n = length(t)
V = [ t[i]^j for i in 1:n, j in 0:n-1 ]
c = V \ temp
```

```{index} Julia; plotting functions
```

The coefficients in vector `c` are used to create a polynomial. Then we create a function that evaluates the polynomial after changing the time variable as we did for the Vandermonde matrix.
```{tip}
:class: dropdown
If you `plot` a function, then the points are chosen automatically to make a smooth curve.
```

```{code-cell}
using Polynomials, Plots
p = Polynomial(c)
f = yr -> p((yr - 1950) / 10)
plot!(f, 1955, 2000, label="interpolant")
```

As you can see, the interpolant does represent the data, in a sense. However it's a crazy-looking curve for the application. Trying too hard to reproduce all the data exactly is known as _overfitting_.


::::

```{index} data fitting
```

In many cases we can get better results by relaxing the interpolation requirement. In the polynomial case this allows us to lower the degree of the polynomial, which limits the number of local max and min points. Let $(t_i,y_i)$ for $i=1,\ldots,m$ be the given points. We will represent the data by the polynomial

```{math}
:label: vanderls
y \approx f(t) = c_1 + c_2t + \cdots + c_{n-1} t^{n-2} + c_n t^{n-1},
```

with $n<m$. Just as in {eq}`vandersystem`, we can express a vector of $f$-values by a matrix-vector multiplication. In other words, we seek an approximation

```{math}
:label: vandersystemrect
\begin{bmatrix} y_1 \\ y_2 \\ y_3 \\ \vdots \\ y_m \end{bmatrix} \approx
\begin{bmatrix}
f(t_1)                               \\
f(t_2)                               \\
f(t_3)                               \\
\vdots                               \\
f(t_m)
\end{bmatrix} =
\begin{bmatrix}
1      & t_1    & \cdots & t_1^{n-1} \\
1      & t_2    & \cdots & t_2^{n-1} \\
1      & t_3    & \cdots & t_3^{n-1} \\
\vdots & \vdots &        & \vdots    \\
1      & t_m    & \cdots & t_m^{n-1} \\
\end{bmatrix}
\begin{bmatrix}
c_1                                  \\
c_2                                  \\
\vdots                               \\
c_n
\end{bmatrix}
= \mathbf{V} \mathbf{c}.
```

```{index} Vandermonde matrix
```

Note that $\mathbf{V}$ has the same structure as the Vandermonde matrix in {eq}`vandersystem` but is $m\times n$, thus taller than it is wide. It's impossible in general to satisfy $m$ conditions with $n<m$ variables, and we say the system is **overdetermined**. Rather than solving the system exactly, we have to find the best approximation. Below we specify precisely what is meant by this, but first we note that Julia uses the same backslash notation to solve the problem in both the square and overdetermined cases.

::::{prf:example} Fitting temperature data
:label: demo-fitting-tempfit

Here are the 5-year temperature averages again.

```{code-cell}
year = 1955:5:2000
t = @. (year - 1950) / 10
temp = [ -0.0480, -0.0180, -0.0360, -0.0120, -0.0040,
          0.1180, 0.2100, 0.3320, 0.3340, 0.4560 ]
```

```{index} Julia; \\
```

The standard best-fit line results from using a linear polynomial that meets the least-squares criterion.
```{tip}
:class: dropdown
Backslash solves overdetermined linear systems in a least-squares sense.
```

```{code-cell}
V = [ t.^0 t ]    # Vandermonde-ish matrix
@show size(V)
c = V \ temp
p = Polynomial(c)
```

```{code-cell}
f = yr -> p((yr - 1955) / 10)
scatter(year, temp, label="data",
    xlabel="year", ylabel="anomaly (degrees C)", leg=:bottomright)
plot!(f, 1955, 2000, label="linear fit")
```

If we use a global cubic polynomial, the points are fit more closely.

```{code-cell}
V = [ t[i]^j for i in 1:length(t), j in 0:3 ]   
@show size(V);
```

Now we solve the new least-squares problem to redefine the fitting polynomial.
```{tip}
:class: dropdown
The definition of `f` above is in terms of `p`. When `p` is changed, then `f` calls the new version.
```

```{code-cell}
p = Polynomial( V \ temp )
plot!(f, 1955, 2000, label="cubic fit")
```

If we were to continue increasing the degree of the polynomial, the residual at the data points would get smaller, but overfitting would increase.

::::

## The least-squares formulation

In the most general terms, our fitting functions take the form

```{math}
:label: fitbasis
f(t) = c_1 f_1(t) + \cdots + c_n f_n(t)
```

where $f_1,\ldots,f_n$ are all known functions with no undetermined parameters. This leaves only $c_1,\ldots,c_n$ to be determined. The essential feature of a linear least-squares problem is that the fit depends only *linearly* on the unknown parameters. For instance, a function of the form $f(t)=c_1 + c_2 e^{c_3 t}$ is not of this type.

At each observation $(t_i,y_i)$, we define a residual, $y_i - f(t_i)$. A sensible formulation of the fitting criterion is to minimize

```{math}
  R(c_1,\ldots,c_n) = \sum_{i=1}^m\, [ y_i - f(t_i) ]^2,
```

over all possible choices of parameters $c_1,\ldots,c_n$. We can apply linear algebra to write the problem in the form $R=\mathbf{r}^T \mathbf{r}$, where

```{math}
\mathbf{r} =
\begin{bmatrix}
y_1 \\ y_2 \\ \vdots \\y_{m-1} \\ y_m
\end{bmatrix} -
\begin{bmatrix}
f_1(t_1) & f_2(t_1) & \cdots & f_n(t_1) \\[1mm]
f_1(t_2) & f_2(t_2) & \cdots & f_n(t_2) \\[1mm]
& \vdots \\
f_1(t_{m-1}) & f_2(t_{m-1}) & \cdots & f_n(t_{m-1}) \\[1mm]
f_1(t_m) & f_2(t_m) & \cdots & f_n(t_m) \\[1mm]
\end{bmatrix}
\begin{bmatrix}
c_1 \\ c_2 \\ \vdots \\ c_n
\end{bmatrix}.
```

Recalling that $\mathbf{r}^T\mathbf{r}=\| \mathbf{r} \|_2^2$, and renaming the variables to standardize the statement, we arrive at the general {term}`linear least-squares problem`.

```{index} ! argmin
```

```{index} ! linear least-squares problem
```

````{prf:definition} Linear least-squares problem
:label: definition-leastsq
Given $\mathbf{A}\in\mathbb{R}^{m \times n}$ and $\mathbf{b}\in\mathbb{R}^m$, with $m>n$, find

```{math}
:label: linls
\argmin_{{\mathbf{x}\in \mathbb{R}^n}} \, \twonorm{\mathbf{b}-\mathbf{A} \mathbf{x}}^2.
```
````

:::{note}
The notation **argmin** in @linls means to find an $\mathbf{x}$ that produces the minimum value.
:::

```{note}
Another way to describe @linls is as an {term}`overdetermined` linear system. Because the system $\mathbf{A}\mathbf{x}=\mathbf{b}$, has more rows than columns, it is imposing too many constraints on $\mathbf{x}$ to permit satisfying all of them. 
```

```{index} norm; vector
```

While we could choose to minimize in any vector norm, the 2-norm is the most common and convenient choice. *For the rest of this chapter we exclusively use the 2-norm.* In the edge case $m=n$ for a nonsingular $\mathbf{A}$, the definitions of the linear least-squares and linear systems problems coincide: the solution of $\mathbf{A}\mathbf{x}=\mathbf{b}$ implies $\mathbf{r}=\boldsymbol{0}$, which is a global minimum of $\| \mathbf{r} \|_2^2 \ge 0$.

## Change of variables

```{index} data fitting; by straight line
```

The most familiar and common case of a polynomial least-squares fit is the straight line, $f(t) = c_1 + c_2 t$. Certain other fit functions can be transformed into this situation. For example, suppose we want to fit data using $g(t)= a_1 e^{a_2 t}$. Then

```{math}
:label: exptransfit
\log y \approx \log g(t) = (\log a_1) + a_2 t = c_1 + c_2 t.
```

```{index} data fitting; power law
```

While the fit of the $y_i$ to $g(t)$ is nonlinearly dependent on fitting parameters, the fit of $\log y$ to a straight line is a linear problem. Similarly, the power-law relationship $y\approx f(t)=a_1 t^{a_2}$ is equivalent to

```{math}
:label: powtransfit
\log y \approx (\log a_1) + a_2 (\log t).
```

Thus, the variable $z=\log y$ can be fit linearly in terms of the variable $s=\log t$. In practice these two cases—exponential fit and power law—are easily detected by using log-linear or log-log plots, respectively.

::::{prf:example} Fitting a power law
:label: demo-fitting-pirate
Finding numerical approximations to $\pi$ has fascinated people for millennia. One famous formula is

$$
\displaystyle \frac{\pi^2}{6} = 1 + \frac{1}{2^2} + \frac{1}{3^2} + \cdots.
$$

Say $s_k$ is the sum of the first $k$ terms of the series above, and $p_k = \sqrt{6s_k}$. Here is a fancy way to compute these sequences in a compact code.


```{code-cell}
a = [1/k^2 for k=1:100] 
s = cumsum(a)        # cumulative summation
p = @. sqrt(6*s)

scatter(1:100, p;
    title="Sequence convergence",
    xlabel=L"k",  ylabel=L"p_k")
```

This graph suggests that maybe $p_k\to \pi$, but it's far from clear how close the sequence gets. It's more informative to plot the sequence of errors, $\epsilon_k= |\pi-p_k|$. By plotting the error sequence on a log-log scale, we can see a nearly linear relationship.

```{code-cell}
ϵ = @. abs(π - p)    # error sequence
scatter(1:100, ϵ;
    title="Convergence of errors",
    xaxis=(:log10,L"k"),  yaxis=(:log10,"error"))
```

The straight line on the log-log scale suggests a power-law relationship where $\epsilon_k\approx a k^b$, or $\log \epsilon_k \approx b (\log k) + \log a$.

```{code-cell}
k = 1:100
V = [ k.^0 log.(k) ]     # fitting matrix
c = V \ log.(ϵ)          # coefficients of linear fit
```

In terms of the parameters $a$ and $b$ used above, we have

```{code-cell}
a, b = exp(c[1]), c[2];
@show b;
```

It's tempting to conjecture that the slope $b\to -1$ asymptotically. Here is how the numerical fit compares to the original convergence curve.

```{code-cell}
plot!(k, a * k.^b, l=:dash, label="power-law fit")
```

::::

## Exercises

``````{exercise}
:label: problem-fitting-minimize 
✍ Suppose $f$ is a twice-differentiable, nonnegative real function. Show that if there is an $x^*$ such that $f'(x^*)=0$ and $f''(x^*)>0$, then $x^*$ is a local minimizer of the function $[f(x)]^2$. (This justifies dropping the square root when minimizing in the least-squares problem.)

``````

``````{exercise}
:label: problem-fitting-census
⌨  Here are counts of the U.S. population in millions from the census performed every ten years, beginning with 1790 and ending with 2010.

``` julia
3.929, 5.308, 7.240, 9.638, 12.87, 17.07, 23.19, 31.44, 39.82, 50.19, 62.95, 76.21,
92.22, 106.0, 122.8, 132.2, 150.7, 179.3, 203.3, 226.5, 248.7, 281.4, 308.7
```

**(a)** Find a best-fitting cubic polynomial for the data. Plot the data as points superimposed on a (smooth) graph of the cubic over the full range of time. Label the axes. What does the fit predict for the population in the years 2000, 2010, and 2020? (In MATLAB, )

**(b)** Look up the actual U.S. population in 2000, 2010, and 2020 and compare to the predictions of part (a).
``````

``````{exercise}
:label: problem-fitting-movie
⌨  The following are weekly box office earnings (in dollars) in the U.S. for the 2012 film *The Hunger Games*. (Source: boxofficemojo.com.)

```julia
189932838, 79406327, 46230374, 26830921, 18804290,
13822248, 7474688, 6129424, 4377675, 3764963, 2426574,
1713298, 1426102, 1031985, 694947, 518242, 460578, 317909
```

Fit these values to a function of the form $y(t)\approx a e^{b t}$. Plot the data together with the fit using standard linear scales on the axes, and then plot them again using a log scale on the vertical axis.
``````

``````{exercise}
:label: problem-fitting-periodic
⌨  In this problem you are trying to find an approximation to the periodic function $g(t)=e^{\sin(t-1)}$ over one period, $0 < t \le 2\pi$. As data, define

```{math}
:numbered: false
t_i = \frac{2\pi i}{60}, \quad y_i = g(t_i), \quad i=1,\ldots,60.
```

**(a)** Find the coefficients of the least-squares fit

```{math}
:numbered: false
y(t) \approx c_1 + c_2t + \cdots + c_7 t^6.
```

Superimpose a plot of the data values as points with a curve showing the fit.

**(b)** Find the coefficients of the least-squares fit

```{math}
:numbered: false
y \approx d_1 + d_2\cos(t) + d_3\sin(t) + d_4\cos(2t) + d_5\sin(2t).
```

Unlike part (a), this fitting function is itself periodic. Superimpose a plot of the data values as points with a curve showing the fit. (It should look better than the polynomial fit, especially near the endpoints.)
``````

``````{exercise}
:label: problem-fitting-tanh
⌨ Define a vector with elements $0, 0.5, 1, 1.5, \dots, 10$, and another vector with elements $y_i = \tanh(t_i)$. 

**(a)** Fit the data to a cubic polynomial. Plot the data together with the polynomial fit over the interval $0 \le t \le 10$. (It will not be a good fit.)

**(b)** Fit the data to the function $c_1 + c_2z + c_3z^2 + c_4z^3$, where $z=t^2/(1+t^2)$. Plot the data together with the fit. What feature of the new variable $z$ makes this fit much better than the original cubic?

``````

``````{exercise}
:label: problem-fitting-pi
⌨  One series for finding $\pi$ is

```{math}
:numbered: false
\frac{\pi}{2} = 1 + \frac{1}{3} + \frac{1\cdot 2}{3\cdot5} + \frac{1\cdot 2\cdot 3}{3\cdot 5\cdot 7} + \cdots.
```

Define $s_k$ to be the sum of the first $k$ terms on the right-hand side, and let $e_k=|s_k-\pi/2|$. 

**(a)** Calculate $e_k$ for $k=1,\ldots,20$, and plot the sequence on a log-linear scale.

**(b)** Determine $a$ and $b$ in a least-squares fit $e_k \approx a \cdot b^k$, and superimpose the fit on the plot from part (a).
``````

``````{exercise}
:label: problem-fitting-planets
 ⌨  Kepler found that the orbital period $\tau$ of a planet depends on its mean distance $R$ from the sun according to $\tau=c R^{\alpha}$ for a simple rational number $\alpha$. Perform a linear least-squares fit from the following table in order to determine the most likely simple rational value of $\alpha$.

| Planet  | Distance from sun in Mkm |  Orbital period in days              |
|---------|:--------------------------------------------:|:-----------:|
| Mercury | 57.9                                      | 88.0  |
| Venus   | 108.2                                     | 224.7  |
| Earth   | 149.6                                     | 365.2 |
| Mars    | 228.0                                     | 687.0 |
| Jupiter | 778.5                                     | 4331 |
| Saturn  | 1432                                        | 10750   |
| Uranus  | 2867                                     | 30590   |
| Neptune | 4515                                     | 59800   |
``````

``````{exercise}
:label: problem-fitting-inverse
 ✍ Show that finding a fit of the form

```{math}
:numbered: false
y(t) \approx \frac{a}{t+b}
```

can be transformed into a linear fitting problem (with different undetermined coefficients) by rewriting the equation.
``````

``````{exercise}
:label: problem-fitting-lft
 ✍ Show how to find the constants $a$ and $b$ in a data fitting problem of the form $y(t)\approx t/(at+b)$ for $t>1$ by transforming it into a linear least-squares fitting problem.
``````
