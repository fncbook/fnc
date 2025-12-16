---
numbering:
  enumerator: 5.6.%s
---
(section-localapprox-integration)=

# Numerical integration

In calculus, you learn that the elegant way to evaluate a definite integral is to apply the Fundamental Theorem of Calculus and find an antiderivative. The connection is so profound and pervasive that it's easy to overlook that a definite integral is a numerical quantity existing independently of antidifferentiation.  However, most conceivable integrands have no antiderivative in terms of familiar functions.

::::{prf:example} Numerical integration
:label: demo-int-antideriv

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-int-antideriv-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-int-antideriv-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-int-antideriv-python
:::
```` 
`````

::::

```{index} ! numerical integration
```

```{index} see: quadrature; numerical integration
```

Numerical integration, which also goes by the older name *quadrature*, is performed by combining values of the integrand sampled at nodes. In this section we will assume equally spaced nodes using the definitions

```{math}
:label: nc-nodes
  t_i = a +i h, \quad h=\frac{b-a}{n}, \qquad i=0,\ldots,n.
```

::::{prf:definition} Numerical integration formula
:label: definition-integrationformula
A {term}`numerical integration formula` is a list of **weights** $w_0,\ldots,w_n$ chosen so that for all $f$ in some class of functions,

```{math}
:label: quadrature
  \begin{split}
    \int_a^b f(x)\, dx \approx h \sum_{i=0}^n w_if(t_i) =  h \bigl[ w_0f(t_0)+w_1f(t_1)+\cdots w_nf(t_n) \bigr],
  \end{split}
```

with the $t_i$ defined in {eq}`nc-nodes`. The weights are independent of $f$ and $h$.
::::

Numerical integration formulas can be applied to sequences of data values even if no function is explicitly known to generate them. For our presentation and implementations, however, we assume that $f$ is known and can be evaluated anywhere.

```{index} interpolation; by piecewise polynomials
```

A straightforward way to derive integration formulas is to mimic the approach taken for finite differences: find an interpolant and operate exactly on it.

## Trapezoid formula

One of the most important integration formulas results from integration of the piecewise linear interpolant (see {numref}`section-localapprox-pwlin`). Using the cardinal basis form of the interpolant in {eq}`plbasis`, we have

```{math}
\int_a^b f(x) \, dx \approx \int_a^b \sum_{i=0}^n f(t_i) H_i(x)\, dx = \sum_{i=0}^n f(t_i) \left[ \int_a^b H_i(x)\right]\, dx.
```

Thus, we can identify the weights as $w_i = h^{-1} \int_a^b H_i(x)\, dx$. Using areas of triangles, it's trivial to derive that

```{math}
:label: hatintegral
w_i = \begin{cases}
1, & i=1,\ldots,n-1,\\
\frac{1}{2}, & i=0,n.
\end{cases}
```

Putting everything together, the resulting formula is

```{math}
:label: trapezoid
\begin{split}
  \int_a^b f(x)\, dx \approx T_f(n) &= h\left[
    \frac{1}{2}f(t_0) + f(t_1) + f(t_2) + \cdots + f(t_{n-1}) +
    \frac{1}{2}f(t_n) \right].
\end{split}
```

```{index} ! trapezoid formula; for integration
```

::::{prf:definition} Trapezoid formula
:label: definition-trapezoidint
The {term}`trapezoid formula` is a numerical integration formula in the form {eq}`quadrature`, with

$$
w_i = \begin{cases}
  \frac{1}{2},& i=0 \text{ or } i=n, \\
  1, & 0 < i < n.
  \end{cases}
$$
::::

Geometrically, as illustrated in {numref}`fig-trapezoid`, the trapezoid formula sums of the areas of trapezoids approximating the region under the curve $y=f(x)$.[^comp]

[^comp]: Some texts distinguish between a formula for a single subinterval $[t_{k-1},t_k]$ and a *composite* formula that adds them up over the whole interval to get {eq}`trapezoid`.

The trapezoid formula is the Swiss Army knife of integration formulas. A short implementation is given as {numref}`Function {number} <function-trapezoid>`.

```{figure} figures/trapezoid.svg
:name: fig-trapezoid
Trapezoid formula for integration. The piecewise linear interpolant defines trapezoids that approximate the region under the curve.
```

``````{prf:algorithm} trapezoid
:label: function-trapezoid

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #function-trapezoid-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-trapezoid-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #function-trapezoid-python
:::
````
`````
``````

Like finite-difference formulas, numerical integration formulas have a truncation error.

```{index} ! truncation error; of a numerical integration formula, ! order of accuracy; of numerical integration
```

::::{prf:definition} Truncation error of a numerical integration formula
:label: definition-truncationint
For the numerical integration formula {eq}`quadrature`, the {term}`truncation error` is

```{math}
:label: truncint
\tau_f(h) = \int_a^b f(x) \, dx - h \sum_{i=0}^{n} w_i f(t_i).
```

The {term}`order of accuracy` is as defined in {numref}`Definition {number} <definition-orderfd>`.
::::

In @theorem-pwlin-converge we stated that the pointwise error in a piecewise linear interpolant with equal node spacing $h$ is bounded by $O(h^2)$ as $h\rightarrow 0$. Using $I$ to stand for the exact integral of $f$ and $p$ to stand for the piecewise linear interpolant, we obtain

```{math}
\begin{split}
  I - T_f(n) = I - \int_a^b p(x)\, dx &= \int_a^b \bigl[f(x)-p(x)\bigr] \, dx \\
  &\le (b-a) \max_{x\in[a,b]} |f(x)-p(x)| = O(h^2).
\end{split}
```

```{index} ! Euler–Maclaurin formula
```

A more thorough statement of the truncation error is known as the **Euler–Maclaurin formula**,

```{math}
:label: eulermaclaurin
\begin{split}
\int_a^b f(x)\, dx &= T_f(n) - \frac{h^2}{12} \left[ f'(b)-f'(a) \right] + \frac{h^4}{740} \left[ f'''(b)-f'''(a) \right] + O(h^6) \\
    &= T_f(n) - \sum_{k=1}^\infty \frac{B_{2k}h^{2k}}{(2k)!}  \left[ f^{(2k-1)}(b)-f^{(2k-1)}(a) \right],
\end{split}
```

where the $B_{2k}$ are constants known as *Bernoulli numbers*. Unless we happen to be fortunate enough to have a function with $f'(b)=f'(a)$, we should expect truncation error at second order and no better.

::::{prf:observation}
The trapezoid integration formula is second-order accurate.
::::

::::{prf:example} Trapezoid integration
:label: demo-int-trap

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-int-trap-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-int-trap-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-int-trap-python
:::
```` 
`````

::::

(section-integration-extrapolation)=

## Extrapolation

If evaluations of $f$ are computationally expensive, we want to get as much accuracy as possible from them by using a higher-order formula. There are many routes for doing so; for example, we could integrate a not-a-knot cubic spline interpolant. However, splines are difficult to compute by hand, and as a result different methods were developed before computers came on the scene.

```{index} ! extrapolation
```

Knowing the structure of the error allows the use of {term}`extrapolation` to improve accuracy. Suppose a quantity $A_0$ is approximated by an algorithm $A(h)$ with an
error expansion

```{math}
:label: extraperror
  A_0 = A(h) + c_1 h + c_2 h^2 + c_3 h^3 + \cdots.
```

Crucially, it is not necessary to know the values of the error constants $c_k$, merely that they exist and are independent of $h$.

Using $I$ for the exact integral of $f$, the trapezoid formula has

```{math}
  I = T_f(n) + c_2 h^2 + c_4 h^{4} + \cdots,
```

as proved by the Euler–Maclaurin formula {eq}`eulermaclaurin`. The error constants depend on $f$ and can't be evaluated in general, but we know that this expansion holds. For convenience, we recast the error expansion in terms of $n=O(h^{-1})$:

```{math}
:label: traperrorexpansion
  I = T_f(n) + c_2 n^{-2} + c_4 n^{-4} + \cdots.
```

We now make the simple observation that

```{math}
:label: traperrorexpansion2n
  I = T_f(2n) + \tfrac{1}{4} c_2 n^{-2} + \tfrac{1}{16} c_4 n^{-4} + \cdots.
```

It follows that if we combine {eq}`traperrorexpansion` and {eq}`traperrorexpansion2n` correctly, we can cancel out the second-order term in the error. Specifically, define

```{math}
:label: nc-simpson
  S_f(2n) = \frac{1}{3} \Bigl[ 4 T_f(2n) - T_f(n) \Bigr].
```

(We associate $2n$ rather than $n$ with the extrapolated result because of the total number of nodes needed.) Then

```{math}
:label: extraplevel1
  I = S_f(2n) + O(n^{-4}) =  b_4 n^{-4} + b_6 n^{-6} + \cdots.
```

```{index} ! Simpson's formula
```

The formula {eq}`nc-simpson` is called **Simpson's formula**, or *Simpson's rule*. A different presentation and derivation are considered in @problem-integration-simpson.

Equation {eq}`extraplevel1` is another particular error expansion in the form {eq}`extraperror`, so we can extrapolate again! The details change only a little. Considering that

```{math}
  I = S_f(4n) = \tfrac{1}{16} b_4 n^{-4} + \tfrac{1}{64} b_6 n^{-6} + \cdots,
```

the proper combination this time is

```{math}
:label: nc-sixth
  R_f(4n) = \frac{1}{15} \Bigl[ 16 S_f(4n) - S_f(2n) \Bigr],
```

which is sixth-order accurate. Clearly the process can be repeated to get eighth-order accuracy and beyond. Doing so goes by the name of **Romberg integration**, which we will not present in full generality.

## Node doubling

Note in {eq}`nc-sixth` that $R_f(4n)$ depends on $S_f(2n)$ and $S_f(4n)$, which in turn depend on $T_f(n)$, $T_f(2n)$, and $T_f(4n)$.  There is a useful benefit realized by doubling of the nodes in each application of the trapezoid formula. As shown in {numref}`figure-node-doubling`, when doubling $n$, only about half of the nodes are new ones, and previously computed function values at the other nodes can be reused.

```{figure} figures/node-doubling.svg
:name: figure-node-doubling
Dividing the node spacing by half introduces new nodes only at midpoints, allowing the function values at existing nodes to be reused for extrapolation.
```

Specifically, we have

```{math}
:label: nc-doubling
\begin{split}
  T_f(2m) & = \frac{1}{2m} \left[  \frac{1}{2} f(a) + \frac{1}{2} f(b) + \sum_{i=1}^{2m-1}  f\Bigl( a + \frac{i}{2m} \Bigr) \right]\\[1mm]
  & = \frac{1}{2m} \left[  \frac{1}{2} f(a) + \frac{1}{2} f(b)\right] + \frac{1}{2m} \sum_{k=1}^{m-1}  f\Bigl( a+\frac{2k}{2m} \Bigr)  + \frac{1}{2m} \sum_{k=1}^{m} f\Bigl( a+\frac{2k-1}{2m} \Bigr) \\[1mm]
  &=  \frac{1}{2m} \left[  \frac{1}{2} f(a) + \frac{1}{2} f(b) + \sum_{k=1}^{m-1} f\Bigl( a+\frac{k}{m} \Bigr) \right] + \frac{1}{2m} \sum_{k=1}^{m}  f\Bigl( a+\frac{2k-1}{2m} \Bigr)  \\[1mm]
  &= \frac{1}{2} T_f(m) + \frac{1}{2m} \sum_{k=1}^{m-1}  f\left(t_{2k-1} \right),
\end{split}
```

where the nodes referenced in the last line are relative to $n=2m$. Hence in passing from $n=m$ to $n=2m$, new integrand evaluations are needed only at the odd-numbered nodes of the finer grid.

::::{prf:example} Integration by extrapolation
:label: demo-int-extrap

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-int-extrap-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-int-extrap-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-int-extrap-python
:::
```` 
`````

::::

## Exercises

% must be kept as #1

``````{exercise}
:label: problem-integration-tests
⌨ For each integral below, use {numref}`Function {number} <function-trapezoid>` to estimate the integral for $n=10\cdot 2^k$ nodes for $k=1,2,\ldots,10$. Make a log-log plot of the errors and confirm or refute second-order accuracy. (These integrals were taken from {cite}`baileyComparisonThree2005`.)

**(a)** $\displaystyle \int_0^1 x\log(1+x)\, dx = \frac{1}{4}$

**(b)** $\displaystyle \int_0^1 x^2 \tan^{-1}x\, dx = \frac{\pi-2+2\log 2}{12}$

**(c)** $\displaystyle \int_0^{\pi/2}e^x \cos x\, dx = \frac{e^{\pi/2}-1}{2}$

**(d)** $\displaystyle \int_0^1 \sqrt{x} \log(x) \, dx = -\frac{4}{9}$ (Note: Although the integrand has the limiting value zero as $x\to 0$, it cannot be evaluated naively at $x=0$. You can start the integral at $x=\macheps$ instead.)

**(e)** $\displaystyle \int_0^1 \sqrt{1-x^2}\,\, dx = \frac{\pi}{4}$
``````

```{index} ! Gregory integration formula
```

``````{exercise}
:label: problem-integration-gregory
✍ The Euler–Maclaurin error expansion {eq}`eulermaclaurin` for the trapezoid formula implies that if we could cancel out the term due to $f'(b)-f'(a)$, we would obtain fourth-order accuracy. We should not assume that $f'$ is available, but approximating it with finite differences can achieve the same goal. Suppose the forward difference formula {eq}`forwardFD21` is used for $f'(a)$, and its reflected backward difference is used for $f'(b)$. Show that the resulting modified trapezoid formula is

```{math}
:label: gregory
G_f(h) = T_f(h) - \frac{h}{24} \left[ 3\Bigl( f(t_n)+f(t_0) \Bigr) -4\Bigr( f(t_{n-1}) + f(t_1) \Bigr) + \Bigl( f(t_{n-2})+f(t_2)   \Bigr) \right],
```

which is known as a **Gregory integration formula**.
``````

``````{exercise}
:label: problem-integration-gregoryuse
⌨ For the integrals in (a)–(c) of @problem-integration-tests, use Gregory integration {eq}`gregory` instead of the trapezoid formula. Compare the observed errors to fourth-order convergence.

``````

``````{exercise}
:label: problem-integration-simpson
✍  Simpson's formula can be derived without appealing to extrapolation.

**(a)** Show that

```{math}
:numbered: false
p(x) = \beta + \frac{\gamma-\alpha}{2h}\, x + \frac{\alpha-2\beta+\gamma}{2h^2}\, x^2
```

interpolates the three points $(-h,\alpha)$, $(0,\beta)$, and $(h,\gamma)$.

**(b)** Find

```{math}
:numbered: false
\int_{-h}^h p(s)\, ds,
```

where $p$ is the quadratic polynomial from part (a), in terms of $h$, $\alpha$, $\beta$, and $\gamma$.

**(c)** Assume equally spaced nodes in the form $t_i=a+ih$, for $h=(b-a)/n$ and $i=0,\ldots,n$. Suppose $f$ is approximated by $p(x)$ over the subinterval $[t_{i-1},t_{i+1}]$. Apply the result from part (b) to find

```{math}
\int_{t_{i-1}}^{t_{i+1}} f(x)\, dx \approx \frac{h}{3} \bigl[ f(t_{i-1}) + 4f(t_i) + f(t_{i+1}) \bigr].
```

(Use the change of variable $s=x-t_i$.)

**(d)** Now also assume that $n=2m$ for an integer $m$. Derive Simpson's formula,

```{math}
:label: simpson
\begin{split}
\int_a^b f(x)\, dx \approx  \frac{h}{3}\bigl[ &f(t_0) + 4f(t_1) + 2f(t_2) + 4f(t_3) + 2f(t_4) + \cdots\\
&+ 2f(t_{n-2}) + 4f(t_{n-1}) + f(t_n) \bigr].
\end{split}
```
``````

``````{exercise}
:label: problem-integration-simpsonextrap
✍ Show that the Simpson formula {eq}`simpson` is equivalent to $S_f(n),$ given the definition of $S_f$ in {eq}`nc-simpson`.
``````

``````{exercise}
:label: problem-integration-simpsonconverge
⌨ For the integrals in (a)–(c) of @problem-integration-tests, apply the Simpson formula {eq}`simpson` and compare the errors to fourth-order convergence.
``````

``````{exercise}
:label: problem-integration-sublinear
⌨ For $n=10,20,30,\ldots,200$, compute the trapezoidal approximation to

```{math}
\int_{0}^1 \frac{1}{2.01+\sin (6\pi x)-\cos(2\pi x)} \,d x \approx 0.9300357672424684.
```

Make two separate plots of the absolute error as a function of $n$, one using a log-log scale and the other using log-linear. The graphs suggest that the error asymptotically behaves as $C \alpha^n$ for some $C>0$ and some $0<\alpha<1$. How does this result relate to {eq}`eulermaclaurin`?

``````

``````{exercise}
:label: problem-integration-extrap6
⌨ For the integrals in (a)–(c) of @problem-integration-tests, extrapolate the trapezoidal results for $n=32,48,72,108,162$ two levels to get sixth-order accurate results, and make table of the errors.
``````

``````{exercise}
:label: problem-integration-extrap8
✍ Find a formula like {eq}`nc-sixth` that extrapolates two values of $R_f$ to obtain an eighth-order accurate one.
``````
