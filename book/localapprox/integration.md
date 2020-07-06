# Numerical integration

```{index} numerical integration
```

````{sidebar} Demo
:class: demo
{doc}`demos/int-antideriv`
````

In calculus you learn that the elegant way to evaluate a definite integral is to apply the Fundamental Theorem of Calculus and find an antiderivative. The connection is so profound and pervasive that it's easy to overlook that a definite integral is a numerical quantity existing independently of antidifferentiation.  In fact, most conceivable integrands have no antiderivative in terms of familiar functions.

Numerical integration[^quad] is done by combining values of the integrand sampled at nodes, much like finite differences. In this section we will assume equally spaced nodes using the definitions

[^quad]: Numerical integration also goes by the older name **quadrature**.

```{math}
  :label: nc-nodes
  t_i = a +i h, \quad h=\frac{b-a}{n}, \qquad i=0,\ldots,n.
```

The integration formulas are expressed as

```{math}
  :label: quadrature
  \begin{split}
    I = \int_a^b f(x)\, dx \approx Q &= h \sum_{i=0}^n w_if(t_i) \\
    &=  h \bigl[ w_0f(t_0)+w_1f(t_1)+\cdots w_nf(t_n) \bigr].
  \end{split}
```

```{margin}
The weights of numerical integration formulas are chosen independently of the function being integrated.
```

The constants $w_i$ appearing in the formula are called **weights**.  As with finite difference formulas, the weights of numerical integration formulas are chosen independently of the function being integrated, and they determine the formula completely. We can apply quadrature formulas to sequences of data values even if no function is explicitly known to generate them, but for presentation and implementations we assume that we can evaluate $f(x)$ anywhere.

```{index} interpolation; by piecewise polynomials
```

```{index} Newton--Cotes formula
```

A straightforward way to derive integration formulas is to mimic the approach taken for finite differences: find an interpolant and operate exactly on it. If the interpolant is a piecewise polynomial, the result is a {term}`Newton--Cotes formula`.

## Trapezoid formula

One of the most important Newton--Cotes formulas results from integration of the piecewise linear interpolant (see {doc}`pwlin`). Using the cardinal basis form of the interpolant in {eq}`plbasis`, we have

```{math}
I \approx \int_a^b \sum_{i=0}^n f(t_i) H_i(x)\, dx = \sum_{i=0}^n f(t_i) \left[ \int_a^b H_i(x)\right]\, dx.
```

Thus we can identify the weights as $w_i = h^{-1} \int_a^b H_i(x)\, dx$. Using areas of triangles, it's trivial to derive that

```{math}
    :label: hatintegral
    w_i = \begin{cases}
    1, & i=1,\ldots,n-1,\\
    \frac{1}{2}, & i=0,n.
    \end{cases}
```

Putting everything together, the resulting quadrature formula is

```{math}
:label: trapezoid
\begin{split}
  I = \int_a^b f(x)\, dx \approx T_f(n) &= h\left[
    \frac{1}{2}f(t_0) + f(t_1) + f(t_2) + \cdots + f(t_{n-1}) +
    \frac{1}{2}f(t_n) \right].
\end{split}
```

```{margin}
The trapezoid formula results from integration of the piecewise linear interpolant.
```

```{index} trapezoid formula; for integration
```

This is called the {term}`trapezoid formula` or trapezoid rule.[^comp] The trapezoid formula results from integration of the piecewise linear interpolant, or equivalently, as illustrated in {numref}`fig-trapezoid`, from using the area of approximating trapezoids to estimate the area under a curve. The trapezoid formula is the Swiss Army knife of integration formulas. A short implementation is given as {ref}`function-trapezoid`.

[^comp]: Some texts distinguish between a formula for a single subinterval $[t_{k-1},t_k]$ and a "composite" formula that adds them up over the whole interval to get something like our {eq}`trapezoid`.

```{figure} figures/trapezoid.svg
:name: fig-trapezoid
Trapezoid formula for integration.
```

(function-trapezoid)=

````{proof:function} trapezoid
**Trapezoid formula for numerical integration.**

```{code-block} julia
:lineno-start: 1
"""
trapezoid(f,a,b,n)

Apply the trapezoid integration formula for integrand `f` over
interval [`a`,`b`], broken up into `n` equal pieces. Returns
estimate, vector of nodes, and vector of integrand values at the
nodes.
"""
function trapezoid(f,a,b,n)
    h = (b-a)/n
    t = LinRange(a,b,n+1)
    y = f.(t)
    T = h * ( sum(y[2:n]) + 0.5*(y[1] + y[n+1]) )

    return T,t,y
end
```
````

In the [PL convergence theorem](theorem-placcuracy) we stated that the pointwise error in a piecewise linear interpolant with equal node spacing $h$ is bounded by $O(h^2)$ as $h\rightarrow 0$. Using $p$ to stand for the piecewise linear interpolant, we obtain

```{math}
\begin{split}
  I - T_f(n) = I - \int_a^b p(x)\, dx &= \int_a^b \bigl[f(x)-p(x)\bigr] \, dx \\
  &\le (b-a) \max_{x\in[a,b]} |f(x)-p(x)| = O(h^2).
\end{split}
```

```{margin}
Tne trapezoid formula has order of accuracy equal to two.
```

Hence the trapezoid formula has second-order error. This fact is embedded rigorously in "one of the most remarkable formulas in mathematics," the **Euler--Maclaurin formula**,  which may be stated as

```{math}
    :label: eulermaclaurin
    \begin{split}
    I = \int_a^b f(x)\, dx &= T_f(n) - \frac{h^2}{12} \left[ f'(b)-f'(a) \right] + \frac{h^4}{740} \left[ f'''(b)-f'''(a) \right] + O(h^6) \\
        &= T_f(n) - \sum_{k=1}^\infty \frac{B_{2k}h^{2k}}{(2k)!}  \left[ f^{(2k-1)}(b)-f^{(2k-1)}(a) \right],
      \end{split}
```

````{sidebar} Demo
:class: demo
{doc}`demos/int-trap`
````

where the $B_{2k}$ are constants known as **Bernoulli numbers**. Unless we happen to be fortunate enough to have a function with $f'(b)=f'(a)$, we should expect truncation error at second order and no better.

## Extrapolation

If evaluations of $f$ are computationally expensive, we want to get as much accuracy as possible from them by using a higher-order formula. There are many routes for doing so; for example, we could integrate a not-a-knot cubic spline interpolant. However, splines are difficult to compute by hand, and as a result different methods were developed before computers came on the scene.

```{index} extrapolation
```

Knowing the structure of the error allows the use of {term}`extrapolation` to improve accuracy. Suppose a quantity $A_0$ is approximated by an algorithm $A(h)$ with an
error expansion

```{math}
  :label: extraperror
  A_0 = A(h) + c_1 h + c_2 h^2 + c_3 h^3 + \cdots.
```

Crucially, it is not necessary to know the values of the error constants $c_k$, merely that they exist and are independent of $h$. In the case of the trapezoid formula, we have

```{math}
  I = T_f(n) + c_2 h^2 + c_4 h^{4} + \cdots,
```

as proved by the Euler--Maclaurin formula {eq}`eulermaclaurin`. The error constants depend on $f$ and can't be evaluated in general, but we know that this expansion holds.

For convenience we recast the error expansion in terms of $n=O(h^{-1})$:

```{math}
  :label: traperrorexpansion
  I = T_f(n) + c_2 n^{-2} + c_4 n^{-4} + \cdots,
```

We make the simple observation that

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

```{index} Simpson's formula
```

The formula {eq}`nc-simpson` is called **Simpson's formula**. A different presentation and derivation are considered in [an exercise](problem-simpson).

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

Note in {eq}`nc-sixth` that $R_f(4n)$ depends on $S_f(2n)$ and $S_f(4n)$, which in turn depend on $T_f(n)$, $T_f(2n)$, and $T_f(4n)$.  There is a useful benefit realized by doubling of the nodes in each application of the trapezoid formula. For simplicity, suppose that $[a,b]=[0,1]$ and that $n=2m$ for some positive integer $m$.  The nodes are

```{math}
\Bigl[\, 0, \;\frac{1}{2m}, \;\frac{2}{2m}, \;\frac{3}{2m}, \;\frac{4}{2m}, \;\ldots \;\frac{2m-3}{2m}, \;\frac{2m-2}{2m}, \;\frac{2m-1}{2m}, \; 1 \,\Bigr].
```

Suppose we delete every other node:

```{math}
\Bigl[\, 0, \;\frac{2}{2m}, \;\frac{4}{2m}, \;\ldots \;\frac{2m-2}{2m}, \; 1 \,\Bigr]
\quad = \quad \Bigl[\, 0, \;\frac{1}{m}, \;\frac{2}{m}, \;\ldots \;\frac{m-1}{m}, \; 1 \,\Bigr] .
```

What remains are the nodes with $n=m$. That is, if we have computed $T_f(m)$ and want to compute $T_f(2m)$, we begin with half of the evaluations of $f$ already in our pocket. More specifically,

```{math}
:label: nc-doubling
\begin{split}
  T_f(2m) & = \frac{1}{2m} \left[  \tfrac{1}{2} f(0) + \tfrac{1}{2} f(1) + \sum_{k=1}^{m-1}  f\Bigl( \tfrac{2k-1}{2m} \Bigr) +  f\Bigl( \tfrac{2k}{2m} \Bigr) \right] \\
  &=  \frac{1}{2m} \left[  \tfrac{1}{2} f(0) + \tfrac{1}{2} f(1) + \sum_{k=1}^{m-1} f\Bigl( \tfrac{k}{m} \Bigr) \right] + \frac{1}{2m} \sum_{k=1}^{m}  f\Bigl( \tfrac{2k-1}{2m} \Bigr)  \\
  &= \frac{1}{2} T_f(m) + \frac{1}{2m} \sum_{k=1}^{m-1}  f\left( t_{2k-1} \right),
\end{split}
```

````{sidebar} Demo
:class: demo
{doc}`demos/int-extrap`
````

where the nodes referenced in the last line are relative to $n=2m$. To summarize: when $n$ is doubled, new integrand evaluations are needed only at the odd-numbered nodes of the finer grid. Although we derived this result in the particular interval $[0,1]$, it is valid for any interval.

## Exercises

(problem-quadraturetests)=
% must be kept as #1

1. ⌨ For each integral below, use {ref}`function-trapezoid` to estimate the integral for $n=10\cdot 2^k$ nodes for $k=1,2,\ldots,10$. Make a log--log plot of the errors and confirm or refute second-order accuracy. (These integrals were taken from {cite}`baileyComparisonThree2005`.)

    **(a)** $\displaystyle \int_0^1 x\log(1+x)\, dx = \frac{1}{4}$

    **(b)** $\displaystyle \int_0^1 x^2 \tan^{-1}x\, dx = \frac{\pi-2+2\log 2}{12}$

    **(c)** $\displaystyle \int_0^{\pi/2}e^x \cos x\, dx = \frac{e^{\pi/2}-1}{2}$

    **(d)** $\displaystyle \int_0^1 \sqrt{x} \log(x) \, dx = -\frac{4}{9}$ (Note: Although the integrand has the limiting value zero as $x\to 0$, it cannot be evaluated naively at $x=0$. You can start the integral at $x=\macheps$ instead.)

    **(e)** $\displaystyle \int_0^1 \sqrt{1-x^2}\,\, dx = \frac{\pi}{4}$
  
    ````{only} solutions
    f = {
        @(x) x.*log(1+x);
        @(x) x.^2.*atan(x);
        @(x) exp(x).*cos(x);
        @(x) (x>0).*sqrt(x).*log(x);
        @(x) sqrt(1-x.^2);
        };

    a = [ 0;0;0;eps;0];
    b = [1;1;pi/2;1;1];
    I = [.25;(pi-2+2*log(2))/12;(exp(pi/2)-1)/2;-4/9;pi/4];

    n_ = 10*2.^(1:10)';
    err_ = [];

    for n = n_'
        err = [];
        for i = 1:length(f)
            err = [err, I(i) - trapezoid(f{i},a(i),b(i),n)];
        end
        err_ = [err_; err];
    end

    loglog(n_,abs(err_),'o-')
    ````

2. ✍ The Euler--Maclaurin error expansion {eq}`eulermaclaurin` for the trapezoid formula implies that if we could cancel out the term due to $f'(b)-f'(a)$, we would obtain fourth-order accuracy. We should not assume that $f'$ is available, but approximating it with finite differences can achieve the same goal. Suppose the forward difference formula {eq}`forwardFD21` is used for $f'(a)$, and its reflected backward difference is used for $f'(b)$. Show that the resulting modified trapezoid formula is

    ```{math}
      :label: gregory
        G_f(h) = T_f(h) - \frac{h}{24} \left[ 3\Bigl( f(x_n)+f(x_0) \Bigr) -4\Bigr( f(x_{n-1}) + f(x_1) \Bigr) + \Bigl( f(x_{n-2})+f(x_2)   \Bigr) \right],
    ```

    ```{index} Gregory integration formula
    ```

    which is known as a **Gregory integration formula**.

    ````{only} solutions
    ````

3. ⌨ Repeat each integral in exercise 1 above using Gregory integration {eq}`gregory` instead of the trapezoid formula. Compare the observed errors to fourth-order convergence.

    ````{only} solutions
    f = {
        @(x) x.*log(1+x);
        @(x) x.^2.*atan(x);
        @(x) exp(x).*cos(x);
        @(x) sqrt(x).*log(x);
        @(x) sqrt(1-x.^2);
        };

    a = [ 0;0;0;eps;0];
    b = [1;1;pi/2;1;1];
    I = [.25;(pi-2+2*log(2))/12;(exp(pi/2)-1)/2;-4/9;pi/4];

    n_ = 10*2.^(1:10)';
    err_ = [];

    for n = n_'
        err = [];
        for i = 1:length(f)
            [T,t,y] = trapezoid(f{i},a(i),b(i),n);
            h = (b(i)-a(i))/n;
            delta = [3;-4;1] .* (y([n+1 n n-1]) + y([1 2 3]));
            G = T - h/24 * sum(delta);

            err = [err, I(i) - G];
        end
        err_ = [err_; err];
    end

    loglog(n_,abs(err_),'o-')
    ````

    (problem-simpson)=
4. ✍  Simpson's formula can be derived without appealing to extrapolation.

    **(a)** Show that

    ```{math}
    p(x) = \beta + \frac{\gamma-\alpha}{2h}\, x + \frac{\alpha-2\beta+\gamma}{2h^2}\, x^2
    ```

    interpolates the three points $(-h,\alpha)$, $(0,\beta)$, and $(h,\gamma)$.

    **(b)** Find

    ```{math}
      \int_{-h}^h p(s)\, ds,
    ```

    where $p$ is the quadratic polynomial from part~(a), in terms of $h$, $\alpha$, $\beta$, and $\gamma$.
  
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

    ````{only} solutions
    ````

5. ✍ Show that the Simpson formula {eq}`simpson` is equivalent to $S_f(n/2)$, given the definition of $S_f$ in {eq}`nc-simpson`.

    ````{only} solutions
    ````

6. ⌨ For each integral in exercise 1 above, apply the Simpson formula {eq}`simpson` and compare the errors to fourth-order convergence.

    ````{only} solutions
    ````

7. ⌨ For $n=10,20,30,\ldots,200$, compute the trapezoidal approximation to

    ```{math}
    \int_{0}^1 \frac{1}{2.01+\sin (6\pi x)-\cos(2\pi x)} \,d x \approx 0.9300357672424684.
    ```

    Make two separate plots of the absolute error as a function of $n$, one using log--log scales and the other using log only for the $y$-axis. The graphs suggest that the error asymptotically behaves as $C \alpha^n$ for some $C>0$ and some $0<\alpha<1$. How does this result relate to {eq}`eulermaclaurin`?

    ````{only} solutions
    % We apply the trapezoidal quadrature to the periodic function
    f = @(x) 1./(2.01 + sin(6*pi*x) - cos(2*pi*x));
    exact = 0.9300357672424684;

    n_ = 10:20:200;
    err_ = [];
    for n = n_
      h = 1/n;
      x = h*(0:n)';
      c = h*[ 1/2, ones(1,n-1), 1/2 ]';
      q = c'*f(x);  % using an inner product
      err_ = [ err_; abs(exact-q) ];
    end

    %%
    clf
    subplot(1,2,1), loglog(n_,err_,'o-')
    xlabel('n'), ylabel('error in trapezoidal')
    subplot(1,2,2), semilogy(n_,err_,'o-')
    xlabel('n'), ylabel('error in trapezoidal')

    %%
    % The straight line observed on the semilog scale indicates that the
    % convergence is geometric, like $\alpha^n$.
    ````

8. ⌨ For each integral in exercise 1 above, extrapolate the trapezoidal results two levels to get sixth-order accurate results, and compare the expected convergence rate to the observed errors.

    ````{only} solutions
    ````

9. ✍ Find a formula like {eq}`nc-sixth` that extrapolates two values of $R_f$ to obtain an 8th-order accurate one.

    ````{only} solutions
    ````
