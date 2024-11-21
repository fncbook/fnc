---
jupytext:
  cell_metadata_filter: -all
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: Julia 1.7.1
  language: julia
  name: julia-fast
---

```{code-cell}
:tags: [remove-cell]
using FundamentalsNumericalComputation
FNC.init_format()
```

(section-localapprox-finitediffs)=
# Finite differences

Now we turn to one of the most common and important applications of interpolants: finding derivatives of functions. Because differentiation is a linear operation, we will constrain ourselves to formulas that are linear in the nodal values. 

```{index} ! finite differences
```

::::{proof:definition} Finite-difference formula
A **finite-difference** formula is a list of values $a_{-p},\ldots,a_q$, called **weights**, such that for all $f$ in some class of functions,

```{math}
  :label: fdformula
  f'(x) \approx \frac{1}{h} \sum_{k=-p}^{q} a_k f(x + kh).
```

The weights are independent of $f$ and $h$. The formula is said to be **convergent** if the approximation becomes equality in the limit $h\to 0$ for a suitable class of functions.
::::

Note that while {eq}`fdformula` is about finding the derivative at a single point $x$, the same formula can be applied for different $x$. The usual situation is a regularly spaced grid of nodes, $a,a+h,a+2h,\ldots,b$, and then the value of $f$ at each node takes part in multiple applications of the formula. This will be demonstrated in {numref}`Example {number} <example-finitediffs-fd1bd1>` below.

## Common examples

```{index} ! forward difference, ! backward difference, ! centered difference
```

There are three appealing special cases of {eq}`fdformula` that get special attention.

::::{proof:definition} Forward, backward, and centered FD formulas
A **forward difference formula** is characterized by {eq}`fdformula` with $p=0$, a **backward difference formula** has $q=0$, and a **centered difference formula** has $p=q$.
::::

The simplest example of a forward difference formula is inspired by the familiar limit definition of a derivative:

```{math}
  :label: forwardFD11
  f'(x) \approx \frac{f(x+h)-f(x)}{h},
```

which is {eq}`fdformula` with $p=0$, $q=1$, $a_0=-1$, and $a_1=1$. Analogously, we have the backward difference

```{math}
  :label: backwardFD11
  f'(x) \approx \frac{f(x)-f(x-h)}{h},
```

in which $p=1$, $q=0$.

(example-finitediffs-fd1bd1)=
::::{proof:example}
Suppose $f(x)=x^2$, and we take $h=\frac{1}{4}$ over the interval $[0,1]$. This results in the nodes $0,\frac{1}{4},\frac{1}{2},\frac{3}{4},1$. We evaluate $f$ at the nodes to get

$$
f(0) = 0, \; f\left(\tfrac{1}{4}\right) = \frac{1}{16},\; f\left(\tfrac{1}{2}\right)=\frac{1}{4},\; f\left(\tfrac{3}{4}\right)=\frac{9}{16}, \; f(1)=1.
$$

This gives four forward difference estimates,

$$
f'(0) & \approx 4\left(\frac{1}{16}-0\right), &\quad 
f'\left(\tfrac{1}{4}\right)& \approx 4\left(\frac{1}{4}-\frac{1}{16}\right), \\
f'\left(\tfrac{1}{2}\right)& \approx 4\left(\frac{9}{16}-\frac{1}{4}\right), &\quad 
f'\left(\tfrac{3}{4}\right) &\approx 4\left(1-\frac{9}{16}\right).
$$

We also get four backward difference estimates,

$$
f'\left(\tfrac{1}{4}\right) &\approx 4\left(\frac{1}{16}-0\right), &\quad 
f'\left(\tfrac{1}{2}\right) &\approx 4\left(\frac{1}{4}-\frac{1}{16}\right), \\ 
f'\left(\tfrac{3}{4}\right) &\approx 4\left(\frac{9}{16}-\frac{1}{4}\right), &\quad 
f'\left(1\right) &\approx 4\left(1-\frac{9}{16}\right).
$$

Notice that it's the same four differences each time, but we're interpreting them as derivative estimates at different nodes.
::::

As pointed out in {numref}`Example {number} <example-finitediffs-fd1bd1>`, the only real distinction between {eq}`forwardFD11` and {eq}`backwardFD11` is whether we think that $f'$ is being evaluated at the left node or the right one. Symmetry would suggest that we should evaluate it halfway between. That is the motivation behind centered difference formulas.


Let's derive the shortest centered formula using $p=q=1$. For simplicity, we will set $x=0$ without affecting the result. This means that $f(-h)$, $f(0)$, and $f(h)$ are all available in {eq}`fdformula`. 

Note that {eq}`forwardFD11` is simply the slope of the line through the points $\bigl(0,f(0)\bigr)$ and $\bigl(h,f(h)\bigr)$. One route to using all three function values is to differentiate the quadratic polynomial that interpolates $\bigl(-h,f(-h)\bigr)$ as well (see [Exercise 1](problem-quadraticFD)):
 
```{math}
:label: fdinterp2
Q(x) = \frac{x(x-h)}{2h^2} f(-h) - \frac{x^2-h^2}{h^2} f(0) + \frac{x(x+h)}{2h^2} f(h).
```

This leads to

```{math}
:label: centerFD12
f'(0) \approx Q'(0) = \frac{f(h)-f(-h)}{2h}.
```

This result is equivalent to {eq}`fdformula` with $p=q=1$ and weights $a_{-1}=-\frac{1}{2}$, $a_0=0$, and $a_1=\frac{1}{2}$. Observe that while the value of $f(0)$ was available during the derivation, its weight ends up being zero.

Besides the aesthetic appeal of symmetry, in {numref}`section-localapprox-fd-converge` we will see another important advantage of {eq}`centerFD12` compared to the one-sided formulas. 

```{index} interpolation; by polynomials
```

We can in principle derive any finite-difference formula from the same process: Interpolate the given function values, then differentiate the interpolant exactly. Some results of the process are given in {numref}`table-FDcenter` for centered differences, and in {numref}`table-FDforward` for forward differences. Both show the weights for estimating the derivative at $x=0$. To get backward differences, you change the signs and reverse the order of the coefficients in any row of {numref}`table-FDforward`; see [Exercise 2](problem-backwardFD).

```{tabularcolumns} |r|ccccccccc
```

(table-FDcenter)=
```{table} Weights for centered finite-difference formulas.

| order            | $-4h$           | $-3h$            | $-2h$          | $-h$           | $0$ | $h$           | $2h$            | $3h$            | $4h$             |
|:----------------:|:---------------:|:----------------:|:--------------:|:--------------:|:---:|:-------------:|:---------------:|:---------------:|:----------------:|
| 2                |                 |                  |                | $-\frac{1}{2}$ | $0$ | $\frac{1}{2}$ |                 |                 |                  |
| 4                |                 |                  | $\frac{1}{12}$ | $-\frac{2}{3}$ | $0$ | $\frac{2}{3}$ | $-\frac{1}{12}$ |                 |                  |
| 6                |                 | $-\frac{1}{60}$  | $\frac{3}{20}$ | $-\frac{3}{4}$ | $0$ | $\frac{3}{4}$ | $-\frac{3}{20}$ | $\frac{1}{60}$  |                  |
| 8                | $\frac{1}{280}$ | $-\frac{4}{105}$ | $\frac{1}{5}$  | $-\frac{4}{5}$ | $0$ | $\frac{4}{5}$ | $-\frac{1}{5}$  | $\frac{4}{105}$ | $-\frac{1}{280}$ |
```

```{tabularcolumns} |r|ccccc
```

(table-FDforward)=
```{table} Weights for forward finite-difference formulas. To get backward differences, change the signs and reverse the order of the coefficients.

| order         | $0$              | $h$ | $2h$           | $3h$          | $4h$           |
|:-------------:|:----------------:|:---:|:--------------:|:-------------:|:--------------:|
|  1            | $-1$             | $1$ |                |               |                |
|  2            | $-\frac{3}{2}$   | 2   | $-\frac{1}{2}$ |               |                |
|  3            | $-\frac{11}{6}$  | 3   | $-\frac{3}{2}$ | $\frac{1}{3}$ |                |
|  4            | $-\frac{25}{12}$ | $4$ | $-3$           | $\frac{4}{3}$ | $-\frac{1}{4}$ |
```

 The main motivation for using more function values in a formula is to improve the accuracy. This is measured by **order of accuracy**, which is shown in the tables and explored in {numref}`section-localapprox-fd-converge`. 

````{proof:example}
According to the tables, here are three specific finite-difference formulas:

```{math}
\begin{split}
f'(0) &\approx \tfrac{1}{h} \left[ \tfrac{1}{12} f(-2h)
- \tfrac{2}{3} f(-h) + \tfrac{2}{3} f(h) - \tfrac{1}{12} f(2h) \right], \\[1mm]
f'(0) &\approx \tfrac{1}{h} \left[ -\tfrac{3}{2} f(0) + 2 f(h) -\tfrac{1}{2} f(2h) \right], \\[1mm]
f'(0) &\approx \tfrac{1}{h} \left[ \tfrac{1}{2} f(-2h) - 2 f(-h) + \tfrac{3}{2} f(0) \right].
\end{split}
```
````

```{proof:demo}
```

```{raw} html
<div class='demo'>
```

```{raw} latex
%%start demo%%
```

If $f(x)=e^{\,\sin(x)}$, then $f'(0)=1$. 

```{code-cell}
f = x -> exp(sin(x));
```

Here are the first two centered differences from {numref}`table-FDcenter`.

```{code-cell}
h = 0.05
CD2 = (       - f(-h)  + f(h)         ) /  2h
CD4 = (f(-2h) - 8f(-h) + 8f(h) - f(2h)) / 12h
@show (CD2,CD4);
```

Here are the first two forward differences from {numref}`table-FDforward`.

```{code-cell}
FD1 = ( -f(0) +  f(h)        ) /  h
FD2 = (-3f(0) + 4f(h) - f(2h)) / 2h 
@show (FD1,FD2);
```

Finally, here are the backward differences that come from reverse-negating the forward differences. 

```{code-cell}
BD1 = (          -f(-h) +  f(0)) /  h
BD2 = ( f(-2h) - 4f(-h) + 3f(0)) / 2h 
@show (BD1,BD2);
```
```{raw} html
</div>
```

```{raw} latex
%%end demo%%
```

## Higher derivatives

Many applications require the second derivative of a function. It's tempting to use the finite difference of a finite difference. For example, applying {eq}`centerFD12` to $f'$ gives 

```{math}
f''(0) \approx  \frac{ f'(h) - f'(h) }{2h}.
```

Then applying {eq}`centerFD12` to approximate the appearances of $f'$ leads to
```{math}
:label: fd-twocenter
f''(0) \approx  \frac{ f(-2h) - 2 f(0) + f(2h) }{4h^2}.
```

This is a valid formula, but it uses values at $\pm 2h$ rather than the closer values at $\pm h$. A better and more generalizable tactic is to return to the quadratic $Q(x)$ in {eq}`fdinterp2` and use $Q''(0)$ to approximate $f''(0)$. Doing so yields

```{math}
  :label: centerFD22
  f''(0) \approx  \frac{ f(-h) - 2 f(0) + f(h) }{h^2},
```

which is the simplest **centered second-difference formula**. As with the first derivative, we can choose larger values of $p$ and $q$ in {eq}`fdformula` to get new formulas, such as

```{math}
:label: forwardFD21
f''(0) \approx \frac{ f(0) - 2 f(h) + f(2h) }{h^2},
```

and

```{math}
:label: forwardFD22
f''(0) \approx \frac{ 2f(0) - 5 f(h) + 4 f(2h) -f(3h) }{h^2}.
```

For the second derivative, converting a forward difference to a backward difference requires reversing the order of the weights, while *not* changing their signs.

```{proof:demo}
```

```{raw} html
<div class='demo'>
```

```{raw} latex
%%start demo%%
```

If $f(x)=e^{\,\sin(x)}$, then $f''(0)=1$. 

```{code-cell}
f = x -> exp(sin(x));
```

Here is a centered estimate given by {eq}`centerFD22`.

```{code-cell}
h = 0.05
CD2 = (f(-h) - 2f(0) + f(h)) / h^2
@show CD2;
```

For the same $h$, here are forward estimates given by {eq}`forwardFD21` and {eq}`forwardFD22`.

```{code-cell}
FD1 = ( f(0) - 2f(h) +  f(2h)        ) / h^2
FD2 = (2f(0) - 5f(h) + 4f(2h) - f(3h)) / h^2
@show (FD1,FD2);
```

Finally, here are the backward estimates that come from reversing {eq}`forwardFD21` and {eq}`forwardFD22`.

```{code-cell}
BD1 = (           f(-2h) - 2f(-h) +  f(0)) / h^2
BD2 = (-f(-3h) + 4f(-2h) - 5f(-h) + 2f(0)) / h^2
@show (BD1,BD2);
```
```{raw} html
</div>
```

```{raw} latex
%%end demo%%
```

## Arbitrary nodes

Although function values at equally spaced nodes are a common and convenient situation, the node locations may be arbitrary. The general form of a finite-difference formula is

```{math}
:label: fdgenformula
  f^{(m)}(0) \approx \sum_{k=0}^{r} c_{k,m} \,f(t_k).
```

We no longer assume equally spaced nodes, so there is no "$h$" to be used in the formula. As before, the weights may be applied after any translation of the independent variable. The weights again follow from the interpolate/differentiate recipe, but the algebra becomes complicated. Fortunately there is an elegant recursion known as **Fornberg's algorithm** that can calculate these weights for any desired formula. We present it without derivation as {numref}`Function {number} <function-fdweights>`.

(function-fdweights)=

````{proof:function} fdweights
**Fornberg's algorithm for finite-difference weights**

```{code-block} julia
:lineno-start: 1
"""
    fdweights(t,m)

Compute weights for the `m`th derivative of a function at zero using
values at the nodes in vector `t`.
"""
function fdweights(t,m)
# This is a compact implementation, not an efficient one.
    # Recursion for one weight. 
    function weight(t,m,r,k)
        # Inputs
        #   t: vector of nodes 
        #   m: order of derivative sought 
        #   r: number of nodes to use from t 
        #   k: index of node whose weight is found

        if (m<0) || (m>r)        # undefined coeffs must be zero
            c = 0
        elseif (m==0) && (r==0)  # base case of one-point interpolation
            c = 1
        else                     # generic recursion
            if k<r
                c = (t[r+1]*weight(t,m,r-1,k) -
                    m*weight(t,m-1,r-1,k))/(t[r+1]-t[k+1])
            else
                numer = r > 1 ? prod(t[r]-x for x in t[1:r-1]) : 1
                denom = r > 0 ? prod(t[r+1]-x for x in t[1:r]) : 1
                β = numer/denom
                c = β*(m*weight(t,m-1,r-1,r-1) - t[r]*weight(t,m,r-1,r-1))
            end
        end
        return c
    end
    r = length(t)-1
    w = zeros(size(t))
    return [ weight(t,m,r,k) for k=0:r ]
end
```
````

(demo-finitediffs-fd-weights)=

```{proof:demo}
```

```{raw} html
<div class='demo'>
```

```{raw} latex
%%start demo%%
```

We will estimate the derivative of $\cos(x^2)$ at $x=0.5$ using five nodes.

```{code-cell}
t = [ 0.35,0.5,0.57,0.6,0.75 ]   # nodes
f = x -> cos(x^2)
dfdx = x -> -2*x*sin(x^2)
exact_value = dfdx(0.5)
```

We have to shift the nodes so that the point of estimation for the derivative is at $x=0$. (To subtract a scalar from a vector, we must use the `.-` operator.)

```{code-cell}
w = FNC.fdweights(t.-0.5,1)
```

The finite-difference formula is a dot product (i.e., inner product) between the vector of weights and the vector of function values at the nodes.

```{code-cell}
fd_value = dot(w,f.(t))
```

We can reproduce the weights in the finite-difference tables by using equally spaced nodes with $h=1$. For example, here is a one-sided formula at four nodes.

```{code-cell}
FNC.fdweights(0:3,1)
```

```{index} ! Julia; Rational
```

By giving nodes of type `Rational`, we can get exact values instead.

```{code-cell}
FNC.fdweights(Rational.(0:3),1)
```
```{raw} html
</div>
```

```{raw} latex
%%end demo%%
```

## Exercises

(problem-quadraticFD)=

1. ✍ This problem refers to $Q(x)$ defined by {eq}`fdinterp2`.  

    **(a)** Show that $Q(x)$ interpolates the three values of $f$ at $x=-h$, $x=0$, and $x=h$.

    **(b)** Show that $Q'(0)$ gives the finite-difference formula defined by {eq}`centerFD12`.

    ````{only} solutions
    ````

    (problem-backwardFD)=
2. **(a)** ✍ {numref}`table-FDforward` lists forward difference formulas in which $p=0$ in {eq}`fdformula`. Show that the change of variable $g(x) = f(-x)$ transforms these formulas into backward difference formulas with $q=0$, and write out the table analogous to {numref}`table-FDforward` for backward differences.

    **(b)** ⌨ Suppose you are given the nodes $t_0=0.9$, $t_1=1$, and $t_2=1.1$, and $f(x) = \sin(2x)$. Using formulas from {numref}`table-FDcenter` and {numref}`table-FDforward`, compute second-order accurate approximations to $f'$ at each of the three nodes.

    ````{only} solutions
    ````

3. ⌨ Let $f(x)=e^{-x}$, $x=0.5$, and $h=0.2$. Using {numref}`Function {number} <function-fdweights>` to get the necessary weights on five nodes centered at $x$, find finite-difference approximations to the first, second, third, and fourth derivatives of $f$. Make a table showing the derivative values and the errors in each case.

4. ⌨ In the manner of {numref}`Demo {number} <demo-finitediffs-fd-weights>`, use {numref}`Function {number} <function-fdweights>` on centered node vectors of length  3, 5, 7, and 9 to produce a table analogous to {numref}`table-FDcenter` for the second derivative $f''(0)$. (You do not need to show the orders of accuracy, just the weights.)

5. ⌨ For this problem, let $f(x)=\tan(2x)$.

    **(a)** ⌨ Apply {numref}`Function {number} <function-fdweights>` to find a finite-difference approximation to $f''(0.3)$ using the five nodes $t_j=0.3+jh$ for $j=-2,\ldots,2$ and $h=0.05$. Compare to the exact value of $f''(0.3)$.

    **(b)** ⌨  Repeat part (a) for $f''(0.75)$ on the nodes $t_j=0.75+jh$. Why is the finite-difference result so inaccurate? (Hint: A plot of $f$ might be informative.)

6. ✍ Find the finite-difference formula for $f''(0)$ that results from applying {eq}`forwardFD11` on $f'$ and then {eq}`backwardFD11` on $f'$ within that result.  

7. **(a)** ✍ Show using L'Hôpital's Rule that the centered formula approximation {eq}`centerFD12` converges to an equality as $h\to 0$.

    **(b)** ✍ Derive two conditions on the finite-difference weights in {eq}`fdformula` that arise from requiring convergence as $h\to 0$. (Hint: Consider what is required in order to apply L'Hôpital's Rule, as well as the result of applying it.)
