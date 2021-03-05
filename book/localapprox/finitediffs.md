# Finite differences

Much more can be said about interpolation. But now we turn to one of the most common and important applications of interpolants: finding derivatives. For the moment, we will continue to use $x$ as the independent variable name and $a=t_0,\ldots,t_n=b$ as the interpolation nodes. We also continue to use an equispaced grid, so that $t_i=a+i h$ for $i=0,\ldots,n$, where $h=(b-a)/n$.

Considering the most common definition of a derivative,

```{math}
f'(x)=\lim_{h\to0} \frac{f(x+h)-f(x)}{h},
```

```{index} finite differences
```

it stands to reason that an approximation to $f'$ at a node $t_i$ should depend on the values of $f$ at the nodes closest to $t_i$. Also, because differentiation is a linear operation, we will constrain ourselves to formulas that are linear in these nodal values. Consequently, the sort of formula we seek is the {term}`finite difference` formula

```{math}
  :label: fdformula
  f'(t_i) \approx \frac{1}{h} \sum_{k=-p}^{q} a_k f(t_i+kh),
```

```{margin}
Finite difference weights are independent of the function being differentiated.
```

where $p$, $q$ are integers, and the $a_k$'s are constants known as the **weights** of the formula. Crucially, the finite difference weights are independent of $f$, although they do depend on the nodes. The factor of $h^{-1}$ is present to make the expression more convenient in what follows.

Before deriving some finite difference formulas, we make an important observation about them. Define the new variable $s=x-t_i$ and let $\tilde{f}(s) = f(x-t_i)$. Then it is elementary that

```{math}
  \left. \frac{d f}{d x} \right|_{x=t_i} = \left. \dd{\tilde{f}}{s} \right|_{s=0}.
```

Applying the change of variables to {eq}`fdformula` yields

```{math}
   f'(t_i) = \tilde{f}'(0) \approx \frac{1}{h} \sum_{k=i-p}^{i+q} a_k \tilde{f}(kh),
```

with the same constants as before. These manipulations express a property that is simpler than it may appear: we can always derive and write the finite difference formula {eq}`fdformula` with $t_i=0$ without losing generality. In fact, $t_i$ is just a "dummy" variable that we can replace by $x$, as in

```{math}
:label: fdxformula
f'(x) \approx \frac{1}{h} \sum_{k=-p}^{q} a_k f(x+kh).
```

This property is *translation invariance*. The formula combines values of the function at points always placed the same way relative to $x$.

An obvious candidate for a finite difference formula is based on the limit definition above:

```{math}
  :label: forwardFD11
  f'(x) \approx \frac{f(x+h)-f(x)}{h},
```

which is {eq}`fdxformula` with $p=0$, $q=1$, $a_0=-1$, and $a_1=1$. This is referred to as a **forward difference formula**, characterized by $p=0$, because $f$ is evaluated only at points "forward" from $x$. Analogously, we could use the **backward difference formula**

```{math}
  :label: backwardFD11
  f'(x) \approx \frac{f(x)-f(x-h)}{h},
```

in which $q=0$.

Both the forward and backward difference formulas surely become equalities in the limit $h\to 0$, provided $f$ is differentiable at $x$. However, they are not the only possibilities. One aesthetic objection is the lack of symmetry about the point $x$. In response, we will derive a formula that uses $p=q=1$, i.e., in which $f(-h)$, $f(0)$, and $f(h)$ are all available.

```{index} interpolation; by polynomials
```

The formula {eq}`forwardFD11` (with $x=0$) is simply the slope of the line through the points $\bigl(0,f(0)\bigr)$ and $\bigl(h,f(h)\bigr)$. A similar observation holds for the backward difference formula. Thus one route to using three function values is to differentiate the quadratic polynomial that interpolates them (see [this exercise](problem-quadraticFD)):

```{math}
:label: fdinterp2
Q(x) = \frac{x(x-h)}{2h^2} f(-h) - \frac{x^2-h^2}{h^2} f(0) + \frac{x(x+h)}{2h^2} f(h).
```

We now use $f'(0)\approx Q'(0)$ to get the **centered finite difference formula**

```{math}
:label: centerFD12
f'(0) \approx \frac{f(h)-f(-h)}{2h}.
```

This result is equivalent to {eq}`fdxformula` with $p=q=1$ and weights $a_{-1}=-1/2$, $a_0=0$, and $a_1=1/2$. Observe that while the value of $f(0)$ was available during the derivation, its weight ends up being zero.

We can verify using L'Hôpital's Rule that the approximation in {eq}`centerFD12` becomes an equality as $h\to 0$. Such an analysis does not, however, reveal a significant accuracy advantage of the centered variant, one that we will take up in {doc}`fd-converge`.

```{index} interpolation; by polynomials
```

```{margin}
Finite difference formulas are derived by interpolating function values, followed by differentiation of the interpolant.
```

We can in principle derive any finite difference formula from the same process: Interpolate the given function values, then differentiate the interpolant exactly. Some results are given here for two important special cases. {numref}`table-FDcenter` is for $p=q$, or centered differences, while {numref}`table-FDforward` is for $p=0$, or forward differences. Both show the weights for estimating the derivative at zero, but they can be uniformly translated to any other point.

```{tabularcolumns} |r|ccccccccc
```

(table-FDcenter)=

```{table} Weights for centered finite difference formulas.

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

```{table} Weights for forward finite difference formulas.

| order         | $0$              | $h$ | $2h$           | $3h$          | $4h$           |
|:-------------:|:----------------:|:---:|:--------------:|:-------------:|:--------------:|
|  1            | $-1$             | $1$ |                |               |                |
|  2            | $-\frac{3}{2}$   | 2   | $-\frac{1}{2}$ |               |                |
|  3            | $-\frac{11}{6}$  | 3   | $-\frac{3}{2}$ | $\frac{1}{3}$ |                |
|  4            | $-\frac{25}{12}$ | $4$ | $-3$           | $\frac{4}{3}$ | $-\frac{1}{4}$ |
```

 The main motivation for using more function values in a formula is to improve the accuracy. This is measured by **order of accuracy**, which is show in the tables and explained in {doc}`fd-converge`. To get backward differences with $q=0$, you can use the change of variable $\hat{f}(x)=f(-x)$, which changes the sign and reverses the order of the coefficients in {numref}`table-FDforward`; see [this exercise](problem-backwardFD).

````{proof:example}
According to the tables, here are two finite difference formulas:

```{math}
\begin{split}
f'(0) &\approx h^{-1} \left[ \tfrac{1}{12} f(-2h)
- \tfrac{2}{3} f(-h) + \tfrac{2}{3} f(h) - \tfrac{1}{12} f(2h) \right], \\
f'(0) &\approx h^{-1} \left[ \tfrac{1}{2} f(-2h) - 2 f(-h) + \tfrac{3}{2} f(0) \right]. \\	
\end{split}
```
````

## Higher derivatives

Many applications require the second derivative of a function. It's tempting to use the finite difference of a finite difference. For example, applying {eq}`centerFD12` twice leads to

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

## Arbitrary nodes

Although function values at equally spaced nodes are a common and convenient situation, the node locations may be arbitrary. The general form of a finite difference formula is

```{math}
:label: fdgenformula
  f^{(m)}(0) \approx \sum_{k=0}^{r} c_{k,m} \,f(t_k).
```

````{proof:example} Julia demo
:class: demo
{doc}`demos/fd-weights`
````

We no longer assume equally spaced nodes, so there is no "$h$" to be used in the formula. As before, the weights may be applied after any translation of the independent variable. The weights again follow from the interpolate/differentiate recipe, but the algebra becomes complicated. Fortunately there is an elegant recursion known as **Fornberg's algorithm** that can calculate these weights for any desired formula. We present it without derivation as {ref}`function-fdweights`.

(function-fdweights)=

````{proof:function} fdweights
**Fornberg's algorithm for finite difference weights.**

```{code-block} julia
:lineno-start: 1
"""
fdweights(t,m)

Return weights for the `m`th derivative of a function at zero using
values at the nodes in vector `t`.
"""
function fdweights(t,m)
    # This is a compact implementation, not an efficient one.

    function weight(t,m,r,k)
        # Recursion for one weight. Input: t   nodes (vector) m
        # order of derivative sought r   number of nodes to use from
        # t (<= length(t)) k   index of node whose weight is found

        if (m<0) || (m>r)        # undefined coeffs must be zero
            c = 0
        elseif (m==0) && (r==0)  # base case of one-point interpolation
            c = 1
        else                     # generic recursion
            if k<r
                c = (t[r+1]*weight(t,m,r-1,k) -
                    m*weight(t,m-1,r-1,k))/(t[r+1]-t[k+1])
            else
                numer = r > 1 ? prod(t[r]-x for x in t[1:r-1]) : 1.0
                denom = r > 0 ? prod(t[r+1]-x for x in t[1:r]) : 1.0
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

## Exercises

(problem-quadraticFD)=

1. ✍ This problem refers to $Q(x)$ defined by {eq}`fdinterp2`.  

    **(a)** Show that $Q(x)$ interpolates the three values of $f$ at $x=-h$, $x=0$, and $x=h$.

    **(b)** Show that $Q'(0)$ gives the finite difference formula defined by {eq}`centerFD12`.

    ````{only} solutions
    ````

    (problem-backwardFD)=
2. **(a)** ✍ {ref}`table-FDforward` lists forward-difference formulas in which $p=0$ in {eq}`fdxformula`. Show that the change of variable $g(x) = f(-x)$ transforms these formulas into backward difference formulas with $q=0$, and write out the table analogous to {ref}`table-FDforward` for backward differences.

    **(b)** ⌨ Suppose you are given the nodes $t_0=0.9$, $t_1=1$, and $t_2=1.1$, and $f(x) = \sin(2x)$. Using formulas from {numref}`table-FDcenter` and {numref}`table-FDforward`, compute second-order accurate approximations to $f'$ at each of the three nodes.

    ````{only} solutions
    ````

3. ⌨ Using {ref}`function-fdweights` to get the necessary weights, find finite-difference approximations to the first, second, third, and fourth derivatives of $f(x)=e^{-x}$ at $x=0.5$. In each case use a centered stencil of minimum possible width. Make a table showing the values and the errors in each case.

    ````{only} solutions
    ````

4. ⌨ Use {ref}`function-fdweights` to write out a table analogous to {ref}`table-FDcenter` that lists centered finite difference weights for the second derivative $f''(0)$. (Hint: The "rat" command will let you express the results as exact rational numbers.)

    ````{only} solutions
    ````

5. ⌨ For this problem, let $f(x)=\tan(2x)$.

    **(a)** ⌨ Apply {ref}`function-fdweights` to find a finite difference approximation to $f''(0.3)$ using the five nodes $t_j=0.3+jh$ for $j=-2,\ldots,2$ and $h=0.05$. Compare to the exact value of $f''(0.3)$.

    **(b)** ⌨  Repeat part~(a) for $f''(0.75)$ and $t_j=0.75+jh$. Why is the finite difference result so inaccurate?

    ````{only} solutions
    ````

6. **(a)** ✍ Derive {eq}`fd-twocenter` by applying applying {eq}`centerFD12` twice (i.e., apply once to get $f'$ values and then apply again to those values).

    **(b)** ✍ Find the formula for $f''(0)$ that results from applying {eq}`forwardFD11` and then {eq}`backwardFD11`.  

    ````{only} solutions
    ````

7. **(a)** ✍ Show using L'Hôpital's Rule that the centered formula approximation {eq}`centerFD12` converges to an equality as $h\to 0$.

    **(b)** ✍ Derive two conditions on the finite difference weights in {eq}`fdxformula` that arise from requiring convergence as $h\to 0$.
