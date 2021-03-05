# Interpolation-based methods

From a practical standpoint, one of the biggest drawbacks of Newton's method is the requirement to supply $f'$ in [`newton`](function-newton). It is both a programming inconvenience and a step that requires computational time. We can avoid using $f'$ however, by making a simple but easily overlooked observation:

```{proof:observation}
 When a step produces an approximate result, you are free to carry it out approximately. 
```

Let's call this the "principle of approximate approximation."

In the Newton context, the principle of approximate approximation begins with the observation that the use of $f'$ is linked to the construction of a linear approximation $q(x)$ equivalent to a tangent line. The root of $q(x)$ is used to define the next iterate in the sequence. We can avoid calculating the value of $f'$ by choosing a different linear approximation.

```{index} secant method
```

````{prf:example} Julia demo
:class: demo
{doc}`demos/secant-line`
````

The example in {doc}`demos/secant-line` demonstrates the {term}`secant method`. In the secant method, one finds the root of the linear approximation through the two most recent root estimates. That is, given previous approximations $x_1,\ldots,x_k$, define the linear model function as the line through $\bigl(x_{k-1},f(x_{k-1})\bigr)$ and $\bigl(x_k,f(x_k)\bigr)$:

```{math}
:label: secantmodel
q(x) = f(x_k) + \frac{f(x_k)-f(x_{k-1})}{x_k-x_{k-1}}(x-x_k).
```

```{margin}
In the secant method, one finds the root of the line through the two most recent root estimates.
```

Solving $q(x_{k+1})=0$ for $x_{k+1}$ gives the formula

```{math}
:label: secant
x_{k+1} = x_k - \frac{f(x_k)(x_k-x_{k-1})}{f(x_k)-f(x_{k-1})}, \quad n=1,2,\ldots.
```

Our implementation of the method based on this formula is given in {ref}`function-secant`.

(function-secant)=

```{proof:function} secant
**Secant method for scalar rootfinding.**

```{code-block} julia
:lineno-start: 1
"""
secant(f,x1,x2)

Use the secant method to find a root of `f` starting from `x1` and
`x2`. Returns a vector of root estimates.
"""
function secant(f,x1,x2)
    # Operating parameters.
    funtol = 100*eps();  xtol = 100*eps();  maxiter = 40;

    x = [x1,x2]
    y1 = f(x1); y2 = 100;
    dx = Inf   # for initial pass below
    k = 2

    while (abs(dx) > xtol) && (abs(y2) > funtol) && (k < maxiter)
        y2 = f(x[k])
        dx = -y2 * (x[k]-x[k-1]) / (y2-y1)   # secant step
        push!(x,x[k]+dx)        # append new estimate

        k = k+1
        y1 = y2    # current f-value becomes the old one next time
    end

    if k==maxiter
        @warn "Maximum number of iterations reached."
    end

    return x
end
```

## Convergence

The secant method uses a different linear model in each iteration than the one Newton's method would use. Is it as good? As before, let $\epsilon_k = r-x_k$ be the errors in the successive root approximations. If the initial errors are small, then a tedious but straightforward Taylor expansion shows that, to lowest order,

```{math}
:label: secanterr
\epsilon_{k+1} \approx -\frac{1}{2}\frac{f''(r)}{f'(r)} \epsilon_k \epsilon_{k-1}.
```

If we make an educated guess that

```{math}
\epsilon_{k+1} = c (\epsilon_k)^\alpha, \qquad \epsilon_k = c (\epsilon_{k-1})^\alpha, \qquad \alpha>0,
```

then {eq}`secanterr` becomes

```{math}
:label: secantexponents
\left[ \epsilon_{k-1}^{\alpha} \right]^{\alpha} \approx C \epsilon_{k-1}^{\alpha+1},
```

for an unknown constant $C$. Treating the implied approximation as an equality, this becomes solvable if and only if the exponents match, i.e., $\alpha^2 = \alpha+1$. The only positive root of this equation is the golden ratio,

```{math}
  \alpha = \frac{1+\sqrt{5}}{2} \approx 1.618.
```

```{index} superlinear convergence
```

```{index} convergence rate; superlinear
```

Hence the errors in the secant method converge like $\epsilon_{k+1} = c (\epsilon_k)^\alpha$  for $1<\alpha<2$, a situation called {term}`superlinear convergence`.

````{prf:example} Julia demo
:class: demo
{doc}`demos/secant-converge`
````

In terms of error as a function of the iteration number $k$, the secant method converges at a rate between linear and quadratic, which is slower than Newton's method. In that sense we must conclude that the secant line is inferior to the tangent line for approximating $f$ near a root. But error versus iteration count may not be the best means of comparison.

Often we analyze rootfinding methods by assuming that the bulk of computing time is spent evaluating the user-defined functions $f$ and $f'$. (Our simple examples and exercises mostly don't support this assumption, but many practical applications do.) In this light we see that Newton's method requires two evaluations, $f(x_k)$ and $f'(x_k)$, for each iteration. The secant method, on the other hand, while it *uses* the two function values $f(x_k)$ and $f(x_{k-1})$ at each iteration, only needs to *compute* a single new one. Note that {ref}`function-secant` keeps track of one previous function value rather than recomputing it.

Now suppose that $|\epsilon_k|=\epsilon$. Roughly speaking, two units of work (i.e., function evaluations) in Newton's method brings us to an error of $\epsilon^2$. If one spreads out the improvement in the error evenly across the two geometric steps, using

```{math}
\epsilon^2 = \bigl( \epsilon^{\sqrt{2}} \bigr)^{\!\sqrt{2}},
```

```{margin}
The secant method is both easier to apply and faster than Newton's method, on an equal-work basis.
```

it seems reasonable to say that the rate of convergence *per function evaluation* is $\sqrt{2}\approx 1.41$. This is actually less than the comparable rate of about $1.62$ for the secant method. Not only is the secant method easier to apply than Newton's method in practice, it's also more efficient when counting function evaluations—a rare double victory!

## Inverse interpolation

```{index} interpolation; by polynomials
```

At each iteration, the secant method constructs a linear model function that interpolates the two most recently found points on the graph of $f$. Two points determine a straight line, so this seems like a sensible choice. But as the iteration progresses, why use only the *two* most recent points? What would it mean to use more of them?

If we interpolate through three points by a polynomial, we get a unique quadratic function. Unfortunately, a parabola may have zero, one, or two crossings of the $x$-axis, leaving some doubt as to how to define the next root estimate. On the other hand, if we turn a parabola on its side, we get a graph that intersects the $x$-axis exactly once, which is ideal for defining the next root estimate.

This leads to the idea of defining $q(y)$ as the quadratic interpolant to the points $(y_{k-2},x_{k-2})$, $(y_{k-1},x_{k-1})$, and $(y_k,x_k)$, where $y_i=f(x_i)$ for all $i$, and setting $x_{k+1}=q(0)$. The process defined in this way (given three initial estimates) is called **inverse quadratic interpolation**. Rather than deriving lengthy formulas for it here, we demonstrate how to perform inverse quadratic interpolation using `Polynomials.fit` to perform the interpolation step.

````{prf:example} Julia demo
:class: demo
{doc}`demos/secant-iqi`
````

## Bracketing

Like Newton's method, the secant and inverse quadratic interpolation methods cannot guarantee convergence. One final new idea is needed to make a foolproof algorithm.

If $f$ is continuous on the interval $[a,b]$ and $f(a)f(b)<0$—that is, $f$ changes sign on the interval—then $f$ must have (at least) one root in the interval, due to the Intermediate Value Theorem from calculus. If we come up with a new root estimate $c\in(a,b)$, then whatever sign $f(c)$ is, it is different from the sign at one of the endpoints. (Of course, if $f(c)$ is zero, we are done!) So either $[a,c]$ or $[c,b]$ is guaranteed to have a root too, and in this way we can maintain not just individual estimates but an interval that always contains a root.

The best algorithms blend the use of fast-converging methods with the guarantee provided by a bracket. For example, say that an iteration begins with a bracketing interval. Make a list of the inverse quadratic estimate, the secant estimate, and the midpoint of the current interval and pick the first member of the list that lies within the current interval. Replace the interval with the bracketing subinterval, and start a new iteration. This is the idea behind **Brent's method**, which is a very successful rootfinding algorithm.

## Exercises

For each of problems 1--3, do the following steps.
  
**(a)** ✍ Rewrite the equation into the standard form for rootfinding, $f(x) = 0$. **(b)** ⌨ Make a plot of $f$ over the given interval and determine how many roots lie in the interval. **(c)** ⌨ Use `nlsolve` to find an "exact" value for each root. **(d)** ⌨ Determine a bracketing interval for each root. Then use {ref}`function-secant`, starting with the endpoints of the interval, to find each root. **(e)** ⌨ For one of the roots, define `e` as a vector of the errors in the secant sequence. Determine numerically whether the convergence is superlinear.

1. $x^2=e^{-x}$, over $[-2,2]$

    ````{only} solutions

    ``` matlab
    %% Problem 4.3.1 (a)
    % The following function will work (or its opposite):
    %
    % $$ f(x) = x^2 - e^{-x} = 0 $$
    %
    % There is only one root in $$ x \in [0,1] $$.

    f1 = @(x) x.^2 - exp(-x);
    x = linspace(-2,2);
    plot(x,f1(x),'-b',[-2 2],[0 0],'--g','LineWidth',2)
    xlabel('x'), ylabel('f(x)')

    %% Problem 4.3.1(b)
    %
    % Use fzero to get an "exact" value of the root
    ropts = optimset('TolX',1e-14,'TolFun',1e-14);
    rexact = fzero(f1,0.67, ropts)

    %% Problem 4.3.1(c)
    %
    rsec = secant(f1,0,1)

    %% Problem 4.3.1(d)
    %
    % the convergence is superlinear, since the green dashed line is linear
    % convergence, the iterates (blue) converge to zero at a faster rate.

    e = abs(rexact-rsec)
    loglog(e(1:end-1),e(2:end),'-b',e(1:end-1),e(1:end-1),'--g','LineWidth',2)
    xlabel('e_n'), ylabel('e_{n+1}')

    ```
    ````

2. $2x = \tan x$, over $[-0.2,1.4]$

3. $e^{x+1}=2+x$, over $[-2,2]$

4. ⌨ Use a plot to approximately locate all the roots of $f(x)=x^{-2}-\sin(x)$ in the interval $[0.5,4\pi]$. Then find a pair of initial points for each root such that {ref}`function-secant` converges to that root.

5. ✍ Show analytically that the secant method converges in one step for a linear function, regardless of the initialization.

6. ✍ In general, the secant method formula {eq}`secant` cannot be applied if $x_{k}=x_{k-1}$. However, suppose that $f(x)=ax^2+bx+c$ for constants $a$, $b$, and $c$. Show that in this case the formula can be simplified to one that is well defined when $x_{k}=x_{k-1}$. Then show that the resulting $x_{k+1}$ is the same as the result of one step of Newton's method applied to $f$ at $x_k$.

7. ✍ Let $f(x)=x^2$. Show that if $(1/x_1)$ and $(1/x_2)$ are positive integers, and the secant iteration is applied, then the sequence $1/x_1,1/x_2,1/x_3,\ldots$ is a Fibonacci sequence.

8. ✍ Provide the details that show how to derive {eq}`secanterr` from {eq}`secant`.

9. ⌨ Write a function `iqi(f,x1,x2,x3)` that performs inverse quadratic interpolation for finding a root of $f$, given three initial estimates. To find the quadratic polynomial $q(y)$ passing through the three most recent points, use `fit` from the `Polynomials` package. Test your function on the first problem
from this section.
