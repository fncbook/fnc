---
numbering:
  enumerator: 4.4.%s
---
(section-nonlineqn-secant)=

# Interpolation-based methods

From a practical standpoint, one of the biggest drawbacks of Newton's method is the requirement to supply $f'$ in {numref}`Function {number} <function-newton>`. It is both a programming inconvenience and a step that requires computational time. We can avoid using $f'$, however, by making a simple but easily overlooked observation:

```{prf:observation}
 When a step produces an approximate result, you are free to carry it out approximately. 
```

Let's call this the *principle of approximate approximation.*

In the Newton context, the principle of approximate approximation begins with the observation that the use of $f'$ is linked to the construction of a linear approximation $q(x)$ equivalent to a tangent line. The root of $q(x)$ is used to define the next iterate in the sequence. We can avoid calculating the value of $f'$ by choosing a different linear approximation.

::::{prf:example} Graphical interpretation of the secant method
:label: demo-secant-line

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-secant-line-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-secant-line-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-secant-line-python
:::
```` 
`````

::::

The example in @demo-secant-line demonstrates the {term}`secant method`. In the secant method, one finds the root of the linear approximation through the two most recent root estimates. That is, given previous approximations $x_1,\ldots,x_k$, define the linear model function as the line through $\bigl(x_{k-1},f(x_{k-1})\bigr)$ and $\bigl(x_k,f(x_k)\bigr)$:

```{math}
:label: secantmodel
q(x) = f(x_k) + \frac{f(x_k)-f(x_{k-1})}{x_k-x_{k-1}}(x-x_k).
```

Solving $q(x_{k+1})=0$ for $x_{k+1}$ gives the following iteration formula.

```{index} ! secant method
```

:::{prf:definition} Secant method
:label: definition-secant
Given function $f$ and two initial values $x_1$ and $x_2$, define

```{math}
:label: secant
x_{k+1} = x_k - \frac{f(x_k)(x_k-x_{k-1})}{f(x_k)-f(x_{k-1})}, \quad k=2,3,\ldots.
```

Return the sequence $\{x_k\}$.
:::

Our implementation of the secant method is given in {numref}`Function {number} <function-secant>`.

``````{prf:algorithm} secant
:label: function-secant

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #function-secant-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-secant-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #function-secant-python
:::
````
`````
``````

## Convergence

Graphically, a secant line usually looks like a less accurate model of $f$ than the tangent line. How will that affect the convergence?

As before, let $\epsilon_k = x_k-r$ be the errors in the successive root approximations, and assume that $r$ is a simple root. If the initial errors are small, then a tedious but straightforward Taylor expansion shows that, to lowest order,

```{math}
:label: secanterr
\epsilon_{k+1} \approx \frac{1}{2}\frac{f''(r)}{f'(r)} \epsilon_k \epsilon_{k-1}.
```

If we make an educated guess that

```{math}
\epsilon_{k+1} = c (\epsilon_k)^\alpha, \quad \epsilon_k = c (\epsilon_{k-1})^\alpha, \ldots, \qquad \alpha>0,
```

then {eq}`secanterr` becomes

```{math}
:label: secantexponents
\left[ \epsilon_{k-1}^{\alpha} \right]^{\,\alpha} \approx C \epsilon_{k-1}^{\alpha+1}
```

for an unknown constant $C$. Treating the approximation as an equality, this becomes solvable if and only if the exponents match, i.e., $\alpha^2 = \alpha+1$. The only positive root of this equation is the golden ratio,

```{math}
  \alpha = \frac{1+\sqrt{5}}{2} \approx 1.618.
```

```{index} ! superlinear convergence
```

```{index} ! convergence rate; superlinear
```

Hence, the errors in the secant method converge like $\epsilon_{k+1} = c (\epsilon_k)^\alpha$ for $1<\alpha<2$.

::::{prf:definition} Superlinear convergence
:label: definition-superlinearconvergence
Suppose a sequence $x_k$ approaches limit $x^*$. If the error sequence $\epsilon_k=x_k - x^*$ satisfies

```{math}
:label: superlinear-convergence
  \lim_{k\to\infty} \frac{|\epsilon_{k+1}|}{|\epsilon_k|^\alpha} = L
```

for constants $\alpha >1$ and $L>0$, then the sequence has {term}`superlinear convergence` with rate $\alpha$.
::::

Quadratic convergence is a particular case of superlinear convergence. Roughly speaking, we expect

```{math}
\begin{align*}
\label{superlinear-rate}
\log |\epsilon_{k+1}| & \approx \alpha (\log |\epsilon_k|) + \log L, \\
\frac{\log |\epsilon_{k+1}|}{\log |\epsilon_k|} & \approx \alpha + \frac{\log L}{\log |\epsilon_k|} \to \alpha,
\end{align*}
```

as $k\to\infty$.

::::{prf:example} Convergence of the secant method
:label: demo-secant-converge

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-secant-converge-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-secant-converge-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-secant-converge-python
:::
```` 
`````

::::

In terms of the error as a function of the iteration number $k$, the secant method converges at a rate strictly between linear and quadratic, which is slower than Newton's method. But error versus iteration count may not be the best means of comparison.

Often we analyze rootfinding methods by assuming that the bulk of computing time is spent evaluating the user-defined functions $f$ and $f'$. (Our simple examples and exercises mostly don't support this assumption, but many practical applications do.) In this light, we see that Newton's method requires two evaluations, $f(x_k)$ and $f'(x_k)$, for each iteration. The secant method, on the other hand, while it *uses* the two function values $f(x_k)$ and $f(x_{k-1})$ at each iteration, only needs to *compute* a single new one. Note that {numref}`Function {number} <function-secant>` keeps track of one previous function value rather than recomputing it.

Now suppose that $|\epsilon_k|=\delta$. Roughly speaking, two units of work (i.e., function evaluations) in Newton's method brings us to an error of $\delta^2$. If one spreads out the improvement in the error evenly across the two steps, using

```{math}
\delta^2 = \bigl( \delta^{\sqrt{2}} \bigr)^{\!\sqrt{2}},
```

it seems reasonable to say that the rate of convergence in Newton *per function evaluation* is $\sqrt{2}\approx 1.41$. This is actually less than the comparable rate of about $1.62$ for the secant method.

```{prf:observation}
If function evaluations are used to measure computational work, the secant iteration converges more rapidly than Newton's method.
```

Not only is the secant method easier to apply than Newton's method in practice, it's also more efficient—a rare win-win!

## Inverse interpolation

```{index} interpolation; by polynomials
```

At each iteration, the secant method constructs a linear model function that interpolates the two most recently found points on the graph of $f$. Two points determine a straight line, so this seems like a sensible choice. But as the iteration progresses, why use only the *two* most recent points? What would it mean to use more of them?

If we interpolate through three points by a polynomial, we get a unique quadratic function. Unfortunately, a parabola may have zero, one, or two crossings of the $x$-axis, potentially leaving some doubt as to how to define the next root estimate. On the other hand, if we turn a parabola on its side, we get a graph that intersects the $x$-axis exactly once, which is ideal for defining the next root estimate.

This leads to the idea of defining $q(y)$ as the quadratic interpolant to the points $(y_{k-2},x_{k-2})$, $(y_{k-1},x_{k-1})$, and $(y_k,x_k)$, where $y_i=f(x_i)$ for all $i$, and setting $x_{k+1}=q(0)$. The process defined in this way (given three initial estimates) is called **inverse quadratic interpolation**. Rather than deriving lengthy formulas for it here, we demonstrate how to perform inverse quadratic interpolation using `fit` to perform the interpolation step.

::::{prf:example} Inverse quadratic interpolation
:label: demo-secant-iqi

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-secant-iqi-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-secant-iqi-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-secant-iqi-python
:::
```` 
`````

::::

## Bracketing

Like Newton's method, the secant and inverse quadratic interpolation methods cannot guarantee convergence. One final new idea is needed to make a (nearly) foolproof algorithm.

If $f$ is continuous on the interval $[a,b]$ and $f(a)f(b)<0$—that is, $f$ changes sign on the interval—then $f$ must have at least one root in the interval, due to the Intermediate Value Theorem from calculus. If we come up with a new root estimate $c\in(a,b)$, then whatever sign $f(c)$ is, it is different from the sign at one of the endpoints. (Of course, if $f(c)$ is zero, we are done!) So either $[a,c]$ or $[c,b]$ is guaranteed to have a root too, and in this way we can maintain not just individual estimates but an interval that always contains a root.

The best algorithms blend the use of fast-converging methods with the guarantee provided by a bracket. For example, say that an iteration begins with a bracketing interval. Make a list of the inverse quadratic estimate, the secant estimate, and the midpoint of the current interval, and pick the first member of the list that lies within the current interval. Replace the interval with the bracketing subinterval, and start a new iteration. This is the idea behind **Brent's method**, which is a very successful rootfinding algorithm.

## Exercises

For each of Exercises 1–3, do the following steps.
  
**(a)** ✍ Rewrite the equation into the standard form for rootfinding, $f(x) = 0$.

**(b)** ⌨ Make a plot of $f$ over the given interval and determine how many roots lie in the interval.

**(c)** ⌨ Find a reference value for each root.

**(d)** ⌨ Determine a bracketing interval for each root. Then use {numref}`Function {number} <function-secant>`, starting with the endpoints of the bracketing interval, to find each root.

**(e)** ⌨ For one of the roots, use the errors in the Newton sequence to determine numerically whether the convergence is apparently between linear and quadratic.

``````{exercise}
:label: problem-secant-basic1
$x^2=e^{-x}$, over $[-2,2]$
``````

``````{exercise}
:label: problem-secant-basic2
$2x = \tan x$, over $[-0.2,1.4]$
``````

``````{exercise}
:label: problem-secant-basic3
$\sin(\pi x) = 2x$, over $[0, 2]$
``````

---

``````{exercise}
:label: problem-secant-initial
⌨ Use a plot to approximately locate all the roots of $f(x)=x^{-2}-\sin(x)$ in the interval $[0.5,10]$. Then find a pair of initial points for each root such that {numref}`Function {number} <function-secant>` converges to that root.
``````

``````{exercise}
:label: problem-secant-linear
✍ Show analytically that the secant method converges in one step for a linear function, regardless of the initialization.
``````

``````{exercise}
:label: problem-secant-stuck
✍ In general, the secant method formula {eq}`secant` cannot be applied if $x_{k}=x_{k-1}$. However, suppose that $f(x)=ax^2+bx+c$ for constants $a$, $b$, and $c$. Show that in this case the formula can be simplified to one that is well defined when $x_{k}=x_{k-1}$. Then show that the resulting $x_{k+1}$ is the same as the result of one step of Newton's method applied to $f$ at $x_k$.
``````

``````{exercise}
:label: problem-secant-fibonacci
✍ Let $f(x)=x^2$. Show that if $(1/x_1)$ and $(1/x_2)$ are positive integers, and the secant iteration is applied, then the sequence $1/x_1,1/x_2,1/x_3,\ldots$ is a Fibonacci sequence, i.e., satisfying $x_{k+1}=x_k+x_{k-1}$.
``````

``````{exercise}
:label: problem-secant-derivation
✍ Provide the details that show how to derive {eq}`secanterr` from {eq}`secant`.
``````

``````{exercise}
:label: problem-secant-iqi
 ⌨ Write a function `iqi(f,x₁,x₂,x₃)` that performs inverse quadratic interpolation for finding a root of $f$, given three initial estimates. To find the quadratic polynomial $q(y)$ passing through the three most recent points, use `fit`. Test your function on $x^2=e^{-x}$, over $[-2,2]$.
``````
