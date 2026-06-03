---
numbering:
  enumerator: 4.4.%s
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


```{code-cell}
f = lambda x: x * exp(x) - 2
xx = linspace(0.25, 1.25, 400)

fig, ax = subplots()
ax.plot(xx, f(xx), label="function")
ax.set_xlabel("$x$")
ax.set_ylabel("$f(x)$")
ax.grid();
```

From the graph, it's clear that there is a root near $x=1$. To be more precise, there is a root in the interval $[0.5,1]$. So let us take the endpoints of that interval as _two_ initial approximations.

```{code-cell}
x1 = 1
y1 = f(x1)
x2 = 0.5
y2 = f(x2)
ax.plot([x1, x2], [y1, y2], "ko", label="initial points")
ax.legend()
fig
```

Instead of constructing the tangent line by evaluating the derivative, we can construct a linear model function by drawing the line between the two points $\bigl(x_1,f(x_1)\bigr)$ and $\bigl(x_2,f(x_2)\bigr)$. This is called a _secant line_.

```{code-cell}
slope2 = (y2 - y1) / (x2 - x1)
secant2 = lambda x: y2 + slope2 * (x - x2)
ax.plot(xx, secant2(xx), "--", label="secant line")
ax.legend()
fig
```

As before, the next root estimate in the iteration is the root of this linear model.

```{code-cell}
x3 = x2 - y2 / slope2
ax.plot(x3, 0, "o", label="root of secant")
y3 = f(x3)
print(y3)
ax.legend()
fig
```

For the next linear model, we use the line through the two most recent points. The next iterate is the root of that secant line, and so on.

```{code-cell}
slope3 = (y3 - y2) / (x3 - x2)
x4 = x3 - y3 / slope3
print(f(x4))
```

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

```{literalinclude} chapter04.py
:filename: secant.py
:start-line: 36
:end-line: 66
:language: python
:linenos: true
```
```{admonition} About the code
:class: dropdown
Because we want to observe the convergence of the method, {numref}`Function {number} <function-secant>` stores and returns the entire sequence of root estimates. However, only the most recent two are needed by the iterative formula. This is demonstrated by the use of `y₁` and `y₂` for the two most recent values of $f$.
```
``````

## Convergence

Graphically, a secant line usually looks like a less accurate model of $f$ than the tangent line. How will that affect the convergence?

As before, let $\epsilon_k = x_k-r$ be the errors in the successive root approximations, and assume that $r$ is a simple root. If the initial errors are small, then a lengthy but straightforward Taylor expansion shows that, to lowest order,

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

We check the convergence of the secant method from @demo-secant-line.

```{code-cell}
f = lambda x: x * exp(x) - 2
x = FNC.secant(f, 1, 0.5)
print(x)
```

We don't know the exact root, so we use `root_scalar` to get a substitute.

```{code-cell}
from scipy.optimize import root_scalar
r = root_scalar(f, bracket=[0.5, 1]).root
print(r)
```

Here is the sequence of errors.

```{code-cell}
err = r - x
print(err)
```

It's not easy to see the convergence rate by staring at these numbers. We can use {eq}`superlinear-rate` to try to expose the superlinear convergence rate.

```{code-cell}
logerr = log(abs(err))
for i in range(len(err) - 2):
    print(logerr[i+1] / logerr[i])
```

As expected, this settles in at around 1.618.

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

Here we look for a root of $x+\cos(10x)$ that is close to 1.

```{code-cell}
f = lambda x: x + cos(10 * x)
xx = linspace(0.5, 1.5, 400)
fig, ax = subplots()
ax.plot(xx, f(xx), label="function")
ax.grid()
xlabel("$x$"), ylabel("$y$")
fig
```

We choose three values to get the iteration started.

```{code-cell}
x = array([0.8, 1.2, 1])
y = f(x)
ax.plot(x, y, "ko", label="initial points")
ax.legend()
fig
```

If we were using forward interpolation, we would ask for the polynomial interpolant of $y$ as a function of $x$. But that parabola has no real roots.

```{code-cell}
q = poly1d(polyfit(x, y, 2))  # interpolating polynomial
ax.plot(xx, q(xx), "--", label="interpolant")
ax.set_ylim(-0.1, 3), ax.legend()
fig
```

To do inverse interpolation, we swap the roles of $x$ and $y$ in the interpolation.

```{code-cell}
plot(xx, f(xx), label="function")
plot(x, y, "ko", label="initial points")

q = poly1d(polyfit(y, x, 2))  # inverse interpolating polynomial
yy = linspace(-0.1, 2.6, 400)
plot(q(yy), yy, "--", label="inverse interpolant")

grid(), xlabel("$x$"), ylabel("$y$")
legend();
```

We seek the value of $x$ that makes $y$ zero. This means evaluating $q$ at zero.

```{code-cell}
x = hstack([x, q(0)])
y = hstack([y, f(x[-1])])
print("x:", x, "\ny:", y)
```

We repeat the process a few more times.

```{code-cell}
for k in range(6):
    q = poly1d(polyfit(y[-3:], x[-3:], 2))
    x = hstack([x, q(0)])
    y = hstack([y, f(x[-1])])
print(f"final residual is {y[-1]:.2e}")
```

Here is the sequence of errors.

```{code-cell}
from scipy.optimize import root_scalar
r = root_scalar(f, bracket=[0.9, 1]).root
err = x - r
print(err)
```

The error seems to be superlinear, but subquadratic:

```{code-cell}
logerr = log(abs(err))
for i in range(len(err) - 1):
    print(logerr[i+1] / logerr[i])
```

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

**(e)** ⌨ For one of the roots, use the errors in the resulting sequence to determine numerically whether the convergence is apparently between linear and quadratic.

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
$\sin(\pi x) = 2x$, over $[0.1, 2]$
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
✍ In general, the secant method formula {eq}`secant` cannot be applied if $x_{k}=x_{k-1}$. However, something special happens if $f$ is a quadratic polynomial, i.e., $f(x)=ax^2+bx+c$.

**(a)** Show that in this case, the secant update formula can be simplified to one that is well-defined when $x_{k}=x_{k-1}.$ 

**(b)** Show that the resulting $x_{k+1}$ from part (a) is identical to the result of one step of Newton's method applied to $f$ at $x_k.$
``````

``````{exercise}
:label: problem-secant-fibonacci
✍ Let $f(x)=x^2$. Show that if $(1/x_1)$ and $(1/x_2)$ are positive integers, and the secant iteration is applied, then the sequence $1/x_1,1/x_2,1/x_3,\ldots$ is a Fibonacci sequence, i.e., they satisfy $x_{k+1}=x_k+x_{k-1}$.
``````

``````{exercise}
:label: problem-secant-derivation
✍ Provide the details that show how to derive {eq}`secanterr` from {eq}`secant`.
``````

``````{exercise}
:label: problem-secant-iqi
⌨ Write a function `iqi(f, x1, x2, x3)` that performs inverse quadratic interpolation for finding a root of $f$ given three initial estimates. Test your function on $x^2=e^{-x}$ to find a root in the interval $[0,2]$.
``````
