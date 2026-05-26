---
numbering:
  enumerator: 9.2.%s
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

(section-globalapprox-barycentric)=

# The barycentric formula

The Lagrange formula {eq}`lagrangeinterp` is useful theoretically but not ideal for computation. For each new value of $x$, all the cardinal functions $\ell_k$ must be evaluated at $x$, which requires a product of $n$ terms. Thus, the total work is $O(n^2)$ for every value of $x$. Moreover, the formula is numerically unstable. An alternative version of the formula improves on both issues.

## Derivation

We again will use the error indicator function $\Phi$ from {numref}`Definition {number} <definition-polynomial-indicator>`,

```{math}
  \Phi(x) = \prod_{j=0}^n (x-t_j),
```

as well as a set of values derived from the nodes.

```{index} ! barycentric weights
```

::::{prf:definition} Barycentric weights
:label: definition-barycentricweights
The **barycentric weights** for nodes $x_0,\dots,x_n$ are defined as

```{math}
:label: baryweight
w_k = \frac{1}{\displaystyle \prod_{\substack{j=0\\j\neq k}}^n (t_k - t_j)} = \frac{1}{\Phi'(t_k)}, \qquad
k = 0,\ldots,n.
```

::::

The following formula is the key to efficient and stable evaluation of a polynomial interpolant.

```{index} ! barycentric interpolation formula
```

::::{prf:theorem} Barycentric interpolation formula
:label: theorem-barycentric-formula
Given points $(t_k,y_k)$ for $k=0,\ldots,n$ with all the $t_k$ distinct, the unique polynomial of degree $n$ or less that interpolates the points is

```{math}
:label: bary2
  p(x) = \frac{\displaystyle \sum_{k=0}^n \, \dfrac{w_k y_k}{x-t_k}  }{\displaystyle\sum_{k=0}^n \, \dfrac{w_k}{x-t_k}}.
```

::::

::::{prf:proof}
:enumerated: false

The Lagrange cardinal polynomial {eq}`lagrange` can be written as

```{math}
:label: lagrangealt
  \ell_k(x) = \Phi(x) \frac{w_k}{x-t_k},
```

and thus the interpolating polynomial in {eq}`lagrangeinterp` is

```{math}
:label: bary1
p(x) = \Phi(x) \sum_{k=0}^n \frac{w_k}{x-t_k} y_k.
```

Obviously, the constant function $p(x)\equiv 1$ is its own polynomial interpolant on any set of nodes. The uniqueness of the interpolating polynomial, as proved in @theorem-polyinterp, allows us to plug $y_k=1$ for all $k$ into {eq}`bary1` to obtain

$$
1 = \Phi(x) \sum_{k=0}^n \frac{w_k}{x-t_k}.
$$

This is solved for $\Phi(x)$ and put back into {eq}`bary1` to get {eq}`bary2`.
::::

Equation {eq}`bary2` is certainly an odd-looking way to write a polynomial! Indeed, it is technically undefined when $x$ equals one of the nodes, but in fact, $\lim_{x\to t_k} p(x) = y_k$, so a continuous extension to the nodes is justified. (See @problem-barycentric-limit.)

::::{prf:example}
:label: example-writeoutbary2
Let us write out the barycentric formula for the interpolating polynomial for the quadratic case ($n=2$) for {numref}`Example %s <example-ClassicalLagrange>`.  The weights are computed from {eq}`baryweight`:
  
```{math}
w_0 = \frac{1}{(t_0-t_1)(t_0-t_2)} = \frac{1}{\left(0-\frac{\pi}{6}\right)
\left(0-\frac{\pi}{3}\right)} = \frac{18}{\pi^2},
```

and, similarly, $w_1 = -36/\pi^2$ and $w_2=18/\pi^2$.

Note that in {eq}`bary2`, any common factor in the weights cancels out without affecting the results. Hence, it's a lot easier to use $w_0=w_2=1$ and $w_1=-2$. Then

```{math}
:numbered: false
\begin{align*}
    p(x) & = \frac{\rule[-1.2em]{0pt}{1em} \dfrac{w_0}{x-t_0} y_0  + \dfrac{w_1}{x-t_1} y_1 + \dfrac{w_2}{x-t_2} y_2 }{ \rule{0pt}{1.5em} \dfrac{w_0}{x-t_0} + \dfrac{w_1}{x-t_1} + \dfrac{w_2}{x-t_2}}\\[2ex]
    & =\frac{ \rule[-1.2em]{0pt}{1em}\left( \dfrac{1}{x} \right) 0 -  \left( \dfrac{2}{x-\pi/6} \right) \dfrac{1}{\sqrt{3}} + \left( \dfrac{1}{x-\pi/3} \right) \sqrt{3} }{
        \rule{0pt}{1.6em} \dfrac{1}{x} - \dfrac{2}{x-\pi/6} + \dfrac{1}{x-\pi/3}  }.
\end{align*}
```
  
Further algebraic manipulation could return this expression to the classical Lagrange form derived in {numref}`Example %s <example-ClassicalLagrange>`.
::::

## Implementation

For certain important node distributions, simple formulas for the weights $w_k$ are known. Otherwise, the first task of an implementation is to compute the weights $w_k$, or more conveniently, $w_k^{-1}$.  

We begin with the singleton node set $\{t_0\}$, for which one gets the single weight $w_0=1$. The idea is to grow this singleton into the set of all nodes through a recursive formula. Define $\omega_{k,m-1}$ (for $k< m$) as the inverse of the weight for node $k$ using the set $\{t_0,\ldots,t_{m-1}\}$. Then

$$
\omega_{k,m} = \displaystyle \prod_{\substack{j=0\\j\neq k}}^{m} (t_k - t_j)
     = \omega_{k,m-1} \cdot (t_k-t_{m}), \qquad k=0,1,\ldots,m-1.
$$

A direct application of {eq}`baryweight` can be used to find $\omega_{m,m}$. This process is iterated over $m=1,\ldots,n$ to find $w_k=\omega_{k,n}^{-1}$.

In {numref}`Function {number} <function-polyinterp>` we give an implementation of the barycentric formula for polynomial interpolation.

```{index} ! Julia; isinf
```

``````{prf:algorithm} polyinterp
:label: function-polyinterp

```{literalinclude} chapter09.py
:filename: polyinterp.py
:start-at: def polyinterp
:end-at: return np.vectorize(p)
:language: python
:linenos: true
```

````{admonition} About the code
:class: dropdown
As noted in {numref}`Example %s <example-writeoutbary2>`, a common scaling factor in the weights does not affect the barycentric formula {eq}`bary2`. In lines 9--10 this fact is used to rescale the nodes in order to avoid eventual tiny or enormous numbers that could go outside the bounds of double precision.

The return value is a function that evaluates the polynomial interpolant. Within this function, `isinf` is used to detect either `Inf` or `-Inf`, which occurs when $x$ exactly equals one of the nodes. In this event, the corresponding data value is returned.
````
``````

Computing all $n+1$ weights in {numref}`Function {number} <function-polyinterp>` takes $O(n^2)$ operations. Fortunately, the weights depend only on the nodes, not the data, and once they are known, computing $p(x)$ at a particular value of $x$ takes just $O(n)$ operations.

::::{prf:example} Barycentric interpolation
:label: demo-barycentric-example
We show the barycentric formula in action for values from the function $\sin(e^{2x})$ at equally spaced nodes in $[0,1]$ with $n=3$ and $n=6$.

```{code-cell}
f = lambda x: sin(exp(2 * x))
x = linspace(0, 1, 500)
fig, ax = subplots()
ax.plot(x, f(x), label="function")
```

```{code-cell}
t = linspace(0, 1, 4)
y = f(t)
p = FNC.polyinterp(t, y)

ax.plot(x, p(x), label="interpolant")
ax.plot(t, y, "ko", label="nodes")
ax.legend()
ax.set_title("Interpolation on 4 nodes")
fig
```

The curves must intersect at the interpolation nodes. For $n=6$ the interpolant is noticeably better.

```{code-cell}
plot(x, f(x), label="function")
t = linspace(0, 1, 7)
y = f(t)
p = FNC.polyinterp(t, y)
plot(x, p(x), label="interpolant")
plot(t, y, "ko", label="nodes")
legend(),  title("Interpolation on 7 nodes");
```

::::

## Stability

You might suspect that as the evaluation point $x$ approaches a node $t_k$, subtractive cancellation error will creep into the barycentric formula because of the term $1/(x-t_k)$. While such errors do occur, they turn out not to cause trouble, because the *same* cancellation happens in the numerator and denominator. In fact, the stability of the barycentric formula has been proved, though we do not give the details.

## Exercises

``````{exercise}
:label: problem-barycentric-weights
✍ **(a)** Find the barycentric weights for the nodes $t_0=0$, $t_1=1$, $t_2=3$.

**(b)** Compute the interpolant at $x=2$ for the nodes in part (a) and the data $y_0=-2$, $y_1=2$, $y_2=1$.
``````

``````{exercise}
:label: problem-barycentric-explicit
✍ For each case of @problem-polynomial-lagrange, write out the barycentric form of the interpolating polynomial.
``````

``````{exercise}
:label: problem-barycentric-limit
✍  Show using L'Hôpital's rule on {eq}`bary2` that $p(t_i)=y_i$ for all $i=0,\ldots,n$.
``````

``````{exercise}
:label: problem-barycentric-polyinterp
⌨ In each case, use {numref}`Function {number}<function-polyinterp>` to interpolate the given function using $n+1$ evenly spaced nodes in the given interval. In a 3-by-1 grid of plots, show each interpolant together with the exact function and the data points (i.e., values at the nodes).

**(a)** $f(x) = \ln (x), \quad n = 2,3,4, \quad x\in [1,10]$

**(b)** $f(x) = \tanh (x), \quad n = 2,3,4, \quad x \in [-3,2]$

**(c)** $f(x) = \cosh (x), \quad n = 2,3,4, \quad x \in [-1,3]$

**(d)** $f(x) = |x|, \quad n = 3,5,7, \quad x \in [-2,1]$
``````

``````{exercise}
:label: problem-barycentric-weightsequi
⌨ Using code from {numref}`Function {number}<function-polyinterp>`, compute the barycentric weights numerically using $n+1$ equally spaced nodes in $[-1,1]$ for $n=30$, $n=60$, and $n=90$. On a single graph, plot $|w_i|$ as a function of $t_i$ on a log-linear scale. (The resulting graphs are an indication of the trouble with equally spaced nodes that is explored in {numref}`section-globalapprox-stability`.)
``````

``````{exercise}
:label: problem-barycentric-derive
✍ Derive this fact stated implicitly in {eq}`baryweight`:

$$
\Phi'(t_k) = \prod_{\substack{j=0\\j\neq k}}^n (t_k - t_j).
$$
``````

``````{exercise}
:label: problem-barycentric-derivative
✍ Use {eq}`lagrangealt` to show that if $j\neq k$, 

$$
\ell_k'(t_j) = \frac{w_k}{w_j(t_j-t_k)}.
$$
``````
