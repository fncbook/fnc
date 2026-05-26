---
numbering:
  enumerator: 9.3.%s
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

(section-globalapprox-stability)=

# Stability of polynomial interpolation

```{index} stability; of polynomial interpolation
```

With barycentric interpolation available in the form of {numref}`Function {number} <function-polyinterp>`, we can explore polynomial interpolation using a numerically stable algorithm. Any remaining sensitivity to error is due to the conditioning of the interpolation process itself.

::::{prf:example} Ill-conditioning in polynomial interpolation
:label: demo-stability-equispaced

We choose a function over the interval $[0,1]$. 

```{code-cell} 
f = lambda x: sin(exp(2 * x))
```

Here is a graph of $f$ and its polynomial interpolant using seven equally spaced nodes.

```{code-cell} 
:tags: hide-input
x = linspace(0, 1, 500)
plot(x, f(x), label="function")
t = linspace(0, 1, 7)
y = f(t)
p = FNC.polyinterp(t, y)
plot(x, p(x), label="interpolant")
plot(t, y, 'ko', label="nodes")
legend(),  title("Equispaced interpolant, n=6");
```

This looks pretty good. We want to track the behavior of the error as $n$ increases. We will estimate the error in the continuous interpolant by sampling it at a large number of points and taking the max-norm.

```{code-cell} 
:tags: hide-input
N = arange(5, 65, 5)
err = zeros(N.size)
x = linspace(0, 1, 1001)         # for measuring error
for k, n in enumerate(N):
    t = linspace(0, 1, n + 1)    # equally spaced nodes
    y = f(t)  # interpolation data
    p = FNC.polyinterp(t, y)
    err[k] = max(abs(f(x) - p(x)))

semilogy(N, err, "-o")
xlabel("$n$"),  ylabel("max error")
title(("Polynomial interpolation error"));
```

The error initially decreases as one would expect but then begins to grow. Both phases occur at rates that are exponential in $n$, i.e., $O(K^n$) for a constant $K$, appearing linear on a semi-log plot.

::::

(section-stability-runge)=

```{index} ! Runge phenomenon
```

## Runge phenomenon

The disappointing loss of convergence in @demo-stability-equispaced is a sign of ill conditioning due to the use of equally spaced nodes. We will examine this effect using the error formula {eq}`interperror` as a guide:

$$
f(x) - p(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!} \Phi(x), \qquad \Phi(x) = \prod_{i=0}^n (x-t_i).
$$

While the dependence on $f$ is messy here, the error indicator $\Phi(x)$ can be studied as a function of the nodes only.

::::{prf:example} Error indicator function for equispaced nodes
:label: demo-stability-errfun

We plot $|\Phi(x)|$ over the interval $[-1,1]$ with equispaced nodes for different values of $n$. 

```{code-cell} 
:tags: hide-input
x = linspace(-1, 1, 1601)
for n in range(10, 60, 10):
    t = linspace(-1, 1, n + 1)
    Phi = array([prod(xk - t) for xk in x])
    semilogy(x, abs(Phi), ".", markersize=2)
xlabel("$x$")
ylabel("$|\Phi(x)|$")
ylim([1e-25, 1])
title(("Effect of equispaced nodes"));
```

Each time $\Phi$ passes through zero at an interpolation node, the value on the log scale should go to $-\infty$, which explains the numerous cusps on the curves.

::::

Two observations from the result of @demo-stability-errfun are important. First, $|\Phi|$ decreases exponentially at each fixed location in the interval (note that the spacing between curves is constant for constant increments of $n$). Second, $|\Phi|$ is larger at the ends of the interval than in the middle, by an exponentially growing factor. This gap is what can ruin the convergence of polynomial interpolation.

::::{prf:example} Runge phenomenon
:label: demo-stability-runge

This function has infinitely many continuous derivatives on the entire real line and looks easy to approximate over $[-1,1]$.

```{code-cell} 
f = lambda x: 1 / (x**2 + 16)
x = linspace(-1, 1, 1601)
plot(x, f(x))
title("Test function");
```

We start by doing equispaced polynomial interpolation for some small values of $n$.

```{code-cell} 
:tags: hide-input
N = arange(4, 16, 4)
label = []
for k, n in enumerate(N):
    t = linspace(-1, 1, n + 1)  # equally spaced nodes
    y = f(t)  # interpolation data
    p = FNC.polyinterp(t, y)
    err = abs(f(x) - p(x))
    semilogy(x, err, ".", markersize=2)
    label.append(f"degree {n}")

xlabel("$x$"),  ylabel("$|f(x)-p(x)|$")
ylim([1e-20, 1])
legend(label),  title("Error for low degrees");
```

The convergence so far appears rather good, though not uniformly so. However, notice what happens as we continue to increase the degree.

```{code-cell} 
:tags: hide-input
N = 12 + 15 * arange(1, 4)
labels = []
for k, n in enumerate(N):
    t = linspace(-1, 1, n + 1)  # equally spaced nodes
    y = f(t)  # interpolation data
    p = FNC.polyinterp(t, y)
    err = abs(f(x) - p(x))
    semilogy(x, err, ".", markersize=2)
    labels.append(f"degree {n}")
xlabel("$x$"),  ylabel("$|f(x)-p(x)|$"),  ylim([1e-20, 1])
legend(labels),  title("Error for higher degrees");
```

The convergence in the middle can't get any better than machine precision relative to the function values. So maintaining the growing gap between the center and the ends pushes the error curves upward exponentially fast at the ends, wrecking the convergence.

::::

```{index} ! Runge phenomenon
```

The observation of instability in @demo-stability-runge is known as the **Runge phenomenon**. The Runge phenomenon is an instability manifested when the nodes of the interpolant are equally spaced and the degree of the polynomial increases. We reiterate that the phenomenon is rooted in the interpolation convergence theory and not a consequence of the algorithm chosen to implement polynomial interpolation.

Significantly, the convergence observed in @demo-stability-runge is stable within a middle portion of the interval. By redistributing the interpolation nodes, we will next sacrifice a little of the convergence in the middle portion in order to improve it near the ends and rescue the process globally.

## Chebyshev nodes

```{index} ! Chebyshev points; second kind
```

The observations above hint that we might find success by having more nodes near the ends of the interval than in the middle. Though we will not give the details, it turns out that there is a precise asymptotic sense in which this must be done to make polynomial interpolation work over the entire interval. One especially important node family that gives stable convergence for polynomial interpolation is the **Chebyshev points of the second kind** (or *Chebyshev extreme points*) defined by

```{math}
:label: chebextreme
 t_k = - \cos\left(\frac{k \pi}{n}\right), \qquad k=0,\ldots,n.
```

These are the projections onto the $x$-axis of $n$ equally spaced points on a unit circle. They are densely clustered near the ends of $[-1,1]$, and this feature turns out to overcome the Runge phenomenon.

::::{prf:example} Error indicator function for Chebyshev nodes
:label: demo-stability-errcheb

Now we look at the error indicator function $\Phi$ for Chebyshev node sets.

```{code-cell} 
:tags: hide-input
x = linspace(-1, 1, 1601)
labels = []
for n in range(10, 60, 10):
    theta = pi * arange(n + 1) / n
    t = -cos(theta)
    Phi = array([prod(xk - t) for xk in x])
    semilogy(x, abs(Phi), ".")
    labels.append(f"degree {n}")

xlabel("$x$"),  ylabel("$|\\Phi(x)|$"),  ylim([1e-18, 1e-2])
legend(labels),  title("Error indicator for Chebyshev nodes");
```

In contrast to the equispaced case, $|\Phi|$ decreases exponentially with $n$ almost uniformly across the interval.

::::

::::{prf:example} Stability of interpolation with Chebyshev nodes
:label: demo-stability-rungefix

Here again is the function from @demo-stability-runge that provoked the Runge phenomenon when using equispaced nodes.

```{code-cell} 
f = lambda x: 1 / (x**2 + 16)
```

```{code-cell} 
:tags: hide-input
x = linspace(-1, 1, 1601)
labels = []
for k, n in enumerate([4, 10, 16, 40]):
    t = -cos(pi * arange(n + 1) / n)         # Chebyshev nodes
    y = f(t)                                 # interpolation data
    p = FNC.polyinterp(t, y)
    err = abs(f(x) - p(x))
    semilogy(x, err, ".", markersize=2)
    labels.append(f"degree {n}")

xlabel("$x$"),  ylabel("$|f(x)-p(x)|$"),  ylim([1e-20, 1])
legend(labels),  title("Error for Chebyshev interpolants");
```

By degree 16 the error is uniformly within machine epsilon, and, importantly, it stays there as $n$ increases. Note that as predicted by the error indicator function, the error is uniform over the interval at each value of $n$.

::::

As a bonus, for Chebyshev nodes the barycentric weights are simple:

```{math}
:label: weightcheb
  w_k = (-1)^k d_k, \qquad d_k =
  \begin{cases}
    1/2 & \text{if $k=0$ or $k=n$},\\
    1 & \text{otherwise}.
  \end{cases}
```

## Spectral convergence

If we take $n\rightarrow \infty$ and use polynomial interpolation on Chebyshev nodes, the convergence rate is exponential in $n$. The following is typical of the results that can be proved.

::::{prf:theorem} Spectral convergence
:label: theorem-spectral
Suppose $f(x)$ is analytic in an open real interval containing $[-1,1]$. Then there exist constants $C>0$ and $K>1$ such that
  
```{math}
:label: spectral
\max_{x\in[-1,1]} | f(x) - p(x) | \le C K^{-n},
```

where $p$ is the unique polynomial of degree $n$ or less defined by interpolation on $n+1$ Chebyshev second-kind points.
::::

The condition "$f$ is analytic" means that the Taylor series of $f$ converges to $f(x)$ in an open interval containing $[-1,1]$.[^analytic] A necessary condition of analyticity is that $f$ is infinitely differentiable.

[^analytic]: Alternatively, analyticity means that the function is extensible to one that is differentiable in the complex plane.

```{index} ! convergence rate; spectral
```

In other contexts we refer to {eq}`spectral` as linear convergence, but here it is usual to say that the rate is exponential or that one has {term}`spectral convergence`. It achieves constant reduction factors in the error by constant increments of $n$. By contrast, algebraic convergence in the form $O(n^{-p})$ for some $p>0$ requires *multiplying* $n$ by a constant factor in order to reduce error by a constant factor. Graphically, spectral error is a straight line on a log-linear scale, while algebraic convergence is a straight line on a log-log scale.

::::{prf:example} Spectral convergence
:label: demo-stability-spectral

On the left, we use a log-log scale, which makes fourth-order algebraic convergence $O(n^{-4})$ a straight line. On the right, we use a log-linear scale, which makes spectral convergence $O(K^{-n})$ a straight line.


```{code-cell} 
:tags: hide-input
n = arange(20, 420, 20)
algebraic = 100 / n**4
spectral = 10 * 0.85**n

subplot(2, 1, 1)
loglog(n, algebraic, 'o-', label="algebraic")
loglog(n, spectral, 'o-', label="spectral")
xlabel('n'),  ylabel('error') 
title('log–log'), ylim([1e-16, 1]);
legend() 

subplot(2, 1, 2)
semilogy(n, algebraic, 'o-', label="algebraic")
semilogy(n, spectral, 'o-', label="spectral")
xlabel('n'),  ylabel('error')
title('log–linear') ,  ylim([1e-16, 1]);  
legend();
```

::::

## Exercises

``````{exercise}
:label: problem-stability-rate
⌨ Revisit @demo-stability-equispaced and determine an approximate value for the convergent phase of the constant $K$ mentioned in the comments there. 
``````

``````{exercise}
:label: problem-stability-chebyshev
⌨ For each case, compute the polynomial interpolant using $n$ second-kind Chebyshev nodes in $[-1,1]$ for $n=4,8,12,\ldots,60$. At each value of $n$, compute the infinity-norm error (that is, $\max |p(x)-f(x)|$ evaluated for at least 4000 values of $x$). Using a log-linear scale, plot the error as a function of $n$, then determine a good approximation to the constant $K$ in {eq}`spectral`. 

**(a)** $f(x) = 1/(25x^2+1)\qquad$ 
**(b)** $f(x) = \tanh(5 x+2)$

**(c)** $f(x) = \cosh(\sin x)\qquad$
**(d)** $f(x) = \sin(\cosh x)$
``````

``````{exercise}
:label: problem-stability-chebinterp
⌨ Write a function `chebinterp(f,n)` that returns a function representing the polynomial interpolant of the input function `f` using $n+1$ Chebyshev second kind nodes over $[-1,1]$. You should use {eq}`weightcheb` to compute the barycentric weights directly, rather than using the method in {numref}`Function {number} <function-polyinterp>`. Test your function by revisiting @demo-stability-runge to use Chebyshev rather than equally spaced nodes. 
``````

``````{exercise}
:label: problem-stability-nonsmooth
@theorem-spectral assumes that the function being approximated has infinitely many derivatives over $[-1,1]$. But now consider the family of functions $f_m(x)=|x|^m$. 

**(a)** ✍ How many continuous derivatives over $[-1,1]$ does $f_m$ possess?

**(b)** ⌨ Compute the polynomial interpolant using $n$ second-kind Chebyshev nodes in $[-1,1]$ for $n=10,20,30,\ldots,100$. At each value of $n$, compute the max-norm error (that is, $\max |p(x)-f_m(x)|$ evaluated for at least 5000 values of $x$). Using a single log-log graph, plot the error as a function of $n$ for all six values $m=1,3,5,7,9$. Each convergence curve should be nearly a straight line.

**(c)** ✍ Based on the results of parts (a) and (b), form a hypothesis about the asymptotic behavior of the error for fixed $m$ as $n\rightarrow \infty$. (Hint: Find the slopes of the lines in part (b).)
``````

``````{exercise}
:label: problem-stability-changeinterval
The Chebyshev points can be used when the interval of interpolation is $[a,b]$ rather than $[-1,1]$ by means of the change of variable

```{math}
:label: changeinterval
z = \psi(x) = a + (b-a)\frac{(x+1)}{2}.
```

**(a)** ✍  Show that $\psi(-1) = a$, $\psi(1) = b$, and $\psi$ is strictly increasing on $[-1,1]$.

**(b)** ✍ Invert the relation {eq}`changeinterval` to solve for $x$ in terms of $\psi^{-1}(z)$. 

**(c)** ✍ Let $t_0,\ldots,t_n$ be standard second-kind Chebyshev points. Then a polynomial in $x$ can be used to interpolate the function values $f\bigl(\psi(t_i)\bigr)$. This in turn implies an interpolating function $\tilde{P}(z) =P\bigl(\psi^{-1}(z)\bigr)$. Show that $\tilde{P}$ is a polynomial in $z$. 

**(d)** ⌨ Implement the idea of part (c) to plot a polynomial interpolant of $f(z) =\cosh(\sin z)$ over $[0,2\pi]$ using $n+1$ Chebyshev nodes with $n=40$. 
``````

``````{exercise}
:label: problem-stability-infinite
The Chebyshev points can be used for interpolation of functions defined on the entire real line by using the change of variable

```{math}
:label: changeintervalinf
z = \phi(x) = \frac{2x}{1-x^2},
```

which maps the interval $(-1,1)$ in one-to-one fashion to the entire real line.

**(a)** ✍ Find $\displaystyle \lim_{x\to 1^-} \phi(x)$ and $\displaystyle \lim_{x\to -1^+} \phi(x)$. 

**(b)** ✍ Invert {eq}`changeintervalinf` to express $x=\phi^{-1}(z)$. (Be sure to enforce $-1\le x \le 1$.)

**(c)** ⌨ Let $t_0,\ldots,t_n$ be standard second-kind Chebyshev points. These map to the $z$ variable as $\zeta_i=\phi(t_i)$ for all $i$. Suppose that $f(z)$ is a given function whose domain is the entire real line. Then the function values $y_i=f(\zeta_i)$ can be associated with the Chebyshev nodes $t_i$, leading to a polynomial interpolant $p(x)$. This in turn implies an interpolating function on the real line, defined as

$$
q(z)=p\bigl(\phi^{-1}(z)\bigr) = p(x).
$$

Implement this idea to plot an interpolant of $f(z)=(z^2-2z+2)^{-1}$ using $n=30$. Your plot should show $q(z)$ evaluated at 1000 evenly spaced points in $[-6,6]$, with markers at the nodal values (those lying within the $[-6,6]$ window).
``````
