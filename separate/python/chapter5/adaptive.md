---
numbering:
  enumerator: 5.7.%s
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

(section-localapprox-adaptive)=

# Adaptive integration

```{index} numerical integration
```

To this point, we have used only equally spaced nodes to compute integrals. Yet there are problems in which non-uniformly distributed nodes would clearly be more appropriate, as demonstrated in @demo-adapt-motive.

::::{prf:example} Motivation for adaptive integration
:label: demo-adapt-motive

This function gets increasingly oscillatory as $x$ increases.

```{code-cell}
f = lambda x: (x + 1) ** 2 * cos((2 * x + 1) / (x - 4.3))
x = linspace(0, 4, 600)
plot(x, f(x))
xlabel("$x$")
ylabel("$f(x)$");
```

Accordingly, the trapezoid rule is more accurate on the left half of this interval than on the right half.

```{code-cell}
n_ = 50 * 2 ** arange(4)
Tleft = zeros(4)
Tright = zeros(4)
for i, n in enumerate(n_):
    Tleft[i] = FNC.trapezoid(f, 0, 2, n)[0]
    Tright[i] = FNC.trapezoid(f, 2, 4, n)[0]
print("left half:", Tleft)
print("right half:", Tright)
```

```{code-cell}
from scipy.integrate import quad
left_val, err = quad(f, 0, 2, epsabs=1e-13, epsrel=1e-13)
right_val, err = quad(f, 2, 4, epsabs=1e-13, epsrel=1e-13)

print("    n     left error   right error")
for k in range(n_.size):
    print(f"  {n_[k]:4}    {Tleft[k]-left_val:8.3e}    {Tright[k]-right_val:8.3e}")
```

Both the picture and the numerical results suggest that more nodes should be used on the right half of the interval than on the left half.

::::

We would like an algorithm that automatically detects and reacts to a situation like that in @demo-adapt-motive, a trait known as **adaptivity**.

## Error estimation

```{index} Simpson's formula
```

Ideally, we would like to make adaptation decisions based on the error of the integration result. Knowing the error exactly would be equivalent to knowing the exact answer, but we can estimate it using the extrapolation technique of {numref}`section-localapprox-integration`. Consider the Simpson formula {eq}`extraplevel1` resulting from one level of extrapolation from trapezoid estimates:

```{math}
:label: extraplevel1repeat
  S_f(2n) = \frac{1}{3} \Bigl[ 4 T_f(2n) - T_f(n) \Bigr].
```

We expect this method to be fourth-order accurate, i.e.,

```{math}
  \int_a^b f(x)\, dx = S_f(2n) + O(n^{-4}),
```

We can further extrapolate to sixth-order accuracy using {eq}`nc-sixth`:

```{math}
:label: extraplevel2repeat
  R_f(4n) = \frac{1}{15} \Bigl[ 16 S_f(4n) - S_f(2n) \Bigr].
```

By virtue of higher order of accuracy, $R_f(4n)$ should be more accurate than $S_f(4n)$. Hence, a decent estimate of the error in the better of the two Simpson values is

```{math}
:label: adapterr
  E = R_f(4n) - S_f(4n) = \frac{S_f(4n) - S_f(2n)}{15}.
```

```{index} adaptivity; in integration, divide and conquer
```

## Divide and conquer

If $|E|$ is judged to be acceptably small, we are done. This judgment takes some care. For instance, suppose the exact integral is $10^{20}$.  Requiring $|E| < \delta\ll 1$ would be fruitless in double precision, since it would require more than 20 accurate digits. Hence, checking the absolute size of the error alone is not appropriate. Conversely, consider the integral

```{math}
  \int_{10^{-6}}^{2\pi} 2 \sin x\, dx \approx -10^{-12}.
```

We are likely to sample values of the integrand that are larger than, say, $1/2$ in absolute value, so obtaining this very small result has to rely on subtractive cancellation. We cannot hope for more than 4-5 accurate digits, so a strict test of the relative error is also not recommended. In other words, we can seek an error that is small relative to the data (the integrand), which is $O(1)$, but not relative to the answer itself.

Typically, we use both relative and absolute error, stopping when either one is considered small enough. Algebraically, the test is

```{math}
:label: absreltolerance
  |E| < \delta_a + \delta_r |S_f(n)|,
```

where $\delta_a$ and $\delta_r$ are given absolute and relative error tolerances, respectively.

When $|E|$ fails to meet {eq}`absreltolerance`, we bisect the interval $[a,b]$ to exploit the identity

```{math}
  \int_a^b f(x)\, dx = \int_a^{(a+b)/2} f(x)\, dx + \int_{(a+b)/2}^b f(x)\, dx,
```

and independently compute estimates to each of the half-length integrals. Each of these half-sized computations recursively applies Simpson's formula and the error estimation criterion, making further bisections as necessary. Such an approach is called **divide and conquer** in computer science: recursively split the problem into easier pieces and glue the results together.

## Implementation

It is typical to use just the minimal formula $S_f(4)$ and its error estimate $E$ to make decisions about adaptivity. A computation of $S_f(4)$ requires three trapezoid estimates $T_f(1)$, $T_f(2)$, and $T_f(4)$. As observed in {eq}`nc-doubling` and @demo-int-extrap, the five integrand evaluations in $T_f(4)$ are sufficient to compute all of these values.

{numref}`Function {number} <function-intadapt>` shows an implementation. It uses five function values to compute three trapezoid estimates with $n=1$, $n=2$, and $n=4$, applying the updating formula {eq}`nc-doubling` twice. It goes on to find the two Simpson approximations and to estimate the error by {eq}`adapterr`.

If the error estimate passes the test {eq}`absreltolerance`, the better Simpson value is returned as the integral over the given interval. Otherwise, the interval is bisected, integrals over the two pieces are computed using recursive calls, and those results are added to give the complete integral.

```{index} Julia; keyword function arguments
```

(function-intadapt)=

``````{prf:algorithm} intadapt

```{literalinclude} chapter05.py
:filename: intadapt.py
:start-at: def intadapt
:language: python
:linenos: true
```
:::{admonition} About the code
:class: dropdown
The intended way for a user to call {numref}`Function {number} <function-intadapt>` is with only `f`, `a`, `b`, and `tol` provided. We then use default values on the other parameters to compute the function values at the endpoints, the interval's midpoint, and the function value at the midpoint. Recursive calls from within the function itself will provide all of that information, since it was already calculated along the way.
:::
``````

::::{prf:example} Using adaptive integration
:label: demo-adapt-usage

We'll integrate the function from @demo-adapt-motive.

```{code-cell}
from scipy.integrate import quad
f = lambda x: (x + 1) ** 2 * cos((2 * x + 1) / (x - 4.3))
I, errest = quad(f, 0, 4, epsabs=1e-12, epsrel=1e-12)
print("integral:", I)    # 'exact' value
```

We perform the integration and show the nodes selected underneath the curve.

```{code-cell}
Q, t = FNC.intadapt(f, 0, 4, 0.001)
print("number of nodes:", t.size)

x = linspace(0, 4, 600)
plot(x, f(x), "k")
stem(t, f(t))
xlabel("$x$"); ylabel("$f(x)$");
```

The error turns out to be a bit more than we requested. It's only an estimate, not a guarantee.

```{code-cell}
print("error:", I - Q)
```

Let's see how the number of integrand evaluations and the error vary with the requested tolerance.

```{code-cell}
tol_ = 10.0 ** arange(-4, -12, -1)
err_ = zeros(tol_.size)
num_ = zeros(tol_.size, dtype=int)
print("    tol         error     # f-evals")
for i, tol in enumerate(tol_):
    Q, t = FNC.intadapt(f, 0, 4, tol)
    err_[i] = I - Q
    num_[i] = t.size
    print(f"  {tol:6.1e}    {err_[i]:10.3e}    {num_[i]:6d}")
```

As you can see, even though the errors are not smaller than the estimates, the two columns decrease in tandem. If we consider now the convergence not in $h$, which is poorly defined now, but in the number of nodes actually chosen, we come close to the fourth-order accuracy of the underlying Simpson scheme.

```{code-cell}
loglog(num_, abs(err_), "-o", label="results")
order4 = 0.01 * (num_ / num_[0]) ** (-4)
loglog(num_, order4, "--", label="$O(n^{-4})$")
xlabel("number of nodes"), ylabel("error")
legend()
title("Convergence of adaptive quadrature");
```

::::

Although adaptivity and the error estimation that goes with it can be very powerful, they come at some cost. The error estimation cannot be universally perfect, so sometimes the answer will not be as accurate as requested, and sometimes the function will be evaluated more times than necessary. Subtle problems may arise when the integral is a step within a larger computation (see @problem-adaptive-nonsmooth).

## Exercises

% must be kept as #1

``````{exercise}
:label: problem-adaptive-tests
⌨ For each integral below, use {numref}`Function {number} <function-intadapt>` with error tolerance $10^{-2},10^{-3},\ldots,10^{-12}$. Make a table of errors and the number of integrand evaluation nodes used, and use a convergence plot as in @demo-adapt-usage to compare to fourth-order accuracy. (These integrals were taken from {cite}`baileyComparisonThree2005`.)

**(a)** $\displaystyle \int_0^1 x\log(1+x)\, dx = \frac{1}{4}$

**(b)** $\displaystyle \int_0^1 x^2 \tan^{-1}x\, dx = \frac{\pi-2+2\log 2}{12}$

**(c)** $\displaystyle \int_0^{\pi/2}e^x \cos x\, dx = \frac{e^{\pi/2}-1}{2}$

**(d)** $\displaystyle \int_{0}^1 \sqrt{x} \log(x) \, dx = -\frac{4}{9}$ (Note: Although the integrand has the limiting value zero as $x\to 0$, you have to implement the function carefully to return zero as the value of $f(0)$, or start the integral at $x=\macheps$.)

**(e)** $\displaystyle \int_0^1 \sqrt{1-x^2}\, dx = \frac{\pi}{4}$
``````

``````{exercise}
:label: problem-adaptive-estimate
⌨ For each integral below: (i) imitate @demo-adapt-usage to find the true value to at least 12 digits; (ii) use {numref}`Function {number} <function-intadapt>` to evaluate the integral to a tolerance of $10^{-8}$; (iii) compute the absolute error and the number of nodes used; (iv) use the $O(h^2)$ term in the Euler–Maclaurin formula {eq}`eulermaclaurin` to estimate how many nodes are required by the fixed-stepsize trapezoidal formula to reach an absolute error of $10^{-8}$.

**(a)** $\displaystyle \int_{0.1}^3 \operatorname{sech}(\sin(1/x))\, d x$

**(b)** $\rule[2em]{0pt}{0pt} \displaystyle\int_{-0.9}^9 \ln((x+1)^3))\, d x$

**(c)** $\rule[2em]{0pt}{0pt} \displaystyle\int_{-\pi}^\pi \cos(x^3)\, d x$

``````

```{index} improper integral
```

``````{exercise}
:label: problem-adaptive-improper
⌨ An integral such as $\displaystyle \int_0^1 x^{-\gamma}\, dx$ for $\gamma>0$, in which the integrand blows up at one or both ends, is known as an *improper* integral. It has a finite value if $\gamma<1$, despite the singularity. One way to deal with the problem of the infinite value for $f(t_0)$ is to replace the lower limit with a small number $\epsilon$. (A more robust way to handle improper integrals is discussed in Chapter 9.)

Using {numref}`Function {number} <function-intadapt>` with tolerance $10^{-14}$, make a log-log plot of the error as a function of $\epsilon$ when $\gamma=2/3$, for $\epsilon=10^{-15},10^{-16},\ldots,10^{-45}$. 
``````

``````{exercise}
:label: problem-adaptive-extrapolation
⌨ A curious consequence of our logic in {numref}`Function {number} <function-intadapt>` is that the algorithm uses what we believe to be a more accurate, sixth-order answer only for estimating error; the returned value is the supposedly less accurate $S_f(2n)$. The practice of returning the extrapolated $R_f(4n)$ instead is called *local extrapolation*. 

Modify {numref}`Function {number} <function-intadapt>` to use local extrapolation and repeat parts (a) and (e) of @problem-adaptive-tests, comparing the observed convergence to both fourth order and sixth order.
``````

```{index} sine integral function
```

``````{exercise}
:label: problem-adaptive-si
⌨ The *sine integral function* is defined by

```{math}
:numbered: false
\operatorname{Si}(x) = \int_0^x \frac{\sin z}{z}\, dz.
```

Use {numref}`Function {number} <function-intadapt>` to plot Si over the interval $[1,10]$. Note: You will need to replace the lower bound of integration by $\macheps$.
``````

``````{exercise}
:label: problem-adaptive-nonsmooth
⌨  Adaptive integration can have subtle drawbacks. This exercise is based on the *error function*, a smooth function defined as

```{math}
:numbered: false
\operatorname{erf}(x) = \frac{2}{\sqrt{\pi}}\int_0^x e^{-s^2}\,ds.
```

**(a)** Define a function $g$ that approximates erf by applying {numref}`Function {number} <function-trapezoid>` with $n=300$. Make a plot of the error $g(x)-\operatorname{erf}(x)$ at 500 points in the interval $[0,3]$.

**(b)** Define another approximation $h$ that applies {numref}`Function {number} <function-intadapt>` with error tolerance $10^{-7}$. Plot the error in $h$ as in part (a). Why does it look so different from the previous case?

**(c)** Suppose you wished to find $x$ such that $\operatorname{erf}(x) = .95$ by using rootfinding on one of your two approximations. Why is the version from part (a) preferable?
``````
