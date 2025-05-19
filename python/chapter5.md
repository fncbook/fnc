---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
# Chapter 5

## Functions 

(function-hatfun-python)=
``````{dropdown} Hat function
:open:
```{literalinclude} fncbook/fncbook/chapter05.py
:filename: hatfun.py
:start-at: def hatfun
:end-at: return evaluate
:language: python
:linenos: true
```
``````

(function-plinterp-python)=
``````{dropdown} Piecewise linear interpolation
:open:
```{literalinclude} fncbook/fncbook/chapter05.py
:filename: plinterp.py
:start-at: def plinterp
:end-at: return evaluate
:language: python
:linenos: true
```
``````

(function-spinterp-python)=
``````{dropdown} Cubic spline interpolation
:open:
```{literalinclude} fncbook/fncbook/chapter05.py
:filename: spinterp.py
:start-at: def spinterp
:end-at: return evaluate
:language: python
:linenos: true
```
``````

(function-fdweights-python)=
``````{dropdown} Fornberg's algorithm for finite difference weights
:open:
```{literalinclude} fncbook/fncbook/chapter05.py
:filename: fdweights.py
:start-at: def fdweights
:end-at: return np.array
:language: python
:linenos: true
```
``````

(function-trapezoid-python)=
``````{dropdown} Trapezoid formula for numerical integration
:open:
```{literalinclude} fncbook/fncbook/chapter05.py
:filename: trapezoid.py
:start-at: def trapezoid
:end-at: return T, t, y
:language: python
:linenos: true
```
``````

(function-intadapt-python)=
``````{dropdown} Adaptive integration
:open:
```{literalinclude} fncbook/fncbook/chapter05.py
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

## Examples

```{code-cell} ipython3
:tags: remove-cell
exec(open("FNC_init.py").read());
```

### 5.1 @section-localapprox-interpolation

(demo-interpolation-global-python)=
``````{dropdown} @demo-interpolation-global
:open:
Here are some points that we could consider to be observations of an unknown function on $[-1,1]$.

```{code-cell}
n = 5
t = linspace(-1, 1, n + 1)
y = t**2 + t + 0.05 * sin(20 * t)
fig, ax = subplots()
plot(t, y, "o", label="data")
xlabel("$x$"),  ylabel("$y$");
```

```{index} ! Python; polyfit
```

The polynomial interpolant, as computed using `fit`, looks very sensible. It's the kind of function you'd take home to meet your parents.

```{code-cell}
p = poly1d(polyfit(t, y, n))  # interpolating polynomial
tt = linspace(-1, 1, 400)
ax.plot(tt, p(tt), label="interpolant")
ax.legend()
fig
```

But now consider a different set of points generated in almost exactly the same way.

```{code-cell}
n = 18
t = linspace(-1, 1, n + 1)
y = t**2 + t + 0.05 * sin(20 * t)
fig, ax = subplots()
plot(t, y, "o", label="data")
xlabel("$x$"),  ylabel("$y$");
```

The points themselves are unremarkable. But take a look at what happens to the polynomial interpolant.

```{code-cell}
p = poly1d(polyfit(t, y, n))
ax.plot(tt, p(tt), label="interpolant")
ax.legend()
fig
```

Surely there must be functions that are more intuitively representative of those points!
``````

(demo-interpolation-pwise-python)=
``````{dropdown} @demo-interpolation-pwise
:open:
Let us recall the data from @demo-interpolation-global.

```{code-cell}
clf
n = 12
t = linspace(-1, 1, n + 1)
y = t**2 + t + 0.5 * sin(20 * t)
fig, ax = subplots()
scatter(t, y, label="data")
xlabel("$x$"),  ylabel("$y$");
```

Here is an interpolant that is linear between each consecutive pair of nodes, using `plinterp` from {numref}`section-localapprox-pwlin`.

```{code-cell}
from scipy.interpolate import interp1d
tt = linspace(-1, 1, 400)
p = interp1d(t, y, kind="linear")
ax.plot(tt, p(tt), label="piecewise linear")
ax.legend()
fig
```

```{index} ! Python; interp1d
```

We may prefer a smoother interpolant that is piecewise cubic:

```{code-cell}
scatter(t, y, label="data")
p = interp1d(t, y, kind="cubic")
tt = linspace(-1, 1, 400)
plot(tt, p(tt), label="cubic spline")
xlabel("$x$"),  ylabel("$y$")
legend();
```
``````

(demo-interp-cond-python)=
``````{dropdown} @demo-interp-cond
:open:
In @demo-interpolation-pwise we saw a big difference between polynomial interpolation and piecewise polynomial interpolation of some arbitrarily chosen data. The same effects can be seen clearly in the cardinal functions, which are closely tied to the condition numbers.

```{code-cell}
clf
n = 18
t = linspace(-1, 1, n + 1)
y = zeros(n + 1)
y[9] = 1.0
p = interp1d(t, y, kind="cubic")

scatter(t, y, label="data")
tt = linspace(-1, 1, 400)
plot(tt, p(tt), label="cardinal function")
title("Cubic spline cardinal function")
legend();
```

The piecewise cubic cardinal function is nowhere greater than one in absolute value. This happens to be true for all the cardinal functions, ensuring a good condition number for any interpolation with these functions. But the story for global polynomials is very different.

```{code-cell}
p = poly1d(polyfit(t, y, n))
scatter(t, y, label="data")
plot(tt, p(tt), label="cardinal function")
xlabel("$x$")
ylabel("$y$")
title("Polynomial cardinal function")
legend();
```

From the figure we can see that the condition number for polynomial interpolation on these nodes is at least 500.
``````

### 5.2 @section-localapprox-pwlin

(demo-pwlin-hat-python)=
``````{dropdown} @demo-pwlin-hat
:open:
Let's define a set of four nodes (i.e., $n=3$ in our formulas).

```{code-cell}
t = array([0, 0.075, 0.25, 0.55, 0.7, 1])
```

We plot the hat functions $H_0,\ldots,H_3$.

```{code-cell}
x = linspace(0, 1, 300)
for k in range(6):
    plot(x, FNC.hatfun(t, k)(x))
xlabel("$x$"),  ylabel("$H_k(x)$")
title("Hat functions");
```
``````

(demo-pwlin-usage-python)=
``````{dropdown} @demo-pwlin-usage
:open:
We generate a piecewise linear interpolant of $f(x)=e^{\sin 7x}$.

```{code-cell}
f = lambda x: exp(sin(7 * x))
x = linspace(0, 1, 400)
fig, ax = subplots()
plot(x, f(x), label="function")
xlabel("$x$")
ylabel("$f(x)$");
```

First we sample the function to create the data.

```{code-cell}
t = array([0, 0.075, 0.25, 0.55, 0.7, 1])  # nodes
y = f(t)  # function values

ax.plot(t, y, "o", label="nodes")
ax.legend()
fig
```

Now we create a callable function that will evaluate the piecewise linear interpolant at any $x$, and then plot it.

```{code-cell}
p = FNC.plinterp(t, y)
ax.plot(x, p(x), label="interpolant")
ax.legend()
fig
```
``````

(demo-pwlin-converge-python)=
``````{dropdown} @demo-pwlin-converge
:open:
We measure the convergence rate for piecewise linear interpolation of $e^{\sin 7x}$ over $x \in [0,1]$.

```{code-cell}
f = lambda x: exp(sin(7 * x))
x = linspace(0, 1, 10000)  # sample the difference at many points
N = 2 ** arange(3, 11)
err = zeros(N.size)
for i, n in enumerate(N):
    t = linspace(0, 1, n + 1)  # interpolation nodes
    p = FNC.plinterp(t, f(t))
    err[i] = max(abs(f(x) - p(x)))
print(err)
```

As predicted, a factor of 10 in $n$ produces a factor of 100 in the error. In a convergence plot, it is traditional to have $h$ *decrease* from left to right, so we expect a straight line of slope $-2$ on a log-log plot.

```{code-cell}
order2 = 0.1 * (N / N[0]) ** (-2)
loglog(N, err, "-o", label="observed error")
loglog(N, order2, "--", label="2nd order")
xlabel("$n$")
ylabel("$\|f-p\|_\infty$")
legend();
```
``````

### 5.3 @section-localapprox-splines 

(demo-splines-splines-python)=
``````{dropdown} @demo-splines-splines
:open:
For illustration, here is a spline interpolant using just a few nodes.

```{code-cell}
f = lambda x: exp(sin(7 * x))

x = linspace(0, 1, 500)
fig, ax = subplots()
ax.plot(x, f(x), label="function")

t = array([0, 0.075, 0.25, 0.55, 0.7, 1])  # nodes
y = f(t)  # values at nodes

xlabel("$x$")
ylabel("$y$")
ax.scatter(t, y, label="nodes")
```

```{code-cell}
S = FNC.spinterp(t, y)
ax.plot(x, S(x), label="spline")
ax.legend()
fig
```

Now we look at the convergence rate as the number of nodes increases.

```{code-cell}
N = floor(2 ** linspace(3, 8, 17)).astype(int)
err = zeros(N.size)
for i, n in enumerate(N):
    t = linspace(0, 1, n + 1)  # interpolation nodes
    p = FNC.spinterp(t, f(t))
    err[i] = max(abs(f(x) - p(x)))
print(err)
```

Since we expect convergence that is $O(h^4)=O(n^{-4})$, we use a log-log graph of error and expect a straight line of slope $-4$.

```{code-cell}
order4 = (N / N[0]) ** (-4)
loglog(N, err, "-o", label="observed error")
loglog(N, order4, "--", label="4th order")
xlabel("$n$")
ylabel("$\|f-S\|_\infty$")
legend();
```
``````

### 5.4 @section-localapprox-finitediffs

(demo-finitediffs-fd1-python)=
``````{dropdown} @demo-finitediffs-fd1
:open:
If $f(x)=e^{\,\sin(x)}$, then $f'(0)=1$.

```{code-cell}
f = lambda x: exp(sin(x))
```

Here are the first two centered differences from {numref}`table-FDcenter`.

```{code-cell}
h = 0.05
CD2 = (-f(-h) + f(h)) / (2*h)
CD4 = (f(-2*h) - 8*f(-h) + 8*f(h) - f(2*h)) / (12*h)
print(f"CD2 is {CD2:.9f} and CD4 is {CD4:.9f}")
```

Here are the first two forward differences from {numref}`table-FDforward`.

```{code-cell}
FD1 = (-f(0) + f(h)) / h
FD2 = (-3*f(0) + 4*f(h) - f(2*h)) / (2*h)
print(f"FD1 is {FD1:.9f} and FD2 is {FD2:.9f}")
```

Finally, here are the backward differences that come from reverse-negating the forward differences.

```{code-cell}
BD1 = (-f(-h) + f(0)) / h
BD2 = (f(-2*h) - 4*f(-h) + 3*f(0)) / (2*h)
print(f"BD1 is {BD1:.9f} and BD2 is {BD2:.9f}")
```
``````

(demo-finitediffs-fd2-python)=
``````{dropdown} @demo-finitediffs-fd2
:open:
If $f(x)=e^{\,\sin(x)}$, then $f''(0)=1$.

```{code-cell}
f = lambda x: exp(sin(x))
```

Here is a centered estimate given by {eq}`centerFD22`.

```{code-cell}
h = 0.05
CD2 = (f(-h) - 2*f(0) + f(h)) / h**2
print(f"CD2 is {CD2:.9f}")
```

For the same $h$, here are forward estimates given by {eq}`forwardFD21` and {eq}`forwardFD22`.

```{code-cell}
FD1 = (f(0) - 2*f(h) + f(2*h)) / h**2
FD2 = (2*f(0) - 5*f(h) + 4*f(2*h) - f(3*h)) / h**2
print(f"FD1 is {FD1:.9f} and FD2 is {FD2:.9f}")
```

Finally, here are the backward estimates that come from reversing {eq}`forwardFD21` and {eq}`forwardFD22`.

```{code-cell}
BD1 = (f(-2*h) - 2*f(-h) + f(0)) / h**2
BD2 = (-f(-3*h) + 4*f(-2*h) - 5*f(-h) + 2*f(0)) / h**2
print(f"BD1 is {BD1:.9f} and BD2 is {BD2:.9f}")
```
``````

(demo-finitediffs-fd-weights-python)=
``````{dropdown} @demo-finitediffs-fd-weights
:open:
We will estimate the derivative of $\cos(x^2)$ at $x=0.5$ using five nodes.

```{code-cell}
t = array([0.35, 0.5, 0.57, 0.6, 0.75])   # nodes
f = lambda x: cos(x**2)
dfdx = lambda x: -2 * x * sin(x**2)
exact_value = dfdx(0.5)
```

We have to shift the nodes so that the point of estimation for the derivative is at $x=0$. (To subtract a scalar from a vector, we must use the `.-` operator.)

```{code-cell}
w = FNC.fdweights(t - 0.5, 1)
```

The finite-difference formula is a dot product (i.e., inner product) between the vector of weights and the vector of function values at the nodes.

```{code-cell}
fd_value = dot(w, f(t))
```

We can reproduce the weights in the finite-difference tables by using equally spaced nodes with $h=1$. For example, here is a one-sided formula at four nodes.

```{code-cell}
print(FNC.fdweights(linspace(0, 3, 4), 1))
```
``````

### 5.5 @section-localapprox-fd-converge

(demo-fdconverge-order12-python)=
``````{dropdown} @demo-fdconverge-order12
:open:
Let's observe the convergence of the formulas in {numref}`Example {number} <example-fd-converge-FD11>` and {numref}`Example {number} <example-fd-converge-FD12>`, applied to the function $\sin(e^{x+1})$ at $x=0$.

```{code-cell}
f = lambda x: sin(exp(x + 1))
exact_value = exp(1) * cos(exp(1))
```

We'll compute the formulas in parallel for a sequence of $h$ values.

```{code-cell}
h_ = array([5 / 10**(n+1) for n in range(6)])
FD = zeros((len(h_), 2))
for (i, h) in enumerate(h_):
    FD[i, 0] = (f(h) - f(0)) / h 
    FD[i, 1] = (f(h) - f(-h)) / (2*h)
results = PrettyTable()
results.add_column("h", h_)
results.add_column("FD1", FD[:, 0])
results.add_column("FD2", FD[:, 1])
print(results)
```

All that's easy to see from this table is that FD2 appears to converge to the same result as FD1, but more rapidly. A table of errors is more informative.

```{code-cell}
errors = FD - exact_value
results = PrettyTable()
results.add_column("h", h_)
results.add_column("error in FD1", errors[:, 0])
results.add_column("error in FD2", errors[:, 1])
print(results)
```

In each row, $h$ is decreased by a factor of 10, so that the error is reduced by a factor of 10 in the first-order method and 100 in the second-order method.

A graphical comparison can be useful as well. On a log-log scale, the error should (as $h\to 0$) be a straight line whose slope is the order of accuracy. However, it's conventional in convergence plots to show $h$ _decreasing_ from left to right, which negates the slopes.

```{code-cell}
plot(h_, abs(errors), "o-", label=["FD1", "FD2"])
gca().invert_xaxis()
# Add lines for perfect 1st and 2nd order.
loglog(h_, h_, "--", label="$O(h)$")
loglog(h_, h_**2, "--", label="$O(h^2)$")
xlabel("$h$")
ylabel("error")
legend();
```
``````

(demo-fdconverge-round-python)=
``````{dropdown} @demo-fdconverge-round
:open:
Let $f(x)=e^{-1.3x}$. We apply finite-difference formulas of first, second, and fourth order to estimate $f'(0)=-1.3$.

```{code-cell}
f = lambda x: exp(-1.3 * x)
exact = -1.3

h_ = array([1 / 10**(n+1) for n in range(12)])
FD = zeros((len(h_), 3))
for (i, h) in enumerate(h_):
    nodes = h * linspace(-2, 2, 5)
    vals = f(nodes)
    FD[i, 0] = dot(array([0, 0, -1, 1, 0]) / h, vals)
    FD[i, 1] = dot(array([0, -1/2, 0, 1/2, 0]) / h, vals)
    FD[i, 2] = dot(array([1/12, -2/3, 0, 2/3, -1/12]) / h, vals)

results = PrettyTable()
results.add_column("h", h_)
results.add_column("FD1", FD[:, 0])
results.add_column("FD2", FD[:, 1])
results.add_column("FD4", FD[:, 2])
print(results)
```

They all seem to be converging to $-1.3$. The convergence plot reveals some interesting structure to the errors, though.

```{code-cell}
loglog(h_, abs(FD[:, 0] + 1.3), "-o", label="FD1")
loglog(h_, abs(FD[:, 1] + 1.3), "-o", label="FD2")
loglog(h_, abs(FD[:, 2] + 1.3), "-o", label="FD4")
gca().invert_xaxis()
plot(h_, 0.1 * 2 ** (-52) / h_, "--", color="k", label="$O(h^{-1})$")
xlabel("$h$")
ylabel("total error")
title("FD error with roundoff")
legend();
```

Again the graph is made so that $h$ decreases from left to right. The errors are dominated at first by truncation error, which decreases most rapidly for the fourth-order formula. However, increasing roundoff error eventually equals and then dominates the truncation error as $h$ continues to decrease. As the order of accuracy increases, the crossover point moves to the left (greater efficiency) and down (greater accuracy).
``````

### 5.6 @section-localapprox-integration
(demo-int-antideriv-python)=
``````{dropdown} @demo-int-antideriv
:open:
The antiderivative of $e^x$ is, of course, itself. That makes evaluation of $\int_0^1 e^x\,dx$ by the Fundamental Theorem trivial.

```{code-cell}
exact = exp(1) - 1
```

```{index} ! Python; quad
```

The module `scipy.integrate` has multiple functions that estimate the value of an integral numerically without finding the antiderivative first. As you can see here, it's often just as accurate.


```{code-cell}
from scipy.integrate import quad
Q, errest = quad(exp, 0, 1, epsabs=1e-13, epsrel=1e-13)
print(Q)

```

The numerical approach is also far more robust. For example, $e^{\,\sin x}$ has no useful antiderivative. But numerically, it's no more difficult.

```{code-cell}
Q, errest = quad(lambda x: exp(sin(x)), 0, 1, epsabs=1e-13, epsrel=1e-13)
print(Q)
```

When you look at the graphs of these functions, what's remarkable is that one of these areas is basic calculus while the other is almost impenetrable analytically. From a numerical standpoint, they are practically the same problem.

```{code-cell}
x = linspace(0, 1, 300)
subplot(1, 2, 1)
plot(x, exp(x))
ylim([0, 2.7]), title("exp(x)")
subplot(1, 2, 2)
plot(x, exp(sin(x)))
ylim([0, 2.7]), title("exp(sin(x))");
```
``````

(demo-int-trap-python)=
``````{dropdown} @demo-int-trap
:open:
We will approximate the integral of the function $f(x)=e^{\sin 7x}$ over the interval $[0,2]$.

```{code-cell}
f = lambda x: exp(sin(7 * x))
a, b = 0, 2
```

In lieu of the exact value, we will use the `quad` function to find an accurate result.

```{code-cell}
from scipy.integrate import quad
I, errest = quad(f, a, b, epsabs=1e-13, epsrel=1e-13)
print(f"Integral = {I:.14f}")
```

Here is the trapezoid result at $n=40$, and its error.

```{code-cell}
T, t, y = FNC.trapezoid(f, a, b, 40)
print(f"Trapezoid estimate is {T:.14f} with error {I - T:.2e}")
```

In order to check the order of accuracy, we increase $n$ by orders of magnitude and observe how the error decreases.

```{code-cell}
n_ = 40 * 2 ** arange(6)
err = zeros(size(n_))
print("     n     error")
for k, n in enumerate(n_):
    T, t, y = FNC.trapezoid(f, a, b, n)
    err[k] = I - T
    print(f"{n:6d}   {err[k]:8.3e} ")
```

Each increase by a factor of 10 in $n$ cuts the error by a factor of about 100, which is consistent with second-order convergence. Another check is that a log-log graph should give a line of slope $-2$ as $n\to\infty$.

```{code-cell}
loglog(n_, abs(err), "-o", label="results")
loglog(n_, 3e-3 * (n_ / n_[0]) ** (-2), "--", label="2nd order")
gca().invert_xaxis()
xlabel("$n$")
ylabel("error")
legend()
title("Convergence of trapezoidal integration");
``````

(demo-int-extrap-python)=
``````{dropdown} @demo-int-extrap
:open:
We estimate $\displaystyle\int_0^2 x^2 e^{-2x}\, dx$ using extrapolation. First we use `quadgk` to get an accurate value.

```{code-cell}
from scipy.integrate import quad
f = lambda x: x**2 * exp(-2 * x)
a = 0
b = 2
I, errest = quad(f, a, b, epsabs=1e-13, epsrel=1e-13)
print(f"Integral = {I:.14f}")
```

We start with the trapezoid formula on $n=N$ nodes.

```{code-cell}
N = 20    # the coarsest formula
n = N
h = (b - a) / n
t = h * arange(n + 1)
y = f(t)
```

We can now apply weights to get the estimate $T_f(N)$.

```{code-cell}
T = zeros(3)
T[0] = h * (sum(y[1:-1]) + y[0] / 2 + y[-1] / 2)
print(f"error (2nd order): {I - T[0]:.2e}")
```

Now we double to $n=2N$, but we only need to evaluate $f$ at every other interior node and apply {eq}`nc-doubling`.

```{code-cell}
n = 2 * n
h = h / 2
t = h * arange(n + 1)
T[1] = T[0] / 2 + h * sum(f(t[1:-1:2]))
print("error (2nd order):", I - T[:2])
```

As expected for a second-order estimate, the error went down by a factor of about 4. We can repeat the same code to double $n$ again.

```{code-cell}
n = 2 * n
h = h / 2
t = h * arange(n + 1)
T[2] = T[1] / 2 + h * sum(f(t[1:-1:2]))
print("error (2nd order):", I - T[:3])
```

Let us now do the first level of extrapolation to get results from Simpson's formula. We combine the elements `T[i]` and `T[i+1]` the same way for $i=1$ and $i=2$.

```{code-cell}
S = array([(4 * T[i + 1] - T[i]) / 3 for i in range(2)])
print("error (4th order):", I - S)
```

With the two Simpson values $S_f(N)$ and $S_f(2N)$ in hand, we can do one more level of extrapolation to get a sixth-order accurate result.

```{code-cell}
R = (16 * S[1] - S[0]) / 15
print("error (6th order):", I - R)
```

We can make a triangular table of the errors:

```{code-cell}
err = nan * ones((3, 3))
err[0, :] = I - T
err[1, 1:] = I - S
err[2, 2] = I - R
results = PrettyTable(["2nd order", "4th order", "6th order"])
results.add_rows(err.T)
print(results)
```

If we consider the computational time to be dominated by evaluations of $f$, then we have obtained a result with about twice as many accurate digits as the best trapezoid result, at virtually no extra cost.
``````

### 5.7 @section-localapprox-adaptive

(demo-adapt-motive-python)=
``````{dropdown} @demo-adapt-motive
:open:
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
``````

(demo-adapt-usage-python)=
``````{dropdown} @demo-adapt-usage
:open:
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
``````