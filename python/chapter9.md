---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
# Chapter 9

## Functions

(function-polyinterp-python)=
``````{dropdown} Barycentric polynomial interpolation
:open:
```{literalinclude} fncbook/fncbook/chapter09.py
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

(function-triginterp-python)=
``````{dropdown} Trigonometric interpolation
:open:
```{literalinclude} fncbook/fncbook/chapter09.py
:filename: triginterp.py
:start-at: def triginterp
:end-at: return np.vectorize(p)
:language: python
:linenos: true
```

````{admonition} About the code
:class: dropdown
The construct on line 13 is known as a *ternary operator*. It is a shorthand for an `if`–`else` statement, giving two alternative results for the true/false cases. Line 19 uses `eachindex(y)`, which generalizes `1:length(y)` to cases where a vector might have a more exotic form of indexing.
````
``````

(function-ccint-python)=
``````{dropdown} Clenshaw–Curtis integration
:open:
```{literalinclude} fncbook/fncbook/chapter09.py
:filename: ccint.py
:start-at: def ccint
:end-at: return I, x
:language: python
:linenos: true
```
``````

(function-glint-python)=
``````{dropdown} Gauss–Legendre integration
:open:
```{literalinclude} fncbook/fncbook/chapter09.py
:filename: glint.py
:start-at: def glint
:end-at: return I, x
:language: python
:linenos: true
```
``````

(function-intinf-python)=
``````{dropdown} Integration over $(-\infty,\infty)$
:open:
```{literalinclude} fncbook/fncbook/chapter09.py
:filename: intinf.py
:start-at: def intinf
:end-at: return I, x
:language: python
:linenos: true
```

::::{admonition} About the code
:class: dropdown
The test `isinf(x(M))` in line 17 checks whether $x(M)$ is larger than the maximum double-precision value, causing it to *overflow* to `Inf`.
::::
``````

(function-intsing-python)=
``````{dropdown} Integration with endpoint singularities
:open:
```{literalinclude} fncbook/fncbook/chapter09.py
:filename: intsing.py
:start-at: def intsing
:end-at: return I, x
:language: python
:linenos: true
```

::::{admonition} About the code
:class: dropdown
The test `iszero(x(M))` in line 17 checks whether $x(M)$ is less than the smallest positive double-precision value, causing it to *underflow* to zero.
::::
``````

## Examples

```{code-cell}
:tags: remove-cell
exec(open("FNC_init.py").read())
```

### 9.1 @section-globalapprox-polynomial

(demo-polynomial-lagrange-python)=
``````{dropdown} @demo-polynomial-lagrange
Here is a vector of nodes.

```{code-cell}
t = array([1, 1.5, 2, 2.25, 2.75, 3])
n = 5
```

Let's apply the definition of the cardinal Lagrange polynomial for $k=2$. First we define a polynomial $q$ that is zero at all the nodes except $i=k$. Then $\ell_2$ is found by normalizing $q$ by $q(t_k)$.

```{code-cell}
k = 2
q = lambda x: prod([x - t[i] for i in range(n + 1) if i != k])
ell_k = lambda x: q(x) / q(t[k])
```

A plot confirms the cardinal property of the result.

```{code-cell}
x = linspace(1, 3, 500)
plot(x, [ell_k(xx) for xx in x])
y = zeros(n+1)
y[k] = 1
plot(t, y, "ko")
xlabel("$x$"),  ylabel("$\\ell_2(x)$")
title(("Lagrange cardinal function"));
```

Observe that $\ell_k$ is _not_ between zero and one everywhere, unlike a hat function.
``````

(demo-polynomial-error-python)=
``````{dropdown} @demo-polynomial-error

```{code-cell}
from scipy.interpolate import BarycentricInterpolator as interp
t = array([1, 1.6, 1.9, 2.7, 3])
p = interp(t, log(t))
```

```{code-cell}
from scipy.interpolate import BarycentricInterpolator as interp
t = array([1, 1.6, 1.9, 2.7, 3])
p = interp(t, log(t))
x = linspace(1, 3, 500)
Phi = lambda x: prod([x - ti for ti in t])
plot(x, [Phi(xj) / 5 for xj in x], label="$\\frac{1}{5}|\\Phi(x)|$")
plot(x, abs(log(x) - p(x)), label="$|f(x)-p(x)|$")
plot(t, zeros(t.size), "ko", label="nodes")
xlabel("$x$"),  ylabel("error")
title("Interpolation error and upper bound"),  legend();
```

The error is zero at the nodes, by the definition of interpolation. The error bound, as well as the error itself, has one local maximum between each consecutive pair of nodes.
``````

### 9.2 @section-globalapprox-barycentric

(demo-barycentric-example-python)=
``````{dropdown} @demo-barycentric-example
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
``````

### 9.3 @section-globalapprox-stability

(demo-stability-equispaced-python)=
``````{dropdown} @demo-stability-equispaced
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
``````

(demo-stability-errfun-python)=
``````{dropdown} @demo-stability-errfun
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
``````

(demo-stability-runge-python)=
``````{dropdown} @demo-stability-runge
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
``````

(demo-stability-errcheb-python)=
``````{dropdown} @demo-stability-errcheb
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
``````

(demo-stability-rungefix-python)=
``````{dropdown} @demo-stability-rungefix
Here again is the function from {numref}`Demo {number} <demo-stability-runge>` that provoked the Runge phenomenon when using equispaced nodes.

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
``````

(demo-stability-spectral-python)=
``````{dropdown} @demo-stability-spectral
On the left, we use a log-log scale, which makes second-order algebraic convergence $O(n^{-4})$ a straight line. On the right, we use a log-linear scale, which makes spectral convergence $O(K^{-n})$ linear.

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
``````

### 9.4 @section-globalapprox-orthogonal

(demo-orthogonal-approx-python)=
``````{dropdown} @demo-orthogonal-approx
Let's approximate $e^x$ over the interval $[−1,1]$. We can sample it at, say, 15 points, and find the best-fitting straight line to that data.

```{code-cell}
from numpy.linalg import lstsq
t = linspace(-1, 1, 15)
y = exp(t)
plot(t, y, label="function")

V = [[ti**j for j in range(2)] for ti in t]
c = lstsq(V, y, rcond=None)[0]
print("fit coeffs:", c)

x = linspace(-1, 1, 600)
plot(x, c[1] + c[0] * x, label="fit")
xlabel("x"),  ylabel("value")
legend(),  title("Least squares fit of exp(x)");
```

There's nothing special about 15 points. Choosing more doesn't change the result much.

```{code-cell}
t = linspace(-1, 1, 150)
y = exp(t)
plot(t, y, label="function")

V = [[ti**j for j in range(2)] for ti in t]
c = lstsq(V, y, rcond=None)[0]
print("fit coeffs:", c)

x = linspace(-1, 1, 600)
plot(x, c[1] + c[0] * x, label="fit")
xlabel("x"),  ylabel("value")
legend(),  title("Least squares fit of exp(x)");
```

This situation is unlike interpolation, where the degree of the interpolant increases with the number of nodes. Here, the linear fit is apparently approaching a limit that we may think of as a continuous least-squares fit.

```{code-cell}
n = arange(40, 420, 60)
results = PrettyTable(["n", "intercept", "slope"])
slope = zeros(n.size)
intercept = zeros(n.size)

for k in range(n.size):
    t = linspace(-1, 1, n[k])
    y = exp(t)
    V = [[ti**j for j in range(2)] for ti in t]
    c = lstsq(V, y, rcond=None)[0]
    results.add_row([n[k], c[1], c[0]])

print(results)
```
``````

### 9.5 @section-globalapprox-trig

(demo-trig-interp-python)=
``````{dropdown} @demo-trig-interp
We will get a cardinal function without using an explicit formula, just by passing data that is 1 at one node and 0 at the others.
```{tip}
:class: dropdown
The operator `÷`, typed as `\div` then <kbd>Tab</kbd>, returns the quotient without remainder of two integers.
```

```{code-cell}
N = 7
n = int((N - 1) / 2)
t = 2 * arange(-n, n + 1) / N
y = zeros(N)
y[n] = 1

p = FNC.triginterp(t, y)
x = linspace(-1, 1, 600)
plot(x, p(x))
plot(t, y, "ko")

xlabel("x"),  ylabel("tau(x)")
title("Trig cardinal function");
```

Here is a 2-periodic function and one of its interpolants.

```{code-cell}
f = lambda x: exp(sin(pi * x) - 2 * cos(pi * x))

plot(x, f(x), label="periodic function")
y = f(t)

p = FNC.triginterp(t, y)
plot(x, p(x), label="trig interpolant")
plot(t, y, "ko", label="nodes")

xlabel("$x$"),  ylabel("$p(x)$")
legend(),  title("Trig interpolation");
```

The convergence of the interpolant is spectral. We let $N$ go needlessly large here in order to demonstrate that unlike polynomials, trigonometric interpolation is stable on equally spaced nodes. Note that when $N$ is even, the value of $n$ is not an integer but works fine for defining the nodes.

```{code-cell}
:tags: hide-input
N = arange(2, 62, 2)
err = zeros(N.size)

x = linspace(-1, 1, 1601)    # for measuring error
for k in range(N.size):
    n = (N[k] - 1) / 2
    t = 2 * arange(-n, n + 1) / N[k]
    p = FNC.triginterp(t, f(t))
    err[k] = max(abs(f(x) - p(x)))

semilogy(N, err, "-o")
xlabel("N"),  ylabel("max error")
title("Convergence of trig interpolation");
```
``````

(demo-trig-fft-python)=
``````{dropdown} @demo-trig-fft
This function has frequency content at $2\pi$, $-2\pi$, and $\pi$. 

```{code-cell}
f = lambda x: 3 * cos(2 * pi * x) - exp(1j * pi * x)
```

To use `fft`, we set up nodes in the interval $[0,2)$. 

```{code-cell}
n = 4
N = 2 * n + 1
t = 2 * arange(0, N) / N    # nodes in [0,2)
y = f(t)
```

We perform Fourier analysis using `fft` and then examine the resulting coefficients.

```{code-cell}
from scipy.fftpack import fft, ifft, fftshift
c = fft(y) / N
freq = hstack([arange(n+1), arange(-n,0)])
results = PrettyTable()
results.add_column("freq", freq) 
results.add_column("coefficient", c)
results
```

Note that $1.5 e^{2i\pi x}+1.5 e^{-2i\pi x} = 3 \cos(2\pi x)$, so this result is sensible.

Fourier's greatest contribution to mathematics was to point out that *every* periodic function is just a combination of frequencies—infinitely many of them in general, but truncated for computational use. Here we look at the magnitudes of the coefficients for $f(x) = \exp( \sin(\pi x) )$.

```{code-cell}
:tags: hide-input
f = lambda x: exp(sin(pi * x))    # content at all frequencies
n = 9;  N = 2*n + 1;
t = 2 * arange(0, N) / N    # nodes in [0,2)
c = fft(f(t)) / N

semilogy(range(-n, n+1), abs(fftshift(c)), "o")
xlabel("$k$"),  ylabel("$|c_k|$")
title("Fourier coefficients");
```

The Fourier coefficients of smooth functions decay exponentially in magnitude as a function of the frequency. This decay rate is determines the convergence of the interpolation error.
``````

### 9.6 @section-globalapprox-integration

(demo-integration-ellipse-python)=
``````{dropdown} @demo-integration-ellipse
```{code-cell}
f = lambda t: pi * sqrt(cos(pi * t) ** 2 + sin(pi * t) ** 2 / 4)
N = arange(4, 48, 6)
results = PrettyTable(["N", "perimeter estimate"])
for k in range(N.size):
    h = 2 / N[k]
    t = h * arange(N[k]) - 1
    results.add_row([N[k], h * sum(f(t))])

results
```
The approximations gain about one digit of accuracy for each constant increment of $n$, which is consistent with spectral convergence.
``````

(demo-integration-compare-python)=
``````{dropdown} @demo-integration-compare
First consider the integral 

$$
\int_{-1}^1 \frac{1}{1+4x^2} \, dx = \arctan(2).
$$

```{code-cell}
f = lambda x: 1 / (1 + 4 * x**2)
exact = arctan(2)
```

We compare the two spectral integration methods for a range of $n$ values.

```{code-cell}
:tags: hide-input
N = range(8, 100, 4)
errCC = zeros(len(N))
errGL = zeros(len(N))
for k, n in enumerate(N):
    errCC[k] = exact - FNC.ccint(f, n)[0]
    errGL[k] = exact - FNC.glint(f, n)[0]

semilogy(N, abs(errCC), "-o", label="Clenshaw–Curtis")
semilogy(N, abs(errGL), "-o", label="Gauss–Legendre")
xlabel("number of nodes"),  ylabel("error"),  ylim(1e-16, 0.01)
legend(),  title("Spectral integration");
```

(The missing dots are where the error is exactly zero.) Gauss–Legendre does converge faster here, but at something less than twice the rate.

Now we try a more sharply peaked integrand:
 
 $$\int_{-1}^1 \frac{1}{1+16x^2} \, dx = \frac{1}{2}\arctan(4).$$ 

```{code-cell}
f = lambda x: 1 / (1 + 16 * x**2)
exact = atan(4) / 2
```

```{code-cell}
:tags: hide-input
N = range(8, 100, 4)
errCC = zeros(len(N))
errGL = zeros(len(N))
for k, n in enumerate(N):
    errCC[k] = exact - FNC.ccint(f, n)[0]
    errGL[k] = exact - FNC.glint(f, n)[0]

semilogy(N, abs(errCC), "-o", label="Clenshaw–Curtis")
semilogy(N, abs(errGL), "-o", label="Gauss–Legendre")
xlabel("number of nodes"),  ylabel("error"),  ylim(1e-16, 0.1)
legend(),  title("Spectral integration");
```

The two are very close until about $n=40$, when the Clenshaw–Curtis method slows down.

Now let's compare the spectral performance to that of our earlier adaptive method in `intadapt`. We will specify varying error tolerances and record the error as well as the total number of evaluations of $f$.

```{code-cell}
:tags: hide-input
loglog(N, abs(errCC), "-o", label="ccint")
loglog(N, abs(errGL), "-o", label="glint")

tol_ = 1 / 10 ** arange(2, 15)
n = zeros(tol_.size)
errAdapt = zeros(tol_.size)
for k, tol in enumerate(tol_):
    Q, t = FNC.intadapt(f, -1, 1, tol)
    errAdapt[k] = exact - Q
    n[k] = t.size

loglog(n, abs(errAdapt), "-o", label="intadapt")
loglog(n, 1 / (n**4), "--", label="4th order")
xlabel("number of nodes"),  ylabel("error"),  ylim(1e-16, 1)
legend(),  title("Spectral vs 4th order");
```

At the core of `intadapt` is a fourth-order formula, and the results track that rate closely. For all but the most relaxed error tolerances, both spectral methods are far more efficient than the low-order counterpart. For other integrands, particularly those that vary nonuniformly across the interval, the adaptive method might be more competitive.

``````

### 9.7 @section-globalapprox-improper

(demo-improper-decay-python)=
``````{dropdown} @demo-improper-decay
```{code-cell}
f = lambda x: 1 / (1 + x**2)
x = linspace(-4, 4, 500)
subplot(2, 1, 1)
plot(x, f(x)),  yscale('log')
xlabel('x'),  ylabel('f(x)'),  ylim([1e-16, 1])  
title('Original integrand')   

xi = lambda t: sinh( pi * sinh(t) / 2 )
dxi_dt = lambda t: pi/2 * cosh(t) * cosh( pi * sinh(t) / 2 )
integrand = lambda t: f(xi(t)) * dxi_dt(t)
subplot(2, 1, 2)
plot(x, integrand(x)),  yscale('log')
xlabel('t'),  ylabel('f(x(t))'),  ylim([1e-16, 1])  
title('Transformed integrand')   
```

This graph suggests that we capture all of the integrand values that are larger than machine epsilon by integrating in $t$ from $-4$ to $4$.
``````

(demo-improper-intinf-python)=
``````{dropdown} @demo-improper-intinf
```{code-cell}
:tags: hide-input
f = lambda x: 1 / (1 + x**2)
exact = pi
tol = array([1 / 10**d for d in arange(5, 14, 0.5)])
err = zeros((tol.size, 2))
length = zeros((tol.size, 2))
for k in range(tol.size):
    I1, x1 = FNC.intadapt(f, -2/tol[k], 2/tol[k], tol[k])
    I2, x2 = FNC.intinf(f, tol[k])
    err[k] = abs(exact - array([I1, I2]))
    length[k] = [x1.size, x2.size]
loglog(length, err, "-o")
# plot(len,err,m=:o,label=["direct" "double exponential"])
n = array([100, 10000])
loglog(n, 1000 / n**4, 'k--')
xlabel("number of nodes"),  ylabel("error")
title("Comparison of integration methods")
legend(["direct", "double exponential", "4th-order"], loc="lower left");
```

Both methods are roughly fourth-order due to Simpson's formula in the underlying adaptive integration method. At equal numbers of evaluation nodes, however, the double exponential method is consistently 2--3 orders of magnitude more accurate.
``````

(demo-improper-intsing-python)=
``````{dropdown} @demo-improper-intsing

```{code-cell}
:tags: hide-input
f = lambda x: 1 / (10 * sqrt(x))
exact = 0.2
tol = array([1 / 10**d for d in arange(5, 14, 0.5)])
err = zeros((tol.size, 2))
length = zeros((tol.size, 2))
for k in range(tol.size):
    I1, x1 = FNC.intadapt(f, (tol[k]/20)**2, 1, tol[k])
    I2, x2 = FNC.intsing(f, tol[k])
    err[k] = abs(exact - array([I1, I2]))
    length[k] = [x1.size, x2.size]
loglog(length, err, "-o")
# plot(len,err,m=:o,label=["direct" "double exponential"])
n = array([100, 10000])
loglog(n, 30 / n**4, 'k--')
xlabel("number of nodes"),  ylabel("error")
title("Comparison of integration methods")
legend(["direct", "double exponential", "4th-order"], loc="lower left");
```

As in {numref}`Demo {number} <demo-improper-intinf>`, the double exponential method is more accurate than direct integration by a few orders of magnitude. Equivalently, the same accuracy can be reached with many fewer nodes.
``````
