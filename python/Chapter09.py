# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: myst
#     language: python
#     name: python3
# ---

# +
from scipy import *
from numpy import *
from matplotlib.pyplot import *
from scipy.linalg import *
from numpy.linalg import *
from scipy.fftpack import fft, ifft, fftshift
import FNC

# -

# # Example 9.1.1

# We plot a cardinal Lagrange polynomial for $n=5$ and $k=2$.

# +
t = array([1, 1.5, 2, 2.25, 2.75, 3])
n = 5
k = 2

phi = lambda x: prod([x - t[i] for i in range(n + 1) if i != k])
ell_k = lambda x: phi(x) / phi(t[k])
# -

x = linspace(1, 3, 500)
plot(x, [ell_k(xx) for xx in x])
y = zeros(n + 1)
y[k] = 1
plot(t, y, "o", color="black")
xlabel("$x$")
ylabel("$\ell_2(x)$")
title("Lagrange cardinal function")

# Observe that $\ell_k$ is _not_ between zero and one everywhere between the nodes.

# # Example 9.1.3

# Consider the problem of interpolating $\log(x)$ at these nodes:

t = array([1, 1.6, 1.9, 2.7, 3])

# Here $n=4$ and $f^{(5)}(\xi) = 4!/\xi^5$. For $\xi\in[1,3]$ we can say that $|f^{(5)}(\xi)| \le 4!$. Hence
#
# $$ |f(x)-p(x)| \le \frac{1}{5} \Phi(x).$$

Phi = lambda x: prod([x - ti for ti in t])
plot(x, [Phi(xj) / 5 for xj in x])
plot(t, zeros(t.size), "o", color="black")
xlabel("$x$")
ylabel("$\Phi(x)/5$")
title("Interpolation error function")

# The error bound has one local extreme point between each consecutive pair of nodes.

# # Example 9.2.2

# We use $n=3$ and $n=6$ with equally spaced nodes for the function $\sin(e^{2x})$ over $[0,1]$.

f = lambda x: sin(exp(2 * x))
x = linspace(0, 1, 500)
fig, ax = subplots()
ax.plot(x, f(x), label="function")

# +
t = arange(4) / 3
y = f(t)
p = FNC.polyinterp(t, y)

ax.plot(t, y, "o", color="black", label="nodes")
ax.plot(x, p(x), label="interpolant")
ax.legend()
ax.set_title("Interpolation on 4 nodes")
fig
# -

plot(x, f(x), label="function")
t = arange(7) / 6
y = f(t)
p = FNC.polyinterp(t, y)
plot(t, y, "o", color="black", label="nodes")
plot(x, p(x), label="interpolant")
legend()
title("Interpolation on 7 nodes")

# The curves always intersect at the interpolation nodes.

# # Example 9.3.1

# We choose a function over the interval $[0,1]$.

f = lambda x: sin(exp(2 * x))
x = linspace(0, 1, 500)
plot(x, f(x))

# We interpolate it at equally spaced nodes for increasing values of $n$. We will sample the interpolant at a large number of points in order to estimate the interpolation error.

# +
N = arange(5, 65, 5)
err = zeros(N.size)
x = linspace(0, 1, 1001)  # for measuring error
for k, n in enumerate(N):
    t = linspace(0, 1, n + 1)  # equally spaced nodes
    y = f(t)  # interpolation data
    p = FNC.polyinterp(t, y)
    err[k] = max(abs(f(x) - p(x)))

semilogy(N, err, "-o")
xlabel("$n$")
ylabel("max error")
title("Polynomial interpolation error")
# -

# Initially the error decreases exponentially, i.e. as $O(K^{-n})$ for some $K>1$. However, around $n=30$ the error starts to _grow_ exponentially.

# # Example 9.3.2

# We plot $|\Phi(x)|$ over the interval $[-1,1]$ with equispaced nodes for different values of $n$.

x = linspace(-1, 1, 1601)
for n in range(10, 60, 10):
    t = linspace(-1, 1, n + 1)
    Phi = array([prod(xk - t) for xk in x])
    semilogy(x, abs(Phi), ".", markersize=2)
xlabel("$x$")
ylabel("$|\Phi(x)|$")
ylim([1e-25, 1])
title("Effect of equispaced nodes")

# (Each time $\Phi$ passes through zero at an interpolation node, the value on the log scale should go to $-\infty$, which explains the numerous cusps on the curves.) Two observations are important: First, the size of $|\Phi|$ decreases exponentially at each fixed location in the interval (because the spacing between curves is constant for constant increments of $n$). Second, $|\Phi|$ is larger at the ends of the interval than in the middle, by an exponentially growing factor.

# # Example 9.3.3

# This function has infinitely many continuous derivatives on the entire real line and looks very easy to approximate over $[-1,1]$.

f = lambda x: 1 / (x**2 + 16)
x = linspace(-1, 1, 1601)
plot(x, f(x))
title("Test function")

# We start by doing polynomial interpolation for some rather small values of $n$.

# +
N = arange(4, 16, 4)
label = []
for k, n in enumerate(N):
    t = linspace(-1, 1, n + 1)  # equally spaced nodes
    y = f(t)  # interpolation data
    p = FNC.polyinterp(t, y)
    err = abs(f(x) - p(x))
    semilogy(x, err, ".", markersize=2)
    label.append(f"degree {n}")

xlabel("$x$")
ylabel("$|f(x)-p(x)|$")
ylim([1e-20, 1])
legend(label)
title("Error for low degrees")
# -

# The convergence so far appears rather good, though not uniformly so. Now watch what happens as we continue to increase the degree.

N = 12 + 15 * arange(1, 4)
label = []
for k, n in enumerate(N):
    t = linspace(-1, 1, n + 1)  # equally spaced nodes
    y = f(t)  # interpolation data
    p = FNC.polyinterp(t, y)
    err = abs(f(x) - p(x))
    semilogy(x, err, ".", markersize=2)
    label.append(f"degree {n}")
xlabel("$x$")
ylabel("$|f(x)-p(x)|$")
ylim([1e-20, 1])
legend(label)
title("Error for higher degrees")

# The convergence in the middle can't get any better than machine precision. So maintaining the growing gap between the center and the ends pushes the error curves upward exponentially fast at the ends, wrecking the convergence.

# # Example 9.3.4

# +
x = linspace(-1, 1, 1601)
for n in range(10, 60, 10):
    theta = pi * arange(n + 1) / n
    t = -cos(theta)
    Phi = array([prod(xk - t) for xk in x])
    semilogy(x, abs(Phi), ".")

xlabel("$x$")
ylabel("$|\Phi(x)|$")
ylim([1e-18, 1e-2])
title("Effect of Chebyshev nodes")
# -

# The convergence is a bit slower in the middle than with equally spaced points, but far more uniform over the entire interval, which is the key to global convergence.

# # Example 9.3.5

# This function has infinitely many continuous derivatives on the entire real line and looks very easy to approximate over $[-1,1]$.

# +
f = lambda x: 1 / (x**2 + 16)
x = linspace(-1, 1, 1601)

label = []
for k, n in enumerate([4, 10, 16, 40]):
    t = -cos(pi * arange(n + 1) / n)
    y = f(t)  # interpolation data
    p = FNC.polyinterp(t, y)
    err = abs(f(x) - p(x))
    semilogy(x, err, ".", markersize=2)
    label.append(f"degree {n}")

xlabel("$x$")
ylabel("$|f(x)-p(x)|$")
ylim([1e-20, 1])
legend(label)
title("Error for Chebyshev interpolants")
# -

# By degree 16 the error is uniformly within machine epsilon. Note that even as the degree continues to increase, the error near the ends does not grow as with the Runge phenomenon for equally spaced nodes.

# # Example 9.4.1

# Let's approximate $e^x$ over the interval $[-1,1]$. We can sample it at, say, 40 points, and find the best-fitting straight line to that data.

# +
t = linspace(-1, 1, 40)
y = exp(t)
plot(t, y, "-o", label="function")

V = [[ti**j for j in range(2)] for ti in t]
c = lstsq(V, y, rcond=None)[0]
print("fit coeffs:", c)

x = linspace(-1, 1, 600)
plot(x, c[1] + c[0] * x, label="fit")
xlabel("x")
ylabel("value")
legend()
title("Least squares fit of exp(x)")
# -

# There's nothing special about 40 points. By choosing more we get closer to the true function $e^x$.

# +
t = linspace(-1, 1, 200)
y = exp(t)
plot(t, y, label="function")

V = [[ti**j for j in range(2)] for ti in t]
c = lstsq(V, y, rcond=None)[0]
print("fit coeffs:", c)

x = linspace(-1, 1, 600)
plot(x, c[1] + c[0] * x, label="fit")
xlabel("x")
ylabel("value")
legend()
title("Least squares fit of exp(x)")

# +
t = linspace(-1, 1, 1000)
y = exp(t)
plot(t, y, label="function")

V = [[ti**j for j in range(2)] for ti in t]
c = lstsq(V, y, rcond=None)[0]
print("fit coeffs:", c)

x = linspace(-1, 1, 600)
plot(x, c[1] + c[0] * x, label="fit")
xlabel("x")
ylabel("value")
legend()
title("Least squares fit of exp(x)")
# -

# It's quite plausible that the coefficients of the best-fit line are approaching a limit as the number of nodes goes to infinity.

# # Example 9.5.1

# We get a cardinal function if we use data that is one at a node and zero at the others.

# +
N = 7
n = int((N - 1) / 2)
t = 2 * arange(-n, n + 1) / N
y = zeros(N)
y[n] = 1

p = FNC.triginterp(t, y)
x = linspace(-1, 1, 600)
plot(x, p(x))
plot(t, y, "o")

xlabel("x")
ylabel("tau(x)")
title("Trig cardinal function")
# -

# Here is a 2-periodic function and one of its interpolants.

# +
f = lambda x: exp(sin(pi * x) - 2 * cos(pi * x))

plot(x, f(x), label="function")
y = f(t)
plot(t, y, "o", label="nodes")

p = FNC.triginterp(t, y)
plot(x, p(x), label="interpolant")

xlabel("$x$")
ylabel("$p(x)$")
legend()
title("Trig interpolation")
# -

# The convergence of the interpolant is exponential (spectral). We let $N$ go needlessly large here in order to demonstrate that unlike polynomials, trigonometric interpolation is stable on equally spaced nodes.

# +
N_ = range(3, 93, 3)
err = zeros(len(N_))

x = linspace(-1, 1, 1601)  # for measuring error
for k, N in enumerate(N_):
    n = (N - 1) / 2
    t = 2 * arange(-n, n + 1) / N
    p = FNC.triginterp(t, f(t))
    err[k] = max(abs(f(x) - p(x)))

semilogy(N_, err, "-o")
xlabel("N")
ylabel("max error")
title("Error in trig interpolation")
# -

# # Example 9.5.2

# This function has two distinct frequencies.

f = lambda x: 3 * cos(5 * pi * x) - exp(2j * pi * x)

# We set up to use `fft` from `scipy.fftpack`. Note how the definition of the nodes has changed.

n = 10
N = 2 * n + 1
t = 2 * arange(0, N) / N  # nodes in [0,2)
y = f(t)

# We perform Fourier analysis using `fft` and then examine the coefficients.

c = fft(y) / N
print(c)

cs = fftshift(c)  # put zero frequency in center
plot(range(-n, n + 1), real(cs), "o")
xlabel("$k$")
ylabel("$c_k$")
title("Interpolant coefficients")

# Note that $1.5 e^{5i\pi x}+1.5 e^{-5i\pi x} = 3 \cos(5\pi x)$ by Euler's formula, so this result is sensible.

# Fourier's greatest contribution to mathematics was to point out that every periodic function is just a combination of frequencies---infinitely many of them in general, but truncated for computational use.

# +
f = lambda x: exp(sin(pi * x))
c = fft(f(t)) / N

semilogy(range(-n, n + 1), abs(fftshift(c)), "o")
xlabel("$k$")
ylabel("$|c_k|$")
title("Fourier coefficients")
# -

# The Fourier coefficients of smooth functions decay exponentially in magnitude as a function of the frequency. This decay rate is directly linked to the convergence of the interpolation error.

# # Example 9.6.1

# We use the trapezoidal integration formula to compute the perimeter of an ellipse with semi-axes 1 and 1/2. Parameterizing the ellipse as $x=\cos \pi t$, $y=\frac{1}{2}\sin \pi t$ leads to the integral
#
# $$\int_{-1}^1 \pi\sqrt{ \cos^2(\pi t) + \tfrac{1}{4}\sin^2(\pi t)}\,dt.$$

# +
f = lambda t: pi * sqrt(cos(pi * t) ** 2 + sin(pi * t) ** 2 / 4)
N_ = arange(4, 64, 4)
C = zeros(N_.size)

for i, N in enumerate(N_):
    h = 2 / N
    t = h * arange(N) - 1
    C[i] = h * sum(f(t))
print("perimeter:", C[-1])
# -

err = abs(C - C[-1])
semilogy(N_, err, "-o")
xlabel("number of nodes")
ylabel("error")
ylim(1e-15, 1)
title("Convergence of perimeter calculation")

# The approximations gain about one digit of accuracy for each constant increase in $N$, consistent with geometric (linear) convergence.

# # Example 9.6.3

# First consider the integral  $$\int_{-1}^1 \frac{1}{1+4x^2} \, dx = \arctan(2).$$

f = lambda x: 1 / (1 + 4 * x**2)
exact = arctan(2)

# We compare the two spectral integration methods for a range of $n$ values.

N = range(8, 100, 4)
errCC = zeros(len(N))
errGL = zeros(len(N))
for k, n in enumerate(N):
    errCC[k] = exact - FNC.ccint(f, n)[0]
    errGL[k] = exact - FNC.glint(f, n)[0]

semilogy(N, abs(errCC), "o", label="CC")
semilogy(N, abs(errGL), "o", label="GL")
xlabel("number of nodes")
ylabel("error")
ylim(1e-16, 1)
legend()
title("Spectral integration")

# (The missing dots are where the error is exactly zero.) Gauss--Legendre does converge faster here, but at something less than twice the rate. Now we try a more sharply peaked integrand:

#  $$\int_{-1}^1 \frac{1}{1+16x^2} \, dx = \frac{1}{2}\arctan(4).$$

f = lambda x: 1 / (1 + 16 * x**2)
exact = arctan(4) / 2

# +
N = range(8, 100, 4)
errCC = zeros(len(N))
errGL = zeros(len(N))
for k, n in enumerate(N):
    errCC[k] = exact - FNC.ccint(f, n)[0]
    errGL[k] = exact - FNC.glint(f, n)[0]

semilogy(N, abs(errCC), "o", label="CC")
semilogy(N, abs(errGL), "o", label="GL")
xlabel("number of nodes")
ylabel("error")
ylim(1e-16, 1)
legend()
title("Spectral integration")
# -

# The two are very close until about $n=40$, when the Clenshaw--Curtis method slows down.

# Now let's compare the spectral performance to that of our earlier adaptive method in `adaptquad`. We will specify varying error tolerances and record the error as well as the total number of evaluations of $f$.

# +
loglog(N, abs(errCC), "o", label="CC")
loglog(N, abs(errGL), "o", label="GL")

tol_ = 1 / 10 ** arange(2, 15)
n = zeros(tol_.size)
errAdapt = zeros(tol_.size)
for k, tol in enumerate(tol_):
    Q, t = FNC.intadapt(f, -1, 1, tol)
    errAdapt[k] = exact - Q
    n[k] = t.size

loglog(n, abs(errAdapt), "-o", label="intadapt")
loglog(n, 1 / (n**4), "--", label="4th order")
xlabel("number of nodes")
ylabel("error")
ylim(1e-16, 1)
legend()
title("Spectral vs 4th order")
# -

# At the core of `intadapt` is a fourth-order formula, and the results track that rate closely. For all but the most relaxed error tolerances, both spectral methods are far more efficient than the low-order counterpart.

# # Example 9.7.2

# +
f = lambda x: 1 / (1 + x**2)
x = lambda t: sinh(pi * sinh(t) / 2)
chain = lambda t: pi / 2 * cosh(t) * cosh(pi * sinh(t) / 2)
integrand = lambda t: f(x(t)) * chain(t)

t = linspace(-4, 4, 500)
semilogy(t, integrand(t))
xlabel("$t$")
ylabel("$f(x(t))$")
title("Doubly exponential integrand decay")
# -

# This graph suggests that we may integrate $t$ from $-4$ to $4$ and capture all of the integrand values that are larger than machine epsilon.

# # Example 9.7.3

# +
f = lambda x: 1 / (1 + x**2)

h = 1 / 10 ** linspace(3, 1, 60)
M = 10 ** linspace(3, 16, 60)

I = array([[FNC.intde(f, hi, Mj)[0] for Mj in M] for hi in h])
err = abs(I - pi)
# -

contour(h, M, -log10(err + 1e-20))
xlabel("$h$")
xscale("log")
xlim(0.1, 0.001)
ylabel("$M$")
yscale("log")
title("Number of accurate digits")
colorbar()

# As predicted, the error can't much beat $1/M$. Even for very large $M$, however, not many nodes are needed, as seen in the very weak dependence of the accuracy on $h$. For instance,

# +
I, x = FNC.intde(f, 0.1, 1e10)
err = abs(pi - I)
number_of_nodes = x.size

xpos = x > 0
semilogx(x[xpos], f(x[xpos]), "-o")
xlabel("$x$")
ylabel("$f(x)$")
ylim(0, 1)
title("Positive nodes used for integration")
# -

# # Example 9.7.4

I, x = FNC.intadapt(sqrt, 0, 1, 1e-10)
print("error:", I - 2 / 3)
print("number of nodes:", x.size)

# The adaptive integrator was reasonably successful. But if we integrate $1/\sqrt{x}$, which in unbounded at the origin, the number of nodes goes up dramatically.

I, x = FNC.intadapt(lambda x: 1 / sqrt(x), 2**-52, 1, 1e-10)
print("error:", I - 2)
print("number of nodes:", x.size)

# The nodes are packed in very closely at the origin; in fact they are placed with exponential spacing.

loglog(x[:-1], diff(x), "o")
xlabel("x")
xlim(1e-15, 1)
ylabel("distance")
ylim(1e-15, 1)
title("Distance between nodes")

# # Example 9.7.5

# We return to the problem of computing $\int_0^1 \sqrt{x}\,dx$. In order to apply `intsing`, we first have to transform the interval of integration to $[-1,1]$. We can do this through $z=2x-1$. Note that

# $$\int_0^1 \sqrt{x}\,dx = \int_{-1}^1 \sqrt{\tfrac{1}{2}(z+1)}\cdot \tfrac{1}{2}dz.$$

f = lambda z: 0.5 * sqrt((1 + z) / 2)
I, z = FNC.intsing(f, 0.1, 1e-12)
print("error:", I - 2 / 3)
print("number of nodes:", z.size)

# The integration required very few nodes. For the more difficult integral of $1/\sqrt{x}$, the results are limited by how accurately we can represent $-1+\delta$.

f = lambda z: 0.5 / sqrt((1 + z) / 2)
I, z = FNC.intsing(f, 0.1, 1e-12)
print("error:", I - 2)
print("number of nodes:", z.size)

# If we make $\delta$ any smaller, the outermost trapezoid nodes will be indistinguishable from $z=\pm 1$, i.e. the exact endpoints of the interval. We would need to use special code to evaluate $f$ indirectly in the limits $t\to \pm \infty$.
