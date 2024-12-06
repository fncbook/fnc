---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.4
kernelspec:
  display_name: myst
  language: python
  name: python3
---

```{code-cell}
from scipy import *
from numpy import *
from matplotlib.pyplot import *
from scipy.linalg import *
from numpy.linalg import *
from scipy.optimize import root_scalar
import scipy.special as special
import FNC
```

```{code-cell}
# This (optional) block is for improving the display of plots.
from IPython.display import set_matplotlib_formats
```

```{code-cell}
set_matplotlib_formats("svg", "pdf")
rcParams["figure.figsize"] = [7, 4]
rcParams["lines.linewidth"] = 2
rcParams["lines.markersize"] = 4
rcParams["animation.html"] = "jshtml"  # or try "html5"
```

 # Example 4.1.1

+++

 In the theory of vibrations of a circular drum, the displacement of the drumhead can be expressed in terms of pure harmonic modes,

 $$J_m(\omega_{k,m} r) \cos(m\theta) \cos(c \omega_{k,m} t),$$

 where $(r,\theta)$ are polar coordinates, $0\le r\le 1$, $t$ is time, $m$ is a positive integer, $c$ is a material parameter, and $J_m$ is a _Bessel function of the first kind_. The quantity $\omega_{k,m}$ is a resonant frequency and is a positive root of the equation

 $$J_m(\omega_{k,m}) = 0,$$

 which states that the drumhead is clamped around the rim. Tabulating approximations to the zeros of Bessel functions has occupied countless mathematician-hours throughout the centuries.

```{code-cell}
def J3(x):
    return special.jv(3.0, x)
```

```{code-cell}
xx = linspace(0, 20, 500)
fig, ax = subplots()
ax.plot(xx, J3(xx))
ax.grid()
xlabel("$x$")
ylabel("$J_3(x)$")
title("Bessel function")
```

 From the graph we see roots near 6, 10, 13, 16, and 19. We use `root_scalar` from the `scipy.optimize` package to find these roots accurately.

```{code-cell}
omega = []
for guess in [6.0, 10.0, 13.0, 16.0, 19.0]:
    s = root_scalar(J3, bracket=[guess - 0.5, guess + 0.5]).root
    omega.append(s)

print(omega)
```

```{code-cell}
ax.plot(omega, J3(omega), "o")
ax.set_title("Bessel function roots")
fig
```

 # Example 4.1.2

+++

 Consider first the function

```{code-cell}
f = lambda x: (x - 1) * (x - 2)
```

 At the root $r=1$, we have $f'(r)=-1$. If the values of $f$ were perturbed at any point by noise of size, say, $0.05$, we can imagine finding the root of the function as though drawn with a thick line, whose edges we show here.

```{code-cell}
xx = linspace(0.8, 1.2, 400)

plot(xx, f(xx))
plot(xx, f(xx) + 0.02, "k")
plot(xx, f(xx) - 0.02, "k")
axis("equal")
grid(True)
xlabel("x")
ylabel("f(x)")
title("Well-conditioned root")
```

 The possible values for a perturbed root all lie within the interval where the black lines intersect the $x$ axis. The width of that zone is about the same as the vertical distance between the lines.

+++

 By contrast, consider the function

```{code-cell}
f = lambda x: (x - 1) * (x - 1.01)
```

 Now $f'(1)=-0.01$, and the graph of $f$ will be much shallower near $x=1$. Look at the effect this has on our thick rendering:

```{code-cell}
plot(xx, f(xx))
plot(xx, f(xx) + 0.02, "k")
plot(xx, f(xx) - 0.02, "k")
axis("equal")
grid(True)
xlabel("x")
ylabel("f(x)")
title("Poorly conditioned root")
```

The vertical displacements in this picture are exactly the same as before. But the potential _horizontal_ displacement of the root is much wider. In fact, if we perturb the function upward by the amount drawn here, the root disappears entirely!

+++

 # Example 4.2.1

+++

 Let's convert the roots of a quadratic polynomial $f(x)$ to a fixed point problem.

```{code-cell}
f = poly1d([1, -4, 3.5])
r = f.roots
print(r)
```

 We'll define $g(x)=x-f(x)$. Intersections of its graph with the line $y=x$ are fixed points of $g$ and thus roots of $f$. (Only one is shown in the chosen plot range.)

```{code-cell}
fig, ax = subplots()
g = lambda x: x - f(x)
xx = linspace(2, 3, 400)
ax.plot(xx, g(xx), label="y=g(x)")
ax.plot(xx, xx, label="y=x")
axis("equal")
legend()
title("Finding a fixed point")
```

 If we evalaute $g(2.1)$, we get a value of almost 2.6.

```{code-cell}
x = 2.1
y = g(x)
print(y)
```

 So $g(x)$ is considerably closer to a fixed point than $x$ was. The value $y=g(x)$ ought to become our new $x$ value! Changing the $x$ coordinate in this way is the same as following a horizontal line over to the graph of $y=x$.

```{code-cell}
ax.plot([x, y], [y, y], "r:", label="")
fig
```

 Now we can compute a new value for $y$. We leave $x$ alone here, so we travel along a vertical line to the graph of $g$.

```{code-cell}
x = y
y = g(x)
print("y:", y)
ax.plot([x, x], [x, y], "k:")
fig
```

 You see that we are in a position to repeat these steps as often as we like. Let's apply them a few times and see the result.

```{code-cell}
for k in range(5):
    ax.plot([x, y], [y, y], "r:")
    x = y
    # y --> new x
    y = g(x)
    ax.plot([x, x], [x, y], "k:")  # g(x) --> new y
fig
```

 The process spirals in beautifully toward the fixed point we seek. Our last estimate has almost 4 accurate digits.

```{code-cell}
abs(y - max(r)) / max(r)
```

 Now let's try to find the other fixed point $\approx 1.29$ in the same way. We'll use 1.3 as a starting approximation.

```{code-cell}
xx = linspace(1, 2, 400)
fig, ax = subplots()
ax.plot(xx, g(xx), label="y=g(x)")
ax.plot(xx, xx, label="y=x")
ax.set_aspect(1.0)
ax.legend()

x = 1.3
y = g(x)
for k in range(5):
    ax.plot([x, y], [y, y], "r:")
    x = y
    # y --> new x
    y = g(x)
    ax.plot([x, x], [x, y], "k:")  # g(x) --> new y
ylim(1, 2.5)
title("No convergence")
```

 This time, the iteration is pushing us _away_ from the correct answer.

+++

 # Example 4.2.3

```{code-cell}
f = poly1d([1, -4, 3.5])
r = f.roots
print(r)
```

 Here is the fixed point iteration. This time we keep track of the whole sequence of approximations.

```{code-cell}
g = lambda x: x - f(x)
x = zeros(12)
x[0] = 2.1
for k in range(11):
    x[k + 1] = g(x[k])

print(x)
```

 It's easiest to construct and plot the sequence of errors.

```{code-cell}
err = abs(x - max(r))
semilogy(err, "-o")
xlabel("iteration number")
ylabel("error")
title("Convergence of fixed point iteration")
```

 It's quite clear that the convergence quickly settles into a linear rate. We could estimate this rate by doing a least-squares fit to a straight line. Keep in mind that the values for small $k$ should be left out of the computation, as they don't represent the linear trend.

```{code-cell}
p = polyfit(arange(5, 13), log(err[4:]), 1)
print(p)
```

 We can exponentiate the slope to get the convergence constant $\sigma$.

```{code-cell}
print("sigma:", exp(p[0]))
```

 The numerical values of the error should decrease by a factor of $\sigma$ at each iteration. We can check this easily with an elementwise division.

```{code-cell}
err[8:] / err[7:-1]
```

 # Example 4.3.1

```{code-cell}
f = lambda x: x * exp(x) - 2
xx = linspace(0, 1.5, 400)

fig, ax = subplots()
ax.plot(xx, f(xx), label="function")
ax.grid()
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
```

From the graph, it is clear that there is a root near $x=1$. So we call that our initial guess, $x_1$.

```{code-cell}
x1 = 1
f1 = f(x1)
ax.plot(x1, f1, "ko", label="initial point")
ax.legend()
fig
```

Next, we can compute the tangent line at the point $\bigl(x_1,f(x_1)\bigr)$, using the derivative.

```{code-cell}
dfdx = lambda x: exp(x) * (x + 1)
slope1 = dfdx(x1)
tangent1 = lambda x: f1 + slope1 * (x - x1)

ax.plot(xx, tangent1(xx), "--", label="tangent line")
ax.set_ylim(-2, 4)
ax.legend()
fig
```

In lieu of finding the root of $f$ itself, we settle for finding the root of the tangent line approximation, which is trivial. Call this $x_2$, our next approximation to the root.

```{code-cell}
x2 = x1 - f1 / slope1
ax.plot(x2, 0, "ko", label="tangent root")
ax.legend()
fig
```

```{code-cell}
f2 = f(x2)
print(f2)
```

The residual (value of $f$) is smaller than before, but not zero. So we repeat the process with a new tangent line based on the latest point on the curve.

```{code-cell}
xx = linspace(0.83, 0.88, 200)

plot(xx, f(xx))
plot(x2, f2, "ko")
grid()
xlabel("$x$")
ylabel("$y$")

slope2 = dfdx(x2)
tangent2 = lambda x: f2 + slope2 * (x - x2)
plot(xx, tangent2(xx), "--")
x3 = x2 - f2 / slope2
plot(x3, 0, "ko")
title("Second iteration")
```

```{code-cell}
f3 = f(x3)
print(f3)
```

We appear to be getting closer to the true root each time.

+++

# Example 4.3.2

```{code-cell}
f = lambda x: x * exp(x) - 2
dfdx = lambda x: exp(x) * (x + 1)
```

We don't know the exact root, so we use `root_scalar` to determine the "true" value.

```{code-cell}
r = root_scalar(f, bracket=[0.8, 1.0]).root
print(r)
```

We use $x_1=1$ as a starting guess and apply the iteration in a loop, storing the sequence of iterates in a vector.

```{code-cell}
x = ones(7)
for k in range(6):
    x[k + 1] = x[k] - f(x[k]) / dfdx(x[k])

print(x)
```

Here is the sequence of errors.

```{code-cell}
err = x - r
print(err)
```

Glancing at the exponents of the errors, they roughly form a neat doubling sequence until the error is comparable to machine precision. We can see this more precisely by taking logs.

```{code-cell}
print(log(abs(err)))
```

Quadratic convergence isn't as graphically distinctive as linear convergence.

```{code-cell}
semilogy(range(7), abs(err), "o-")
xlabel("$k$")
ylabel("$|x_k-r|$")
title("Quadratic convergence")
```

# Example 4.3.3

+++

Suppose we want to solve $e^x=x+c$ for multiple values of $c$. We can create functions for $f$ and $f'$ in each case.

```{code-cell}
for c in [2, 4, 7.5, 11]:
    f = lambda x: exp(x) - x - c
    dfdx = lambda x: exp(x) - 1
    x = FNC.newton(f, dfdx, 1.0)
    r = x[-1]
    print(f"root with c = {c} is {r}")
```

There's a subtlety about the definition of `f`. It uses whatever value is assigned to `c` at the moment `f` is called. (This is unlike MATLAB, which locks in the value defined for `c` at the moment of definition.) If we later change the value assigned to `c`, the function is changed also.

```{code-cell}
c = 11
f = lambda x: exp(x) - x - c
print(f(0))
```

```{code-cell}
c = 100
print(f(0))
```

# Example 4.4.1

+++

We return to finding a root of the equation $xe^x=2$.

```{code-cell}
f = lambda x: x * exp(x) - 2
xx = linspace(0.25, 1.25, 400)

fig, ax = subplots()
ax.plot(xx, f(xx), label="function")
ax.set_xlabel("$x$")
ax.set_ylabel("$f(x)$")
ax.grid()
```

From the graph, it's clear that there is a root near $x=1$. To be more precise, there is a root in the interval $[0.5,1]$. So let us take the endpoints of that interval as _two_ initial approximations.

```{code-cell}
x1 = 1
f1 = f(x1)
x2 = 0.5
f2 = f(x2)
ax.plot([x1, x2], [f1, f2], "ko", label="initial points")
ax.legend()
fig
```

Instead of constructing the tangent line by evaluating the derivative, we can construct a linear model function by drawing the line between the two points $\bigl(x_1,f(x_1)\bigr)$ and $\bigl(x_2,f(x_2)\bigr)$. This is called a _secant line_.

```{code-cell}
slope2 = (f2 - f1) / (x2 - x1)
secant2 = lambda x: f2 + slope2 * (x - x2)
ax.plot(xx, secant2(xx), "--", label="secant line")
ax.legend()
fig
```

As before, the next value in the iteration is the root of this linear model.

```{code-cell}
x3 = x2 - f2 / slope2
ax.plot(x3, 0, "o", label="root of secant")
f3 = f(x3)
print(f3)
ax.legend()
fig
```

For the next linear model, we use the line through the two most recent points. The next iterate is the root of that secant line, and so on.

```{code-cell}
slope3 = (f3 - f2) / (x3 - x2)
x4 = x3 - f3 / slope3
print(f(x4))
```

# Example 4.4.2

+++

We check the convergence of the secant method from the previous example.

```{code-cell}
f = lambda x: x * exp(x) - 2
x = FNC.secant(f, 1, 0.5)
print(x)
```

We don't know the exact root, so we use `root_scalar` to get a substitute.

```{code-cell}
r = root_scalar(f, bracket=[0.5, 1]).root
print(r)
```

Here is the sequence of errors.

```{code-cell}
err = r - x
print(err)
```

It's not so easy to see the convergence rate by looking at these numbers. But we can check the ratios of the log of successive errors.

```{code-cell}
logerr = log(abs(err))
print(logerr[1:-1] / logerr[:-2])
```

It seems to be heading toward a constant ratio of about 1.6 before it bumps up against machine precision.

+++

# Example 4.4.3

+++

Here we look for a root of $x+\cos(10x)$ that is close to 1.

```{code-cell}
f = lambda x: x + cos(10 * x)
xx = linspace(0.5, 1.5, 400)
fig, ax = subplots()
ax.plot(xx, f(xx), label="function")
ax.grid()
xlabel("$x$")
ylabel("$y$")
```

```{code-cell}
r = root_scalar(f, bracket=[0.9, 1]).root
print(r)
```

We choose three values to get the iteration started.

```{code-cell}
x = array([0.8, 1.2, 1])
y = f(x)
ax.plot(x, y, "ko", label="initial points")
ax.legend()
fig
```

If we were using "forward" interpolation, we would ask for the polynomial interpolant of $y$ as a function of $x$. But that parabola has no real roots.

```{code-cell}
q = poly1d(polyfit(x, y, 2))  # interpolating polynomial
ax.plot(xx, q(xx), "--", label="interpolant")
ax.set_ylim(-0.1, 3)
ax.legend()
fig
```

To do inverse interpolation, we swap the roles of $x$ and $y$ in the interpolation.

```{code-cell}
plot(xx, f(xx), label="function")
plot(x, y, "ko", label="initial points")

q = poly1d(polyfit(y, x, 2))  # inverse interpolating polynomial
yy = linspace(-0.1, 2.6, 400)
plot(q(yy), yy, "--", label="inverse interpolant")

grid()
xlabel("$x$")
ylabel("$y$")
legend()
```

We seek the value of $x$ that makes $y$ zero. This means evaluating $q$ at zero.

```{code-cell}
x = hstack([x, q(0)])
y = hstack([y, f(x[-1])])
print("x:", x, "\ny:", y)
```

We repeat the process a few more times.

```{code-cell}
for k in range(5):
    q = poly1d(polyfit(y[-3:], x[-3:], 2))
    x = hstack([x, q(0)])
    y = hstack([y, f(x[-1])])
```

Here is the sequence of errors.

```{code-cell}
err = x - r
print(err)
```

The error seems to be superlinear, but subquadratic.

```{code-cell}
logerr = log(abs(err))
print("ratios:", logerr[1:] / logerr[:-1])
```

# Example 4.5.2

+++

Let us use Newton's method on the system defined by the function

```{code-cell}
def nlvalue(x):
    return array(
        [exp(x[1] - x[0]) - 2, x[0] * x[1] + x[2], x[1] * x[2] + x[0] ** 2 - x[1]]
    )
```

Here is a function that computes the Jacobian matrix.

```{code-cell}
def nljac(x):
    J = array(
        [
            [-exp(x[1] - x[0]), exp(x[1] - x[0]), 0],
            [x[1], x[0], 1],
            [2 * x[0], x[2] - 1, x[1]],
        ]
    )
    return J
```

(These functions could be embedded within a single function in other implementations.) Our initial guess at a root is the origin.

```{code-cell}
X = zeros([3, 7])
Y = copy(X)
x = X[:, 0]
```

We need the value (residual) of the nonlinear function, and its Jacobian, at this value for $\mathbf{x}$.

```{code-cell}
Y[:, 0] = nlvalue(x)
J = nljac(x)
print("f(x):", Y[:, 0])
```

We solve for the Newton step and compute the new estimate.

```{code-cell}
s = -linalg.solve(J, Y[:, 0])
X[:, 1] = x + s
print(X[:, :2])
```

Here is the new residual.

```{code-cell}
Y[:, 1] = nlvalue(X[:, 1])
print("f(x):", Y[:, 1])
```

We don't seem to be especially close to a root. Let's iterate a few more times.

```{code-cell}
for n in arange(2, 7):
    s = -linalg.solve(nljac(X[:, n - 1]), Y[:, n - 1])
    X[:, n] = X[:, n - 1] + s
    Y[:, n] = nlvalue(X[:, n])
```

We find the infinity norm of the residuals.

```{code-cell}
print("residual norm:", amax(abs(Y), 0))  # max in dimension 1
```

We don't know an exact answer, so we will take the last computed value as its surrogate.

```{code-cell}
r = X[:, -1]
x = X[:, :-1]
```

The following will subtract r from every column of x.

```{code-cell}
e = x - r.reshape([3, 1])
```

Now we take infinity norms and check for quadratic convergence.

```{code-cell}
errs = amax(abs(e), 0)
print("ratios:", log(errs[1:]) / log(errs[:-1]))
```

For a brief time, we see ratios around 2, but then the limitation of double precision makes it impossible for the doubling to continue.

+++

# Example 4.5.3

+++

As before, the system is defined by its residual and Jacobian, but this time we implement them as a single function.

```{code-cell}
def nlsystem(x):
    f = array(
        [exp(x[1] - x[0]) - 2, x[0] * x[1] + x[2], x[1] * x[2] + x[0] ** 2 - x[1]]
    )

    J = array(
        [
            [-exp(x[1] - x[0]), exp(x[1] - x[0]), 0],
            [x[1], x[0], 1],
            [2 * x[0], x[2] - 1, x[1]],
        ]
    )

    return f, J
```

Our initial guess is the origin. The output has one column per iteration.

```{code-cell}
x1 = zeros(3)
x = FNC.newtonsys(nlsystem, x1)
print(x)
```

The last column contains the final Newton estimate. We'll compute the residual there in order to check the quality of the result.

```{code-cell}
r = x[:, -1]
f, J = nlsystem(r)
print("f:", f)
```

Let's use the convergence to the first component of the root as a proxy for the convergence of the vectors.

```{code-cell}
print("digits:", log10(abs(x[1, :-1] - r[1])))
```

The exponents approximately double, as is expected of quadratic convergence.

+++

# Example 4.6.1

+++

To solve a nonlinear system, we need to code only the function defining the system (residual vector), and not its Jacobian.

```{code-cell}
def nlsystem(x):
    return array(
        [exp(x[1] - x[0]) - 2, x[0] * x[1] + x[2], x[1] * x[2] + x[0] ** 2 - x[1]]
    )
```

In all other respects usage is the same as for the `newtonsys` function.

```{code-cell}
x1 = zeros(3)
x = FNC.levenberg(nlsystem, x1)
print(x.T)
```

It's always a good idea to check the accuracy of the root, by measuring the residual (backward error).

```{code-cell}
r = x[:, -1]
print("backward error:", norm(nlsystem(r)))
```

Looking at the convergence of the first component, we find a subquadratic convergence rate, just as with the secant method.

```{code-cell}
print(log10(abs(x[0, :-1] - r[0])))
```

# Example 4.7.1

+++

Inhibited enzyme reactions often follow what are known as _Michaelis–Menten_ kinetics, in which a reaction rate $v$ follows a law of the form

$$v(x) = \frac{V x}{K_m + x},$$

+++

where $x$ is the concentration of a substrate. The real values $V$ and $K_m$ are parameters that are free to fit to data. For this example we cook up some artificial data with $V=2$ and $K_m=1/2$.

```{code-cell}
m = 25
x = linspace(0.05, 6, m)
y = 2 * x / (0.5 + x)  # exactly on the curve
y += 0.15 * cos(2 * exp(x / 16) * x)
# noise added
```

```{code-cell}
fig, ax = subplots()
ax.plot(x, y, "o", label="data")
xlabel("x")
ylabel("v")
```

The idea is to pretend that we know nothing of the origins of this data and use nonlinear least squares on the misfit to find the parameters in the theoretical model function $v(x)$. Note in the Jacobian that the derivatives are _not_ with respect to $x$, but with respect to the two parameters, which are contained in the vector `c`.

```{code-cell}
def misfit(c):
    V, Km = c  # rename components for clarity
    f = V * x / (Km + x) - y
    J = zeros([m, 2])
    J[:, 0] = x / (Km + x)  # d/d(V)
    J[:, 1] = -V * x / (Km + x) ** 2  # d/d(Km)
    return f, J
```

```{code-cell}
c1 = array([1, 0.75])
c = FNC.newtonsys(misfit, c1)
V, Km = c[:, -1]  # final values
print("V:", V, "Km:", Km)
model = lambda x: V * x / (Km + x)
```

The final values are close to the noise-free values of $V=2$, $K_m=0.5$ that we used to generate the data. We can calculate the amount of misfit at the end, although it's not completely clear what a "good" value would be. Graphically, the model looks reasonable.

```{code-cell}
print("final misfit norm:", norm(model(x) - y))
```

```{code-cell}
xx = linspace(0.05, 6, 300)
ax.plot(xx, model(xx), label="MM fit")
ax.legend()
fig
```

For this model, we also have the option of linearizing the fit process. Rewrite the model as $1/v= (a/x)+b$ for the new parameters $\alpha=K_m/V$ and $\beta=1/V$. This corresponds to the misfit function whose entries are

$$f_i(\alpha,\beta) = \alpha \cdot \frac{1}{x_i} + \beta - \frac{1}{y_i},$$

for $i=1,\ldots,m$. Although the misfit is nonlinear in $x$ and $y$, it's linear in the unknown parameters $\alpha$ and $\beta$, and so can be posed and solved as a linear least-squares problem.

```{code-cell}
A = [[1 / x[i], 1.0] for i in range(x.size)]
u = 1 / y
z = lstsq(A, u, rcond=None)[0]
alpha, beta = z
print("alpha:", alpha, "beta:", beta)
```

The two fits are different, because they do not optimize the same quantities.

```{code-cell}
linmodel = lambda x: 1 / (beta + alpha / x)
print("final misfit of linearized:", norm(linmodel(x) - y))
```

```{code-cell}
ax.plot(xx, linmodel(xx), label="linearized fit")
ax.legend()
fig
```
