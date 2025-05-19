---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
# Chapter 4

## Functions

(function-newton-python)=
``````{dropdown} Newton's method
:open:
```{literalinclude} fncbook/fncbook/chapter04.py
:filename: newton.py
:start-line: 4
:end-line: 34
:language: python
:linenos: true
```
```{admonition} About the code
:class: dropdown
{numref}`Function {number} <function-newton>` accepts *keyword arguments*. In the function declaration, these follow the semicolon, and when the function is called, they may be supplied as `keyword=value` in the argument list. Here, these arguments are also given default values by the assignments within the declaration. This arrangement is useful when there are multiple optional arguments, because the ordering of them doesn't matter.

The `break` statement, seen here in line 25, causes an immediate exit from the innermost loop in which it is called. It is often used as a safety valve to escape an iteration that may not be able to terminate otherwise.
```
``````

(function-secant-python)=
``````{dropdown} Secant method
:open:
```{literalinclude} fncbook/fncbook/chapter04.py
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

(function-newtonsys-python)=
``````{dropdown} Newton's method for systems
:open:
```{literalinclude} fncbook/fncbook/chapter04.py
:filename: newtonsys.py
:start-line: 68
:end-line: 96
:language: python
:linenos: true
```
````{admonition} About the code
:class: dropdown
The output of {numref}`Function {number} <function-newtonsys>` is a vector of vectors representing the entire history of root estimates. Since these should be in floating point, the starting value is converted with `float` before the iteration starts.
````
``````

(function-fdjac-python)=
``````{dropdown} Finite differences for Jacobian
:open:
```{literalinclude} fncbook/fncbook/chapter04.py
:filename: fdjac.py
:start-line: 98
:end-line: 112
:language: python
:linenos: true
```
::::{admonition} About the code
:class: dropdown
{numref}`Function {number} <function-fdjac>` is written to accept the case where $\mathbf{f}$ maps $n$ variables to $m$ values with $m\neq n$, in anticipation of {numref}`section-nonlineqn-nlsq`.

Note that a default value is given for the third argument `y₀`, and it refers to earlier arguments in the list. The reason is that in some contexts, the caller of `fdjac` may have already computed `y₀` and can supply it without computational cost, while in other contexts, it must be computed fresh. The configuration here adapts to either situation.
::::
``````

(function-levenberg-python)=
``````{dropdown} Levenberg's method
:open:
```{literalinclude} fncbook/fncbook/chapter04.py
:filename: levenberg.py
:start-line: 114
:end-line: 168
:language: python
:linenos: true
```
``````

## Examples

```{code-cell} ipython3
:tags: remove-cell
exec(open("FNC_init.py").read())
```

### 4.1 @section-nonlineqn-rootproblem
(demo-rootproblem-bessel-python)=
``````{dropdown} @demo-rootproblem-bessel
:open:

```{code-cell}
import scipy.special as special
def J3(x):
    return special.jv(3.0, x)

xx = linspace(0, 20, 500)
fig, ax = subplots()
ax.plot(xx, J3(xx))
ax.grid()
xlabel("$x$"), ylabel("$J_3(x)$")
title("Bessel function");
```
From the graph we see roots near 6, 10, 13, 16, and 19. We use `root_scalar` from the `scipy.optimize` package to find these roots accurately.

```{code-cell}
from scipy.optimize import root_scalar

omega = []
for guess in [6.0, 10.0, 13.0, 16.0, 19.0]:
    s = root_scalar(J3, bracket=[guess - 0.5, guess + 0.5]).root
    omega.append(s)

results = PrettyTable()
results.add_column("root estimate", omega)
results.add_column("function value", [J3(ω) for ω in omega])
print(results)
```

```{code-cell}
ax.scatter(omega, J3(omega))
ax.set_title("Bessel function roots")
fig
```

If instead we seek values at which $J_3(x)=0.2$, then we must find roots of the function $J_3(x)-0.2$.

```{code-cell}
omega = []
for guess in [3., 6., 10., 13.]:
    f = lambda x: J3(x) - 0.2
    s = root_scalar(f, x0=guess).root
    omega.append(s)

ax.scatter(omega, J3(omega))
fig
```
``````

(demo-roots-cond-python)=
``````{dropdown} @demo-roots-cond
:open:
Consider first the function

```{code-cell}
f = lambda x: (x - 1) * (x - 2)
```

At the root $r=1$, we have $f'(r)=-1$. If the values of $f$ were perturbed at every point by a small amount of noise, we can imagine finding the root of the function drawn with a thick ribbon, giving a range of potential roots.

```{code-cell}
xx = linspace(0.8, 1.2, 400)
plot(xx, f(xx))
plot(xx, f(xx) + 0.02, "k")
plot(xx, f(xx) - 0.02, "k")
axis("equal"), grid(True)
xlabel("x"), ylabel("f(x)")
title("Well-conditioned root");
```

The possible values for a perturbed root all lie within the interval where the ribbon intersects the $x$-axis. The width of that zone is about the same as the vertical thickness of the ribbon.

By contrast, consider the function

```{code-cell}
f = lambda x: (x - 1) * (x - 1.01)
```

Now $f'(1)=-0.01$, and the graph of $f$ will be much shallower near $x=1$. Look at the effect this has on our thick rendering:

```{code-cell}
xx = linspace(0.8, 1.2, 400)
plot(xx, f(xx))
plot(xx, f(xx) + 0.02, "k")
plot(xx, f(xx) - 0.02, "k")
axis("equal"), grid(True)
xlabel("x"), ylabel("f(x)")
title("Poorly-conditioned root");
```

The vertical displacements in this picture are exactly the same as before. But the potential _horizontal_ displacement of the root is much wider. In fact, if we perturb the function entirely upward by the amount drawn here, the root disappears!
``````

### 4.2 @section-nonlineqn-fixed-point

(demo-fp-spiral-python)=
``````{dropdown} @demo-fp-spiral
:open:
Let's convert the roots of a quadratic polynomial $f(x)$ to a fixed point problem.

```{code-cell}
f = poly1d([1, -4, 3.5])
r = f.roots
print(r)
```

We define $g(x)=x - f(x)$. 

```{code-cell}
g = lambda x: x - f(x)
```

Intersections of $y=g(x)$ with the line $y=x$ are fixed points of $g$ and thus roots of $f$. (Only one is shown in the chosen plot range.)

```{code-cell}
fig, ax = subplots()
g = lambda x: x - f(x)
xx = linspace(2, 3, 400)
ax.plot(xx, g(xx), label="y=g(x)")
ax.plot(xx, xx, label="y=x")
axis("equal"), legend()
title("Finding a fixed point");
```

If we evaluate $g(2.1)$, we get a value of almost 2.6, so this is not a fixed point.

```{code-cell}
x = 2.1
y = g(x)
print(y)
```

However, $y=g(x)$ is considerably closer to the fixed point at around 2.7 than $x$ is. Suppose then that we adopt $y$ as our new $x$ value. Changing the $x$ coordinate in this way is the same as following a horizontal line over to the graph of $y=x$.

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
    x = y       # y --> new x
    y = g(x)    # g(x) --> new y
    ax.plot([x, x], [x, y], "k:")  
fig
```

The process spirals in beautifully toward the fixed point we seek. Our last estimate has almost 4 accurate digits.

```{code-cell} 
print(abs(y - max(r)) / max(r))
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
    y = g(x)
    ax.plot([x, x], [x, y], "k:")
ylim(1, 2.5)
title("No convergence");
```

This time, the iteration is pushing us _away from_ the correct answer.
``````

(demo-fp-converge-python)=
``````{dropdown} @demo-fp-converge
:open:
We revisit @demo-fp-spiral and investigate the observed convergence more closely. Recall that above we calculated $g'(p)\approx-0.42$ at the convergent fixed point.

```{code-cell}
f = poly1d([1, -4, 3.5])
r = f.roots
print(r)
```

Here is the fixed-point iteration. This time we keep track of the whole sequence of approximations.

```{code-cell}
g = lambda x: x - f(x)
x = zeros(12)
x[0] = 2.1
for k in range(11):
    x[k + 1] = g(x[k])

print(x)
```

It's illuminating to construct and plot the sequence of errors.

```{code-cell}
err = abs(x - max(r))
semilogy(err, "-o")
xlabel("iteration number"), ylabel("error")
title("Convergence of fixed-point iteration");
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

The error should therefore decrease by a factor of $\sigma$ at each iteration. We can check this easily from the observed data.

```{code-cell}
err[8:] / err[7:-1]
```

The methods for finding $\sigma$ agree well.
``````

### 4.3 @section-nonlineqn-newton
(demo-newton-line-python)=
``````{dropdown} @demo-newton-line
:open:

Suppose we want to find a root of this function:

```{code-cell}
f = lambda x: x * exp(x) - 2
xx = linspace(0, 1.5, 400)

fig, ax = subplots()
ax.plot(xx, f(xx), label="function")
ax.grid()
ax.set_xlabel("$x$")
ax.set_ylabel("$y$");
```

From the graph, it is clear that there is a root near $x=1$. So we call that our initial guess, $x_1$.

```{code-cell}
x1 = 1
y1 = f(x1)
ax.plot(x1, y1, "ko", label="initial point")
ax.legend()
fig
```

Next, we can compute the tangent line at the point $\bigl(x_1,f(x_1)\bigr)$, using the derivative.

```{code-cell}
df_dx = lambda x: exp(x) * (x + 1)
slope1 = df_dx(x1)
tangent1 = lambda x: y1 + slope1 * (x - x1)

ax.plot(xx, tangent1(xx), "--", label="tangent line")
ax.set_ylim(-2, 4)
ax.legend()
fig
```

In lieu of finding the root of $f$ itself, we settle for finding the root of the tangent line approximation, which is trivial. Call this $x_2$, our next approximation to the root.

```{code-cell}
x2 = x1 - y1 / slope1
ax.plot(x2, 0, "ko", label="tangent root")
ax.legend()
fig
```

```{code-cell}
y2 = f(x2)
print(y2)
```

The residual (i.e., value of $f$) is smaller than before, but not zero. So we repeat the process with a new tangent line based on the latest point on the curve.

```{code-cell}
xx = linspace(0.83, 0.88, 200)

plot(xx, f(xx))
plot(x2, y2, "ko")
grid(), xlabel("$x$"), ylabel("$y$")

slope2 = df_dx(x2)
tangent2 = lambda x: y2 + slope2 * (x - x2)
plot(xx, tangent2(xx), "--")
x3 = x2 - y2 / slope2
plot(x3, 0, "ko")
title("Second iteration");
```

```{code-cell}
y3 = f(x3)
print(y3)
```

Judging by the residual, we appear to be getting closer to the true root each time.
``````

(demo-newton-converge-python)=
``````{dropdown} @demo-newton-converge
:open:
We again look at finding a solution of $x e^x=2$ near $x=1$. To apply Newton's method, we need to calculate values of both the residual function $f$ and its derivative.

```{code-cell}
f = lambda x: x * exp(x) - 2
df_dx = lambda x: exp(x) * (x + 1)
```

We don't know the exact root, so we use `nlsolve` to determine a proxy for it.

```{code-cell}
r = root_scalar(f, bracket=[0.8, 1.0]).root
print(r)
```

We use $x_1=1$ as a starting guess and apply the iteration in a loop, storing the sequence of iterates in a vector.

```{code-cell}
x = ones(5)
for k in range(4):
    x[k + 1] = x[k] - f(x[k]) / df_dx(x[k])

print(x)
```

Here is the sequence of errors.

```{code-cell}
err = x - r
print(err)
```

The exponents in the scientific notation definitely suggest a squaring sequence. We can check the evolution of the ratio in {eq}`quadratictest`.

```{code-cell}
logerr = log(abs(err))
for i in range(len(err) - 1):
    print(logerr[i+1] / logerr[i])
```

The clear convergence to 2 above constitutes good evidence of quadratic convergence.
``````

(demo-newton-usage-python)=
``````{dropdown} @demo-newton-usage
:open:
Suppose we want to evaluate the inverse of the function $h(x)=e^x-x$. This means solving $y=h(x)$, or $h(x)-y=0$, for $x$ when $y$ is given. That equation has no solution in terms of elementary functions. If a value of $y$ is given numerically, though, we simply have a rootfinding problem for $f(x)=e^x-x-y$.
```{tip}
:class: dropdown
The `enumerate` function produces a pair of values for each iteration: a positional index and the corresponding contents.
```

```{index} ! Python; enumerate
```

```{code-cell}
h = lambda x: exp(x) - x
dh_dx = lambda x: exp(x) - 1
y_ = linspace(h(0), h(2), 200)
x_ = zeros(y_.shape)
for (i, y) in enumerate(y_):
    f = lambda x: h(x) - y
    df_dx = lambda x: dh_dx(x)
    x = FNC.newton(f, df_dx, y)
    x_[i] = x[-1]

plot(x_, y_, label="$y=h(x)$")
plot(y_, x_, label="$y=h^{-1}(x)$")
plot([0, max(y_)], [0, max(y_)], 'k--', label="")
title("Function and its inverse")
xlabel("x"), ylabel("y"), axis("equal")
ax.grid()
legend();
```
``````
### 4.4 @section-nonlineqn-secant
(demo-secant-line-python)=
``````{dropdown} @demo-secant-line
:open:

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
``````

(demo-secant-converge-python)=
``````{dropdown} @demo-secant-converge
:open:
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
``````

(demo-secant-iqi-python)=
``````{dropdown} @demo-secant-iqi
:open:
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
``````

### 4.5 @section-nonlineqn-newtonsys
(demo-newtonsys-converge-python)=
``````{dropdown} @demo-newtonsys-converge
:open:
A system of nonlinear equations is defined by its residual and Jacobian.

```{code-cell}
def func(x):
    return array([
        exp(x[1] - x[0]) - 2, 
        x[0] * x[1] + x[2], 
        x[1] * x[2] + x[0]**2 - x[1]
    ])

def jac(x):
    return array([
            [-exp(x[1] - x[0]), exp(x[1] - x[0]), 0],
            [x[1], x[0], 1],
            [2 * x[0], x[2] - 1, x[1]],
    ])
```

Our initial guess at a root is the origin. 

```{code-cell}
x1 = zeros(3)
x = FNC.newtonsys(func, jac, x1)
print(x)
```

The output has one row per iteration, so the last row contains the final Newton estimate. Let's compute its residual.

```{code-cell}
r = x[-1]
f = func(r)
print("final residual:", f)
```

Let's check the convergence rate:

```{code-cell}
logerr = [log(norm(x[k] - r)) for k in range(x.shape[0] - 1)]
for k in range(len(logerr) - 1):
    print(logerr[k+1] / logerr[k])
```

The ratio is apparently converging toward 2, as expected for quadratic convergence.
``````

### 4.6 @section-nonlineqn-quasinewton
(demo-quasi-levenberg-python)=
``````{dropdown} @demo-quasi-levenberg
:open:
To solve a nonlinear system, we need to code only the function defining the system, and not its Jacobian.

```{code-cell}
def func(x):
    return array([
        exp(x[1] - x[0]) - 2, 
        x[0] * x[1] + x[2], 
        x[1] * x[2] + x[0]**2 - x[1]
    ])
```

In all other respects usage is the same as for the `newtonsys` function.

```{code-cell}
x1 = zeros(3)
x = FNC.levenberg(func, x1)
print(f"Took {len(x) - 1} iterations.")
```

It's always a good idea to check the accuracy of the root, by measuring the residual (backward error).

```{code-cell}
r = x[-1]
print("backward error:", norm(func(r)))
```
Looking at the convergence in norm, we find a convergence rate between linear and quadratic, like with the secant method:

```{code-cell}
logerr = [log(norm(x[k] - r)) for k in range(len(x) - 1)]
for k in range(len(logerr) - 1):
    print(logerr[k+1] / logerr[k])
```
``````
### 4.7 @section-nonlineqn-nlsq
(demo-nlsq-converge-python)=
``````{dropdown} @demo-nlsq-converge
:open:
We will observe the convergence of {numref}`Function {number} <function-levenberg>` for different levels of the minimum least-squares residual. We start with a function mapping from $\real^2$ into $\real^3$, and a point that will be near the optimum.

```{code-cell}
g = lambda x: array([sin(x[0] + x[1]), cos(x[0] - x[1]), exp(x[0] - x[1])])
p = array([1, 1])
```

The function $\mathbf{g}(\mathbf{x}) - \mathbf{g}(\mathbf{p})$ obviously has a zero residual at $\mathbf{p}$. We'll make different perturbations of that function in order to create nonzero residuals.

```{code-cell}
for R in [1e-3, 1e-2, 1e-1]:
    # Define the perturbed function.
    f = lambda x: g(x) - g(p) + R * array([-1, 1, -1]) / sqrt(3)
    x = FNC.levenberg(f, [0, 0])
    r = x[-1]
    err = [norm(x[j] - r) for j in range(len(x) - 1)]
    normres = norm(f(r))
    semilogy(err, label=f"R={normres:.2g}")
title("Convergence of Gauss–Newton")
xlabel("iteration"), ylabel("error")
legend();
```

In the least perturbed case, where the minimized residual is less than $10^{-3}$, the convergence is plausibly quadratic. At the next level up, the convergence starts similarly but suddenly stagnates for a long time. In the most perturbed case, the quadratic phase is nearly gone and the overall shape looks linear.
``````

(demo-nlsq-MM-python)=
``````{dropdown} @demo-nlsq-MM
:open:
```{code-cell}
m = 25
V, Km = 2, 0.5
s = linspace(0.05, 6, m)
model = lambda x: V * x / (Km + x)
w = model(s) + 0.15 * cos(2 * exp(s / 16) * s)    # noise added

fig, ax = subplots()
ax.scatter(s, w, label="data")
ax.plot(s, model(s), 'k--', label="unperturbed model")
xlabel("s"), ylabel("w")
legend();
```

```{index} ! Python; destructuring
```

The idea is to pretend that we know nothing of the origins of this data and use nonlinear least squares to find the parameters in the theoretical model function $v(s)$. In {eq}`nlsq-misfit`, the $s$ variable plays the role of $t$, and $v$ plays the role of $g$.
```{tip}
:class: dropdown
Putting comma-separated values on the left of an assignment will **destructure** the right-hand side, drawing individual assignments from entries of a vector, for example.
```

```{code-cell}
def misfit(c):
    V, Km = c  # rename components for clarity
    f = V * s / (Km + s) - w
    return f
```

In the Jacobian the derivatives are with respect to the parameters in $\mathbf{x}$.

```{code-cell}
def misfitjac(x):
    V, Km = x   # rename components for clarity
    J = zeros([m, 2])
    J[:, 0] = s / (Km + s)          # d/d(V)
    J[:, 1] = -V * s / (Km + s)**2  # d/d(Km)
    return J
```

```{code-cell}
x1 = [1, 0.75]
x = FNC.newtonsys(misfit, misfitjac, x1)
V, Km = x[-1]  # final values
print(f"estimates are V = {V:.3f}, Km = {Km:.3f}")
```

The final values are reasonably close to the values $V=2$, $K_m=0.5$ that we used to generate the noise-free data. Graphically, the model looks close to the original data.

```{code-cell}
# since V and Km have been updated, model() is too
ax.plot(s, model(s), label="nonlinear fit")
```

For this particular model, we also have the option of linearizing the fit process. Rewrite the model as 

```{math}
:enumerated: false
\frac{1}{w} = \frac{\alpha}{s} + \beta = \alpha \cdot s^{-1} + \beta
```

for the new fitting parameters $\alpha=K_m/V$ and $\beta=1/V$. This corresponds to the misfit function whose entries are

$$f_i([\alpha,\beta]) = \left(\alpha \cdot \frac{1}{s_i} + \beta\right) - \frac{1}{w_i}$$
for $i=1,\ldots,m$. Although this misfit is nonlinear in $s$ and $w$, it's linear in the unknown parameters $\alpha$ and $\beta$. This lets us pose and solve it as a linear least-squares problem.

```{code-cell}
from numpy.linalg import lstsq
A = array( [[1 / s[i], 1.0] for i in range(len(s))] )
z = lstsq(A, 1 / w, rcond=None)[0]
alpha, beta = z
print("alpha:", alpha, "beta:", beta)
```

The two fits are different; they do not optimize the same quantities.

```{code-cell}
linmodel = lambda x: 1 / (beta + alpha / x)
ax.plot(s, linmodel(s), label="linear fit")
ax.legend()
fig
```

The truly nonlinear fit is clearly better in this case. It optimizes a residual for the original measured quantity rather than a transformed one we picked for algorithmic convenience.
``````