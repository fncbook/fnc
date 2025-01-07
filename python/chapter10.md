---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
# Chapter 10

## Functions
(function-shoot-python)=
``````{dropdown} Shooting method for a two-point boundary-value problem
:open:
```{literalinclude} pkg/fncbook/chapter10.py
:filename: shoot.py
:start-at: def shoot
:end-at: return x, u, du_dx
:language: python
:linenos: true
```

``````

(function-diffmats2-python)=
``````{dropdown} Second-order differentiation matrices
:open:
```{literalinclude} pkg/fncbook/chapter10.py
:filename: diffmats2.py
:start-at: def diffmat2
:end-at: return x, Dx, Dxx
:language: python
:linenos: true
```
``````

(function-diffcheb-python)=
``````{dropdown} Chebyshev differentiation matrices
:open:
```{literalinclude} pkg/fncbook/chapter10.py
:filename: diffcheb.py
:start-at: def diffcheb
:end-at: return x, Dx, Dxx
:language: python
:linenos: true
```
``````

(function-bvplin-python)=
``````{dropdown} Solution of a linear boundary-value problem
:open:
```{literalinclude} pkg/fncbook/chapter10.py
:filename: bvplin.py
:start-at: def bvplin(
:end-at: return x, u
:language: python
:linenos: true
```

::::{admonition} About the code
:class: dropdown
Note that there is no need to explicitly form the row-deletion matrix $\mathbf{E}$ from {eq}`rowdeletion`. Since it only appears as left-multiplying $\mathbf{L}$ or $\mathbf{r}$, we simply perform the row deletions as needed using indexing.
::::
``````

(function-bvp-python)=
``````{dropdown} Solution of a nonlinear boundary-value problem
:open:
```{literalinclude} pkg/fncbook/chapter10.py
:filename: bvp.py
:start-at: def bvp(
:end-at: return x
:language: python
:linenos: true
```
:::{admonition} About the code
:class: dropdown
The nested function `residual` uses differentiation matrices computed externally to it, rather than computing them anew on each invocation. As in {numref}`Function {number} <function-bvplin>`, there is no need to form the row-deletion matrix $\mathbf{E}$ explicitly. In lines 23--24, we divide the values of $g_1$ and $g_2$ by a factor of $h$. This helps scale the residual components more uniformly and improves the robustness of convergence a bit.
:::
``````

(function-fem-python)=
``````{dropdown} Piecewise linear finite elements for a linear BVP
:open:
```{literalinclude} pkg/fncbook/chapter10.py
:filename: fem.py
:start-at: def fem
:end-at: return x, u
:language: python
:linenos: true
```
``````

## Examples

```{code-cell}
:tags: [remove-cell]
import os
print(os.chdir("/Users/driscoll/Documents/GitHub/fnc/python"))
exec(open("FNC_init.py").read())
```

### 10.1 @section-bvp-tpbvp

(demo-tpbvp-mems-python)=
``````{dropdown} @demo-tpbvp-mems
To solve this problem, we have to define functions for the ODE and boundary conditions. The first returns the computed values of $y_1'$ and $y_2'$.  

```{code-cell}
lamb = 0.6
def ode(r, y):
    return array([
        y[1],
        lamb / y[0]**2 - y[1] / r
    ])
```

To encode the boundary conditions $y_2(0)=0$, $y_1(2)=1$, we define a function for their residual values.

```{code-cell}
def bc(ya, yb):    # given y(a), y(b)
    return array([
        ya[1],
        yb[0] - 1
    ])
```

The domain of the mathematical problem is $r\in [0,1]$. However, there is a division by $r$ in the ODE, so we want to avoid $r=0$ by truncating the domain a bit.

```{code-cell}
a, b = finfo(float).eps, 1
```

We need one last ingredient that is not part of the mathematical setup: an initial estimate for the solution. As we will see, this plays the same role as initialization in Newton's method for rootfinding. Here, we try a constant value for each component.

```{code-cell}
r = linspace(a, b, 50)
y_init = vstack([ones(r.size), zeros(r.size)])
```

Now we can solve the problem using `solve_bvp` from `scipy.integrate`.

```{code-cell}
from scipy.integrate import solve_bvp
sol = solve_bvp(ode, bc, r, y_init)
print(f"Solved at {sol.x.size} nodes.")
plot(sol.x, sol.y[0])
xlabel("$r$"),  ylabel("$w(r)$")
title("Solution of MEMS problem for $\\lambda=0.6$");
```
``````

### 10.2 @section-bvp-shooting
(demo-shooting-naive-python)=
``````{dropdown} @demo-shooting-naive

Let's first examine the shooting approach for the TPBVP from {numref}`Example {number} <example-tpbvp-mems>` with $\lambda=0.6$. 

```{code-cell}
lamb = 0.6
phi = lambda r, w, dw_dr: lamb / w**2 - dw_dr / r
```

We convert the ODE to a first-order system in order to apply a numerical method. We also have to truncate the domain to avoid division by zero.

```{code-cell}
f = lambda r, y: hstack([y[1], phi(r, y[0], y[1])])
a, b = finfo(float).eps, 1
```

The BVP specifies $w'(0)=y_2(0)=0$. We can try multiple values for the unknown $w(0)=y_1(0)$ and plot the solutions.

```{code-cell}
from scipy.integrate import solve_ivp
t = linspace(a, b, 400)
for w0 in arange(0.4, 1.0, 0.1):
    sol = solve_ivp(f, [a, b], [w0, 0], t_eval=t)
    plot(t, sol.y[0], label=f"$w_0$ = {w0:.1f}")

xlabel("$r$"),  ylabel("$w(r)$")
legend(),  grid(True)
title("Solutions for choices of w(0)");
```

On the graph, it's the curve starting at $w(0)=0.8$ that comes closest to the required condition $w(1)=1$, but it's a bit too large.
``````

(demo-shooting-mems-python)=
``````{dropdown} @demo-shooting-mems
We revisit {numref}`Demo {number} <demo-shooting-naive>` but let {numref}`Function {number} <function-shoot>` do the heavy lifting.

```{code-cell}
lamb = 0.6
phi = lambda r, w, dwdr: lamb / w**2 - dwdr / r
a, b = finfo(float).eps, 1
```

We specify the given and unknown endpoint values.

```{code-cell}
ga = lambda w, dw : dw       # w'=0 at left
gb = lambda w, dw : w - 1    # w=1 at right
```

In this setting, we need to provide initial guesses for $w(a)$ and $w'(a)$.

```{code-cell}
init = array([0.8, 0])
r, w, dw_dx = FNC.shoot(phi, a, b, ga, gb, init)
plot(r, w)
title("Shooting solution")
xlabel("$r$"),  ylabel("$w(r)$");
```

The value of $w$ at $r=1$, meant to be exactly one, was computed to be

```{code-cell}
print(f"w at right end is {w[-1]}")
```

The accuracy is consistent with the error tolerance used for the IVP solution by `shoot`. The initial value $w(0)$ that gave this solution is

```{code-cell}
print(f"w at left end is {w[0]}")
```
``````

(demo-shooting-unstable-python)=
``````{dropdown} @demo-shooting-unstable

```{code-cell}
ga = lambda u, du : u + 1    # u=-1 at left
gb = lambda u, du : u        # u= 0 at right
init = array([-1, 0])
for lamb in range(6, 22, 4):
    phi = lambda x, u, du_dx: lamb**2 * u + lamb**2
    x, u, du_dx = FNC.shoot(phi, 0.0, 1.0, ga, gb, init)
    plot(x, u, label=f"$\\lambda$ = {lamb:.1f}")

xlabel("$x$"),  ylabel("$u(x)$"),  ylim(-1.0, 0.25)
grid(True),  legend(loc="upper left")
title("Shooting instability");
```

The numerical solutions evidently don't satisfy the right boundary condition as $\lambda$ increases, which makes them invalid. 

``````

### 10.3 @section-bvp-diffmats
(demo-diffmats-2nd-python)=
``````{dropdown} @demo-diffmats-2nd
We test first-order and second-order differentiation matrices for the function $x + \exp(\sin 4x)$ over $[-1,1]$.

```{code-cell}
f = lambda x: x + exp( sin(4 * x) )
```

For reference, here are the exact first and second derivatives.

```{code-cell}
df_dx = lambda x: 1 + 4 * exp(sin(4 * x)) * cos(4 * x)
d2f_dx2 = lambda x: 4 * exp(sin(4 * x)) * (4 * cos(4 * x)**2 - 4 * sin(4 * x))
```

We discretize on equally spaced nodes and evaluate $f$ at the nodes.

```{code-cell}
t, Dx, Dxx = FNC.diffmat2(12, [-1, 1])
y = f(t)
```

Then the first two derivatives of $f$ each require one matrix-vector multiplication.

```{code-cell}
yx = Dx @ y
yxx = Dxx @ y
```

The results show poor accuracy for this small value of $n$.

```{code-cell}
x = linspace(-1, 1, 500)
subplot(2, 1, 1)
plot(x, df_dx(x))
plot(t, yx, "ko")
xlabel("$x$"),  ylabel("$f'(x)$")

subplot(2, 1, 2)
plot(x, d2f_dx2(x))
plot(t, yxx, "ko")
xlabel("$x$"),  ylabel("$f''(x)$");
```

A convergence experiment confirms the order of accuracy. Because we expect an algebraic convergence rate, we use a log-log plot of the errors.

```{code-cell}
N = array([int(2**k) for k in arange(4, 11.5, 0.5)])
err1 = zeros(len(N))
err2 = zeros(len(N))
for k, n in enumerate(N):
    t, Dx, Dxx = FNC.diffmat2(n, [-1, 1])
    y = f(t)
    err1[k] = norm(df_dx(t) - Dx @ y, inf)
    err2[k] = norm(d2f_dx2(t) - Dxx @ y, inf)

loglog(N, err1, "-o", label="$f'$")
loglog(N, err2, "-o", label="$f''$")
plot(N, 10 * 10 / N**2, "k--", label="2nd order")
xlabel("$n$"),  ylabel("max error")
legend(loc="lower left")
title("Convergence of finite differences");
```
``````

(demo-diffmats-cheb-python)=
``````{dropdown} @demo-diffmats-cheb
Here is a $4\times 4$ Chebyshev differentiation matrix.

```{code-cell}
t, Dx, Dxx = FNC.diffcheb(3, [-1, 1])
print(Dx)
```

We again test the convergence rate.

```{code-cell}
f = lambda x: x + exp(sin(4 * x))
df_dx = lambda x: 1 + 4 * exp(sin(4 * x)) * cos(4 * x)
d2f_dx2 = lambda x: 4 * exp(sin(4 * x)) * (4 * cos(4 * x) ** 2 - 4 * sin(4 * x))

N = range(5, 75, 5)
err1 = zeros(len(N))
err2 = zeros(len(N))
err = zeros((len(N), 2))
for k, n in enumerate(N):
    t, Dx, Dxx = FNC.diffcheb(n, [-1, 1])
    y = f(t)
    err[k, 0] = norm(df_dx(t) - Dx @ y, inf)
    err[k, 1] = norm(d2f_dx2(t) - Dxx @ y, inf)
```

Since we expect a spectral convergence rate, we use a semi-log plot for the error.

```{code-cell}
semilogy(N, err, "-o")
xlabel("$n$"), ylabel("max error")
legend(["$f'$", "$f''$"], loc="lower left")
title("Convergence of Chebyshev derivatives");
```
``````

### 10.4 @section-bvp-linear
(demo-linear-solve-python)=
``````{dropdown} @demo-linear-solve

```{code-cell}
exact = lambda x: exp( sin(x) )
```

The problem is presented above in our standard form, so we can identify the coefficient functions in the ODE. Each should be coded as a function.

```{code-cell}
p = lambda x: -cos(x)
q = sin
r = lambda x: 0 * x    # must be a function 
```

We solve the BVP and compare the result to the exact solution.

```{code-cell}
x, u = FNC.bvplin(p, q, r, [0, pi/2], 1, exp(1), 25)
```

```{code-cell}
subplot(2, 1, 1)
plot(x, u)
ylabel("solution"),  title("Solution of the BVP")

subplot(2, 1, 2)
plot(x, exact(x) - u, "-o")
ylabel("error");
```
``````

(demo-linear-converge-python)=
``````{dropdown} @demo-linear-converge
```{code-cell}
lamb = 10
exact = lambda x: sinh(lamb * x) / sinh(lamb) - 1
```

The following functions define the ODE.

```{code-cell}
p = lambda x: zeros(size(x))
q = lambda x: -(lamb**2) * ones(len(x))
r = lambda x: lamb**2 * ones(len(x))
```

We compare the computed solution to the exact one for increasing $n$.

```{code-cell}
N = array([int(2 * 10**d) for d in arange(1, 3.1, 0.25)])
err = zeros(len(N))
results = PrettyTable(["n", "error"])
for k, n in enumerate(N):
    x, u = FNC.bvplin(p, q, r, [0, 1], -1, 0, n)
    err[k] = norm(exact(x) - u, inf)
    results.add_row([n, err[k]])
print(results)
```

Each factor of 10 in $n$ reduces error by a factor of 100, which is indicative of second-order convergence.

```{code-cell}
loglog(N, err, "-o", label="observed")
loglog(N, 1 / N**2, "--", label="2nd order")
xlabel("$n$"),  ylabel("max error")
legend(),  title("Convergence of finite differences");
```
``````

### 10.5 @section-bvp-nonlinear
(demo-nonlinear-pendulum-python)=
``````{dropdown} @demo-nonlinear-pendulum
The first step is to define the function $\phi$ that equals $\theta''$.

```{code-cell}
phi = lambda t, theta, omega: -0.05 * omega - sin(theta)
```

Next, we define the boundary conditions.

```{code-cell}
ga = lambda u, du: u - 2.5
gb = lambda u, du: u + 2
```

The last ingredient is an initial estimate of the solution. Here we choose $n=100$ and a linear function between the endpoint values. 

```{code-cell}
init = linspace(2.5, -2, 101)
```

We find a solution with negative initial slope, i.e., the pendulum is initially pushed back toward equilibrium.

```{code-cell}
t, theta = FNC.bvp(phi, [0, 5], ga, gb, init)
plot(t, theta)
xlabel("$t$")
ylabel("$\theta(t)$")
title("Pendulum over [0,5]");
```

If we extend the time interval longer for the same boundary values, then the initial slope must adjust.

```{code-cell}

t, theta = FNC.bvp(phi, [0, 8], ga, gb, init)
plot(t, theta)
xlabel("$t$")
ylabel("$\theta(t)$")
title("Pendulum over [0,8]");
```

This time, the pendulum is initially pushed toward the unstable equilibrium in the upright vertical position before gravity pulls it back down.
``````

(demo-nonlinear-mems-python)=
``````{dropdown} @demo-nonlinear-mems
Here is the problem definition. We use a truncated domain to avoid division by zero at $r=0$.

```{code-cell}
lamb = 0.5
phi = lambda r, w, dwdr: lamb / w**2 - dwdr / r
a, b = finfo(float).eps, 1
ga = lambda w, dw: dw
gb = lambda w, dw: w - 1
```

First we try a constant function as the initialization.

```{code-cell}
init = ones(201)
r, w1 = FNC.bvp(phi, [a, b], ga, gb, init)
plot(r, w1)
fig, ax = gcf(), gca()
xlabel("$r$"),  ylabel("$w(r)$")
title("Solution of the MEMS problem");
```

It's not necessary that the initialization satisfy the boundary conditions. In fact, by choosing a different constant function as the initial guess, we arrive at another valid solution.

```{code-cell}
r, w2 = FNC.bvp(phi, [a, b], ga, gb, 0.5 * init)
ax.plot(r, w2)
ax.set_title("Multiple solutions of the MEMS problem");
fig
```
``````

(demo-nonlinear-allencahn-python)=
``````{dropdown} @demo-nonlinear-allencahn

```{code-cell}
phi = lambda x, u, dudx: (u**3 - u) / epsilon
ga = lambda u, du: du
gb = lambda u, du: u - 1
```

Finding a solution is easy at larger values of $\epsilon$.

```{code-cell}
epsilon = 0.05
init = linspace(-1, 1, 141)
x, u1 = FNC.bvp(phi, [0, 1], ga, gb, init)

plot(x, u1, label="$\\epsilon = 0.05$")
fig, ax = gcf(), gca()
xlabel("$x$"),  ylabel("$u(x)$")
legend(),  title("Allen-Cahn solution");
```

Finding a good initialization is not trivial for smaller values of $\epsilon$. But the iteration succeeds if we use the first solution as the initialization at the smaller $\epsilon$.


```{code-cell}
epsilon = 0.002
x, u2 = FNC.bvp(phi, [0, 1], ga, gb, u1)
ax.plot(x, u2, label="$\\epsilon = 0.002$")
ax.legend()
fig
```

In this case we can continue further.

```{code-cell}
ϵ = 0.0005
x, u3 = FNC.bvp(phi, [0, 1], ga, gb, u2)
ax.plot(x, u3, label="$\\epsilon = 0.005$")
ax.legend()
fig
```
``````

### 10.6 @section-bvp-galerkin
(demo-galerkin-fem-python)=
``````{dropdown} @demo-galerkin-fem
Here are the coefficient function definitions. Even though $s$ is a constant, it has to be defined as a function for {numref}`Function {number} <function-fem>` to use it.

```{code-cell}
c = lambda x: x**2
q = lambda x: 4 * ones(len(x))
f = lambda x: sin(pi * x)
```

```{code-cell}
x, u = FNC.fem(c, q, f, 0, 1, 50)
plot(x, u)
xlabel("$x$"),  ylabel("$u$")
title("Solution by finite elements");
```
``````
