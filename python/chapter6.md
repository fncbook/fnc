---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
# Chapter 6

## Functions

(function-euler-python)=
``````{dropdown} Euler's method for an initial-value problem
```{literalinclude} ../python/pkg/FNC/FNC06.py
:filename: euler.py
:start-at: def euler
:end-at: return t
:language: python
:linenos: true
```
``````

(function-ie2-python)=
``````{dropdown} Improved Euler method for an IVP
```{literalinclude} ../python/pkg/FNC/FNC06.py
:filename: ie2.py
:start-at: def ie2
:end-at: return t
:language: python
:linenos: true
```
``````

(function-rk4-python)=
``````{dropdown} Fourth-order Runge-Kutta for an IVP
```{literalinclude} ../python/pkg/FNC/FNC06.py
:filename: rk4.py
:start-at: def rk4
:end-at: return t
:language: python
:linenos: true
```
``````

(function-rk23-python)=
``````{dropdown} Adaptive IVP solver based on embedded RK formulas
```{literalinclude} ../python/pkg/FNC/FNC06.py
:filename: rk23.py
:start-at: def euler
:end-at: return
:language: python
:linenos: true
```
::::{admonition} About the code
:class: dropdown
The check `t[i]+h==t[i]`on line 19 is to detect when $h$ has become so small that it no longer changes the floating-point value of $t_i$. This may be a sign that the underlying exact solution has a singularity near $t=t_i$, but in any case, the solver must halt by using a `break` statement to exit the loop.

On line 30, we use a combination of absolute and relative tolerances to judge the acceptability of a solution value, as in {eq}`absreltolerance`. In lines 41--43 we underestimate the step factor $q$ a bit and prevent a huge increase in the step size, since a rejected step is expensive, and then we make sure that our final step doesn't take us past the end of the domain.

Finally, line 37 exploits a subtle property of the BS23 formula called *first same as last* (FSAL).
While {eq}`bs23` calls for four stages to find the paired second- and third-order estimates, the final stage computed in stepping from $t_i$ to $t_{i+1}$ is identical to the first stage needed to step from $t_{i+1}$ to $t_{i+2}$. By repurposing `s₄` as `s₁` for the next pass, one of the stage evaluations comes for free, and only three evaluations of $f$ are needed per successful step.
::::

``````

(function-ab4-python)=
``````{dropdown} 4th-order Adams–Bashforth formula for an IVP
```{literalinclude} ../python/pkg/FNC/FNC06.py
:filename: ab4.py
:start-at: def ab4
:end-at: return t
:language: python
:linenos: true
```
::::{admonition} About the code
:class: dropdown
Line 15 sets `σ` to be the coefficients of the generating polynomial $\sigma(z)$ of AB4. Lines 19--21 set up the IVP over the time interval $a \le t \le a+3 h$, call `rk4` to solve it using the step size $h$, and use the result to fill the first four values of the solution. Then line 24 computes the vector $[f_2,f_1,f_0]$.

Line 28 computes $f_i$, based on the most recent solution value and time. That goes into the first spot of `f`, followed by the three values that were previously most recent. These are the four values that appear in {eq}`ab4`. Each particular $f_i$ value starts at the front of `f`, moves through each position in the vector over three iterations, and then is forgotten.
::::
``````

(function-am2-python)=
``````{dropdown} 2nd-order Adams–Moulton (trapezoid) formula for an IVP
```{literalinclude} ../python/pkg/FNC/FNC06.py
:filename: am2.py
:start-at: def am2
:end-at: return t
:language: python
:linenos: true
```
::::{admonition} About the code
:class: dropdown
Lines 22-23 define the function $\mathbf{g}$ and call `levenberg` to find the new solution value, using an Euler half-step as its starting value. A robust code would have to intercept the case where `levenberg` fails to converge, but we have ignored this issue for the sake of brevity.
::::
``````

## Examples

```{code-cell} ipython3
exec(open("FNC_init.py").read())
```

### Section 6.1

(demo-basics-first-python)=
``````{dropdown} @demo-basics-first
Let's use `solve_ivp` from `scipy.integrate` to define and solve an initial-value problem for $u'=\sin[(u+t)^2]$ over $t \in [0,4]$, such that $u(0)=-1$.

To create an initial-value problem for $u(t)$, you must supply a function that computes $u'$, an initial value for $u$, and the endpoints of the interval for $t$. The $t$ interval should be defined as `(a,b)`, where at least one of the values is a float.

```{index} ! Julia; ODEProblem, ! Julia; solve
```

```{code-cell}
f = lambda t, u: sin((t + u) ** 2)
tspan = [0.0, 4.0]
u0 = [-1.0]
```

Note above that even though this is a problem for a scalar function $u(t)$, we had to set the initial condition as a "one-dimensional vector."

```{code-cell}
from scipy.integrate import solve_ivp
sol = solve_ivp(f, tspan, u0)
```

The resulting solution object has fields `t` and `y` that contain the values of the independent and dependent variables, respectively; those field names are the same regardless of what we use in our own codes.

```{code-cell}
print("t shape:", sol.t.shape)
print("u shape:", sol.y.shape)
plot(sol.t, sol.y[0, :], "-o")
xlabel("$t$"), ylabel("$u(t)$")
title("Solution of $u' = sin((t+u)^2)$")
```

You can see above that the solution was not computed at enough points to make a smooth graph. There is a way to request output at times of your choosing.

```{code-cell}
sol = solve_ivp(f, tspan, u0, t_eval=linspace(0, 4, 200))
plot(sol.t, sol.y[0, :], "-")
xlabel("$t$"), ylabel("$u(t)$")
title("Solution of $u' = sin((t+u)^2)$")
```

Another option is to enable interpolation to evaluate the solution anywhere after the fact:

```{code-cell}
sol = solve_ivp(f, tspan, u0, dense_output=True)
for t in linspace(0, 4, 6):
    print(f"u({t:.2f}) = {sol.sol(t)[0]:.4f}")
```
``````

(demo-basics-sing-python)=
``````{dropdown} @demo-basics-sing

::::{grid} 1 1 2 2
The equation $u'=(u+t)^2$ gives us some trouble.
:::{card}
It's a good idea to check `sol.success` after calling `solve_ivp`. If it's `False`, the solution may not be reliable. 
:::
::::

```{code-cell}
f = lambda t, u: (t + u) ** 2
sol = solve_ivp(f, [0.0, 1.0], [1.0])
if not sol.success:
    print(sol.message)
```

The warning message we received can mean that there is a bug in the formulation of the problem. But if everything has been done correctly, it suggests that the solution may not exist past the indicated time. This is a possibility in nonlinear ODEs.

```{code-cell}
semilogy(sol.t, sol.y[0, :])
xlabel("$t$")
ylabel("$u(t)$")
title("Blowup in finite time")
```
``````

(demo-basics-cond-python)=
``````{dropdown} @demo-basics-cond
Consider the ODEs $u'=u$ and $u'=-u$. In each case we compute $\partial f/\partial u = \pm 1$, so the condition number bound from {numref}`Theorem %s <theorem-depIC>` is $e^{b-a}$ in both problems. However, they behave quite differently. In the case of exponential growth, $u'=u$, the bound is the actual condition number.

```{code-cell}
t = linspace(0, 3, 200)
u = array([exp(t) * u0 for u0 in [0.7, 1, 1.3]])
plot(t, u.T)
xlabel("$t$")
ylabel("$u(t)$")
title("Exponential divergence of solutions")
```

But with $u'=-u$, solutions actually get closer together with time.

```{code-cell}
t = linspace(0, 3, 200)
u = array([exp(-t) * u0 for u0 in [0.7, 1, 1.3]])
plot(t, u.T)
xlabel("$t$")
ylabel("$u(t)$")
title("Exponential convergence of solutions")
```

In this case the actual condition number is one, because the initial difference between solutions is the largest over all time. Hence, the exponentially growing upper bound $e^{b-a}$ is a gross overestimate.
``````

### Section 6.2
(demo-euler-converge-python)=
``````{dropdown} @demo-euler-converge
We consider the IVP $u'=\sin[(u+t)^2]$ over $0\le t \le 4$, with $u(0)=-1$.

```{code-cell}
f = lambda t, u: sin((t + u) ** 2)
tspan = [0.0, 4.0]
u0 = -1.0
t, u = FNC.euler(f, tspan, u0, 20)

fig, ax = subplots()
ax.plot(t, u[0, :], "-o", label="$n=20$")
xlabel("$t$"), ylabel("$u(t)$")
title("Solution by Euler's method")
legend()
```

We could define a different interpolant to get a smoother picture above, but the derivation of Euler's method assumed a piecewise linear interpolant. We can instead request more steps to make the interpolant look smoother.

```{code-cell}
t, u = FNC.euler(f, tspan, u0, 200)
ax.plot(t, u[0, :], label="$n=200$")
ax.legend()
fig
```

Increasing $n$ changed the solution noticeably. Since we know that interpolants and finite differences become more accurate as $h\to 0$, we should anticipate the same behavior from Euler's method. We don't have an exact solution to compare to, so we will use `solve_ivp` to construct an accurate reference solution.

```{code-cell}
from scipy.integrate import solve_ivp
sol = solve_ivp(f, tspan, [u0], dense_output=True, atol=1e-8, rtol=1e-8)
ax.plot(t, sol.sol(t)[0, :], "--", label="accurate")
ax.legend()
fig
```

Now we can perform a convergence study.

```{code-cell}
n_ = array([int(5 * 10**k) for k in arange(0, 3, 0.5)])
err_ = zeros(6)
results = PrettyTable(["n", "error"])
for j, n in enumerate(n_):
    t, u = FNC.euler(f, tspan, u0, n)
    err_[j] = norm(sol.sol(t)[0, :] - u[0, :], inf)
    results.add_row((n, err_[j]))
print(results)
```

The error is approximately cut by a factor of 10 for each increase in $n$ by the same factor. A log-log plot also confirms first-order convergence. Keep in mind that since $h=(b-a)/n$, it follows that $O(h)=O(n^{-1})$.

```{code-cell}
loglog(n_, err_, "-o", label="results")
plot(n_, 0.5 * (n_ / n_[0])**(-1), "--", label="1st order")
xlabel("$n$"), ylabel("inf-norm error")
title("Convergence of Euler's method")
legend()
```
``````

### Section 6.3
(demo-systems-predator-python)=
``````{dropdown} @demo-systems-predator
We encode the predator–prey equations via a function.

```{code-cell}
def predprey(t, u):
    y, z = u                        # rename for convenience
    s = (y * z) / (1 + beta * y)    # appears in both equations
    return array([y * (1 - alpha * y) - s, -z + s])
```

As before, the ODE function must accept three inputs, `u`, `p`, and `t`, even though in this case there is no explicit dependence on `t`. The second input is used to pass parameters that don't change throughout a single instance of the problem.

To specify the IVP we must also provide the initial condition, which is a 2-vector here, and the interval for the independent variable. These are given in the call to `solve_ivp`.

```{code-cell}
from scipy.integrate import solve_ivp
u0 = array([1, 0.01])
tspan = [0.0, 80.0]
alpha, beta = 0.1, 0.25
sol = solve_ivp(predprey, tspan, u0, dense_output=True)
print(f"solved with {sol.y.shape[1]} time steps")
```

As in scalar problems, the solution object has fields `t` and `y` that contain the values of the independent and dependent variables, respectively. Each row of `y` represents one component of the solution at every time step, and each column of `y` is the entire solution vector at one time step. Since we used `dense_output=True`, there is also a method `sol` that can be used to evaluate the solution at any time. 

```{code-cell}
t = linspace(0, 80, 1200)
u = vstack([sol.sol(t[i]) for i in range(t.size)]).T    # same shape as sol.y
fig, ax = subplots()
ax.plot(t, u[0, :], label="prey")
ax.plot(t, u[1, :], label="predator")
xlabel("$t$"), ylabel("population")
title("Predator-prey solution")
```

We can also use {numref}`Function {number} <function-euler>` to find the solution.

```{code-cell}
t_E, u_E = FNC.euler(predprey, tspan, u0, 800)
ax.scatter(t_E, u_E[0, :], label="prey (Euler)", s=1)
ax.scatter(t_E, u_E[1, :], label="predator (Euler)", s=2)
ax.legend()
fig
```

You can see above that the Euler solution is not very accurate. When the solution has two components, it's common to plot the it in the _phase plane_, i.e., with $u_1$ and $u_2$ along the axes and time as a parameterization of the curve.

```{code-cell}
plot(u[0, :], u[1, :])
xlabel("prey"), ylabel("predator")
title("Predator-prey phase plane")
```
From this plot we can see that the solution approaches a periodic one, which in the phase plane is represented by a closed path.
``````

(demo-systems-coupledpendula-python)=
``````{dropdown} @demo-systems-coupledpendula
Let's implement the coupled pendulums from {numref}`Example {number} <example-systems-coupledpendula>`. The pendulums will be pulled in opposite directions and then released together from rest.

```{code-cell}
def couple(t, u, params):
    gamma, L, k = params
    g = 9.8
    udot = copy(u)
    udot[:2] = u[2:4]
    udot[2] = -gamma * u[2] - (g / L) * sin(u[0]) + k * (u[1] - u[0])
    udot[3] = -gamma * u[3] - (g / L) * sin(u[1]) + k * (u[0] - u[1])
    return udot

u0 = array([1.25, -0.5, 0, 0])
tspan = [0.0, 50.0]
```

::::{grid} 1 1 2 2
First we check the behavior of the system when the pendulums are uncoupled, i.e., when $k=0$.
:::{card}
We use a closure here to pass the fixed parameter values into `couple`.
:::
::::

```{code-cell}
gamma, L, k = 0.01, 0.5, 0.0
du_dt = lambda t, u: couple(t, u, (gamma, L, k))
sol = solve_ivp(du_dt, tspan, u0, t_eval=linspace(0, 50, 1000))
plot(sol.t, sol.y[:2, :].T)    # first two components of solution
xlabel("t"), ylabel("angle")
title("Uncoupled pendulums");
```

You can see that the pendulums swing independently. Because the model is nonlinear and the initial angles are not small, they have slightly different periods of oscillation, and they go in and out of phase.

With coupling activated, a different behavior is seen.

```{code-cell}
k = 0.75    # changes the value in the du_dt closure
sol = solve_ivp(du_dt, tspan, u0, t_eval=linspace(0, 50, 1000))
plot(sol.t, sol.y[:2, :].T)
xlabel("t"), ylabel("angle")
title("Coupled pendulums");
```

The coupling makes the pendulums swap energy back and forth.
``````

### Section 6.4
(demo-rk-converge-python)=
``````{dropdown} @demo-rk-converge
We solve the IVP $u'=\sin[(u+t)^2]$ over $0\le t \le 4$, with $u(0)=-1$. We start by getting a reference solution to validate against.

```{code-cell}
from scipy.integrate import solve_ivp
du_dt = lambda t, u: sin((t + u)**2)
tspan = (0.0, 4.0)
u0 = -1.0
sol = solve_ivp(du_dt, tspan, [u0], dense_output=True, atol=1e-13, rtol=1e-13)
u_ref = sol.sol
```

Now we perform a convergence study of our two Runge–Kutta implementations.

```{code-cell}
n = array([int(2 * 10**k) for k in linspace(0, 3, 7)])
err = {"IE2" : [], "RK4" : []}
results = PrettyTable(["n", "IE2 error", "RK4 error"])
for i in range(len(n)):
    t, u = FNC.ie2(du_dt, tspan, u0, n[i])
    err["IE2"].append( abs(u_ref(4)[0] - u[0][-1]) )
    t, u = FNC.rk4(du_dt, tspan, u0, n[i])
    err["RK4"].append( abs(u_ref(4)[0] - u[0][-1]) )
    results.add_row([n[i], err["IE2"][-1], err["RK4"][-1]])

print(results)
```

The amount of computational work at each time step is assumed to be proportional to the number of stages. Let's compare on an apples-to-apples basis by using the number of $f$-evaluations on the horizontal axis.

```{code-cell}
loglog(2 * n, err["IE2"], "-o", label="IE2")
loglog(4 * n, err["RK4"], "-o", label="RK4")
plot(2 * n, 0.5 * err["IE2"][-1] * (n / n[-1])**(-2), "--", label="2nd order")
plot(4 * n, 0.5 * err["RK4"][-1] * (n / n[-1])**(-4), "--", label="4th order")

xlabel("f-evaluations"),  ylabel("inf-norm error")
legend()
title("Convergence of RK methods");
```

The fourth-order variant is more efficient in this problem over a wide range of accuracy.
``````

### Section 6.5
(demo-adapt-basic-python)=
``````{dropdown} @demo-adapt-basic
Let's run adaptive RK on  $u'=e^{t-u\sin u}$.

```{code-cell}
f = lambda t, u: exp(t - u * sin(u))
t, u = FNC.rk23(f, [0.0, 5.0], [0.0], 1e-5)
scatter(t, u[0, :])
xlabel("$t$"), ylabel("$u(t)$")
title("Adaptive IVP solution")
```

The solution makes a very abrupt change near $t=2.4$. The resulting time steps vary over three orders of magnitude.

```{code-cell}
dt = [t[i + 1] - t[i] for i in range(t.size - 1)]
semilogy(t[:-1], dt)
xlabel("$t$"), ylabel("time step")
title("Adaptive step sizes")
```

If we had to run with a uniform step size to get this accuracy, it would be

```{code-cell}
print(f"min step size was {min(dt):.2e}")
```

On the other hand, the average step size that was actually taken was

```{code-cell}
print(f"mean step size was {mean(dt):.2e}")
```

We took fewer steps by a factor of 1000! Even accounting for the extra stage per step and the occasional rejected step, the savings are clear.

``````

(demo-adapt-sing-python)=
``````{dropdown} @demo-adapt-sing
In {numref}`Demo %s <demo-basics-sing>` we saw an IVP that appears to blow up in a finite amount of time. Because the solution increases so rapidly as it approaches the blowup, adaptive stepping is required even to get close.

```{code-cell}
du_dt = lambda t, u: (t + u)**2
tspan = (0.0, 2.0)
u0 = [1.0]
t, u = FNC.rk23(du_dt, tspan, u0, 1e-5)
```

In fact, the failure of the adaptivity gives a decent idea of when the singularity occurs.

```{code-cell}
semilogy(t, u[0, :])
xlabel("$t$"),  ylabel("$u(t)$")
title("Finite-time blowup")

tf = t[-1]
axvline(x=tf, color='k', linestyle='--', label=f"t = {tf:.6f}")
legend();
```
``````

### Section 6.6
(demo-implicit-ab4-python)=
``````{dropdown} @demo-implicit-ab4
We study the convergence of AB4 using the IVP $u'=\sin[(u+t)^2]$ over $0\le t \le 4$, with $u(0)=-1$. As usual, `solve_ivp` is called to give an accurate reference solution.

```{code-cell}
from scipy.integrate import solve_ivp
du_dt = lambda t, u: sin((t + u)**2)
tspan = (0.0, 4.0)
u0 = [-1.0]
u_ref = solve_ivp(du_dt, tspan, u0, dense_output=True, rtol=1e-13, atol=1e-13).sol
```

Now we perform a convergence study of the AB4 code.

```{code-cell}
n = array([int(4 * 10**k) for k in linspace(0, 3, 7)])
err = []
results = PrettyTable(["n", "AB4 error"])
for i in range(len(n)):
    t, u = FNC.ab4(du_dt, tspan, u0, n[i])
    err.append( abs(u_ref(4)[0] - u[0][-1]) )
    results.add_row([n[i], err[-1]])

print(results)
```

The method should converge as $O(h^4)$, so a log-log scale is appropriate for the errors.

```{code-cell}
loglog(n, err, "-o", label="AB4")
loglog(n, 0.5 * err[-1] * (n / n[-1])**(-4), "--", label="4th order")

xlabel("$n$"),  ylabel("final error")
legend(), title("Convergence of AB4");
```
``````

(demo-implicit-stiff-python)=
``````{dropdown} @demo-implicit-stiff
The following simple ODE uncovers a surprise.

```{code-cell}
f = lambda t, u: u**2 - u**3
u0 = array([0.005])
tspan = [0, 400]
```

We will solve the problem first with the implicit AM2 method using $n=200$ steps.

```{code-cell}
tI, uI = FNC.am2(f, [0.0, 400.0], u0, 200)
fig, ax = subplots()
ax.plot(tI, uI[0], label="AM2")
xlabel("$t$"), ylabel("$y(t)$");
```

So far, so good. Now we repeat the process using the explicit AB4 method.

```{code-cell}
tE, uE = FNC.ab4(f, [0.0, 400.0], u0, 200)
ax.scatter(tE, uE[0], label="AB4")
ax.set_ylim([-4, 2]), ax.legend()
fig
```

Once the solution starts to take off, the AB4 result goes catastrophically wrong.

```{code-cell}
uE[0, 104:111]
```

We hope that AB4 will converge in the limit $h\to 0$, so let's try using more steps.

```{code-cell}
plot(tI, uI[0], color="k", label="AM2")
tE, uE = FNC.ab4(f, [0, 400], u0, 1000)
plot(tE, uE[0], ".-", label="AM4, n=1000")
tE, uE = FNC.ab4(f, [0, 400], u0, 1600)
plot(tE, uE[0], ".-", label="AM4, n=1600")
xlabel("$t$"),  ylabel("$u(t)$")
legend()
```

So AB4, which is supposed to be _more_ accurate than AM2, actually needs something like 8 times as many steps to get a reasonable-looking answer!
``````

### Section 6.7
(demo-zs-LIAF-python)=
``````{dropdown} @demo-zs-LIAF
We'll measure the error at the time $t=1$.

```{code-cell}
du_dt = lambda t, u: u
u_exact = exp
a, b = (0.0, 1.0)

def LIAF(du_dt, tspan, u0, n):
    a, b = tspan
    h = (b - a) / n
    t = linspace(a, b, n+1)
    u = np.tile(np.array(u0), (n+1, 1))
    u[1] = u_exact(t[1])    # use an exact starting value
    f = copy(u)
    f[0] = du_dt(t[0], u[0])
    for i in range(n):
        f[i] = du_dt(t[i], u[i])
        u[i + 1] = -4 * u[i] + 5 * u[i-1] + h * (4 * f[i] + 2 * f[i-1])

    return t, u.T

n = [5, 10, 20, 40, 60]
results = PrettyTable(["n", "error"])
for j in range(5):
    t, u = LIAF(du_dt, [a, b], [1.0], n[j])
    err = abs(u_exact(b) - u[0, -1])
    results.add_row([n[j], err])
print(results)
```

There is no convergence in sight! A graph of the last numerical attempt yields a clue:

```{code-cell}
semilogy(t, abs(u[0]), "-o")
xlabel("$t$"), ylabel("$|u|$")
title("LIAF solution")
```

It's clear that the solution is growing exponentially in time.
``````
