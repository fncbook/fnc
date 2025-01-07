---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
# Chapter 11

## Functions

(function-diffper-python)=
``````{dropdown} Differentiation matrices for periodic end conditions
:open:
```{literalinclude} pkg/FNC/FNC11.py
:filename: diffper.py
:start-at: def diffper
:end-at: return x, Dx, Dxx
:language: python
:linenos: true
```
``````

(function-parabolic-python)=
``````{dropdown} Solution of parabolic PDEs by the method of lines
:open:
```{literalinclude} pkg/FNC/FNC11.py
:filename: parabolic.py
:start-at: def parabolic
:end-at: return x, lambda t
:language: python
:linenos: true
```
``````

## Examples

```{code-cell}
exec(open("FNC_init.py").read())
```

### 11.1 @section-diffusion-blackscholes

(demo-blackscholes-solve-python)=
``````{dropdown} @demo-blackscholes-solve
We consider the Black–Scholes problem for the following parameter values:

```{code-cell}
Smax, T = 8, 6
K = 3
sigma = 0.06
r = 0.08
```

We discretize space and time.

```{code-cell}
m = 200
h = Smax / m
x = h * arange(m + 1)
n = 1000
tau = T / n
t = tau * arange(n + 1)
lamb, mu = tau / h**2, tau / h
```

We set the initial condition and then march forward in time.

```{code-cell}
V = zeros([m + 1, n + 1])
V[:, 0] = maximum(0, x - K)
for j in range(n):
    # Fictitious value from Neumann condition.
    Vfict = 2 * h + V[m-1, j]
    Vj = hstack([V[:, j], Vfict])
    # First row is zero by the Dirichlet condition.
    for i in range(1, m+1):
        diff1 = Vj[i+1] - Vj[i-1]
        diff2 = Vj[i+1] - 2 * Vj[i] + Vj[i-1]
        V[i, j+1] = (
            Vj[i]
            + (lamb * sigma**2 * x[i] ** 2 / 2) * diff2
            + (r * x[i] * mu) / 2 * diff1
            - r * tau * Vj[i]
        )
```

Here is a plot of the solution after every 250 time steps.

```{code-cell}
select_times = 250 * arange(5)
show_times = t[select_times]

for j, col in enumerate(select_times):
    plot(x, V[:, col], label=f"t={show_times[j]:.1f}")

legend()
title("Black-Scholes solution")
xlabel("stock price"),  ylabel("option value");
```

```{index} ! Python; animation
```

Alternatively, here is an animation of the solution.

```{code-cell}
:tags: remove-output
from matplotlib import animation
fig = figure()
ax = fig.add_subplot(autoscale_on=False, xlim=(0, 8), ylim=(0, 6))
ax.grid()
ax.set_title("Black-Scholes solution")

line, = ax.plot([], [], '-', lw=2)
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def animate(j):
    line.set_data(x, V[:, j])
    time_text.set_text(f"t = {t[j]:.2f}")
    return line, time_text

anim = animation.FuncAnimation(
    fig, animate, frames=range(0, n+1, 10), blit=True);
anim.save("figures/black-scholes-6.mp4", fps=30)
close()
```

![Black–Scholes solution](figures/black-scholes-6.mp4)

The results are easy to interpret, recalling that the time variable really means *time until strike*. Say you are close to the option's strike time. If the current stock price is, say, $S=2$, then it's not likely that the stock will end up over the strike price $K=3$, and therefore the option has little value. On the other hand, if presently $S=3$, then there are good odds that the option will be exercised at the strike time, and you will need to pay a substantial portion of the stock price in order to take advantage. As the time to strike increases, there is an expectation that the stock price is more likely to rise somewhat, making the value of the option larger at each fixed $S$. 
``````

(demo-blackscholes-unstable-python)=
``````{dropdown} @demo-blackscholes-unstable
Let's try to do everything the same as in {numref}`Demo {number} <demo-blackscholes-solve>`, but extending the simulation time to $T=8$.

```{code-cell}
T = 8
n = 1000;  tau = T / n
t = tau * arange(n + 1)
lamb, mu = tau / h**2, tau / h

V = zeros([m+1, n+1])
V[:, 0] = maximum(0, x - K)
for j in range(n):
    # Fictitious value from Neumann condition.
    Vfict = 2 * h + V[m - 1, j]
    Vj = hstack([V[:, j], Vfict])
    # First row is zero by the Dirichlet condition.
    for i in range(1, m + 1):
        diff1 = Vj[i + 1] - Vj[i - 1]
        diff2 = Vj[i + 1] - 2 * Vj[i] + Vj[i - 1]
        V[i, j + 1] = (
            Vj[i]
            + (lamb * sigma**2 * x[i] ** 2 / 2) * diff2
            + (r * x[i] * mu) / 2 * diff1
            - r * tau * Vj[i]
        )

select_times = 250 * arange(5)
show_times = t[select_times]

for j, col in enumerate(select_times):
    plot(x, V[:, col], label=f"t={show_times[j]:.1f}")

legend()
title("Black-Scholes solution")
xlabel("stock price")
ylim([0, 6]);  ylabel("option value")
```

```{code-cell}
:tags: remove-output
from matplotlib import animation
fig = figure()
ax = fig.add_subplot(autoscale_on=False, xlim=(0, 8), ylim=(0, 6))
ax.grid()
ax.set_title("Black-Scholes solution...?")

line, = ax.plot([], [], '-', lw=2)
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def animate(j):
    line.set_data(x, V[:, j])
    time_text.set_text(f"t = {t[j]:.2f}")
    return line, time_text

anim = animation.FuncAnimation(
    fig, animate, frames=range(0, n+1, 10), blit=True);
anim.save("figures/black-scholes-8.mp4", fps=30)
close()
```

![Trouble in Black–Scholes solution](figures/black-scholes-8.mp4)

This so-called solution is nonsense!
``````

### 11.2 @section-diffusion-methodlines

(demo-methodlines-heatFE-python)=
``````{dropdown} @demo-methodlines-heatFE
Let's implement the method of {numref}`Example {number} <example-methodlines-heatFE>` with second-order space semidiscretization.

```{code-cell}
m = 100
x, Dx, Dxx = FNC.diffper(m, [0, 1])

tfinal = 0.15  
n = 2400                 # number of time steps  
tau = tfinal / n         # time step
t = tau * arange(n+1)    # time values
```

Next we set an initial condition. It isn't mathematically periodic, but the end values and derivatives are so small that for numerical purposes it may as well be.

```{code-cell}
U = zeros([m, n+1])
U[:, 0] = exp(-60 * (x - 0.5) ** 2)
plot(x, U[:, 0])
xlabel("x");  ylabel("u(x,0)")
title("Initial condition");
```

The Euler time stepping simply multiplies the solution vector by the constant matrix in {eq}`Eulerxx` at each time step. Since that matrix is sparse, we will declare it as such, even though the run-time savings may not be detectable for this small value of $m$.

```{code-cell}
import scipy.sparse as sp
I = sp.eye(m)
A = I + tau * sp.csr_array(Dxx)
for j in range(n):
    U[:, j+1] = A @ U[:, j]

plot(x, U[:, :31:10])
ylim([-0.25, 1])
xlabel("$x$");  ylabel("$u(x,t)$")
legend([f"$t={tj:.2e}$" for tj in t[:60:20]])
title("Heat equation by forward Euler");
```

You see above that things seem to start well, with the initial peak widening and shrinking. But then there is a nonphysical growth in the solution.

```{code-cell}
:tags: hide-input
from matplotlib import animation
fig = figure()
ax = fig.add_subplot(autoscale_on=False, xlim=(0, 1), ylim=(-1, 2))
ax.grid()

line, = ax.plot([], [], '-', lw=2)
ax.set_title("Heat equation by forward Euler")
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def animate(j):
    line.set_data(x, U[:, j])
    time_text.set_text(f"t = {t[j]:.2e}")
    return line, time_text

anim = animation.FuncAnimation(
    fig, animate, frames=range(0, 100), blit=True)
anim.save("figures/diffusionFE.mp4", fps=30)
close()
```

![Instability in Euler solution](figures/diffusionFE.mp4)

The growth in norm is exponential in time.

```{code-cell}
M = abs(U).max(axis=0)  # max in each column
semilogy(t, M)
xlabel("$t$");  ylabel("$\\max_x |u(x,t)|$")
title("Nonphysical growth");
```
``````

(demo-methodlines-heatBE-python)=
``````{dropdown} @demo-methodlines-heatBE
Now we apply backward Euler to the heat equation. Mathematically this means multiplying by the *inverse* of a matrix, but we interpret that numerically as a linear system solution. We will reuse the setup from {numref}`Demo {number} <demo-methodlines-heatFE>`. 

```{code-cell}
from scipy.sparse.linalg import spsolve
B = sp.csr_matrix(I - tau * Dxx)
for j in range(n):
    U[:, j + 1] = spsolve(B, U[:, j])

plot(x, U[:, ::500])
xlabel("$x$")
ylabel("$u(x,t)$")
legend([f"$t={tj:.2g}$" for tj in t[::500]])
title("Heat equation by backward Euler");
```

```{code-cell}
:tags: hide-input
fig = figure()
ax = fig.add_subplot(autoscale_on=False, xlim=(0, 1), ylim=(-0.25, 1))
ax.grid()

line, = ax.plot([], [], '-', lw=2)
ax.set_title("Backward Euler")
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

anim = animation.FuncAnimation(
    fig, animate, frames=range(0, n+1, 20), blit=True)
anim.save("figures/diffusionBE.mp4")
close()
```

![Stable Backward Euler solution](figures/diffusionBE.mp4)

This solution looks physically plausible, as the large concentration in the center diffuses outward until the solution is essentially constant. Observe that the solution remains periodic in space for all time.
``````
(demo-methodlines-auto-python)=
``````{dropdown} @demo-methodlines-auto
We set up the semidiscretization and initial condition in $x$ just as before.

```{code-cell}
m = 100
x, Dx, Dxx = FNC.diffper(m, [0, 1])
u0 = exp(-60 * (x - 0.5) ** 2)
```

Now, however, we apply a standard solver using `solve_ivp` to the initial-value problem $\mathbf{u}'=\mathbf{D}_{xx}\mathbf{u}$.

```{code-cell}
from scipy.integrate import solve_ivp
tfinal = 0.05
f = lambda t, u: Dxx @ u
sol = solve_ivp(f, [0, tfinal], u0, method="RK45", dense_output=True)

t = linspace(0, 0.05, 5)
plot(x, sol.sol(t))
xlabel("$x$"),  ylabel("$u(x,t)$")
legend([f"$t={tj:.4g}$" for tj in t])
title("Heat equation by RK45");
```

The solution appears to be correct. But the number of time steps that were selected automatically is surprisingly large, considering how smoothly the solution changes.

```{code-cell}
print(f"RK45 took {len(sol.t) - 1} steps")
```

Now we apply a different solver called `BDF`.

```{code-cell}
sol = solve_ivp(f, [0, tfinal], u0, method="BDF")
print(f"BDF took {len(sol.t) - 1} steps")
```

The number of steps selected was reduced by a factor of 20!
``````

### 11.3 @section-diffusion-absstab
(demo-absstab-regions-python)=
``````{dropdown} @demo-absstab-regions
Euler and Backward Euler time-stepping methods were used to solve $\mathbf{u}'=\mathbf{D}_{xx}\mathbf{u}$.

```{code-cell}
m = 40
Dxx = FNC.diffper(m, [0, 1])[2]
```

The eigenvalues of this matrix are real and negative:

```{code-cell}
from scipy.linalg import eigvals
lamb = eigvals(Dxx)
plot(real(lamb), imag(lamb), "o")
xlabel("Re $\\lambda$")
ylabel("Im $\\lambda$")
title("Eigenvalues");
```

The Euler method is absolutely stable in the region $|\zeta+1| \le 1$ in the complex plane:

```{code-cell}
:tags: hide-input
phi = 2 * pi * arange(361) / 360
z = exp(1j * phi) - 1  # unit circle shifted to the left by 1

fill(real(z), imag(z), color=(0.8, 0.8, 1))
xlabel("Re $\\lambda$")
ylabel("Im $\\lambda$")
axis("equal")
title("Stability region");
```

In order to get inside this region, we have to find $\tau$ such that $\lambda \tau > -2$ for all eigenvalues $\lambda$. This is an upper bound on $\tau$. 

```{code-cell}
lambda_min = min(real(lamb))
max_tau = -2 / lambda_min
print(f"predicted max time step is {max_tau:.3e}")
```

Here we plot the resulting values of $\zeta=\lambda \tau$. 

```{code-cell}
zeta = lamb * max_tau
fill(real(z), imag(z), color=(0.8, 0.8, 1))
plot(real(zeta), imag(zeta), "o")
xlabel("Re $\\lambda$")
ylabel("Im $\\lambda$")
axis("equal")
title("Stability region and $\zeta$ values");
```

In backward Euler, the region is $|\zeta-1|\ge 1$. Because they are all on the negative real axis, all of the $\zeta$ values will fit no matter what $\tau$ is chosen.

```{code-cell}
:tags: hide-input
fill([-6, 6, 6, -6], [-6, -6, 6, 6], color=(0.8, 0.8, 1))
z = exp(1j * phi) + 1
# unit circle shifted right by 1
fill(real(z), imag(z), color="w")

plot(real(zeta), imag(zeta), "o")
axis([-4, 2, -3, 3])
axis("equal")
xlabel("Re $\\lambda$")
ylabel("Im $\\lambda$")
title("Stability region and $\\zeta$ values");
```
``````


### 11.4 @section-diffusion-stiffness
(demo-stiffness-oregon-python)=
``````{dropdown} @demo-stiffness-oregon
In {numref}`Example {number} <example-stiffness-oregon>` we derived a Jacobian matrix for the Oregonator model. Here is a numerical solution of the ODE.

```{code-cell}
:tags: hide-input
from scipy.integrate import solve_ivp
q, s, w = (8.375e-6, 77.27, 0.161)

def ode(t, u):
    return array(
        [
            s * (u[1] - u[0] * u[1] + u[0] - q * u[0]**2),
            (-u[1] - u[0] * u[1] + u[2]) / s,
            w * (u[0] - u[2]),
        ]
    )

u0 = array([1.0, 2.0, 3.0])
tspan = (0, 500)
start = timer()
sol = solve_ivp(ode, tspan, u0, method="BDF")
semilogy(sol.t, sol.y.T)
xlabel("$t$"),  ylabel("$u(t)$")
title("Oregonator");
```

At each value of the numerical solution, we can compute the eigenvalues of the Jacobian. Here we plot all of those eigenvalues in the complex plane.

```{code-cell}
:tags: hide-input
J = lambda u: array(
    [
        [-s * (u[1] + 1 - 2 * q * u[0]), s * (1 - u[0]), 0],
        [-u[1] / s, (-1 - u[0]) / s, 1 / s],
        [w, 0, -w],
    ]
)

from scipy.linalg import eigvals

lamb = array([eigvals(J(u)) for u in sol.y.T])
ax = figure().add_subplot(projection='3d')
for i in range(3):
    ax.plot(real(lamb[:, i]), imag(lamb[:, i]), sol.t, ".")
ax.set_xlabel("Re $\\lambda$")
ax.set_ylabel("Im $\\lambda$")
ax.set_zlabel("$t$")
ax.set_title("Oregonator eigenvalues")
```

You can see that there is one eigenvalue that ranges over a wide portion of the negative real axis and dominates stability considerations.
``````

(demo-stiffness-explicit-python)=
``````{dropdown} @demo-stiffness-explicit
The `BDF` solver is good for stiff problems and needs few time steps to solve the Oregonator from {numref}`Demo {number} <demo-stiffness-oregon>`.

```{code-cell}
tspan = (0, 25)
start = timer()
sol = solve_ivp(ode, tspan, u0, method="BDF")
print(f"stiff solver took {timer() - start:.3f} seconds with {len(sol.t) - 1} time steps")
```

But if we apply {numref}`Function {number} <function-rk23>` to the problem, the step size will be made small enough to cope with the large negative eigenvalue. 

```{code-cell}
start = timer()
t, u = FNC.rk23(ode, tspan, u0, 1e-6)
print(f"rk23 solver took {timer() - start:.3f} seconds with {len(t) - 1} time steps")
```

Starting from the eigenvalues of the Jacobian matrix, we can find an effective $\zeta(t)$ by multiplying with the local time step size. The values of $\zeta(t)$ for each time level are plotted below and color coded by component of the diagonalized system.

```{code-cell}
:tags: hide-input
zeta = zeros([len(t)- 1, 3]) + 0j    # complex array
for i in range(len(t) - 1):
    dt = t[i+1] - t[i]
    lamb = eigvals(J(u[:, i]))
    zeta[i] = lamb * dt
plot(real(zeta), imag(zeta), ".")
axis("equal")
xlabel("Re $\\zeta$")
ylabel("Im $\\zeta$")
title("Oregonator stability")
```

Roughly speaking, the $\zeta$ values stay within or close to the RK2 stability region in {numref}`figure-stabreg_bd_rk`. Momentary departures from the region are possible, but time stepping repeatedly in that situation would cause instability. 

``````

### 11.5 @section-diffusion-boundaries
(demo-boundaries-heat-python)=
``````{dropdown} @demo-boundaries-heat
First, we define functions for the PDE and each boundary condition.

```{code-cell}
phi = lambda t, x, u, ux, uxx: uxx
ga = lambda u, ux: u
gb = lambda u, ux: u - 2
```

Our next step is to write a function to define the initial condition. This one satisfies the boundary conditions exactly.

```{code-cell}
init = lambda x: 1 + sin(pi * x/2) + 3 * (1 - x**2) * exp(-4*x**2)
```

Now we can use {numref}`Function {number} <function-parabolic>` to solve the problem.

```{code-cell}
x, u = FNC.parabolic(phi, (-1, 1), 60, ga, gb, (0, 0.75), init)

for t in arange(0,0.5,0.1):
    plot(x, u(t), label=f"t={t:.2f}")
xlabel("$x$"),  ylabel("$u(x,t)$")
legend()
title("Heat equation with Dirichlet boundaries");
```

```{code-cell}
:tags: hide-input
from matplotlib import animation
fig = figure()
ax = fig.add_subplot(autoscale_on=False, xlim=(-1, 1), ylim=(0, 4.2))

line, = ax.plot([], [], '-')
ax.set_title("Heat equation with Dirichlet boundaries")
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def animate(t):
    line.set_data(x, u(t))
    time_text.set_text(f"t = {t:.2e}")
    return line, time_text

anim = animation.FuncAnimation(
    fig, animate, frames=linspace(0, 0.75, 201), blit=True)
anim.save("figures/boundaries-heat.mp4", fps=30)
close()
```
![Heat equation with Dirichlet boundaries](figures/boundaries-heat.mp4)

``````


(demo-boundaries-bratu-python)=
``````{dropdown} @demo-boundaries-bratu

```{code-cell}
phi = lambda t, x, u, ux, uxx: u**2 + uxx
ga = lambda u, ux: u
gb = lambda u, ux: ux
init = lambda x: 400 * x**4 * (1 - x)**2
x, u = FNC.parabolic(phi, (0, 1), 60, ga, gb, (0, 0.1), init);
```

```{code-cell}
:tags: hide-input
fig = figure()
ax = fig.add_subplot(autoscale_on=False, xlim=(0, 1), ylim=(0, 10))

line, = ax.plot([], [], '-')
ax.set_title("Heat equation with source")
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

anim = animation.FuncAnimation(
    fig, animate, frames=linspace(0, 0.1, 101), blit=True)
anim.save("figures/boundaries-source.mp4", fps=30)
close()
```

![Heat equation with source](figures/boundaries-source.mp4)
``````

(demo-boundaries-bs-python)=
``````{dropdown} @demo-boundaries-bs

```{code-cell}
K = 3;  sigma = 0.06;  r = 0.08;  Smax = 8;
phi = lambda t, x, u, ux, uxx: sigma**2/2 * (x**2 * uxx) + r*x*ux - r*u
ga = lambda u, ux: u
gb = lambda u, ux: ux - 1
```

```{code-cell}
u0 = lambda x: maximum(0, x - K)
x, u = FNC.parabolic(phi, (0, Smax), 80, ga, gb, (0, 15), u0);
```

```{code-cell}
:tags: hide-input
fig = figure()
ax = fig.add_subplot(autoscale_on=False, xlim=(0, Smax), ylim=(-0.5, 8))

line, = ax.plot([], [], '-')
ax.set_title("Black–Scholes equation with boundaries")
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

anim = animation.FuncAnimation(
    fig, animate, frames=linspace(0, 15, 151), blit=True)
anim.save("figures/boundaries-bs.mp4", fps=30)
close()
```
![Black–Scholes equation with boundaries](figures/boundaries-bs.mp4)

Recall that $u$ is the value of the call option, and time runs backward from the strike time. The longer the horizon, the more value the option has due to anticipated growth in the stock price.
``````