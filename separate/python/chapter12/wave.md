---
numbering:
  enumerator: 12.4.%s
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

(section-advection-wave)=
# The wave equation

A close cousin of the advection equation is the {term}`wave equation`.

```{index} ! wave equation
```

````{prf:definition} Wave equation
```{math}
:label: wavepde
    u_{tt} - c^2 u_{xx} = 0.
```
````

The wave equation describes the propagation of waves, such as for sound or light. This is our first PDE having a second derivative in time. Consequently, we should expect to have two initial conditions:

```{math}
:label: waveIC
\begin{split}
u(x,0) &= f(x),  \\
u_t(x,0) &= g(x).
\end{split}
```

There are also two space derivatives, so if boundaries are present, two boundary conditions are required.

In the absence of boundaries, $u(x,t)=\phi(x-ct)$ is a solution of {eq}`wavepde`, just as it is for the advection equation. Now, though, so is $u(x,t)=\phi(x+c t)$ for any twice-differentiable $\phi$ (see @problem-wave-twodir).

```{note}
The wave equation in one dimension supports advection in both directions simultaneously.
```

<!-- We will use $x \in [0,1]$ as the domain.  -->

## First-order system

In order to apply semidiscretization with the standard IVP solvers that we have encountered, we must recast {eq}`wavepde` as a first-order system in time. Using our typical methodology, we would define $y=u_t$ and derive

```{math}
:label: wavefirst1
\begin{split}
  u_t &= y, \\
  y_t &= c^2 u_{xx}.
\end{split}
```

This is a valid approach. However, the wave equation has another, less obvious option for transforming to a first-order system:

```{math}
:label: wavefirst2
\begin{split}
    u_t &= z_x, \\
    z_t &= c^2 u_{x}.
\end{split}
```

```{index} Maxwell's equations
```

This second form is appealing in part because it's equivalent to [Maxwell's equations](wiki:Maxwell%27s_equations) for electromagnetism. In this form, we typically replace the velocity initial condition in {eq}`waveIC` with a condition on $z$,

```{math}
:label: waveIC2
\begin{split}
u(x,0) &= f(x),  \\
z(x,0) &= g(x).
\end{split}
```

This alternative to {eq}`waveIC` is useful in many electromagnetic applications, where $u$ and $z$ represent the electric and magnetic fields, respectively.

## Semidiscretization

Because waves travel in both directions, there is no preferred upwind direction. This makes a centered finite difference in space appropriate. Before application of the boundary conditions, semidiscretization of {eq}`wavefirst2` leads to

```{math}
:label: waveMOL
  \begin{bmatrix}
    \mathbf{u}'(t) \\[2mm]  \mathbf{z}'(t)
  \end{bmatrix}
  =
  \begin{bmatrix}
    \boldsymbol{0} & \mathbf{D}_x \\[2mm] c^2 \mathbf{D}_x & \boldsymbol{0}
  \end{bmatrix}
  \begin{bmatrix}
    \mathbf{u}(t) \\[2mm] \mathbf{z}(t)
  \end{bmatrix}.
```

We now suppose that the domain is $[0,1]$ and impose the Dirichlet conditions  

```{math}
:label: waveBC
u(0,t) = u(1,t) = 0, \qquad t \ge 0.
```

The boundary conditions {eq}`waveBC` suggest that we should remove both of the end values of $\mathbf{u}$ from the discretization, but retain all the $\mathbf{z}$ values. We use $\mathbf{w}(t)$ to denote the vector of all the unknowns in the semidiscretization. If the nodes are numbered $x_0,\ldots,x_m$, then we have

```{math}
:label: wavew
\mathbf{w}(t) =  \begin{bmatrix}
  u_1(t) \\ \vdots \\ u_{m-1}(t) \\ z_0(t) \\ \vdots \\ z_m(t)
\end{bmatrix} =  \begin{bmatrix}
  \mathbf{v}(t) \\[1mm] \mathbf{z}(t) 
\end{bmatrix} \in \mathbb{R}^{2m}.
```

When computing $\mathbf{w}'(t)$, we extract the $\mathbf{v}$ and $\mathbf{z}$ components, and we define two helper functions: `extend`, which pads the $\mathbf{v}$ component with the zero end values, and `chop`, which deletes them from $\mathbf{u}$ to give $\mathbf{v}$.

::::{prf:example} Wave equation with boundaries
:label: demo-wave-boundaries

We solve the wave equation {eq}`wavefirst2` with speed $c=2$, subject to {eq}`waveBC` and initial conditions {eq}`waveIC2`.


```{code-cell}
m = 200
x, Dx, Dxx = FNC.diffmat2(m, [-1, 1])
```

The boundary values of $u$ are given to be zero, so they are not unknowns in the ODEs. Instead they are added or removed as necessary.

```{code-cell}
chop = lambda u: u[1:-1]
extend = lambda v: hstack([0, v, 0])
```

The following function computes the time derivative of the system at interior points.

```{code-cell}
def dw_dt(t, w):
    u = extend(w[:m-1])
    z = w[m-1:]
    du_dt = Dx @ z
    dz_dt = c**2 * (Dx @ u)
    return hstack([chop(du_dt), dz_dt])
```

Our initial condition is a single hump for $u$.

```{code-cell}
u_init = exp(-100 * x**2)
z_init = -u_init
w_init = hstack([chop(u_init), z_init])
```

Because the wave equation is hyperbolic, we can use a nonstiff explicit solver.

```{code-cell}
from scipy.integrate import solve_ivp
c = 2
sol = solve_ivp(dw_dt, (0, 2), w_init, dense_output=True)
u = lambda t: extend(sol.sol(t)[:m-1])   # extract the u component
```

We plot the results for the original $u$ variable only. Its interior values are at indices `1:m-1` of the composite $\mathbf{w}$ variable.

```{code-cell}
t = linspace(0, 2, 80)
U = [u(tj) for tj in t]
contour(x, t, U, levels=24, cmap="RdBu", vmin=-1, vmax=1)
xlabel("$x$"),  ylabel("$t$")
title("Wave equation with boundaries");
```

```{code-cell}
:tags: hide-input
from matplotlib import animation
fig, ax = subplots()
curve = ax.plot(x, u_init)[0]
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
ax.set_xlabel("$x$")
ax.set_ylabel("$u(x,t)$")
ax.set_ylim(-1.05, 1.05)
ax.set_title("Wave equation with boundaries")

def snapshot(t):
    curve.set_ydata(u(t))
    time_text.set_text(f"t = {t:.2f}")

anim = animation.FuncAnimation(
    fig, snapshot, frames=linspace(0, 2, 161)
    )
anim.save("wave-boundaries.mp4", fps=30)
close()
```

![Wave equation with boundaries](wave-boundaries.mp4)

The original hump breaks into two pieces of different amplitudes, each traveling with speed $c=2$. They pass through one another without interference. When a hump encounters a boundary, it is perfectly reflected, but with inverted shape. At time $t=2$, the solution looks just like the initial condition.


::::

## Variable speed

An interesting situation is when the wave speed $c$ in @wavepde changes discontinuously, as when light passes from one material into another. For this we must replace the term $c^2$ in {eq}`waveMOL` with the matrix $\operatorname{diag}\bigl(c^2(x_0),\ldots,c^2(x_m)\bigr)$.

::::{prf:example} Wave equation with variable speed
:label: demo-wave-speed

We now use a wave speed that is discontinuous at $x=0$; to the left, $c=1$, and to the right, $c=2$.

The variable wave speed is set to be re-used

```{code-cell}
m = 120
x, Dx, Dxx = FNC.diffcheb(m, [-1, 1])
c = 1 + (sign(x) + 1) / 2
u_init = exp(-100 * x**2)
z_init = -u_init
w_init = hstack([chop(u_init), z_init])

sol = solve_ivp(dw_dt, (0, 5), w_init, dense_output=True, method="Radau")
u = lambda t: extend(sol.sol(t)[:m-1])
```

```{code-cell}
t = linspace(0, 5, 150)
U = [u(tj) for tj in t]
contour(x, t, U, levels=24, cmap="RdBu", vmin=-1, vmax=1)
xlabel("$x$"),  ylabel("$t$")
title("Wave equation with variable speed");
```

```{code-cell}
:tags: hide-input
fig, ax = subplots()
curve = ax.plot(x, u_init)[0]
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
ax.set_xlabel("$x$")
ax.set_ylabel("$u(x,t)$")
ax.set_ylim(-1.05, 1.05)
ax.set_title("Wave equation with variable speed")

anim = animation.FuncAnimation(
    fig, snapshot, frames=linspace(0, 5, 251)
    )
anim.save("wave-speed.mp4", fps=30)
close()
```

![Wave equation with variable speed](wave-speed.mp4)

Each pass through the interface at $x=0$ generates a reflected and transmitted wave. By conservation of energy, these are both smaller in amplitude than the incoming bump.

::::

## Exercises

``````{exercise}
:label: problem-wave-maxwell
✍ Consider the Maxwell equations {eq}`wavefirst2` with smooth solution $u(x,t)$ and $z(x,t)$.

**(a)** Show that $u_{tt} = c^2 u_{xx}$.

**(b)** Show that $z_{tt} = c^2 z_{xx}$.
``````

``````{exercise}
:label: problem-wave-twodir
✍ Suppose that $\phi(s)$ is any twice-differentiable function.

**(a)** Show that $u(x,t) = \phi(x-c t)$ is a solution of $u_{tt}=c^2u_{xx}$. (As in the advection equation, this is a traveling wave of velocity $c$.)

**(b)** Show that $u(x,t) = \phi(x+c t)$ is another solution of $u_{tt}=c^2u_{xx}$. (This is a traveling wave of velocity $-c$.)
``````

```{index} D'Alembert's solution
```

``````{exercise}
:label: problem-wave-dalembert
✍ Show that the following is a solution to the wave equation $u_{tt}=c^2u_{xx}$ with initial and boundary conditions {eq}`waveBC` and {eq}`waveIC`:

$$
u(x,t) = \frac{1}{2} \left[ f(x-ct)+f(x+ct)\right] + \frac{1}{2c} \int_{x-ct}^{x+ct} g(\xi) \, d\xi
$$

This is known as [D'Alembert's formula](wiki:D%27Alembert%27s_formula).
``````

``````{exercise}
:label: problem-wave-neumann
⌨ Suppose the wave equation has homogeneous Neumann conditions on $u$ at each boundary instead of Dirichlet conditions. Using the Maxwell formulation {eq}`wavefirst2`, we have $z_t=c^2u_x$, so $z$ is constant in time at each boundary. Therefore, the endpoint values of $\mathbf{z}$ can be taken from the initial condition and removed from the ODE, while the entire $\mathbf{u}$ vector is now part of the ODE. 

Modify @demo-wave-boundaries to solve the PDE there with Neumann instead of Dirichlet conditions, and make an animation or space-time portrait of the solution. In what major way is it different from the Dirichlet case?
``````

```{index} sine–Gordon equation
```

``````{exercise}
:label: problem-wave-sinegordon
The nonlinear [sine–Gordon equation](wiki:Sine-Gordon_equation) $u_{tt}-u_{xx}= - \sin u$ is a model of multiple pendulums hanging from a single string. It has some interesting solutions.

**(a)** ✍ Write the equation as a first-order system in the variables $u$ and $y=u_t$.

**(b)** ⌨ Assume periodic end conditions on $[-10,10]$ and discretize at $m=200$ points. Let 

$$
u(x,0) = 4\operatorname{arctan}\left( \frac{1}{\cosh(x/\sqrt{2})} \right) 
$$ 

and $u_t(x,0) = 0$. Solve the system using a nonstiff solver for $t \in [0,50],$ and make an animation of the solution $u(x,t)$. (This type of solution is called a [breather](wiki:Breather).)
``````

```{index} beam equation
```

``````{exercise}
:label: problem-wave-beam

The vibrations of a stiff beam, such as a ruler, are governed by the PDE $u_{tt}=-u_{xxxx}$, where $u(x,t)$ is the deflection of the beam at position $x$ and time $t$. This is called the [dynamic beam equation](Euler–Bernoulli_beam_theory).

**(a)** ✍ Show that solutions of the first-order system

\begin{align*}
u_t &= v_{xx}, \\
v_t &= -u_{xx}.
\end{align*}

also satisfy the dynamic beam equation. 

**(b)** ⌨ Assuming periodic end conditions on $[-1,1]$, use $m=100$, let $u(x,0) =\exp(-24x^2)$, $v(x,0) = 0$, and simulate the solution of the beam equation for $0\le t \le 1$ using {numref}`Function {number} <function-am2>` with $n=100$ time steps. Make a plot or animation of the solution.
``````
