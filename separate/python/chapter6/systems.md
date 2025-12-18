---
numbering:
  enumerator: 6.3.%s
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

(section-ivp-systems)=

# IVP systems

Few applications involve an initial-value problem with just a single dependent variable. Usually there are multiple unknowns and a system of equations to define them.  

````{prf:example}
:label: example-systems-predprey
Variations of the following model are commonly seen in ecology:

```{math}
:label: predprey
  \begin{split}
    \frac{d y}{d t} &= y(1-\alpha y) - \frac{yz}{1+\beta y}, \\
    \frac{d z}{d t} &= -z + \frac{yz}{1+\beta y},
  \end{split}
```

```{index} predator–prey model
```

where $\alpha$ and $\beta$ are positive constants. This model is a system of two differential equations for the unknown functions $y(t)$, which could represent a prey species or susceptible host, and $z(t)$, which could represent a predator species or infected population.  We refer to this as a **predator–prey model**. Both of the equations involve both of the unknowns, with no clear way to separate them.

We can pack the two dependent variables $y$ and $z$ into a vector-valued function of time, $\mathbf{u}(t)$, writing

```{math}
\begin{split}
  u_1'(t) &= f_1(t,\mathbf{u}) =  u_1(1-au_1) - \frac{u_1 u_2}{1+bu_1},\\
  u_2'(t) &= f_2(t,\mathbf{u}) = -u_2 + \frac{u_1 u_2}{1+bu_1},
\end{split}
```

and identifying $u_1=y$, $u_2=z$.
````

We now upgrade our IVP definition, {numref}`Definition {number} <definition-scalarivp>`.

::::{prf:definition} Vector-valued IVP / IVP system
:label: definition-vectorivp
A vector-valued first-order **initial-value problem** (IVP) is
  
```{math}
:label: IVPsys
  \mathbf{u}'(t) = \mathbf{f}\bigl(t,\mathbf{u}(t)\bigr), \qquad a \le t \le b, \qquad
  \mathbf{u}(a)=\mathbf{u}_0,
```

where $\mathbf{u}(t)$ is $m$-dimensional. If $\mathbf{f}(t,\mathbf{u})=\mathbf{A}(t)\mathbf{u}(t)+ \mathbf{g}(t)$, the differential equation is **linear**; otherwise, it is **nonlinear**.
::::

We use the terms *IVP system* and *vector-valued IVP* interchangeably; a system of scalar IVPs can be put into the form of {eq}`IVPsys` by appropriate definitions of $\mathbf{u}$ and $\mathbf{f}$, as shown in {numref}`Example {number} <example-systems-predprey>`.

## Numerical solutions

```{index} Euler's method
```

The generalization of any scalar IVP solver to handle systems is straightforward. Consider Euler's method, which in system form becomes

```{math}
:label: eulersys
  \begin{split}
    \mathbf{u}_{i+1} &= \mathbf{u}_i + h\,\mathbf{f}(t_i,\mathbf{u}_i), \qquad i=0,\ldots,n-1.
  \end{split}
```

The vector difference equation {eq}`eulersys` is just Euler's formula applied simultaneously to each component of the ODE system. Because operations such as addition and multiplication translate easily from scalars to vectors, {numref}`Function {number} <function-euler>` that we wrote for scalar IVPs works for systems as well. Practically speaking, the only changes that must be made are that the initial condition and the ODE function have to be coded to use vectors.

::::{prf:example} Predator-prey model
:label: demo-systems-predator

We encode the predator–prey equations via a function.

```{code-cell}
def predprey(t, u):
    y, z = u                        # rename for convenience
    s = (y * z) / (1 + beta * y)    # appears in both equations
    return array([y * (1 - alpha * y) - s, -z + s])
```

As before, the ODE function must accept three inputs, `u`, `p`, and `t`, even though in this case there is no explicit dependence on `t`. The second input is used to pass parameters that don't change throughout a single instance of the problem.

To specify the IVP, we must also provide the initial condition, which is a 2-vector here, and the interval for the independent variable. These are given in the call to `solve_ivp`.

```{code-cell}
from scipy.integrate import solve_ivp
u0 = array([1, 0.01])
tspan = [0.0, 60.0]
alpha, beta = 0.1, 0.25
sol = solve_ivp(predprey, tspan, u0, dense_output=True)
print(f"solved with {sol.y.shape[1]} time steps")
```

As in scalar problems, the solution object has fields `t` and `y` that contain the values of the independent and dependent variables, respectively. Each row of `y` represents one component of the solution at every time step, and each column of `y` is the entire solution vector at one time step. Since we used `dense_output=True`, there is also a method `sol` that can be used to evaluate the solution at any time. 

```{code-cell}
t = linspace(0, 60, 1200)
u = vstack([sol.sol(t[i]) for i in range(t.size)]).T    # same shape as sol.y
fig, ax = subplots()
ax.plot(t, u[0, :], label="prey")
ax.plot(t, u[1, :], label="predator")
xlabel("$t$"), ylabel("population")
title(("Predator-prey solution"));
```

We can also use {numref}`Function {number} <function-euler>` to find the solution.

```{code-cell}
t_E, u_E = FNC.euler(predprey, tspan, u0, 1001)
ax.scatter(t_E, u_E[0, :], label="prey (Euler)", s=1)
ax.scatter(t_E, u_E[1, :], label="predator (Euler)", s=2)
ax.legend()
fig
```

Yotice above that the accuracy of the Euler solution deteriorates rapidly.

When the solution has two components, it's common to plot the it in the _phase plane_, i.e., with $u_1$ and $u_2$ along the axes and time as a parameterization of the curve.

```{code-cell}
plot(u[0, :], u[1, :])
xlabel("prey, $y(t)$"), ylabel("predator, $z(xt)$")
title("Predator-prey solution in the phase plane");
```

As time progresses, the point in the phase plane spirals inward toward a limiting closed loop called a *limit cycle* representing a periodic solution:

```{index} Python; animation
```

```{code-cell}
:tags: [hide-input, remove-output]
from matplotlib import animation
fig, ax = subplots()
curve_ = ax.plot(u[0, 0], u[1, 0])[0]
head_ = ax.plot(u[0, 0], u[1, 0], 'ro', markersize=6)[0]

ax.set_xlabel("prey, $y(t)$")
ax.set_ylabel("predator, $z(xt)$")
ax.set_xlim(0, 9)
ax.set_ylim(0, 6)
ax.set_title("Predator-prey solution in the phase plane")

predator = []
prey = []
def snapshot(t):
    y, z = sol.sol(t)
    prey.append(y)
    predator.append(z)
    head_.set_data([y], [z])
    curve_.set_data(prey, predator)

anim = animation.FuncAnimation(
    fig, snapshot, frames=linspace(0, 60, 801)
    )
anim.save("figures/predator-prey.mp4", fps=24)
close()
```

![Predator–prey solution](figures/predator-prey.mp4)

::::

In the rest of this chapter we present methods as though they are for scalar equations, but their application to systems is taken for granted. The generalization of error analysis can be more complicated, but our statements about order of accuracy and other properties are true for systems as well as scalars. The codes are all written to accept systems.

## Transformation of high-order systems

Fortunately, the ability to solve first-order ODE systems implies the ability to solve systems of higher differential order, too. The reason is that there is a systematic way to turn a higher-order problem into a first-order one of higher dimension.

(example-nlosc3)=

````{prf:example}
 Consider the nonlinear initial-value problem
  
```{math}
  y''+(1+y')^3 y = 0, \qquad y(0)= y_0, \quad y'(0) = 0.
```

In order to write this problem as a first-order system we define two scalar unknown functions, $u_1 = y$ and $u_2 = y'$. With these definitions, we have the two differential equations

```{math}
\begin{split}
  u_1' &= u_2, \\
  u_2' &= -(1+u_2)^3 u_1,
\end{split}
```

which is a first-order system in two dimensions. The initial
condition of the system is

```{math}
  u_1(0) = y_0, \quad u_2(0) = 0.
```
````

````{prf:example}
:label: example-systems-coupledpendula
Two identical pendulums suspended from the same rod and swinging in parallel planes can be modeled as the second-order system

```{math}
\begin{split}
  \theta_1''(t) +\gamma \theta_1' + \frac{g}{L} \sin \theta_1 +
  k(\theta_1-\theta_2) &= 0,\\
  \theta_2''(t) +\gamma \theta_2' + \frac{g}{L} \sin \theta_2 +
  k(\theta_2-\theta_1) &= 0,
\end{split}
```

where $\theta_1$ and $\theta_2$ are angles made by the two pendulums, $L$ is the length of each pendulum, $\gamma$ is a frictional parameter, and $k$ is a parameter describing a torque produced by the rod when it is twisted. We can convert this problem into a first-order system using the substitutions

```{math}
  u_1 = \theta_1, \quad u_2 = \theta_2, \quad u_3 = \theta_1', \quad
  u_4 = \theta_2'.
```

With these definitions the system becomes

```{math}
\begin{split}
  u_1' &= u_3, \\
  u_2' &= u_4, \\
  u_3' &= -\gamma u_3 - \frac{g}{L}\sin u_1 + k(u_2-u_1), \\
  u_4' &= -\gamma u_4 - \frac{g}{L}\sin u_2 + k(u_1-u_2),
\end{split}
```

which is a first-order system in four dimensions. To complete the description of the problem, you need to specify values for $\theta_1(0)$, $\theta_1'(0)$, $\theta_2(0)$, and $\theta_2'(0)$.
````

The trick illustrated in the preceding examples is always available. Suppose $y$ is a scalar dependent variable in the system. You should introduce a component of $\mathbf{u}$ for $y$, $y'$, etc., up to but not including the highest derivative appearing anywhere for $y$. This is done for each scalar variable in the original system. There should be one component of $\mathbf{u}$ for each scalar initial condition given. Many equations for the first-order system then come from the trivial relationships among all the lower derivatives. The remaining equations for the system come from the original, high-order equations. In the end, there must be as many scalar component equations as unknown first-order variables.

::::{prf:example} Coupled pendulums
:label: demo-systems-coupledpendula

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

First we check the behavior of the system when the pendulums are uncoupled, i.e., when $k=0$.
```{tip}
:class: dropdown
We use a closure here to pass the fixed parameter values into `couple`.
```

```{code-cell}
gamma, L, k = 0.01, 0.5, 0.0
du_dt = lambda t, u: couple(t, u, (gamma, L, k))
sol = solve_ivp(du_dt, tspan, u0, t_eval=linspace(0, 50, 1001))
plot(sol.t, sol.y[:2, :].T)    # first two components of solution
xlabel("t"), ylabel("angle")
title("Uncoupled pendulums");
```

You can see that the pendulums swing independently:

```{code-cell}
:tags: [hide-input, remove-output]
from matplotlib import animation
fig, ax = subplots(1, 2)
rod_ = [] 
bob_ = []
for k in [0, 1]:
    rod_.append(ax[k].plot([], [], linewidth=3)[0])
    bob_.append(ax[k].scatter([0], [0], color='black', s=50))
    ax[k].set_aspect('equal')
    ax[k].set_xlim([-1.1, 1.1])
    ax[k].set_ylim([-1.1, 0.4])
    ax[k].set_axis_off()
time_ = ax[0].text(-0.95, 0.25, "", fontsize=11)

def snapshot(i):
    q = sol.y[:2, i]
    for k in [0, 1]:
        rod_[k].set_data([0, sin(q[k])], [0, -cos(q[k])])
        bob_[k].set_offsets([sin(q[k]), -cos(q[k])])
    time_.set_text(f"t = {sol.t[i]:.1f} s")

anim = animation.FuncAnimation(fig, snapshot, frames=401)
anim.save("figures/pendulums-weak.mp4", fps=24)
```

![Uncoupled pendulums](figures/pendulums-weak.mp4)

Because the model is nonlinear and the initial angles are not small, they have slightly different periods of oscillation, and they go in and out of phase.

With coupling activated, a different behavior is seen.

```{code-cell}
k = 0.75    # changes the value in the du_dt closure
sol = solve_ivp(du_dt, tspan, u0, t_eval=linspace(0, 50, 1001))
plot(sol.t, sol.y[:2, :].T)
xlabel("t"), ylabel("angle")
title("Coupled pendulums");
```

The coupling makes the pendulums swap energy back and forth:

```{code-cell}
:tags: [hide-input, remove-output]
anim = animation.FuncAnimation(fig, snapshot, frames=1001)
anim.save("figures/pendulums-strong.mp4", fps=24)
```

![Coupled pendulums](figures/pendulums-strong.mp4)

::::

## Exercises

``````{exercise}
:label: problem-systems-rewrite
✍ Rewrite the given higher order problems as first-order systems.

**(a)** $y'''-3y''+3 y' -y = t, \: y(0) = 1, \: y'(0) = 2, \: y''(0) = 3$

**(b)** $y'' + 4 (x^2-1)y' + y = 0, \: y(0) = 2, \: y'(0) = -1$

**(c)** For a given constant $a$,

```{math}
:numbered: false
\begin{split}
x'' + \frac{a x}{(x^2+y^2)^{3/2}} &= 0,\\
y'' + \frac{a y}{(x^2+y^2)^{3/2}} &= 0,
\end{split}
```

with initial values $x(0) = 1$, $x'(0)=y(0) = 0$, $y'(0)=3$

**(d)** $y^{(4)} -y = e^{-t}, \: y(0) = 0, \: y'(0) = 0, \: y''(0) = 1,\: y'''(0) = 0$

**(e)** $y'''-y''+y'-y = t, \: y(0) = 1, \: y'(0) = 2, \: y''(0) = 3$
``````

``````{exercise}
:label: problem-systems-byhand
✍ Write the given IVP as a system. Then do two steps of Euler's method by hand (perhaps with a calculator) with the indicated step size $h$. Using the given exact solution, compute the error after the second step.

**(a)** $y''+ 4y = 4t, \: y(0) = 1,\: y'(0) = 1; \: \hat{y}(t) = t+\cos (2t),\: h=0.1$

**(b)** $y''- 4y = 4t, \: y(0) = 2,\: y'(0) = -1; \: \hat{y}(t) = e^{2t} + e^{-2t}-t,\: h=0.1$

**(c)** $2 x^2 y'' +3xy' - y = 0, \: y(2) = 1, \: y'(2) = -1/2,  \: \hat{y}(x) = 2/x, h = 1/8$

**(d)** $2 x^2 y'' +3xy' - y = 0,\: y(1) = 4, \: y'(1) = -1, \: \hat{y}(x) = 2(x^{1/2} + x^{-1}), h=1/4$
``````

``````{exercise}
:label: problem-systems-euler
⌨ Solve the following IVPs using {numref}`Function {number} <function-euler>` with $n=1000$ steps. Plot the solution and its first derivative together on one plot, and plot the error in each component as functions of time on another.

**(a)** $y''+ 4y = 4t, \: 0< t< 2\pi, \: y(0) = 1,\: y'(0) = 1; \: \hat{y}(t) = t+\cos (2t)$

**(b)** $y''+ 9y = \sin(2t), \: 0< t< 2\pi, \: y(0) = 2,\: y'(0) = 1$; $\quad \hat{y}(t) = (1/5) \sin(3t) + 2 \cos (3t)+  (1/5) \sin (2t)$

**(c)** $y''- 4y = 4t \: 0< t< 1.5, \: y(0) = 2,\: y'(0) = -1; \: \hat{y}(t) = e^{2t} + e^{-2t}-t$

**(d)** $y''+ 4y'+ 4y = t, \: 0< t< 4, \: y(0) = 1,\: y'(0) = 3/4; \: \hat{y}(t) = (3t+5/4)e^{-2t} + (t-1)/4$

**(e)** $x^2 y'' +5xy' + 4y = 0,\: 1<x<e^2, \: y(1) = 0, \: y'(1) = 2, \: \hat{y}(x) = (2/x^2) \ln x$

**(f)** $x^2 y'' +5xy' + 4y = 0,\: 1<x<e^2, \: y(1) = 1, \: y'(1) = -1, \: \hat{y}(x) = x^{-2}( 1 + \ln x)$

**(g)** $2 x^2 y'' +3xy' - y = 0,\: 2<x<20, \: y(2) = 1, \: y'(2) = -1/2, \: \hat{y}(x) = 2/x$

**(h)** $2 x^2 y'' +3xy' - y = 0,\: 1<x<16, \: y(1) = 4, \: y'(1) = -1, \: \hat{y}(x) = 2(x^{1/2} + x^{-1})$

**(i)** $x^2 y'' -xy' + 2y = 0,\: 1<x<e^{\pi}, \: y(1) = 3, \: y'(1) = 4$; $\quad \hat{y}(x) = x \left[ 3 \cos \left( \ln x \right)+\sin \left( \ln x \right) \right]$

``````

``````{exercise}
:label: problem-systems-SIR
⌨ A disease that is endemic to a population can be modeled by tracking the fraction of the population that is susceptible to infection, $v(t)$, and the fraction that is infectious, $w(t)$. (The rest of the population is considered to be recovered and immune.) A typical model is the *SIR model* (see {cite}`brittonEssentialMathematical2003`)

```{math}
:numbered: false
\frac{dv}{dt} = 0.2(1-v) - 3vw, \qquad \frac{dw}{dt} = (3v-1)w.
```

Starting with $v(0) = 0.95$ and $w(0) = 0.05$, imitate the method of @demo-systems-coupledpendula to find the long-term steady values of $v(t)$ and $w(t)$. Plot both components of the solution as functions of time.
``````

``````{exercise}
:label: problem-systems-phaseplane
⌨ In each case below, imitate the method of @demo-systems-coupledpendula to solve the given ODE for $0\le t \le 10$ with the given initial conditions. Plot the results together as curves in the phase plane (that is, with $x$ and $y$ as the axes of the plot).

**(a)** 

```{math}
:numbered: false
\begin{split}
x'(t) & = - 4y + x(1-x^2-y^2),\\
y'(t) & = 4x + y(1-x^2-y^2),
\end{split}
```

with $[x(0),y(0)]=[0.1,0]$ and $[x(0),y(0)]=[0,1.9]$.

**(b)** 

```{math}
:numbered: false
\begin{split}
x'(t) & = - 4y - \tfrac{1}{4}x(1-x^2-y^2)(4-x^2-y^2),\\
y'(t) & = 4x - \tfrac{1}{4}y(1-x^2-y^2)(4-x^2-y^2),
\end{split}
```

with $[x(0),y(0)]=[0.95,0]$, $[0,1.05]$, and $[-2.5,0]$.

``````

``````{exercise}
:label: problem-systems-fitznag
⌨ The [FitzHugh–Nagumo equations](wiki:FitzHugh–Nagumo_model) are a simple model of the repeated firing of a neuron. They are given by

```{math}
:numbered: false
\begin{split}
\frac{d v_1}{dt} &= - v_1(v_1-1)(v_1-a) - v_2 + I, \\
\frac{d v_2}{dt} &= \epsilon ( v_1 - \gamma v_2).
\end{split}
```

Assume $v_1(0) = 0.5$, $v_2(0) = 0.1$, $a = 0.1$, $\epsilon = 0.008$, $\gamma = 1$. For each value of $I$ below, imitate the method of @demo-systems-coupledpendula tofind and plot the solution for $0\le t \le 600$. The solutions are highly sensitive to $I$, and you need to change the requested absolute and relative error tolerances to $10^{-9}$. In each case the solution quickly approaches a periodic oscillation.

**(a)** $I = 0.05527,\quad$
**(b)** $I = 0.05683,\quad$
**(c)** $I = 0.0568385,\quad$
**(d)** $I = 0.05740$.

This exploration was carried out by Baer and Erneux {cite}`baerSingularHopf1986`.
``````
