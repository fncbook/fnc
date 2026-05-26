---
numbering:
  enumerator: 12.1.%s
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

(section-advection-traffic)=

# Traffic flow

Have you ever been driving on a highway when you suddenly came upon a traffic jam? Maybe you had to brake harder than you would like to admit, endured a period of bumper-to-bumper progress, and experienced a much more gradual emergence from dense traffic than the abrupt entry into it. The mathematics of this phenomenon are well understood.

Consider a one-dimensional road extending in the $x$ direction. We represent the vehicles by a continuous density function $\rho(x)$.  The flow rate or [flux](wiki:Flux) of vehicles, expressed as the number of cars per unit time crossing a fixed point on the road, is denoted by $q$. We assume that this flux depends on the local density of cars. It's reasonable to suppose that the flux will be zero when $\rho=0$ (no cars), reach a maximum $q_m$ at some $\rho=\rho_m$, and approach zero again as the density approaches a critical density $\rho_c$. These conditions are met by the model

```{math}
:label: trafficQ
Q_0(\rho) = \frac{4 q_m \rho_m \rho (\rho-\rho_c) (\rho_m-\rho_c)}{[\rho (\rho_c-2 \rho_m)+\rho_c \rho_m]^2}.
```

Observations ({cite}`whithamLinearNonlinear1974`, Chapter 3) suggest that good values for a three-lane highway are $\rho_c = 1080$ vehicles per km, $\rho_m=380$ vehicles per km, and $q_m=4500$ vehicles per hour. In addition, we want to account for the fact that drivers anticipate slowing down or speeding up when they perceive changes in density, and therefore use

$$
q=Q_0(\rho)-\epsilon \rho_x
$$

for a small $\epsilon > 0$.

## Conservation laws

```{index} ! conservation law
```

*Conservation laws* play a major role in science and engineering. They are typically statements that matter, energy, momentum, or some other meaningful quantity cannot be created or destroyed. In one dimension they take the form

```{math}
:label: conservation
u_t + q_x = 0,
```

where $u$ represents the density of a conserved quantity, such as mass, and $q$ represents the flux (flow rate) of $u.$ Typically, the flux is a function of $u$, $q=Q(u)$. By the chain rule,

```{math}
:label: conservation-velocity
u_t + Q'(u) u_x = 0.
```

For the particular $Q_0(\rho)$ for traffic flow, we arrive at the evolutionary PDE

```{math}
:label: trafficpde
  \rho_t + Q_0'(\rho) \rho_x = \epsilon \rho_{xx}.
```

Note that $Q_0'$ has the dimensions of (cars per time) over (cars per length), or length over time. We recognize the first and last terms as indicative of diffusion, but the middle term has a different effect. A similar term appears in the Black–Scholes equation.

## Advection equation

Let's momentarily consider a simpler, more fundamental PDE.

```{index} ! advection equation, ! hyperbolic PDE
```

::::{prf:definition} Advection equation
:label: definition-advectionequation
The {term}`advection equation` in one dimension is

```{math}
:label: advectpde
  u_t + c u_x = 0,
```

where $c$ is constant.
::::

The advection equation is the archetype of a **hyperbolic PDE**, which is a separate class from the parabolic PDEs.

It's surprisingly easy to produce a solution of {eq}`advectpde`. Let $u(t,x)=\psi(x-c t)$, where $\psi$ is any differentiable function of one variable. Then the chain rule says that

$$
u_t + cu_x = \psi'(x-c t)\cdot(-c)  + c [\psi'(x-c t)] = 0,
$$

and the PDE is satisfied. Furthermore, $u(x,0) = \psi(x)$, so the initial value of $u$ determines $\psi$. We have proved the following.

::::{prf:theorem} Advection equation
If $u_0(x)$ is a continuously differentiable function, then the solution of the advection equation {eq}`advectpde` satisfying $u(x,0)=u_0(x)$ is

```{math}
:label: eq-advectsol
u(x,t) = u_0(x-c t),
```

if there are no boundary conditions.
::::

The form of @eq-advectsol tells us that $u$ remains constant along any path on which $x-c t$ is constant, i.e., $x=a + c t$ for a constant $a$. Thus, if $c>0$, a fixed value of $u$ moves rightward with speed $|c|$, and if $c<0$, it moves leftward with speed $|c|$.  

::::{prf:observation}
Solutions to the advection equation propagates with constant speed and fixed shape.
::::

```{note}
Keep in mind that $c$ is a velocity, not a speed: if $c>0$, solutions travel rightward, and if $c<0$, they travel leftward.
```

We can solve {eq}`advectpde` by the method of lines as in [Chapter 11](../diffusion/overview). We need the second-order first-derivative matrix for periodic end conditions,

```{math}
:label: trafficdiffmat
  \mathbf{D}_x =
 \frac{1}{2h}
    \begin{bmatrix}
      0 & 1 & & & -1 \\
      -1 & 0 & 1 & & \\
      & \ddots & \ddots & \ddots & \\
      & & -1 & 0 & 1 \\
      1 & & & -1 & 0
    \end{bmatrix}.
```

This matrix is returned by {numref}`Function {number} <function-diffper>`.

``````{prf:example} Advection equation
:label: demo-traffic-advection

We solve the advection equation on $[-4,4]$ with periodic end conditions using the method of lines.


 In the following definition we allow the velocity $c$ to be specified as a parameter in the `ODEProblem`.

```{code-cell}
x, Dx, Dxx = FNC.diffper(300, [-4, 4])
f = lambda t, u: -c * (Dx @ u)
```

The following initial condition isn't mathematically periodic, but the deviation is less than machine precision. We specify RK4 as the solver.  

```{code-cell}
from scipy.integrate import solve_ivp
u_init = 1 + exp(-3 * x**2)
c = 2
sol = solve_ivp(f, [0, 3.0], u_init, method="Radau", dense_output=True)
```

```{code-cell}
:tags: hide-input
for t in arange(0, 3, 2/3):
    plot(x, sol.sol(t), label=f"t = {t:.1f}")
legend()
xlabel("$x$"),  ylabel("$u(x,t)$")
title("Advection with periodic boundary");
```

An animation shows the solution nicely. The bump moves with speed 2 to the right, reentering on the left as it exits to the right because of the periodic conditions. 

```{code-cell}
:tags: hide-input
from matplotlib import animation
fig, ax = subplots()
curve = ax.plot(x, u_init)[0]
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

ax.set_xlabel("$x$")
ax.set_ylabel("$u(x,t)$")
ax.set_ylim(0.9, 2.1)
ax.set_title("Advection equation with periodic boundary")

def snapshot(t):
    curve.set_ydata(sol.sol(t))
    time_text.set_text(f"t = {t:.2f}")

anim = animation.FuncAnimation(
    fig, snapshot, frames=linspace(0, 3, 201)
    )
anim.save("advection-periodic.mp4", fps=30)
close()
```

![Advection with periodic boundary](advection-periodic.mp4)


``````

If you look carefully at @demo-traffic-advection, you'll notice that we used the time integrator `RK4`, a nonstiff method. As we will see later in this chapter, the pure advection equation is not inherently stiff.

## Solutions for traffic flow

```{index} advection-diffusion equation
```

Returning to the traffic flow PDE {eq}`trafficpde`, we can interpret the conservation law as advection at velocity $Q_0'(\rho)$.  The dependence of velocity on the solution is different from the linear PDE {eq}`advectpde` and is the key to the peculiar behavior of traffic jams. In particular, the velocity is low at high density, and higher at low density.  The term on the right side of {eq}`trafficpde` provides a bit of diffusion, and the parameter $\epsilon$ determines the balance between the two effects—advection effects dominate if $\epsilon$ is small.

Exact solutions of {eq}`trafficpde` are much harder to come by than for the standard advection equation, but the method of lines is still effective.

::::{prf:example} Traffic flow model
:label: demo-traffic-solve

We solve for traffic flow using periodic boundary conditions.

The following are parameters and a function relevant to defining the problem. 

```{code-cell}
rho_c = 1080
rho_m = 380
q_m = 10000
Q0prime = (
    lambda rho: q_m
    * 4
    * rho_c**2
    * (rho_c - rho_m)
    * rho_m
    * (rho_m - rho)
    / (rho * (rho_c - 2 * rho_m) + rho_c * rho_m) ** 3
)
```

Here we create a discretization on $m=800$ points.

```{code-cell}
x, Dx, Dxx = FNC.diffper(800, [0, 4])
```

Next we define the ODE resulting from the method of lines.

```{code-cell}
ode = lambda t, rho: -Q0prime(rho) * (Dx @ rho) + ep * (Dxx @ rho)
```

Our first initial condition has moderate density with a small bump. Because of the diffusion present, we use a stiff solver for the IVP.

```{code-cell}
rho_init = 400 + 10 * exp(-20 * (x - 3) ** 2)
ep = 0.02
sol = solve_ivp(ode, [0, 1.0], rho_init, method="Radau", dense_output=True)
```

```{code-cell}
:tags: hide-input
for t in linspace(0, 1, 6):
    plot(x, sol.sol(t), label=f"t = {t:.1f}")

xlabel("$x$"),  ylabel("car density")
legend(),  title("Traffic flow");
```

The bump slowly moves backward on the roadway, spreading out and gradually fading away due to the presence of diffusion.

```{code-cell}
:tags: hide-input
fig, ax = subplots()
curve = ax.plot(x, rho_init)[0]
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
ax.set_xlabel("$x$")
ax.set_ylabel("density")
ax.set_ylim(400, 410)
ax.set_title("Traffic flow")

anim = animation.FuncAnimation(
    fig, snapshot, frames=linspace(0, 1, 101)
    )
anim.save("traffic-small.mp4", fps=30)
close()
```

![Traffic flow simulation](traffic-small.mp4)

Now we use an initial condition with a larger bump. Note that the scale on the $y$-axis is much different for this solution.

```{code-cell}
rho_init = 400 + 80 * exp(-16 * (x - 3) ** 2)
sol = solve_ivp(ode, [0, 0.5], rho_init, method="Radau", dense_output=True)
```

```{code-cell}
:tags: hide-input
for t in linspace(0, 0.5, 6):
    plot(x, sol.sol(t), label=f"t = {t:.1f}")

xlabel("$x$"),  ylabel("car density")
legend(),  title("Traffic jam");
```

```{code-cell}
:tags: hide-input
fig, ax = subplots()
curve = ax.plot(x, rho_init)[0]
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
ax.set_xlabel("$x$")
ax.set_ylabel("density")
ax.set_ylim(400, 480)
ax.set_title("Traffic jam")

anim = animation.FuncAnimation(
    fig, snapshot, frames=linspace(0, 0.5, 101)
    )
anim.save("traffic-jam.mp4", fps=30)
close()
```

![Traffic jam simulation](traffic-jam.mp4)

In this case the density bump travels backward along the road. It also steepens on the side facing the incoming traffic and decreases much more slowly on the other side. A motorist would experience this as an abrupt increase in density, followed by a much more gradual decrease in density and resulting gradual increase in speed. (You also see some transient, high-frequency oscillations. These are caused by instabilities, as we discuss in simpler situations later in this chapter.)


::::

```{index} shock wave
```

The phenomenon in the second plot of @demo-traffic-solve is called a *shock wave*. The underlying mathematics is not much different from the shock wave that comes off of the wing of a supersonic aircraft in the form of a sonic boom, or from cresting waves under certain conditions in the ocean. In the absence of diffusion ($\epsilon=0$), the shock becomes a jump discontinuity in the solution, which breaks down both the finite differences and the original PDE, requiring different approaches. For many applications, the addition of a small amount of diffusion is appropriate and simple. However, we will first try to come to terms with pure advection in a linear problem.

## Exercises

``````{exercise}
:label: problem-traffic-backward
✍ By the analogy between {eq}`trafficpde` and {eq}`advectpde`, use {eq}`trafficQ` to confirm that constant traffic density moves backward (right to left) for $\rho_m<\rho<\rho_c$, as observed in @demo-traffic-solve. (Note that the derivative of $Q_0$ is given in the code for the example.)
``````

``````{exercise}
:label: problem-traffic-shockspeed

⌨ **(a)** Using as large a discretization and as small a dissipation parameter $\epsilon$  as you can get away with, perform experiments to estimate the speed of the shockwave in @demo-traffic-solve between times $t=0.11$ and $t=0.15$. (Hint: You can use `argmax` to locate the peak of the solution vector at a particular time.)

**(b)** Theory predicts that the speed of the shockwave is the average of $Q_0'$ evaluated at the values of $\rho$ at the top and bottom of the shock. Perform this calculation and compare to the result of part (a).
``````

```{index} Burgers equation
```

``````{exercise}
:label: problem-traffic-burgers

⌨ The simplest model that includes both diffusion and nonlinear advection is the *viscous Burgers equation*, $u_t+u u_x=\epsilon u_{xx}$. Assume periodic end conditions on $-4\le x < 4$, let $\epsilon=0.04$, and suppose $u(x,0)=e^{-2x^2}$. Solve the problem numerically with $m=200$, plotting the solution on one graph at times $t=0,0.5,1,\ldots,3$. (You should see that the bump decays but also steepens as it moves to the right.)
``````

```{index} Kuramoto–Sivashinsky equation
```

``````{exercise}
:label: problem-traffic-ks

⌨ The *Kuramoto–Sivashinsky equation*, $u_t+u u_x=-u_{xx}-\epsilon u_{xxxx}$, exhibits solutions that are considered chaotic. Assume periodic end conditions, $\epsilon=0.05$, and initial condition $u(x,0)=1+e^{-2x^2}$. Using the method of lines, solve the problem numerically with $m=200$ for $-4\le x \le 4$ and $0\le t \le 20$. (You should discretize the fourth derivative as repeated second derivatives.) Make an animation of the plot (preferred), or plot the solution as a surface over $x$ and $t$. 
``````

``````{exercise}
:label: problem-traffic-period

⌨ (Adapted from {cite}`trefethenSpectralMethods2000`.) Consider the problem $u_t + c(x) u_x =0$, for $0 \le x \le 2\pi$, with periodic boundary conditions and variable speed $c(x) = 0.2 + \sin^2(x-1)$. This problem has a solution that is $T$-periodic in time, for some $T\approx 13$.

**(a)** Find this value $T$ accurately to 10 digits by converting the integral $T
= \int_0^T \,dt$ to an integral in $x$ via $\frac{dx}{dt}=c(x)$, then applying numerical integration.

**(b)** Using $u(x,0) =e^{\sin x}$, solve the PDE numerically for $0\le t \le T$ using the answer from part (a). Make an animation of the solution, or a surface plot of the solution as a function of $x$ and $t$.

**(c)**  Adjust the accuracy in time and space until you can verify that $\|u(x,T)-u(x,0)\|_\infty < 0.03$.

``````
