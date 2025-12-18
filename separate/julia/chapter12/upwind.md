---
numbering:
  enumerator: 12.2.%s
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.12
---
```{code-cell}
:tags: [remove-cell]
import Pkg; Pkg.activate("/Users/driscoll/Documents/GitHub/fnc")

using FNCFunctions

using Plots
default(
    titlefont=(11,"Helvetica"),
    guidefont=(11,"Helvetica"),
    linewidth = 2,
    markersize = 3,
    msa = 0,
    size=(500,320),
    label="",
    html_output_format = "svg"
)

using PrettyTables, LaTeXStrings, Printf
using LinearAlgebra

@ptconf backend = Val(:html) tf = tf_html_simple
```

(section-advection-upwind)=

# Upwinding and stability

Advection phenomena are characterized by the transport of information at finite speed. This implies that the solution at a given point in space and time can be affected only by the initial data in a limited region of space. For now, we continue to ignore the effects of boundaries.

```{index} ! domain of dependence
```

::::{prf:definition} Domain of dependence
:label: definition-domaindependence
Let $u(x,t)$ be the solution of an evolutionary PDE with initial condition $u_0(x)$. The **domain of dependence** of the solution at $(x,t)$ is the set of all $s$ such that $u_0(s)$ can possibly affect $u(x,t)$.
::::

```{index} advection equation
```

````{prf:example} Domain of dependence
The linear advection equation $u_t+cu_x = 0$ has solution $u(x,t) = u_0(x-ct)$, where $u_0(x)$ is the initial condition. Therefore, the solution at $(x,t)$ depends on the initial data only at the single point $x-ct$. The domain of dependence is $\{x-ct\}$.
````

````{prf:example}
The heat equation $u_t = u_{xx}$ has a solution that, at any positive time $t$, depends on all the initial data. That is, the domain of dependence everywhere is the entire real line.
````

```{index} ! numerical domain of dependence
```

Any numerical method we choose to solve a PDE has analogous property.

::::{prf:definition} Numerical domain of dependence
:label: definition-domaindependencenum
Let a numerical method be used to solve an evolutionary PDE on a grid, such that $U_{i,j}$ is the approximate solution at $x=x_i$, $t=t_j$. The **numerical domain of dependence** of the method at $(x_i,t_j)$ is the set of all $x_k$ such that the initial data $U_{k,0}$ can possibly affect $U_{i,j}$.
::::

::::{prf:example}
:label: example-upwind-centered
In $u_t+cu_x=0$, suppose we discretize $u_x$ by a centered difference:

```{math}
:label: cflcentral
  u_x(x_i,t_j) \approx \frac{U_{i+1,j}-U_{i-1,j}}{2h}.
```

If we use the Euler time discretization with step size $\tau$, then

```{math}
:label: cflcenteuler
  U_{i,j+1} = U_{i,j} - \frac{c\tau}{2h} (U_{i+1,j}-U_{i-1,j}).
```

Starting with $j=0$, we find that $U_{i,1}$ depends on $U_{i-1,0}$, $U_{i,0}$, and $U_{i+1,0}$. Hence, the numerical domain of dependence at $(x_i,t_1)$ is $\{x_{i-1},x_i,x_{i+1}\}$.

Now we set $j=1$. From @cflcenteuler, we see that $U_{i,2}$ depends on $U_{i-1,1}$, $U_{i,1}$, and $U_{i+1,1}$. In turn, each of these values at time $t_1$ depends on three values at time $t_0$:

$$
U_{i-1,1} &\;\text{ uses }\; U_{i-2,0},\, U_{i-1,0},\, U_{i,0}, \\
U_{i,1} &\;\text{ uses }\; U_{i-1,0},\, U_{i,0},\, U_{i+1,0}, \\
U_{i+1,1} &\;\text{ uses }\; U_{i,0},\, U_{i+1,0},\, U_{i+2,0}.
$$

Considering only the unique values, the numerical domain of dependence at $(x_i,t_2)$ is $\{x_{i-2},x_{i-1},x_i,x_{i+1},x_{i+2}\}$. Continuing this kind of reasoning backward from any $U_{i,j}$, we find that the numerical domain of dependence at $(x_i,t_j)$ is $\{x_{i-j},\ldots,x_{i+j}\}$.
::::

## The CFL condition

We now state an important principle about a necessary relationship between domains of dependence.

```{index} ! CFL condition
```

::::{prf:theorem} Courant–Friedrichs–Lewy (CFL) condition
:label: theorem-upwind-cfl
In order for a numerical method for an evolutionary equation to converge to the correct solution, the numerical domain of dependence in the limit $h \to 0,$ $\tau\to 0$ must contain the exact domain of dependence.
::::

Although we will not provide the rigor behind this theorem, its conclusion is not difficult to justify. If the CFL condition does not hold, the exact solution at $(x,t)$ could be affected by a change in the initial data while having no effect on the numerical solution. Hence, there is no way for the method to get the solution correct for all problems. By contradiction, then, the CFL criterion is necessary for convergence.

:::{caution}
The CFL condition is a *necessary* criterion for convergence, but not a *sufficient* one. For instance, we could define $U_{i,j}$ to be any weighted convergent sum of all values of $U_{i,0}$. While that would make the numerical domain of dependence equal to the entire real line, this method has nothing to do with solving a PDE correctly!
:::

::::{prf:example}
:label: example-upwind-centeredcfl
In @example-upwind-centered we concluded that the numerical domain of dependence for a centered Euler discretization of the advection equation at $(x_i,t_j)$ is $\{x_{i-j},\ldots,x_{i+j}\}$. @figure-cflpicture illustrates what happens as $h$ and $\tau$ go to zero in a manner that leaves the ratio $h/\tau$ constant.

```{figure} figures/cflpicture.svg
:name: figure-cflpicture
:width: 450px
Numerical domain of dependence for the explicit time stepping scheme in {numref}`Example {number} <example-upwind-centered>`. As  $\tau$ and $h$ approach zero, the shaded region is filled in.
```

In order to observe both the exact and numerical domains of dependence at a fixed location $(x,t)$, we must have $x=ih$ and $t=j\tau$. The numerical domain of dependence fills in the interval $[x_{i-j},x_{i+j}]$, which is $[x-jh,x+jh]$, or

$$
[x - \frac{h}{\tau} t, x + \frac{h}{\tau} t].
$$

This interval captures the exact domain of dependence $\{x-ct\}$ only if

```{math}
:label: cfl-speed
  |c| \le \frac{h}{\tau}.
```

::::

Equation {eq}`cfl-speed` is the implication of the CFL condition for the stated discretization. Notice that $h/\tau$ is the speed at which information moves in the numerical method, which leads to the following restatement.

:::{prf:observation}
The CFL condition requires that the maximum propagation speed in the numerical method be at least as large as the maximum speed in the original PDE problem.
:::

We can rearrange {eq}`cfl-speed` to imply a necessary time step restriction $\tau \le h/|c|$. This restriction for advection is much less severe than the $\tau = O(h^2)$ restriction we derived for Euler in the heat equation in {numref}`section-diffusion-stiffness`. This is our first clear indication that advection is less stiff than diffusion.

::::{prf:example} The CFL condition in action
:label: demo-upwind-cfl

We solve linear advection with velocity $c=2$ and periodic end conditions. The initial condition is numerically, though not mathematically, periodic.

For time stepping, we use the adaptive explicit method `RK4`.

```{code-cell}
using OrdinaryDiffEq
m = 400;
x, Dₓ = FNC.diffper(m, [0, 1])
u_init = x -> exp( -80 * (x - 0.5)^2 )
ode = (u, c, t) -> -c * (Dₓ*u)
IVP = ODEProblem(ode, u_init.(x), (0., 2.), 2.)
u = solve(IVP, RK4());
```

```{code-cell}
:tags: hide-input
using Plots
t = range(0, 2, 81);
U = reduce(hcat, u(t) for t in t)
contour(x, t, U'; 
    color=:redsblues,  clims=(-1, 1),
    xaxis=(L"x"),  yaxis=(L"t"),
    title="Linear advection",  right_margin=3Plots.mm)
```

In the space-time plot above, you can see the initial hump traveling rightward at constant speed. It fully traverses the domain once for each integer multiple of $t=1/2$. 

If we cut $h$ by a factor of 2 (i.e., double $m$), then the CFL condition suggests that the time step should be cut by a factor of 2 also.

```{code-cell}
println("Number of time steps for m = 400: $(length(u.t))")

m = 800;
x, Dₓ = FNC.diffper(m, [0, 1])
IVP = ODEProblem(ode, u_init.(x), (0., 2.), 2.)
u = solve(IVP, RK4())
println("Number of time steps for m = 800: $(length(u.t))")
```

::::

## Upwinding

The advection equation allows information to propagate in only one direction. This asymmetry can have important repercussions.

```{index} ! upwind direction
```

````{prf:definition} Upwind direction
:label: definition-upwind
If the domain of dependence at $(x,t)$ always lies to one side of $x$, that side is called the **upwind direction** of the PDE, and the other direction is called the **downwind direction**.
````

Numerical methods can also have a directional preference.

::::{prf:example}
:label: example-upwind-onesided
Suppose that in $u_t+cu_x=0$ we use the backward difference

```{math}
:label: cflbackward 
  u_x(x_i,t_j)  \approx \frac{U_{i,j}-U_{i-1,j}}{h}  
```

together with an Euler scheme in time. It should be clear that $U_{i,j}$ depends only on points to the left of $x_i$, i.e., the upwind direction of the numerical method is to the left. But if $c<0$, the upwind direction of the PDE is to the right. Hence, it is impossible to satisfy the CFL condition under these choices.

Similarly, the forward difference

```{math}
:label: cflforward
  u_x(x_i,t_j)  \approx \frac{U_{i+1,j}-U_{i,j}}{h}
```

with explicit time stepping leads to the inverse conclusion: its upwind direction is to the right, and it must fail if $c>0$.
::::

The reasoning of @example-upwind-onesided is readily generalized: if the numerical method has an upwind direction, the CFL condition requires that it must agree with the upwind direction of the PDE.

:::{prf:observation} Upwinding
If PDE and a numerical method have different upwinding directions, the method cannot converge to the exact solution.
:::

```{note}
It probably seems like one should always use a centered difference scheme, so that upwinding is not an issue. However, at a shock front, this requires differencing across a jump in the solution, which causes its own set of difficulties.
````

## Inflow boundary condition

```{index} ! boundary conditions; inflow
```

Now suppose that the linear advection equation is posed on a finite domain $x \in [a,b]$.
Since the PDE has only a first-order derivative in $x$, we should have only one boundary condition. But should it be specified at the left end, or at the right end?

Upwinding considerations provide the answer. If we impose a condition at the downwind side of the domain, there is no way for that boundary information to propagate into the interior of the domain as time advances. Conversely, at points close to the upwind boundary, the domain of dependence soon wants to move past the boundary, which is impossible, so there has to be information provided there.

```{important}
For a PDE with an upwind direction, the boundary condition must be specified at the inflow (i.e., upwind) end of the domain.
```

For $c>0$ in the advection equation, the inflow is at the left end, and for $c<0$, it is at the right end.

::::{prf:example} Upwind versus downwind
:label: demo-upwind-direction

If we solve advection over $[0,1]$ with velocity $c=-1$, the right boundary is in the upwind/inflow direction. Thus a well-posed boundary condition is $u(1,t)=0$.

We'll pattern a solution after {numref}`Function {number} <function-parabolic>`. Since $u(x_m,t)=0$, we define the ODE interior problem {eq}`mol-interior` for $\mathbf{v}$ without $u_m$. For each evaluation of $\mathbf{v}'$, we must extend the data back to $x_m$ first.

```{code-cell}
m = 100
x, Dₓ = FNC.diffmat2(m, [0, 1])

interior = 1:m
extend = v -> [v; 0]

function ode!(f, v, c, t)
    u = extend(v)
    uₓ = Dₓ * u
    @. f = -c * uₓ[interior]
end;
```

Now we solve for an initial condition that has a single hump.

```{code-cell}
init = @. exp( -80*(x[interior] - 0.5)^2 )
ivp = ODEProblem(ode!, init, (0., 1), -1)
u = solve(ivp);
```

```{code-cell}
t = range(0, 0.75, 80)
U = reduce(hcat, extend(u(t)) for t in t)
contour(x, t, U';
    color=:blues,  clims=(0, 1), 
    xaxis=(L"x"),  yaxis=(L"t"),
    title="Advection with inflow BC")
```

We find that the hump gracefully exits out the downwind end.

```{code-cell}
:tags: hide-input, remove-output
anim = @animate for t in range(0, 1, 161) 
    plot(x, extend(u(t));
        label=@sprintf("t = %.4f", t), 
        xaxis=(L"x"),  yaxis=(L"u(x, t)", (0, 1)), 
        title="Advection equation with inflow BC",  dpi=150)
end
mp4(anim,"figures/advection-inflow.mp4")
```

![Advection with inflow BC](figures/advection-inflow.mp4)


If instead of $u(1,t)=0$ we were to try to impose the downwind condition $u(0,t)=0$, we only need to change the index of the interior nodes and where to append the zero value.

```{code-cell}
interior = 2:m+1
extend = v -> [0; v]

init = @. exp( -80*(x[interior] - 0.5)^2 )
ivp = ODEProblem(ode!, init, (0., 0.25), -1)
u = solve(ivp);
```

```{code-cell}
:tags: hide-input
t = range(0, 0.2, 61)
U = reduce(hcat, extend(u(t)) for t in t)
contour(x, t, U'; 
    color=:redsblues,  clims=(-1, 1),
    xaxis=(L"x"),  yaxis=(L"t"), 
    title="Advection with outflow BC",  right_margin=3Plots.mm)
```

This time, the solution blows up as soon as the hump runs into the boundary because there are conflicting demands there.

```{code-cell}
:tags: hide-input, remove-output
anim = @animate for t in range(0, 0.2, 41) 
    plot(x, extend(u(t));
        label=@sprintf("t = %.4f", t),
        xaxis=(L"x"),  yaxis=(L"u(x,t)", (0, 1)), 
        title="Advection equation with outflow BC",  dpi=150)
end
mp4(anim,"figures/advection-outflow.mp4")
```
![Advection with outflow BC](figures/advection-outflow.mp4)


::::

## Exercises

``````{exercise}
:label: problem-upwind-weather
✍ Suppose you want to model the weather, including winds up to speed $200$ km/hr, using an explicit method with a second-order centered spatial discretization. If the shortest time step you can take is 4 hr, what is the CFL limit on the spatial resolution of the model? Is this a lower bound or an upper bound?
``````

``````{exercise}
:label: problem-upwind-freeway
✍ Suppose you want to model the traffic on a high-speed freeway using an explicit method with a second-order centered spatial discretization. Derive a CFL condition on the allowable time step, stating your assumptions carefully.
``````

``````{exercise}
:label: problem-upwind-limit
✍ For the heat equation, the domain of dependence at any $(x,t)$ with $t>0$ is all of $x \in (-\infty,\infty)$. Show that the CFL condition implies that $\tau/h\to 0$ is required for convergence as $h\to 0$.
``````

``````{exercise}
:label: problem-upwind-inflow
✍ Suppose you wish to solve $u_t = u u_x$ for $x\in[-1,1]$.

**(a)** If $u(x,0) = -2+\sin(\pi x)$,  which end of the domain is the inflow? 

**(b)** Does the answer to part (a) change if $u(x,0) = 1 + e^{-16x^2}$?
``````
