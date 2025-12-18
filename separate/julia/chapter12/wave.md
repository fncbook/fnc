---
numbering:
  enumerator: 12.4.%s
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
x, Dₓ = FNC.diffcheb(m, [-1, 1]);
```

The boundary values of $u$ are given to be zero, so they are not unknowns in the ODEs. Instead they are added or removed as necessary.

```{code-cell}
extend = v -> [0; v; 0]
chop = u -> u[2:m];
```

The following function computes the time derivative of the system at interior points.

```{code-cell}
ode = function(w, c, t)
    u = extend(w[1:m-1])
    z = w[m:2m]
    du_dt = Dₓ * z
    dz_dt = c^2 * (Dₓ * u)
    return [ chop(du_dt); dz_dt ]
end;
```

Our initial condition is a single hump for $u$.

```{code-cell}
u_init = @. exp( -100*(x + 0.5)^2 )
z_init = -u_init
w_init = [ chop(u_init); z_init ];  
```

Because the wave equation is hyperbolic, we can use a nonstiff explicit solver.

```{code-cell}
using OrdinaryDiffEq
IVP = ODEProblem(ode, w_init ,(0., 2.), 2)
w = solve(IVP, RK4());
```

We plot the results for the original $u$ variable only. Its interior values are at indices `1:m-1` of the composite $\mathbf{w}$ variable.

```{code-cell}
using Plots
t = range(0, 2, 80)
U = [extend(w(t)[1:m-1]) for t in t]
contour(x, t, hcat(U...)';
    levels=24,
    color=:redsblues,  clims=(-1, 1),
    xlabel=L"x",  ylabel=L"t",
    title="Wave equation",  right_margin=3Plots.mm)
```

```{code-cell}
:tags: hide-input, remove-output
anim = @animate for t in range(0 ,2, 120)
    plot(x, extend(w(t)[1:m-1]);
        label=@sprintf("t=%.3f",t),
        xaxis=(L"x"),  yaxis=([-1, 1], L"u(x,t)"),
        dpi=150,  title="Wave equation")
end
mp4(anim, "figures/wave-boundaries.mp4")
```

![Wave equation with boundaries](figures/wave-boundaries.mp4)

The original hump breaks into two pieces of different amplitudes, each traveling with speed $c=2$. They pass through one another without interference. When a hump encounters a boundary, it is perfectly reflected, but with inverted shape. At time $t=2$, the solution looks just like the initial condition.


::::

## Variable speed

An interesting situation is when the wave speed $c$ in @wavepde changes discontinuously, as when light passes from one material into another. For this we must replace the term $c^2$ in {eq}`waveMOL` with the matrix $\operatorname{diag}\bigl(c^2(x_0),\ldots,c^2(x_m)\bigr)$.

::::{prf:example} Wave equation with variable speed
:label: demo-wave-speed

We now use a wave speed that is discontinuous at $x=0$; to the left, $c=1$, and to the right, $c=2$.

The ODE implementation has to change slightly.

```{code-cell}
ode = function(w,c,t)
    u = extend(w[1:m-1])
    z = w[m:2m]
    du_dt = Dₓ*z
    dz_dt = c.^2 .* (Dₓ*u)
    return [ chop(du_dt); dz_dt ]
end;
```

The variable wave speed is passed as an extra parameter through the IVP solver.

```{code-cell}
c = @. 1 + (sign(x)+1)/2
IVP = ODEProblem(ode, w_init, (0., 5.), c)
w = solve(IVP, RK4());
```

```{code-cell}
t = range(0, 5, 80)
U = [extend(w(t)[1:m-1]) for t in t]
contour(x, t, hcat(U...)';
    color=:redsblues,  clims=(-1,1),
    levels=24,
    xlabel=L"x",  ylabel=L"t",
    title="Wave equation",
    right_margin=3Plots.mm
    )
```

```{code-cell}
:tags: hide-input, remove-output
anim = @animate for t in range(0,5,181)
    plot(Shape([-1, 0, 0, -1], [-1, -1, 1, 1]), color=RGB(.8, .8, .8), l=0, label="")
    plot!(x, extend(w(t, idxs=1:m-1));
        label=@sprintf("t=%.2f", t), 
        xaxis=(L"x"),  yaxis=([-1, 1], L"u(x,t)"),
        dpi=150,  title="Wave equation, variable speed")
end
mp4(anim, "figures/wave-speed.mp4")
```

![Wave equation with variable speed](figures/wave-speed.mp4)

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
