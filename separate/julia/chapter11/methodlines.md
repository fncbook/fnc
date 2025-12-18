---
numbering:
  enumerator: 11.2.%s
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

(section-diffusion-methodlines)=

# The method of lines

Our strategy in {numref}`section-diffusion-blackscholes` was to discretize both the time and space derivatives using finite differences, then rearrange so that we could march the solution forward through time. It was partially effective, but as @demo-blackscholes-unstable shows, not a sure thing, for reasons we look into starting in the next section.

```{index} ! boundary conditions; periodic
```

First, though, we want to look at a broader version of the discretization approach. To introduce ideas, let's use the simpler heat equation, $u_t = u_{xx}$, as a model. Because boundaries always complicate things, we will start by doing the next best thing to having no boundaries at all: **periodic end conditions**. Specifically, we will solve the PDE over $0\le x < 1$ and require at all times that

$$
u(x+1,t)=u(x,t) \quad \text{for all $x$}.
$$

This is a little different from simply $u(1,t)=u(0,t)$, as {numref}`figure-periodicfun` illustrates.

:::{figure} figures/periodicfun.svg
:label: figure-periodicfun
:alt: periodic function illustration
:align: center

Left: A function whose values are the same at the endpoints of an interval does not necessarily extend to a smooth periodic function. Right: For a truly periodic function, the function values and all derivatives match at the endpoints of one period.
:::

(section-methodlines-semidiscretization)=

## Semidiscretization

As a reminder, we use $\hat{u}$ when we specifically refer to the exact solution of the PDE. In order to avoid carrying along redundant information about the function, we use $x_i = ih$ only for $i=0,\ldots,m-1$, where $h=1/m$, and it's understood that a reference to $x_m$ is silently translated to one at $x_0$. More generally, we have the identity

```{math}
:label: periodicmod
  \hat{u}(x_i,t) = \hat{u}\bigl(x_{(i \bmod{m})},t \bigr)
```

for the exact solution $\hat{u}$ at any value of $i$.

```{index} ! semidiscretization; see method of lines
```

```{index} finite differences; for parabolic PDE, differentiation matrix
```

Next we define a vector $\mathbf{u}$ by

$$
\mathbf{u}(t) = \begin{bmatrix} u_0(t) \\ u_1(t) \\ \vdots \\ u_n(t) \end{bmatrix}.
$$

This step is called **semidiscretization**, since space is discretized but time is not. As in [Chapter 10](../bvp/overview.md), we will replace $u_{xx}$ with multiplication of $\mathbf{u}$ by a differentiation matrix $\mathbf{D}_{xx}$. The canonical choice is the three-point finite-difference formula {eq}`centerFD22`, which in light of the periodicity {eq}`periodicmod` leads to

```{math}
:label: heatFD22
\mathbf{D}_{xx} = \frac{1}{h^2}
  \begin{bmatrix}
  -2 & 1 & & & 1 \\
  1 & -2 & 1 & & \\
  & \ddots & \ddots & \ddots & \\
  & & 1 & -2 & 1 \\
  1 &  & & 1 & -2
  \end{bmatrix}.
```

Note well how the first and last rows have elements that "wrap around" from one end of the domain to the other by periodicity. Because we will be using this matrix quite a lot, we create {numref}`Function {number} <function-diffper>` to compute it, as well as the corresponding second-order first derivative matrix $\mathbf{D}_x$ for periodic end conditions.

``````{prf:algorithm} diffper
:label: function-diffper

```{literalinclude} chapter11.jl
:filename: diffper.jl
:start-after: # begin diffper
:end-before: # end diffper
:language: julia
:linenos: true
```
``````

The PDE $u_t=u_{xx}$ is now approximated by the semidiscrete problem

```{math}
:label: heatMOL
  \frac{d \mathbf{u}(t)}{d t} = \mathbf{D}_{xx} \mathbf{u}(t),
```

which is simply a linear, constant-coefficient system of *ordinary* differential equations. Given the initial values $\mathbf{u}(0)$ obtained from $u(x_i,0)$, we have an initial-value problem that we already know how to solve!

```{index} ! method of lines
```

Semidiscretization is often called the **method of lines**. Despite the name, it is not exactly a single method because both space and time discretizations have to be specified in order to get a concrete algorithm. The key concept is the separation of those two discretizations, and in that way, it's related to separation of variables in analytic methods for the heat equation.

::::{prf:example}
:label: example-methodlines-heatFE
Suppose we solve {eq}`heatMOL` using the Euler IVP integrator {eq}`euler1` from {numref}`section-ivp-euler` (and also AB1 from {numref}`section-ivp-multistep`). We select a time step $\tau$ and discrete times $t_j=j\tau$, $j=0,1,\ldots,n$. We can discretize the vector $\mathbf{u}$ in time as well to get a sequence $\mathbf{u}_j \approx \mathbf{u}(t_j)$ for varying $j$. (Remember the distinction in notation between $\mathbf{u}_j$, which is a vector, and $u_j$, which is a single element of a vector.)

Thus, a fully discrete method for the heat equation is

```{math}
:label: Eulerxx
\mathbf{u}_{j+1} = \mathbf{u}_j + \tau ( \mathbf{D}_{xx} \mathbf{u}_j) = (\mathbf{I} + \tau \mathbf{D}_{xx} ) \mathbf{u}_j.
```

::::

::::{prf:example} Forward Euler for the heat equation
:label: demo-methodlines-heatFE

Let's implement the method of {numref}`Example {number} <example-methodlines-heatFE>` with second-order space semidiscretization.

```{code-cell}
m = 100
x, Dx, Dxx = FNC.diffper(m, [0, 1]);
tfinal = 0.15 
n = 2400           # number of time steps
τ = tfinal / n     # time step    
t = τ * (0:n)      # time values
```

Next we set an initial condition. It isn't mathematically periodic, but the end values and derivatives are so small that for numerical purposes it may as well be.

```{code-cell}
using Plots
U = zeros(m, n+1);
U[:, 1] = @. exp( -60 * (x - 0.5)^2 )
plot(x, U[:, 1];
    xaxis=(L"x"),  yaxis=(L"u(x,0)"),
    title="Initial condition")
```

The Euler time stepping simply multiplies $\mathbf{u}_j$ by the constant matrix in {eq}`Eulerxx` at each time step. Since that matrix is sparse, we will declare it as such, even though the run-time savings may not be detectable for this small value of $m$.

```{code-cell}
using SparseArrays
A = sparse(I + τ * Dxx)
for j in 1:n
    U[:, j+1] = A * U[:, j]
end

plot_idx = 1:10:31
plot_times = round.(t[plot_idx], digits=4)
labels = ["t = $t" for t in plot_times]
plot(x, U[:, plot_idx];
    label=reshape(labels, 1, :),  legend=:topleft,  
    title="Heat equation by forward Euler",
    xaxis=(L"x"),  yaxis=(L"u(x,0)", [-0.25, 1]))
```

Things seem to start well, with the initial peak widening and shrinking. But then there is a nonphysical growth in the solution.

```{code-cell}
:tags: hide-input, remove-output
anim = @animate for j in 1:101
    plot(x, U[:, j];
    label=@sprintf("t=%.5f", t[j]),
    xaxis=(L"x"),  yaxis=(L"u(x,t)", [-1, 2]),
    dpi=150,  title="Heat equation by forward Euler")
end
mp4(anim, "figures/diffusionFE.mp4")
```

![Instability in Euler solution](figures/diffusionFE.mp4)

The growth in norm is exponential in time.

```{code-cell}
M = vec( maximum(abs, U, dims=1) )   
plot(t[1:1000], M[1:1000];
    xaxis=(L"t"),  yaxis=(:log10, L"\max_x |u(x,t)|"),
    title="Nonphysical growth") 
```

::::

The method in {numref}`Example {number} <example-methodlines-heatFE>` and @demo-methodlines-heatFE is essentially the same one we used for the Black–Scholes equation in {numref}`section-diffusion-blackscholes`. By changing the time integrator, we can get much better results.

::::{prf:example}
:label: example-methodlines-heatBE
An alternative time discretization of {eq}`heatMOL` is to use the backward Euler (AM1) method, resulting in

```{math}
:label: BExx
\begin{split}
    \mathbf{u}_{j+1} &= \mathbf{u}_j + \tau (\mathbf{D}_{xx} \mathbf{u}_{j+1})\\
    (\mathbf{I} - \tau \mathbf{D}_{xx}) \mathbf{u}_{j+1} &= \mathbf{u}_j.
\end{split}
```

Because backward Euler is an implicit method, a linear system must be solved for $\mathbf{u}_{j+1}$ at each time step.
::::

::::{prf:example} Backward Euler for the heat equation
:label: demo-methodlines-heatBE

Now we apply backward Euler to the heat equation. We will reuse the setup from @demo-methodlines-heatFE. Since the matrix in {eq}`BExx` never changes during the time stepping, we do the necessary LU factorization only once.

```{code-cell}
using SparseArrays
B = sparse(I - τ * Dxx)
factor = lu(B)
for j in 1:n
    U[:, j+1] = factor \ U[:, j]
end
```

```{code-cell}
:tags: hide-input
using Plots
idx = 1:600:n+1
times = round.(t[idx], digits=4)
label = reshape(["t = $t" for t in times], 1, length(idx))
plot(x,U[:, idx];
    label, legend=:topleft,
    title="Heat equation by backward Euler",
    xaxis=(L"x"),  yaxis=(L"u(x,0)", [0, 1]))
```

```{code-cell}
:tags: hide-input, remove-output
anim = @animate for j in 1:20:n+1
    plot(x, U[:, j];
    label=@sprintf("t=%.5f", t[j]),
    xaxis=(L"x"),  yaxis=(L"u(x,t)", [0, 1]),
    dpi=150,  title="Heat equation by backward Euler")
end
mp4(anim, "figures/diffusionBE.mp4")
```
![Stable Backward Euler solution](figures/diffusionBE.mp4)

This solution looks physically plausible, as the large concentration in the center diffuses outward until the solution is essentially constant. Observe that the solution remains periodic in space for all time.

::::

@demo-methodlines-heatBE suggests that implicit time stepping methods have an important role in diffusion. We will analyze the reason in the next few sections.

## Black-box IVP solvers

Instead of coding one of the Runge–Kutta or multistep formulas directly for a method of lines solution, we could use any of the IVP solvers from Chapter 6, or a solver from the `DifferentialEquations` package, to solve the ODE initial-value problem {eq}`heatMOL`.

::::{prf:example} Adaptive time stepping for the heat equation
:label: demo-methodlines-auto

We set up the semidiscretization and initial condition in $x$ just as before.

```{code-cell}
m = 100
x, Dx, Dxx = FNC.diffper(m, [0, 1])
u0 = @. exp( -60*(x - 0.5)^2 );
```

Now, however, we apply {numref}`Function {number} <function-rk23>` (`rk23`) to the initial-value problem $\mathbf{u}'=\mathbf{D}_{xx}\mathbf{u}$.

```{code-cell}
using OrdinaryDiffEq
tfinal = 0.25
ODE = (u, p, t) -> Dxx * u  
IVP = ODEProblem(ODE, u0, (0, tfinal))
t, u = FNC.rk23(IVP, 1e-5);
```

We check that the resulting solution looks realistic.

```{code-cell}
:tags: hide-input
plt = plot(
    title="Heat equation by rk23",
    legend=:topleft,  
    xaxis=(L"x"),  yaxis=(L"u(x,0)", [0, 1]))
for idx in 1:600:n+1
    plot!(x, u[idx]; label="t = $(round.(t[idx], digits=4))")
end
plt
```

```{code-cell}
:tags: hide-input, remove-output
anim = @animate for j in 1:20:1600
    plot(x, u[j];
    label=@sprintf("t=%.4f", t[j]),
      xaxis=(L"x"),  yaxis=(L"u(x,t)", [0, 1]),
      dpi=150,  title="Heat equation by rk23")
end
mp4(anim, "figures/diffusionRK23.mp4")
```

![RK23 solution](figures/diffusionRK23.mp4)


The solution appears to be correct. But the number of time steps that were selected automatically is surprisingly large, considering how smoothly the solution changes.

```{code-cell}
println("Number of time steps for rk23: $(length(t)-1)")
```

Now we apply a solver from `DifferentialEquations`.

```{code-cell}
u = solve(IVP, Rodas4P());
println("Number of time steps for Rodas4P: $(length(u.t) - 1)")
```

The number of steps selected is reduced by a factor of more than 100!

::::

The adaptive time integrators can all produce solutions. But, as seen in @demo-methodlines-auto, they are not equivalent in every important sense. Whether we choose to implement a method directly with a fixed step size, or automatically with adaptation, there is something crucial to understand about the semidiscrete problem {eq}`heatMOL` that will occupy our attention in the next two sections.

## Exercises

``````{exercise}
:label: problem-methodlines-heatstability
⌨ Revisit @demo-methodlines-heatFE. For each $m=20,30,\dots,120$ points in space, let $n=20,30,40,\dots$ in turn until you reach the smallest $n$ such that the numerical solution remains bounded above by 2 for all time; call this value $N(m)$. Make a log-log plot of $N(m)$ as a function of $m$. If you suppose that $N=O(m^p)$ for a simple rational number $p$, what is a reasonable hypothesis for $p$?
``````

``````{exercise}
:label: problem-methodlines-average
In @demo-methodlines-auto, as $t\to \infty$ the solution $u(x,t)$ approaches a value that is constant in both space and time.

**(a)** ⌨ Set $m=400$ and use a native IVP solver, as shown in @demo-methodlines-auto, to find this constant value to at least eight digits of accuracy.

**(b)** ✍ Prove that $Q = \int_0^1 u(x,t) \,dx$ is constant in time. (Hint: If its derivative is zero, then it is constant. Take the derivative of the integral with respect to $t$ and apply the PDE and periodic end conditions.)

**(c)** ⌨ Use {numref}`Function {number} <function-trapezoid>` to find $Q$ at $t=0$, and compare this number to the result of part (a).
``````

```{index} Crank–Nicolson method
```

``````{exercise}
:label: problem-methodlines-cranknicolson
✍ Apply the trapezoid IVP formula (AM2 in @table-adams) to the semidiscretization {eq}`heatMOL` and derive what is known as the *Crank–Nicolson* method:

```{math}
:label: CNxx
(\mathbf{I} - \tfrac{1}{2}\tau \mathbf{D}_{xx}) \mathbf{u}_{j+1} =  (\mathbf{I} + \tfrac{1}{2}\tau
\mathbf{D}_{xx}) \mathbf{u}_{j}.
```

Note that each side of the method is evaluated at a different time level.
``````

``````{exercise}
:label: problem-methodlines-cnstable
⌨ Repeat @demo-methodlines-heatBE using the Crank–Nicolson method {eq}`CNxx`. Then try for $n=240$ as well, which uses a time step ten times larger than before. Does the solution remain stable? 
``````

``````{exercise}
:label: problem-methodlines-reactdiff
The PDE $u_t = 2u + u_{xx}$ combines growth with diffusion. 

**(a)** ✍ Derive an equation analogous to {eq}`BExx` that combines second-order semidiscretization in space with the backward Euler solver in time.

**(b)** ⌨ Apply your formula from part (a) to solve this PDE with periodic boundary conditions for the same initial condition as in @demo-methodlines-heatBE. Use  $m=200$ points in space and $n=1000$ time levels. Plot the solution on one graph at times $t=0,0.04,0.08,\ldots,0.2$, or animate the solution over $0\le t \le 0.2$.
``````

``````{exercise}
:label: problem-methodlines-eulerxx
✍ In this problem, you will analyze the convergence of the explicit method given by {eq}`Eulerxx`.  Recall that the discrete approximation $u_{i,j}$ approximates the solution at $x_i$ and $t_j$.

**(a)** Write the method in scalar form as

$$
u_{i,j+1} = (1-2\lambda) u_{i,j} + \lambda u_{i+1,j} + \lambda u_{i-1,j},
$$

where $\lambda = \tau/h^2>0$.

**(b)** Taylor series of the exact solution $\hat{u}$ imply that

\begin{align*}
\hat{u}_{i,j+1} &= u_{i,j} + \frac{\partial \hat{u}}{\partial t} (x_i,t_j) \tau + O(\tau^2),\\
% \frac{\partial^2 u}{\partial t^2} (x_i,\bar{t}) \frac{\tau^2}{2}
\hat{u}_{i\pm1,j} &= \hat{u}_{i,j} \pm \frac{\partial \hat{u}}{\partial x} (x_i,t_j) h + \frac{\partial^2 \hat{u}}{\partial x^2} (x_i,t_j)
\frac{h^2}{2} \pm \frac{\partial^3 \hat{u}}{\partial x^3} (x_i,t_j)
\frac{h^3}{6}+  O(h^4).
%\frac{\partial^4 u}{\partial x^4} (\bar{x}_\pm,t_j) \frac{h^4}{24}.
\end{align*}

Use these to show that

\begin{align*}
\hat{u}_{i,j+1} & = \left[ (1-2\lambda) \hat{u}_{i,j} + \lambda \hat{u}_{i+1,j} + \lambda \hat{u}_{i-1,j}\right]
+  O\Bigl(\tau^2+h^2 \Bigr)\\
&= F\left( \lambda,\hat{u}_{i,j}, \hat{u}_{i+1,j} , \hat{u}_{i-1,j}\right) + O\Bigl(\tau^2+h^2\Bigr).
\end{align*}

(The last line should be considered a definition of the function $F$.)

**(c)** The numerical solution satisfies 

$$
u_{i,j+1}=F\bigl( \lambda,u_{i,j}, u_{i+1,j} , u_{i-1,j}\bigr) 
$$ 

exactly. Using this fact, subtract $u_{i,j+1}$ from both sides of the last line in part (b) to show that

$$
e_{i,j+1} = F\left( \lambda,e_{i,j}, e_{i+1,j} ,e_{i-1,j}\right)  + O\Bigl(\tau^2+h^2\Bigr),
$$

where $e_{i,j}=\hat{u}_{i,j}-u_{i,j}$ is the error in the numerical solution for all $i$ and $j$ .

**(d)** Define $E_j$ as the maximum of $|e_{i,j}|$ over all values of $i$, and use the result of part (c) to show that if $\lambda<1/2$ is kept fixed as $h$ and $\tau$ approach zero, then for sufficiently small $\tau$ and $h$,

$$
E_{j+1}  = E_{j} + O\Bigl(\tau^2+h^2\Bigr) \le E_{j} + K_j\bigl(\tau^2+h^2\bigr)
$$

for a positive $K_j$ independent of $\tau$ and $h$.

**(e)** If the initial conditions are exact, then $E_0=0$. Use this to show finally that if the $K_j$ are bounded above and $\lambda<1/2$ is kept fixed, then $E_n = O(\tau)$ as $\tau\to 0$.
``````
