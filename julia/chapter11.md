---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---

# Chapter 11

## Functions

(function-diffper-julia)=
``````{dropdown} Differentiation matrices for periodic end conditions
:open:
```{literalinclude} FNCFunctions/src/chapter11.jl
:filename: diffper.jl
:start-after: # begin diffper
:end-before: # end diffper
:language: julia
:linenos: true
```
``````

(function-parabolic-julia)=
``````{dropdown} Solution of parabolic PDEs by the method of lines
:open:
```{literalinclude} FNCFunctions/src/chapter11.jl
:filename: parabolic.jl
:start-after: # begin parabolic
:end-before: # end parabolic
:language: julia
:linenos: true
```

::::{admonition} About the code
:class: dropdown
Line 29 uses the macro `@.` to assign into the vector `f` elementwise. Without it, the function would allocate space for the result of `phi` and then change `f` to point at that vector, and that would defeat the purpose of using the preallocated `f` for speed.
::::
``````

## Examples

```{code-cell}
:tags: remove-output
include("FNC_init.jl")
```

### 11.1 @section-diffusion-blackscholes

(demo-blackscholes-solve-julia)=
``````{dropdown} @demo-blackscholes-solve
We consider the Black–Scholes problem for the following parameter values:

```{code-cell}
Smax = 8 
T = 6
K, σ, r = (3, 0.06, 0.08);
```

We discretize space and time.

```{code-cell}
m = 200;  h = Smax / m;
x = h * (0:m)
n = 1000;  τ = T / n;
t = τ * (0:n)
λ = τ / h^2
μ = τ / h
```

We set the initial condition and then march forward in time.

```{code-cell}
V = zeros(m+1, n+1)
V[:, 1] = @. max(0, x - K)
for j in 1:n
    # Fictitious value from Neumann condition.
    Vfict = 2*h + V[m,j]
    Vj = [ V[:, j]; Vfict ]
    # First row is zero by the Dirichlet condition.
    for i in 2:m+1 
        diff1 = (Vj[i+1] - Vj[i-1])
        diff2 = (Vj[i+1] - 2Vj[i] + Vj[i-1])
        V[i,j+1] = Vj[i] +
            (λ * σ^2 * x[i]^2 / 2) * diff2 
            + (r * x[i] * μ) / 2 * diff1 
            - (r * τ) * Vj[i]
    end 
end
```

Here is a plot of the solution after every 250 time steps.

```{code-cell}
using Plots
idx = 1:250:n+1
label = reshape(["t = $t" for t in t[idx]], 1, length(idx))
plot(x, V[:, idx]; 
    label, legend=:topleft,
    xaxis=("stock price"),  yaxis=("option value"),
    title="Black–Scholes solution")
```

```{index} ! Julia; @animate
```

Alternatively, here is an animation of the solution.

```{code-cell}
anim = @animate for j in 1:10:n+1
    plot(x, V[:, j];
        xaxis=(L"S"),  yaxis=([0,6],L"v(S,t)"),
        title="Black–Scholes solution",
        dpi=150,    
        label=@sprintf("t = %.2f", t[j]))
end
mp4(anim, "figures/black-scholes-6.mp4")
```

The results are easy to interpret, recalling that the time variable really means *time until strike*. Say you are close to the option's strike time. If the current stock price is, say, $S=2$, then it's not likely that the stock will end up over the strike price $K=3$, and therefore the option has little value. On the other hand, if presently $S=3$, then there are good odds that the option will be exercised at the strike time, and you will need to pay a substantial portion of the stock price in order to take advantage. As the time to strike increases, there is an expectation that the stock price is more likely to rise somewhat, making the value of the option larger at each fixed $S$. 
``````

(demo-blackscholes-unstable-julia)=
``````{dropdown} @demo-blackscholes-unstable
Let's try to do everything the same as in {numref}`Demo {number} <demo-blackscholes-solve>`, but extending the simulation time to $T=8$.

```{code-cell}
T = 8;

m = 200;  h = Smax / m;
x = h*(0:m)
n = 1000;  τ = T / n;
t = τ*(0:n)
λ = τ / h^2;  μ = τ / h;

for j in 1:n
    # Fictitious value from Neumann condition.
    Vfict = 2h + V[m,j]
    Vj = [ V[:, j]; Vfict ]
    # First row is zero by the Dirichlet condition.
    for i in 2:m+1 
        diff1 = (Vj[i+1] - Vj[i-1])
        diff2 = (Vj[i+1] - 2Vj[i] + Vj[i-1])
        V[i,j+1] = Vj[i] +
            (λ * σ^2 * x[i]^2 / 2) * diff2 
            + (r * x[i] * μ) / 2 * diff1 
            - (r * τ) * Vj[i]
    end   
end

idx = 1:250:n+1
label = reshape(["t = $t" for t in t[idx]], 1, length(idx))
plot(x, V[:, idx];
    label, legend=:topleft,
    title="Black–Scholes solution",
    xaxis=("stock price"),  yaxis=("option value",[0, 6]))
```

```{code-cell}
anim = @animate for j in 1:10:n+1 
    plot(x, V[:, j];
        xaxis=(L"S"),  yaxis=([0,6],L"v(S,t)"),
        title="Black–Scholes solution...?",
        dpi=150,  label=@sprintf("t = %.2f",t[j]))
end
mp4(anim, "figures/black-scholes-8.mp4")
```

This so-called solution is nonsense!
``````

### 11.2 @section-diffusion-methodlines

(demo-methodlines-heatFE-julia)=
``````{dropdown} @demo-methodlines-heatFE
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
:tags: hide-input
anim = @animate for j in 1:101
    plot(x, U[:, j];
    label=@sprintf("t=%.5f", t[j]),
    xaxis=(L"x"),  yaxis=(L"u(x,t)", [-1, 2]),
    dpi=150,  title="Heat equation by forward Euler")
end
mp4(anim, "figures/diffusionFE.mp4")
```

The growth in norm is exponential in time.

```{code-cell}
M = vec( maximum(abs, U, dims=1) )   
plot(t[1:1000], M[1:1000];
    xaxis=(L"t"),  yaxis=(:log10, L"\max_x |u(x,t)|"),
    title="Nonphysical growth") 
```
``````

(demo-methodlines-heatBE-julia)=
``````{dropdown} @demo-methodlines-heatBE
Now we apply backward Euler to the heat equation. We will reuse the setup from {numref}`Demo {number} <demo-methodlines-heatFE>`. Since the matrix in {eq}`BExx` never changes during the time stepping, we do the necessary LU factorization only once.

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
:tags: hide-input
anim = @animate for j in 1:20:n+1
    plot(x, U[:, j];
    label=@sprintf("t=%.5f", t[j]),
    xaxis=(L"x"),  yaxis=(L"u(x,t)", [0, 1]),
    dpi=150,  title="Heat equation by backward Euler")
end
mp4(anim, "figures/diffusionBE.mp4")
```

This solution looks physically plausible, as the large concentration in the center diffuses outward until the solution is essentially constant. Observe that the solution remains periodic in space for all time.
``````

(demo-methodlines-auto-julia)=
``````{dropdown} @demo-methodlines-auto
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
:tags: hide-input
anim = @animate for j in 1:20:1600
    plot(x, u[j];
    label=@sprintf("t=%.4f", t[j]),
      xaxis=(L"x"),  yaxis=(L"u(x,t)", [0, 1]),
      dpi=150,  title="Heat equation by rk23")
end
mp4(anim, "figures/diffusionRK23.mp4")
```

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
``````

### 11.3 @section-diffusion-absstab

(demo-absstab-regions-julia)=
``````{dropdown} @demo-absstab-regions
Euler and Backward Euler time-stepping methods were used to solve $\mathbf{u}'=\mathbf{D}_{xx}\mathbf{u}$.

```{code-cell}
m = 40
_, _, Dₓₓ = FNC.diffper(m, [0, 1]);
```

The eigenvalues of this matrix are real and negative:

```{code-cell}
using Plots
λ = eigvals(Dₓₓ)
scatter(real(λ), imag(λ);
    title="Eigenvalues",
    frame=:zerolines,  aspect_ratio=1,
    xaxis=("Re λ"),  yaxis=("Im λ", (-1000, 1000)))
```

The Euler method is absolutely stable in the region $|\zeta+1| \le 1$ in the complex plane:

```{code-cell}
:tags: hide-input
phi = 2π * (0:360) / 360
z = @. cis(phi) - 1;    # unit circle shifted to the left by 1

plot(Shape(real(z), imag(z));
    color=RGB(.8, .8, 1),
    xaxis=("Re ζ"),  yaxis=("Im ζ"),
    aspect_ratio=1,  frame=:zerolines,
    title="Stability region") 
```

In order to get inside this region, we have to find $\tau$ such that $\lambda \tau > -2$ for all eigenvalues $\lambda$. This is an upper bound on $\tau$. 

```{code-cell}
λ_min = minimum(λ)
@show max_τ = -2 / λ_min;
```

Here we plot the resulting values of $\zeta=\lambda \tau$. 

```{code-cell}
ζ = λ * max_τ
scatter!(real(ζ), imag(ζ), title="Stability region and ζ values")
```

In backward Euler, the region is $|\zeta-1|\ge 1$. Because they are all on the negative real axis, all of the $\zeta$ values will fit no matter what $\tau$ is chosen.

```{code-cell}
:tags: hide-input
plot(Shape([-6, 6, 6, -6], [-6, -6, 6, 6]), color=RGB(.8, .8, 1))
z = @. cis(phi) + 1;   # unit circle shifted right by 1
plot!(Shape(real(z), imag(z)), color=:white)

scatter!(real(ζ), imag(ζ);
    xaxis=([-4, 2], "Re ζ"),  yaxis=([-3, 3], "Im ζ"),
    aspect_ratio=1,  frame=:zerolines,
    title="Stability region and ζ values")
```
``````

### 11.4 @section-diffusion-stiffness

(demo-stiffness-oregon-julia)=
``````{dropdown} @demo-stiffness-oregon
In {numref}`Example {number} <example-stiffness-oregon>` we derived a Jacobian matrix for the Oregonator model. Here is a numerical solution of the ODE.

```{code-cell}
:tags: hide-input
using OrdinaryDiffEq, Plots
function ode(u,p,t)
    s,w,q = p
    f = [ 
        s * ( u[2]*(1 - u[1]) + u[1]*(1 - q*u[1]) ),
        (u[3] - u[2] - u[1] * u[2]) / s,   
        w * (u[1] - u[3])
        ]
    return f
end
s, w, q = (77.27, .161, 8.375e-6)
oregon = ODEProblem(ode, [1., 2, 3], (0., 500.), [s, w, q])
sol = solve(oregon)
plot(sol, yscale=:log10, legend=:none, title="Solution of the Oregonator")
```

At each value of the numerical solution, we can compute the eigenvalues of the Jacobian. Here we plot all of those eigenvalues in the complex plane.

```{code-cell}
:tags: hide-input
t,u = sol.t[1:2:end], sol.u[1:2:end]
λ = fill(0.0im, length(t), 3)
for (k, u) in enumerate(u)
    J = [
    s*(1-u[2]-2q*u[1]) s*(1-u[1])        0 
         -u[2]/s       -(1+u[1])/s     1/s 
            w               0           -w
        ]
    λ[k, :] .= eigvals(J)
end

scatter(real(λ), imag(λ), t;
    xaxis=("Re(λ)", 25000*(-5:2:-1)),  ylabel="Im(λ)",  zlabel="t",
    title="Oregonator eigenvalues")
```

You can see that there is one eigenvalue that ranges over a wide portion of the negative real axis and dominates stability considerations.
``````

(demo-stiffness-explicit-julia)=
``````{dropdown} @demo-stiffness-explicit
The `Rodas4P` solver is good for stiff problems, and needs few time steps to solve the Oregonator from {numref}`Demo {number} <demo-stiffness-oregon>`.

```{code-cell}
oregon = remake(oregon, tspan=(0., 25.))
sol = solve(oregon, Rodas4P())
println("Number of time steps for Rodas4P: $(length(sol.t) - 1)")
```

But if we apply {numref}`Function {number} <function-rk23>` to the problem, the step size will be made small enough to cope with the large negative eigenvalue. 

```{code-cell}
t,u = FNC.rk23(oregon,1e-4)
println("Number of time steps for RK23: $(length(t) - 1)")
```

Starting from the eigenvalues of the Jacobian matrix, we can find an effective $\zeta(t)$ by multiplying with the local time step size. The values of $\zeta(t)$ for each time level are plotted below and color coded by component of the diagonalized system.

```{code-cell}
:tags: hide-input
λ = fill(1.0im, length(t),3)
for (k, u) in enumerate(u)
    J = [
    s*(1-u[2]-2q*u[1]) s*(1-u[1])        0 
         -u[2]/s       -(1+u[1])/s     1/s 
            w               0           -w
        ]
    λ[k, :] .= eigvals(J)
end

ζ = diff(t) .* λ[1:end-1,:]
scatter(real(ζ), imag(ζ), m=2,
    xlabel="Re(ζ)",  ylabel="Im(ζ)",
    title="Oregonator stability")
```

Roughly speaking, the $\zeta$ values stay within or close to the RK2 stability region in {numref}`figure-stabreg_bd_rk`. Momentary departures from the region are possible, but time stepping repeatedly in that situation would cause instability. 

``````

### 11.5 @section-diffusion-boundaries

(demo-boundaries-heat-julia)=
``````{dropdown} @demo-boundaries-heat
First, we define functions for the PDE and each boundary condition.

```{code-cell}
ϕ = (t, x, u, uₓ, uₓₓ) -> uₓₓ
g₁ = (u, uₓ) -> u
g₂ = (u, uₓ) -> u - 2;
```

Our next step is to write a function to define the initial condition. This one satisfies the boundary conditions exactly.

```{code-cell}
init = x -> 1 + sinpi(x/2) + 3 * (1-x^2) * exp(-4x^2);
```

Now we can use {numref}`Function {number} <function-parabolic>` to solve the problem.

```{code-cell}
using Plots
x, u = FNC.parabolic(ϕ, (-1, 1), 60, g₁, g₂, (0, 0.75), init)
plt = plot(
    xlabel=L"x",  ylabel=L"u(x,t)",
    legend=:topleft,  title="Solution of the heat equation")
for t in 0:0.1:0.4
    plot!(x, u(t), label=@sprintf("t=%.2f", t))
end
plt
```

```{code-cell}
:tags: hide-input
anim = @animate for t in range(0,0.75,length=201) 
    plot(x, u(t);
        label=@sprintf("t=%.2f", t),  legend=:topleft,
        xaxis=(L"x"),  yaxis=(L"u(x,t)", (0, 4.2)), 
        title="Heat equation",  dpi=150)
end
mp4(anim, "figures/boundaries-heat.mp4", fps=30)
```
``````

(demo-boundaries-bratu-julia)=
``````{dropdown} @demo-boundaries-bratu

```{code-cell}
ϕ = (t, x, u, uₓ, uₓₓ) -> u^2 + uₓₓ
g₁ = (u, uₓ) -> u
g₂ = (u, uₓ) -> uₓ
init = x -> 400x^4 * (1 - x)^2
x, u = FNC.parabolic(ϕ, (0, 1), 60, g₁, g₂, (0, 0.1), init);
```

```{code-cell}
:tags: hide-input
anim = @animate for t in range(0, 0.1, length=101) 
    plot(x, u(t);
        label=@sprintf("t=%.4f", t),  legend=:topleft,
        xaxis=(L"x"),  yaxis=(L"u(x,t)", (0, 10)),
        dpi=150, title="Heat equation with source")
end
mp4(anim, "figures/boundaries-source.mp4", fps=30)
```
``````

(demo-boundaries-bs-julia)=
``````{dropdown} @demo-boundaries-bs

```{code-cell}
K = 3;  σ = 0.06;  r = 0.08;  Smax = 8;
ϕ = (t, x, u, uₓ, uₓₓ) -> σ^2/2 * (x^2 * uₓₓ) + r*x*uₓ - r*u
g₁ = (u, uₓ) -> u
g₂ = (u, uₓ) -> uₓ - 1;
```

```{code-cell}
u₀ = x -> max(0, x - K)
x, u = FNC.parabolic(ϕ, (0, Smax), 80, g₁, g₂, (0, 15), u₀);
```

```{code-cell}
:tags: hide-input
anim = @animate for t in range(0, 15, 151) 
    plot(x, u(t);
        label=@sprintf("t=%.4f", t),  legend=:topleft,
        xaxis=(L"x"),  yaxis=(L"u(x,t)", (-0.5, 8)), 
        dpi=150,  title="Black–Scholes equation")
end
mp4(anim, "figures/boundaries-bs.mp4", fps=30)
```

Recall that $u$ is the value of the call option, and time runs backward from the strike time. The longer the horizon, the more value the option has due to anticipated growth in the stock price.
``````
