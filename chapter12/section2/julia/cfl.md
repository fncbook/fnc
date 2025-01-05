---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
include("../../../julia/FNC_init.jl")
```
[**Demo %s**](#demo-upwind-cfl)

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
:tags: [hide-input]
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
