---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---
```{code-cell}
:tags: [remove-cell]
include("../../../julia/FNC_init.jl")
```
[**Demo %s**](#demo-stiffness-explicit)

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
:tags: [hide-input]
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

