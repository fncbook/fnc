---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---
In {numref}`Example {number} <example-stiffness-oregon>` we derived a Jacobian matrix for the Oregonator model. Here is a numerical solution of the ODE.

```{code-cell}
:tags: [hide-input]
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
:tags: [hide-input]
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
