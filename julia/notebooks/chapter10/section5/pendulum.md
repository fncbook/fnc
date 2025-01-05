---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
The first step is to define the function $\phi$ that equals $\theta''$.

```{code-cell}
ϕ = (t, θ, ω) -> -0.05 * ω - sin(θ);
```

Next, we define the boundary conditions.

```{code-cell}
g₁(u, du) = u - 2.5
g₂(u, du) = u + 2;
```

```{index} ! Julia; collect
```

The last ingredient is an initial estimate of the solution. Here we choose $n=100$ and a linear function between the endpoint values. 
```{tip}
The `collect` function turns a range object into a true vector.
```

```{code-cell}
init = collect(range(2.5, -2, length = 101));
```

We find a solution with negative initial slope, i.e., the pendulum is initially pushed back toward equilibrium.

```{code-cell}
using Plots
t, θ = FNC.bvp(ϕ, [0, 5], g₁, g₂, init)
plot(t, θ;
    xaxis=(L"t"),  yaxis=(L"\theta(t)"),
    title="Pendulum over [0,5]" )
```

If we extend the time interval longer for the same boundary values, then the initial slope must adjust.

```{code-cell}
t, θ = FNC.bvp(ϕ, [0, 8], g₁, g₂, init)
plot(t, θ;
    xaxis=(L"t"),  yaxis=(L"\theta(t)"),
    title="Pendulum over [0,8]" )
```

This time, the pendulum is initially pushed toward the unstable equilibrium in the upright vertical position before gravity pulls it back down.
