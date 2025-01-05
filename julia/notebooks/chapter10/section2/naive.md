---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
Let's first examine the shooting approach for the TPBVP from {numref}`Example {number} <example-tpbvp-mems>` with $\lambda=0.6$. 
```{tip}
The character `ϕ` is typed as `\phi`<kbd>Tab</kbd>. 
```

```{code-cell}
λ = 0.6
ϕ = (r, w, dwdr) -> λ / w^2 - dwdr / r;
```

We convert the ODE to a first-order system in order to apply a numerical method. We also have to truncate the domain to avoid division by zero.

```{code-cell}
f = (y, p, r) -> [y[2]; ϕ(r, y[1], y[2])]
a, b = eps(), 1.0;
```

The BVP specifies $w'(0)=y_2(0)=0$. We can try multiple values for the unknown $w(0)=y_1(0)$ and plot the solutions.

```{code-cell}
using OrdinaryDiffEq, Plots
plt = plot(
    xaxis = (L"x"),  yaxis = (L"w(x)"),
    title = "Different initial values",  legend = :bottomright)

for w0 in 0.4:0.1:0.9
    IVP = ODEProblem(f, [w0, 0], (a, b))
    y = solve(IVP, Tsit5())
    plot!(y, idxs = [1], label = "w(0) = $w0")
end
plt
```

On the graph, it's the curve starting at $w(0)=0.8$ that comes closest to the required condition $w(1)=1$, but it's a bit too large.
