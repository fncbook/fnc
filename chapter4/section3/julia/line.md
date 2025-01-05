---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{code-cell}
include("../../../julia/FNC_init.jl")
```
[**Demo %s**](#demo-newton-line)


Suppose we want to find a root of the function

```{code-cell}
f(x) = x * exp(x) - 2
plot(f, 0, 1.5; 
    label="function",  legend=:topleft,
    grid=:y,  ylim=[-2, 4],  xlabel=L"x",  ylabel=L"y")
```

From the graph, it is clear that there is a root near $x=1$. So we call that our initial guess, $x_1$.

```{code-cell}
x₁ = 1
y₁ = f(x₁)
scatter!([x₁], [y₁], label="initial point")
```

Next, we can compute the tangent line at the point $\bigl(x_1,f(x_1)\bigr)$, using the derivative.

```{code-cell}
df_dx(x) = exp(x) * (x + 1)
m₁ = df_dx(x₁)
tangent = x -> y₁ + m₁ * (x - x₁)

plot!(tangent, 0, 1.5, l=:dash, label="tangent line",
    title="Tangent line approximation")
```

In lieu of finding the root of $f$ itself, we settle for finding the root of the tangent line approximation, which is trivial. Call this $x_2$, our next approximation to the root.

```{code-cell}
@show x₂ = x₁ - y₁ / m₁
scatter!([x₂], [0], label="tangent root", title="First iteration")
```

```{code-cell}
y₂ = f(x₂)
```

The residual (i.e., value of $f$) is smaller than before, but not zero. So we repeat the process with a new tangent line based on the latest point on the curve.

```{code-cell}
plot(f, 0.82, 0.87;
    label="function",  legend=:topleft,
    xlabel=L"x",  ylabel=L"y",
    title="Second iteration")

scatter!([x₂], [y₂], label="starting point")

m₂ = df_dx(x₂)
tangent = x -> y₂ + m₂ * (x - x₂)
plot!(tangent, 0.82, 0.87; l=:dash, label="tangent line")

@show x₃ = x₂ - y₂ / m₂
scatter!([x₃], [0], label="tangent root")
```

```{code-cell}
y₃ = f(x₃)
```

Judging by the residual, we appear to be getting closer to the true root each time.
