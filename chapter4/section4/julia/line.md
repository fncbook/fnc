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
[**Demo %s**](#demo-secant-line)



We return to finding a root of the equation $x e^x=2$.

```{code-cell}
f(x) = x * exp(x) - 2;
plot(f, 0.25, 1.25;
    label="function",  legend=:topleft,
    xlabel=L"x",  ylabel=L"y")
```

From the graph, it's clear that there is a root near $x=1$. To be more precise, there is a root in the interval $[0.5,1]$. So let us take the endpoints of that interval as _two_ initial approximations.

```{code-cell}
x₁ = 1;
y₁ = f(x₁);
x₂ = 0.5;
y₂ = f(x₂);
scatter!([x₁, x₂], [y₁, y₂];
    label="initial points",
    title="Two initial values")
```

Instead of constructing the tangent line by evaluating the derivative, we can construct a linear model function by drawing the line between the two points $\bigl(x_1,f(x_1)\bigr)$ and $\bigl(x_2,f(x_2)\bigr)$. This is called a _secant line_.

```{code-cell}
m₂ = (y₂ - y₁) / (x₂ - x₁)
secant = x -> y₂ + m₂ * (x - x₂)
plot!(secant, 0.25, 1.25, label="secant line", l=:dash, color=:black,
    title="Secant line")
```

As before, the next root estimate in the iteration is the root of this linear model.

```{code-cell}
x₃ = x₂ - y₂ / m₂
@show y₃ = f(x₃)
scatter!([x₃], [0], label="root of secant", title="First iteration")
```

For the next linear model, we use the line through the two most recent points. The next iterate is the root of that secant line, and so on.

```{code-cell}
m₃ = (y₃ - y₂) / (x₃ - x₂)
x₄ = x₃ - y₃ / m₃
```
