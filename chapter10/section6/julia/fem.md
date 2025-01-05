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
[**Demo %s**](#demo-galerkin-fem)


Here are the coefficient function definitions. Even though $s$ is a constant, it has to be defined as a function for {numref}`Function {number} <function-fem>` to use it.

```{code-cell}
c = x -> x^2;
q = x -> 4;
f = x -> sin(π * x);
```

```{code-cell}
x, u = FNC.fem(c, q, f, 0, 1, 50)
plot(x, u;
    xaxis=(L"x"),  yaxis = (L"u"),
    title = "Solution by finite elements", legend=:none)
```
