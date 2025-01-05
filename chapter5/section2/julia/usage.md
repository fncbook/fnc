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
[**Demo %s**](#demo-pwlin-usage)

We generate a piecewise linear interpolant of $f(x)=e^{\sin 7x}$.

```{code-cell}
f = x -> exp(sin(7x))

plot(f, 0, 1, label="function", xlabel=L"x", ylabel=L"y")
```

First we sample the function to create the data.

```{code-cell}
t = [0, 0.075, 0.25, 0.55, 0.7, 1]    # nodes
y = f.(t)                             # function values
scatter!(t, y, label="values at nodes")
```

Now we create a callable function that will evaluate the piecewise linear interpolant at any $x$, and then plot it.

```{code-cell}
p = FNC.plinterp(t, y)
plot!(p, 0, 1, label="interpolant", title="PL interpolation")
```
