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
[**Demo %s**](#demo-quasi-levenberg)

To solve a nonlinear system, we need to code only the function defining the system, and not its Jacobian.

```{code-cell}
f(x) = 
    [
        exp(x[2] - x[1]) - 2,
        x[1] * x[2] + x[3],
        x[2] * x[3] + x[1]^2 - x[2]
    ]
```

In all other respects usage is the same as for the `newtonsys` function.

```{code-cell}
x₁ = [0.0, 0.0, 0.0]
x = FNC.levenberg(f, x₁)
```

It's always a good idea to check the accuracy of the root, by measuring the residual (backward error).

```{code-cell}
r = x[end]
println("backward error = $(norm(f(r)))")
```

Looking at the convergence in norm, we find a convergence rate between linear and quadratic, like with the secant method.

```{code-cell}
logerr = [log(norm(r - x[k])) for k in 1:length(x)-1]
ratios = [NaN; [logerr[i+1] / logerr[i] for i in 1:length(logerr)-1]]
@pt :header=["iteration", "log error", "ratio"] [eachindex(logerr) logerr ratios]
```
