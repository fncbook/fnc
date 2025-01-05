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
[**Demo %s**](#demo-qr-stable)

We'll repeat the experiment of {numref}`Demo {number} <demo-normaleqns-instab>`, which exposed instability in the normal equations. 

```{code-cell}
t = range(0, 3, 400)
f = [ x -> sin(x)^2, x -> cos((1 + 1e-7) * x)^2, x -> 1. ]
A = [ f(t) for t in t, f in f ]
x = [1., 2, 1]
b = A * x;
```

The error in the solution by {numref}`Function {number} <function-lsqrfact>` is similar to the bound predicted by the condition number.

```{code-cell}
observed_error = norm(FNC.lsqrfact(A, b) - x) / norm(x);
@show observed_error;
@show error_bound = cond(A) * eps();
```
