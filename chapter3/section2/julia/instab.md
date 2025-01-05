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
[**Demo %s**](#demo-normaleqns-instab)


Because the functions $\sin^2(t)$, $\cos^2(t)$, and $1$ are linearly dependent, we should find that the following matrix is somewhat ill-conditioned.
```{tip}
:class: dropdown
The local variable scoping rule for loops applies to comprehensions as well.
```

```{code-cell}
t = range(0, 3, 400)
f = [ x -> sin(x)^2, x -> cos((1 + 1e-7) * x)^2, x -> 1. ]
A = [ f(t) for t in t, f in f ]
@show κ = cond(A);
```

Now we set up an artificial linear least-squares problem with a known exact solution that actually makes the residual zero.

```{code-cell}
x = [1., 2, 1]
b = A * x;
```

Using backslash to find the least-squares solution, we get a relative error that is well below $\kappa$ times machine epsilon.

```{code-cell}
x_BS = A \ b
@show observed_error = norm(x_BS - x) / norm(x);
@show error_bound = κ * eps();
```

If we formulate and solve via the normal equations, we get a much larger relative error. With $\kappa^2\approx 10^{14}$, we may not be left with more than about 2 accurate digits.

```{code-cell}
N = A' * A
x_NE = N \ (A'*b)
@show observed_err = norm(x_NE - x) / norm(x);
@show digits = -log10(observed_err);
```
