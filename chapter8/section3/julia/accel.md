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
[**Demo %s**](#demo-inviter-accel)

```{code-cell}
λ = [1, -0.75, 0.6, -0.4, 0]
# Make a triangular matrix with eigenvalues on the diagonal.
A = triu(ones(5, 5), 1) + diagm(λ)
```

We begin with a shift $s=0.7$, which is closest to the eigenvalue 0.6.

```{code-cell}
s = 0.7
x = ones(5)
y = (A - s * I) \ x
β = x[1] / y[1] + s
```

Note that the result is not yet any closer to the targeted 0.6. But we proceed (without being too picky about normalization here).

```{code-cell}
s = β
x = y / y[1]
y = (A - s * I) \ x
β = x[1] / y[1] + s
```

Still not much apparent progress. However, in just a few more iterations the results are dramatically better.

```{code-cell}
for k in 1:4
    s = β
    x = y / y[1]
    y = (A - s * I) \ x
    @show β = x[1] / y[1] + s
end
```
