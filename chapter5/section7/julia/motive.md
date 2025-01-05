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
[**Demo %s**](#demo-adapt-motive)

This function gets increasingly oscillatory as $x$ increases.

```{code-cell}
f = x -> (x + 1)^2 * cos((2x + 1) / (x - 4.3))
plot(f, 0, 4, xlabel=L"x", ylabel=L"f(x)")
```

Accordingly, the trapezoid rule is more accurate on the left half of this interval than on the right half.

```{code-cell}
using QuadGK
left_val, _ = quadgk(f, 0, 2, atol=1e-14, rtol=1e-14)
right_val, _ = quadgk(f, 2, 4, atol=1e-14, rtol=1e-14)

n = [50 * 2^k for k in 0:3]
err = zeros(length(n), 2)
for (k, n) in enumerate(n)
    T, _ = FNC.trapezoid(f, 0, 2, n)
    err[k, 1] = T - left_val

    T, _ = FNC.trapezoid(f, 2, 4, n)
    err[k, 2] = T - right_val
end

@pt :header=["n", "left error", "right error"] [n err]
```

Both the picture and the numerical results suggest that more nodes should be used on the right half of the interval than on the left half.
