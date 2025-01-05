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
[**Demo %s**](#demo-diffadv-vec)


```{code-cell}
m = 2;
n = 3;
V = rand(1:9, m, n);
v = vec(V)
```

The `unvec` operation is the inverse of vec.

```{code-cell}
unvec = z -> reshape(z, m, n)
unvec(v)
```
