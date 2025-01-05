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
[**Demo %s**](#demo-finitediffs-fd1)

If $f(x)=e^{\,\sin(x)}$, then $f'(0)=1$.

```{code-cell}
f = x -> exp(sin(x));
```

Here are the first two centered differences from {numref}`table-FDcenter`.

```{code-cell}
h = 0.05
CD2 = (-f(-h) + f(h)) / 2h
CD4 = (f(-2h) - 8f(-h) + 8f(h) - f(2h)) / 12h
@show (CD2, CD4);
```

Here are the first two forward differences from {numref}`table-FDforward`.

```{code-cell}
FD1 = (-f(0) + f(h)) / h
FD2 = (-3f(0) + 4f(h) - f(2h)) / 2h
@show (FD1, FD2);
```

Finally, here are the backward differences that come from reverse-negating the forward differences.

```{code-cell}
BD1 = (-f(-h) + f(0)) / h
BD2 = (f(-2h) - 4f(-h) + 3f(0)) / 2h
@show (BD1, BD2);
```
