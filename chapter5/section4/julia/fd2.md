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
[**Demo %s**](#demo-finitediffs-fd2)

If $f(x)=e^{\,\sin(x)}$, then $f''(0)=1$.

```{code-cell}
f = x -> exp(sin(x));
```

Here is a centered estimate given by {eq}`centerFD22`.

```{code-cell}
h = 0.05
CD2 = (f(-h) - 2f(0) + f(h)) / h^2
@show CD2;
```

For the same $h$, here are forward estimates given by {eq}`forwardFD21` and {eq}`forwardFD22`.

```{code-cell}
FD1 = (f(0) - 2f(h) + f(2h)) / h^2
FD2 = (2f(0) - 5f(h) + 4f(2h) - f(3h)) / h^2
@show (FD1, FD2);
```

Finally, here are the backward estimates that come from reversing {eq}`forwardFD21` and {eq}`forwardFD22`.

```{code-cell}
BD1 = (f(-2h) - 2f(-h) + f(0)) / h^2
BD2 = (-f(-3h) + 4f(-2h) - 5f(-h) + 2f(0)) / h^2
@show (BD1, BD2);
```
