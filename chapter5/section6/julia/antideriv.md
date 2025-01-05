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
[**Demo %s**](#demo-int-antideriv)

The antiderivative of $e^x$ is, of course, itself. That makes evaluation of $\int_0^1 e^x\,dx$ by the Fundamental Theorem trivial.

```{code-cell}
exact = exp(1) - 1
```

```{index} ! Julia; quadgk
```

The Julia package `QuadGK` has an all-purpose numerical integrator that estimates the value without finding the antiderivative first. As you can see here, it's often just as accurate.

```{code-cell}
using QuadGK
Q, errest = quadgk(x -> exp(x), 0, 1)
@show Q;
```

The numerical approach is also far more robust. For example, $e^{\,\sin x}$ has no useful antiderivative. But numerically, it's no more difficult.

```{code-cell}
Q, errest = quadgk(x -> exp(sin(x)), 0, 1)
@show Q;
```

When you look at the graphs of these functions, what's remarkable is that one of these areas is basic calculus while the other is almost impenetrable analytically. From a numerical standpoint, they are practically the same problem.

```{code-cell}
plot([exp, x -> exp(sin(x))], 0, 1, fill=0, layout=(2, 1),
    xlabel=L"x", ylabel=[L"e^x" L"e^{\sin(x)}"], ylim=[0, 2.7])
```
