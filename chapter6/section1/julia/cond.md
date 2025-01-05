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
[**Demo %s**](#demo-basics-cond)

Consider the ODEs $u'=u$ and $u'=-u$. In each case we compute $\partial f/\partial u = \pm 1$, so the condition number bound from {numref}`Theorem %s <theorem-depIC>` is $e^{b-a}$ in both problems. However, they behave quite differently. In the case of exponential growth, $u'=u$, the bound is the actual condition number.

```{code-cell}
:tags: [remove-input]
t = range(0, 3, length=800)
u = @. exp(t) * 1
lower, upper = @. exp(t) * 0.7, @. exp(t) * 1.3
plot(t, u;
    l=:black, ribbon=(lower, upper),
    leg=:none,  xlabel=L"t",  ylabel=L"u(t)",
    title="Exponential divergence of solutions")
```

But with $u'=-u$, solutions actually get closer together with time.

```{code-cell}
:tags: [remove-input]
u = @. exp(-t) * 1
lower, upper = @. exp(-t) * 0.7, @. exp(-t) * 1.3
plot(t, u;
    l=:black,  ribbon=(lower, upper),
    leg=:none,  xlabel=L"t",  ylabel=L"u(t)",
    title="Exponential convergence of solutions")
```

In this case the actual condition number is one, because the initial difference between solutions is the largest over all time. Hence the exponentially growing bound $e^{b-a}$ is a gross overestimate.
