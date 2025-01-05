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
[**Demo %s**](#demo-polynomial-error)

Consider the problem of interpolating $\log(x)$ at these nodes:

```{code-cell}
t =  [ 1, 1.6, 1.9, 2.7, 3 ]
n = length(t) - 1;
```

Here $n=4$ and $f^{(5)}(\xi) = 4!/\xi^5$. For $\xi\in[1,3]$ we can say that $|f^{(5)}(\xi)| \le 4!$. Hence 

$$ |f(x)-p(x)| \le \frac{1}{5} \left| \Phi(x) \right|.$$
```{tip}
:class: dropdown
Character Φ is typed as `\Phi`<kbd>Tab</kbd>. (Note the capitalization.)
```

```{code-cell}
using Polynomials
Φ(x) = prod(x - tᵢ for tᵢ in t)
plot(x -> 0.2 * abs(Φ(x)), 1, 3, label=L"\frac{1}{5}|\Phi(t)|")
p = Polynomials.fit(t, log.(t))
plot!(t -> abs(log(t) - p(t)), 1, 3, label=L"|f(x)-p(x)|")
scatter!(t, zeros(size(t)), color=:black,
    xaxis=(L"x"), title="Interpolation error and upper bound")
```

The error is zero at the nodes, by the definition of interpolation. The error bound, as well as the error itself, has one local maximum between each consecutive pair of nodes.
