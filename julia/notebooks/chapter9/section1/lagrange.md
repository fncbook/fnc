---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
Here is a vector of nodes.

```{code-cell}
t = [ 1, 1.5, 2, 2.25, 2.75, 3 ]
n = length(t) - 1;
```

Let's apply the definition of the cardinal Lagrange polynomial for $k=2$. First we define a polynomial $q$ that is zero at all the nodes except $i=k$. Then $\ell_2$ is found by normalizing $q$ by $q(t_k)$.
```{tip}
Character ℓ is typed as `\ell`<kbd>Tab</kbd>.
```

```{code-cell}
k = 2
q(x) = prod(x - t[i] for i in [0:k-1; k+1:n] .+ 1)
ℓₖ(x) = q(x) / q(t[k+1]);
```

A plot confirms the cardinal property of the result.

```{code-cell}
using Plots
plot(ℓₖ, 1, 3)
y = zeros(n+1);  y[k+1] = 1
scatter!(t, y, color=:black,
    xaxis=(L"x"),  yaxis=(L"\ell_2(x)"),
    title="Lagrange cardinal function")
```

Observe that $\ell_k$ is _not_ between zero and one everywhere, unlike a hat function.
