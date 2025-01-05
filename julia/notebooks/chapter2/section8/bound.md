---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---

```{index} ! Julia; cond
```
Julia has a function `cond` to compute matrix condition numbers. By default, the 2-norm is used. As an example, the family of *Hilbert matrices* is famously badly conditioned. Here is the $6\times 6$  case.
```{tip}
Type `\kappa` followed by <kbd>Tab</kbd> to get the Greek letter $\kappa$.
```

```{code-cell}
A = [ 1 / (i + j) for i in 1:6, j in 1:6 ]
κ = cond(A)
```

Because $\kappa\approx 10^8$, it's possible to lose nearly 8 digits of accuracy in the process of passing from $\mathbf{A}$ and $\mathbf{b}$ to $\mathbf{x}$. That fact is independent of the algorithm; it's inevitable once the data are expressed in finite precision. 

Let's engineer a linear system problem to observe the effect of a perturbation. We will make sure we know the exact answer.

```{code-cell}
x = 1:6
b = A * x
```

Now we perturb the system matrix and vector randomly by $10^{-10}$ in norm.

```{code-cell} 
# type \Delta then Tab to get Δ
ΔA = randn(size(A));  ΔA = 1e-10 * (ΔA / opnorm(ΔA));
Δb = randn(size(b));  Δb = 1e-10 * normalize(Δb);
```

We solve the perturbed problem using pivoted LU and see how the solution was changed.

```{code-cell}
new_x = ((A + ΔA) \ (b + Δb))
Δx = new_x - x
```

Here is the relative error in the solution.

```{code-cell}
@show relative_error = norm(Δx) / norm(x);
```

And here are upper bounds predicted using the condition number of the original matrix.

```{code-cell}
println("Upper bound due to b: $(κ * norm(Δb) / norm(b))")
println("Upper bound due to A: $(κ * opnorm(ΔA) / opnorm(A))")
```

Even if we didn't make any manual perturbations to the data, machine roundoff does so at the relative level of $\macheps$.

```{code-cell}
Δx = A\b - x
@show relative_error = norm(Δx) / norm(x);
@show rounding_bound = κ * eps();
```

Larger Hilbert matrices are even more poorly conditioned:

```{code-cell}
A = [ 1 / (i + j) for i=1:14, j=1:14 ];
κ = cond(A)
```

Note that $\kappa$ exceeds $1/\macheps$. In principle we therefore may end up with an answer that has relative error greater than 100%.

```{code-cell}
rounding_bound = κ*eps()
```

Let's put that prediction to the test.

```{code-cell}
x = 1:14
b = A * x  
Δx = A\b - x
@show relative_error = norm(Δx) / norm(x);
```

As anticipated, the solution has zero accurate digits in the 2-norm.
