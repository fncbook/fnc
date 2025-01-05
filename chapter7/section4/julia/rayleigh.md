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
[**Demo %s**](#demo-symm-eig-rayleigh)

We will use a symmetric matrix with a known EVD and eigenvalues equal to the integers from 1 to 20.

```{code-cell}
n = 20;
λ = 1:n
D = diagm(λ)
V, _ = qr(randn(n, n))   # get a random orthogonal V
A = V * D * V';
```

The Rayleigh quotient is a scalar-valued function of a vector.

```{code-cell}
R = x -> (x' * A * x) / (x' * x);
```

The Rayleigh quotient evaluated at an eigenvector gives the corresponding eigenvalue.

```{code-cell}
R(V[:, 7])
```

If the input to he Rayleigh quotient is within a small $\delta$ of an eigenvector, its output is within $O(\delta^2)$ of the corresponding eigenvalue. In this experiment, we observe that each additional digit of accuracy in an approximate eigenvector gives two more digits to the eigenvalue estimate coming from the Rayleigh quotient.

```{code-cell}
δ = @. 1 ./ 10^(1:5)
eval_diff = zeros(size(δ))
for (k, delta) in enumerate(δ)
    e = randn(n)
    e = delta * e / norm(e)
    x = V[:, 7] + e
    eval_diff[k] = R(x) - 7
end
labels = ["perturbation δ", "δ²", "R(x) - λ"]
@pt :header=labels [δ δ .^ 2 eval_diff]
```
