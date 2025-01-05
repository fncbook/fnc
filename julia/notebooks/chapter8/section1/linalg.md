---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
The following generates a random sparse matrix with prescribed eigenvalues.

```{code-cell}
n = 4000
density = 4e-4
λ = @. 1 + 1 / (1:n)   # exact eigenvalues
A = FNC.sprandsym(n, density, λ);
```

```{index} ! Julia; eigs
```

The `eigs` function from `Arpack` finds a small number eigenvalues meeting some criterion. First, we ask for the 5 of largest (complex) magnitude using `which=:LM`.

```{code-cell}
using Arpack
λmax, V = eigs(A, nev=5, which=:LM)    # Largest Magnitude
fmt = ft_printf("%20.15f")
pretty_table([λmax λ[1:5]], header=["found", "exact"], formatters=fmt)
```

Now we find the 5 closest to the value 1 in the complex plane, via `sigma=1`.

```{code-cell}
λ1, V = eigs(A, nev=5, sigma=1)    # closest to sigma
data = [λ1 λ[end:-1:end-4]]
pretty_table(data, header=["found", "exact"], formatters=fmt)
```

```{index} Julia; \\
```

The time needed to solve a sparse linear system is not easy to predict unless you have some more information about the matrix. But it will typically be orders of magnitude faster than the dense version of the same problem.

```{code-cell}
x = @. 1 / (1:n);
b = A * x;
```

```{code-cell}
norm(x - A \ b);  # force compilation
t = @elapsed sparse_err = norm(x - A \ b)
println("Time for sparse solve: $t")
```

```{code-cell}
D = Matrix(A)  # convert to regular matrix
norm(x - D \ b);
t = @elapsed dense_err = norm(x - D \ b)
println("Time for dense solve: $t")
```

```{code-cell}
@show sparse_err;
@show dense_err;
```
