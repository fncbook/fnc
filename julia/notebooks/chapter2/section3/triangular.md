---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{index} ! Julia; tril, ! Julia; triu
```

It's easy to get just the lower triangular part of any matrix using the `tril` function.
```{tip}
Use `tril` to return a matrix that zeros out everything above the main diagonal. The `triu` function zeros out below the diagonal.
```

```{code-cell}
A = rand(1.:9., 5, 5)
L = tril(A)
```

We'll set up and solve a linear system with this matrix.

```{code-cell}
b = ones(5)
x = FNC.forwardsub(L,b)
```

It's not clear how accurate this answer is. However, the residual should be zero or comparable to $\macheps$.

```{code-cell}
b - L * x
```

```{index} ! Julia; Pair, Julia; diagm
```

Next we'll engineer a problem to which we know the exact answer. Use `\alpha` <kbd>Tab</kbd> and `\beta` <kbd>Tab</kbd> to get the Greek letters.
```{tip}
The notation `0=>ones(5)` creates a `Pair`. In `diagm`, pairs indicate the position of a diagonal and the elements that are to be placed on it.
```

```{code-cell}
α = 0.3;
β = 2.2;
U = diagm( 0=>ones(5), 1=>[-1, -1, -1, -1] )
U[1, [4, 5]] = [ α - β, β ]
U
```

```{code-cell}
x_exact = ones(5)
b = [α, 0, 0, 0, 1]
```

Now we use backward substitution to solve for $\mathbf{x}$, and compare to the exact solution we know already.

```{code-cell}
x = FNC.backsub(U,b)
err = x - x_exact
```

Everything seems OK here. But another example, with a different value for $\beta$, is more troubling.

```{code-cell}
α = 0.3;
β = 1e12;
U = diagm( 0=>ones(5), 1=>[-1, -1, -1, -1] )
U[1, [4, 5]] = [ α - β, β ]
b = [α, 0, 0, 0, 1]

x = FNC.backsub(U,b)
err = x - x_exact
```

It's not so good to get 4 digits of accuracy after starting with 16! The source of the error is not hard to track down. Solving for $x_1$ performs $(\alpha-\beta)+\beta$ in the first row. Since $|\alpha|$ is so much smaller than $|\beta|$, this a recipe for losing digits to subtractive cancellation.
