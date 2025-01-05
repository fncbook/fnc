---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{code-cell}
A = [ 2 0; 1 -1 ]
```

In Julia, one uses `norm` for vector norms and for the Frobenius norm of a matrix, which is like stacking the matrix into a single vector before taking the 2-norm. 

```{code-cell}
Fronorm = norm(A)
```

```{index} ! Julia; opnorm
```

Most of the time we want to use `opnorm`, which is an induced matrix norm. The default is the 2-norm.

```{code-cell}
twonorm = opnorm(A)
```

You can get the 1-norm as well.

```{code-cell}
onenorm = opnorm(A, 1)
```

```{index} ! Julia; maximum, ! Julia; minimum, ! Julia; sum
```

According to {eq}`mxonenorm`, the matrix 1-norm is equivalent to the maximum of the sums down the columns (in absolute value).
```{tip}
Use `sum` to sum along a dimension of a matrix. You can also sum over the entire matrix by omitting the `dims` argument.
The `maximum` and `minimum` functions also work along one dimension or over an entire matrix. To get both values at once, use `extrema`.
```

```{code-cell}
# Sum down the rows (1st matrix dimension):
maximum( sum(abs.(A), dims=1) )   
```

Similarly, we can get the $\infty$-norm and check our formula for it.

```{code-cell}
infnorm = opnorm(A, Inf)
```

```{code-cell}
 # Sum across columns (2nd matrix dimension):
maximum( sum(abs.(A), dims=2) )  
```
Next we illustrate a geometric interpretation of the 2-norm. First, we will sample a lot of vectors on the unit circle in $\mathbb{R}^2$.
```{tip}
You can use functions as values, e.g., as elements of a vector. 
```

```{code-cell}
# Construct 601 unit column vectors.
θ = 2π * (0:1/600:1)   # type \theta then Tab
x = [ fun(t) for fun in [cos, sin], t in θ ];
```

```{index} ! Julia; subplots
```

To create an array of plots, start with a `plot` that has a `layout` argument, then do subsequent `plot!` calls with a `subplot` argument.

```{code-cell}
plot(aspect_ratio=1, layout=(1, 2),
    xlabel=L"x_1",  ylabel=L"x_2")
plot!(x[1, :], x[2, :], subplot=1, title="Unit circle") 
```

The linear function $\mathbf{f}(\mathbf{x}) = \mathbf{A}\mathbf{x}$ defines a mapping from $\mathbb{R}^2$ to $\mathbb{R}^2$. We can apply `A` to every column of `x` by using a single matrix multiplication.

```{code-cell}
Ax = A * x;
```

The image of the transformed vectors is an ellipse. 

```{code-cell}
plot!(Ax[1, :], Ax[2, :], 
    subplot=2, title="Image under x → Ax")
```

That ellipse just touches the circle of radius $\|\mathbf{A}\|_2$.

```{code-cell}
plot!(twonorm*x[1, :], twonorm*x[2, :], subplot=2, l=:dash)
```
