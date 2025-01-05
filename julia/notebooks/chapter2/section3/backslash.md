---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
For a square matrix $\mathbf{A}$, the syntax `A \ b` is mathematically equivalent to $\mathbf{A}^{-1} \mathbf{b}$. 

```{code-cell}
A = [1 0 -1; 2 2 1; -1 -3 0]
```

```{code-cell}
b = [1, 2, 3]
```

```{code-cell}
x = A \ b
```

```{index} residual
```

One way to check the answer is to compute a quantity known as the **residual**. It is (ideally) close to machine precision (relative to the elements in the data).

```{code-cell}
residual = b - A * x
```

If the matrix $\mathbf{A}$ is singular, you may get an error.

```{code-cell} julia
:tags: [raises-exception]
A = [0 1; 0 0]
b = [1, -1]
x = A \ b    # throws an error
```

```{index} ! Julia; rank
```

In this case we can check that the rank of $\mathbf{A}$ is less than its number of columns, indicating singularity.
```{tip}
The function `rank` computes the rank of a matrix. However, it is numerically unstable for matrices that are nearly singular, in a sense to be defined in a later section.
```

```{code-cell}
rank(A)
```

A linear system with a singular matrix might have no solution or infinitely many solutions, but in either case, backslash will fail. Moreover, detecting singularity is a lot like checking whether two floating-point numbers are *exactly* equal: because of roundoff, it could be missed. In {numref}`section-linsys-condition-number` we'll find a robust way to fully describe this situation.
