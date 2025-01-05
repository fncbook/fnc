---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
We make a crude discretization for illustrative purposes.

```{code-cell}
m, n = 6, 5
x, Dx, Dxx = FNC.diffmat2(m, [0, 3])
y, Dy, Dyy = FNC.diffmat2(n, [-1, 1])
mtx, X, Y, unvec, is_boundary = FNC.tensorgrid(x, y)
```

Next, we evaluate $\phi$ on the grid to get the forcing vector of the linear system.

```{code-cell}
ϕ = (x, y) -> x^2 - y + 2
b = vec(mtx(ϕ));
```

Here are the coefficients for the PDE collocation, before any modifications are made for the boundary conditions. The combination of Kronecker products and finite differences produces a characteristic sparsity pattern.

```{code-cell}
using SparseArrays, Plots
A = kron(I(n+1), sparse(Dxx)) + kron(sparse(Dyy), I(m+1))
spy(A;
    color=:blues,  m=3,
    title="System matrix before boundary conditions")
```

The number of equations is equal to $(m+1)(n+1)$, which is the total number of points on the grid.

```{code-cell}
N = length(b)
```

The combination of Kronecker products and finite differences produces a characteristic sparsity pattern.

```{code-cell}

```

We now use the Boolean array that indicates where the boundary points lie in the grid.

```{code-cell}
spy(sparse(is_boundary);
    m=3,  color=:darkblue, 
    legend=:none,  title="Boundary points",
    xaxis=("column index", [0, n+2]), 
    yaxis=("row index", [0, m+2]))
```

In order to impose Dirichlet boundary conditions, we replace the boundary rows of the system by rows of the identity.

```{code-cell}
I_N = I(N)
idx = vec(is_boundary)
A[idx, :] .= I_N[idx, :];     # Dirichlet conditions
```

```{code-cell}
:tags: [hide-input]
spy(A;
    color=:blues,  m=3,
    title="System matrix with boundary conditions")
```

Finally, we must replace the rows in the vector $\mathbf{b}$ by the boundary values being assigned to the boundary points. Here, we let the boundary values be zero everywhere.

```{code-cell}
b[idx] .= 0;                 # Dirichlet values
```

Now we can solve for $\mathbf{u}$ and reinterpret it as the matrix-shaped $\mathbf{U}$, the solution on our grid.

```{code-cell}
u = A \ b
U = unvec(u)
```
