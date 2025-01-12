---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
cd  /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init
```
[**Demo %s**](#demo-laplace-fd)

We make a crude discretization for illustrative purposes.

```{code-cell}
m = 6;  n = 5;
[x, Dx, Dxx] = diffmat2(m, [0, 3]);
[y, Dy, Dyy] = diffmat2(n, [-1, 1]);
[mtx, X, Y, vec, unvec, is_boundary] = tensorgrid(x, y);
```

Next, we define $\phi$ and evaluate it on the grid to get the forcing vector of the linear system.

```{code-cell}
phi = @(x, y) x.^2 - y + 2;
b = vec(mtx(phi));
```

Here are the coefficients for the PDE collocation, before any modifications are made for the boundary conditions. The combination of Kronecker products and finite differences produces a characteristic sparsity pattern.

```{code-cell}
A = kron(speye(n+1), sparse(Dxx)) + kron(sparse(Dyy), speye(m+1));
clf,  spy(A)
title("System matrix before boundary conditions")
```

The number of equations is equal to $(m+1)(n+1)$, which is the total number of points on the grid.

```{code-cell}
N = length(b)
```

We now use the Boolean array that indicates where the boundary points lie in the grid.

```{code-cell}
spy(is_boundary)
title("Boundary points")
```

In order to impose Dirichlet boundary conditions, we replace the boundary rows of the system by rows of the identity.
```{tip}
:class: dropdown
Changing rows of a sparse array requires that the operands be in a particular sparse representation called `lil`. The conversion isn't done automatically because it can be slow and you are encouraged to avoid it when possible. We're just trying to keep things conceptually simple here.
```

```{code-cell}
I = speye(N);
idx = vec(is_boundary);
A(idx, :) = I(idx, :);

spy(A)
title("System matrix with boundary conditions")
```

Finally, we must replace the rows in the vector $\mathbf{b}$ by the boundary values being assigned to the boundary points. Here, we let the boundary values be zero everywhere.

```{code-cell}
b(idx) = 0;
```

Now we can solve for $\mathbf{u}$ and reinterpret it as the matrix-shaped $\mathbf{U}$, the solution on our grid.

```{code-cell}
u = A \ b;
U = unvec(u)
```
