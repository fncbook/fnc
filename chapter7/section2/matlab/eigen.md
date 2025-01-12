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
[**Demo %s**](#demo-evd-eigen)


```{index} ! MATLAB; eig
```

The `eig` function with one output argument returns a vector of the eigenvalues of a matrix.

```{code-cell}
A = pi * ones(2, 2);
lambda = eig(A)
```

With two output arguments given, `eig` returns a matrix eigenvectors and a diagonal matrix with the eigenvalues.

```{code-cell}
[V, D] = eig(A)
```

We can check the fact that this is an EVD.

```{code-cell}
norm( A - V*D/V )   % / V is like * inv(V)
```

If the matrix is not diagonalizable, no message is given, but `V` will be singular. The robust way to detect that circumstance is via $\kappa(\mathbf{V})$.

```{index} condition number; of a matrix
```

```{code-cell}
A = [-1 1; 0 -1];
[V, D] = eig(A)
```

```{code-cell}
cond(V)
```

Even in the nondiagonalizable case, $\mathbf{A}\mathbf{V} = \mathbf{V}\mathbf{D}$ holds.

```{code-cell}
norm(A * V - V * D)
```
