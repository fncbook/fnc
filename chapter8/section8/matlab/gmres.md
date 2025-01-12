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
[**Demo %s**](#demo-precond-gmres)

Here is a random nonsymmetric matrix.

```{code-cell}
n = 8000;
A = speye(n) + sprand(n, n, 0.00035);
```

Without a preconditioner, restarted GMRES makes slow progress.

```{code-cell}
b = rand(n, 1);
[x, ~, ~, ~, resid_plain] = gmres(A, b, 50, 1e-10, 3);  % restart at 50
format short e
resid_plain(1:30:end)
```

```{index} ! MATLAB; ilu
```

This version of incomplete LU factorization simply prohibits fill-in for the factors, freezing the sparsity pattern of the approximate factors to match the original matrix.

```{code-cell}
[L, U] = ilu(A);
clf
subplot(121), spy(L)
title('L')
subplot(122), spy(U)
title('U')
disp(sprintf("There are %d nonzeros in A", nnz(A)))
```

It does _not_ produce a true factorization of $\mathbf{A}$.

```{code-cell}
norm( full(A - L * U) )  
```

The actual preconditioning matrix is $\mathbf{M}=\mathbf{L}\mathbf{U}$. However, the `gmres` function allows setting the preconditioner by giving the factors independently.

```{code-cell}
[x, ~, ~, ~, resid_prec] = gmres(A, b, [], 1e-10, 300, L, U);
```

The preconditioner makes a significant difference in the number of iterations needed.

```{code-cell}
clf, semilogy(resid_plain)
hold on, semilogy(resid_prec)
xlabel('iteration number'), ylabel('residual norm')
title('Precondtioned GMRES ')
legend('no preconditioner', 'with preconditioner');
```
```
