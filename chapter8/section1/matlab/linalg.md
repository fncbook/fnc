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
[**Demo %s**](#demo-structure-linalg)


The following generates a random sparse matrix with prescribed eigenvalues.

```{code-cell}
n = 4000;
density = 4e-4;
lambda = 1 ./ (1:n);
A = sprandsym(n, density, lambda);
clf,  spy(A)
title('Sparse symmetric matrix')  
```

```{index} ! MATLAB; eigs
```

The `eigs` function finds a small number eigenvalues meeting some criterion. First, we ask for the 5 of largest (complex) magnitude.

```{code-cell}
[V, D] = eigs(A, 5);    % largest magnitude
1 ./ diag(D)            % should be 1, 2, 3, 4, 5
```

Now we find the 4 closest to the value 0.03 in the complex plane.

```{code-cell}
[V, D] = eigs(A, 4, 0.03);    % closest to 0.03
diag(D)
```

```{index} MATLAB; \\
```

The time needed to solve a sparse linear system is not easy to predict unless you have some more information about the matrix. But it will typically be orders of magnitude faster than the dense version of the same problem.

```{code-cell}
x = 1 ./ (1:n)';  
b = A * x;
tic, sparse_err = norm(x - A\b), sparse_time = toc
```

```{code-cell}
F = full(A);
tic, dense_err = norm(x - F\b), dense_time = toc
```
