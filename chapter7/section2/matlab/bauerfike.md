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
[**Demo %s**](#demo-evd-bauerfike)


```{index} MATLAB; adjoint, MATLAB; \'
```

We first define a hermitian matrix. Note that the `'` operation is the adjoint and includes complex conjugation.

```{code-cell}
n = 7;
A = randn(n, n) + 1i * randn(n, n);
A = (A + A') / 2;
```

We confirm that the matrix $\mathbf{A}$ is normal by checking that $\kappa(\mathbf{V}) = 1$ (to within roundoff).

```{code-cell}
[V, D] = eig(A);
lambda = diag(D);
cond(V)
```

Now we perturb $\mathbf{A}$ and measure the effect on the eigenvalues. The Bauer–Fike theorem uses absolute differences, not relative ones. Note: since the ordering of eigenvalues can change, we look at all pairwise differences and take the minima.

```{code-cell}
E = randn(n, n) + 1i * randn(n, n);
E = 1e-8 * E / norm(E);
dd = eig(A + E);
dist = [];
for j = 1:n
    dist = [dist; min(abs(dd - lambda(j)))];
end
dist
```

As promised, the perturbations in the eigenvalues do not exceed the normwise perturbation to the original matrix.

Now we see what happens for a triangular matrix.

```{code-cell}
n = 20;
x = (1:n)';
A = triu(x * ones(1, n));
A(1:5, 1:5)
```

This matrix is not at all close to normal.

```{code-cell}
[V, D] = eig(A);
lambda = diag(D);
cond(V)
```

As a result, the eigenvalues can change by a good deal more.

```{code-cell}
E = randn(n, n) + 1i * randn(n, n);
E = 1e-8 * E / norm(E);
dd = eig(A + E);
dist = -Inf;
for j = 1:n
    dist = max(dist, min(abs(dd - lambda(j))));
end
fprintf("max change in eigenvalues: %.2e", dist)
fprintf("Bauer-Fike upper bound: %.2e", cond(V) * norm(E))
```

If we plot the eigenvalues of many perturbations, we get a cloud of points that roughly represents all the possible eigenvalues when representing this matrix with single-precision accuracy.

```{code-cell}
clf
plot(lambda, 0*lambda, 'o')
axis equal; hold on
for k = 1:60
    E = randn(n, n) + 1i * randn(n, n);
    E = eps(single(1)) * E / norm(E);
    dd = eig(A + E);
    plot(real(dd), imag(dd), 'k.', markersize=2)
end
```

The plot shows that some eigenvalues are much more affected than others. This situation is not unusual, but it is not explained by the Bauer–Fike theorem.
