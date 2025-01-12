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
[**Demo %s**](#demo-evd-francisqr)

Let's start with a known set of eigenvalues and an orthogonal eigenvector basis.

```{code-cell}
D = diag([-6, -1, 2, 4, 5]);
[V, R]= qr(randn(5, 5));    % V is unitary
A = V * D * V';
```

```{code-cell}
sort(eig(A))
```

Now we will take the QR factorization and just reverse the factors.

```{code-cell}
[Q, R] = qr(A);
A = R * Q;
```

It turns out that this is a similarity transformation, so the eigenvalues are unchanged.

```{code-cell}
sort(eig(A))
```

What's remarkable, and not elementary, is that if we repeat this transformation many times, the resulting matrix converges to $\mathbf{D}$.

```{code-cell}
for k = 1:40
    [Q, R] = qr(A);
    A = R * Q;
end
format short e
A
```
