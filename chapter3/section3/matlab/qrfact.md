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
[**Demo %s**](#demo-qr-qrfact)


MATLAB provides access to both the thin and full forms of the QR factorization.

```{code-cell}
A = magic(5);
A = A(:, 1:4);
[m, n] = size(A)
```

Here is the full form:

```{code-cell}
[Q, R] = qr(A);
szQ = size(Q), szR = size(R)
```

We can test that $\mathbf{Q}$ is an orthogonal matrix:

```{code-cell}
QTQ = Q' * Q
norm(QTQ - eye(m))
```

With a second input argument given to `qr`, the thin form is returned. (This is usually the one we want in practice.)

```{code-cell}
[Q_hat, R_hat] = qr(A, 0);
szQ_hat = size(Q_hat), szR_hat = size(R_hat)
```

Now $\hat{\mathbf{Q}}$ cannot be an orthogonal matrix, because it is not square, but it is still ONC. Mathematically, $\hat{\mathbf{Q}}^T \hat{\mathbf{Q}}$ is a $4\times 4$ identity matrix.

```{code-cell}
Q_hat' * Q_hat - eye(n)
```
