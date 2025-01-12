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
[**Demo %s**](#demo-structure-cholesky)

A randomly chosen matrix is extremely unlikely to be symmetric. However, there is a simple way to symmetrize one.

```{code-cell}
A = magic(4) + eye(4);
B = A + A'
```

```{index} ! MATLAB; chol
```

Similarly, a random symmetric matrix is unlikely to be positive definite. The Cholesky algorithm always detects a non-PD matrix by quitting with an error.
```{tip}
:class: dropdown
The `chol` function computes a Cholesky factorization if possible, or throws an error for a non-positive-definite matrix. 
```

```{warning} 
The `chol` function does *not* check for symmetry. It may give a nonsensical result if the input is not symmetric.
```

```{code-cell}
:tags: [raises-exception]
chol(B)    % throws an error
```

It's not hard to manufacture an SPD matrix to try out the Cholesky factorization.

```{code-cell}
B = A' * A;
R = chol(B)
```

Here we validate the factorization:

```{code-cell}
norm(R' * R - B) / norm(B)
```
