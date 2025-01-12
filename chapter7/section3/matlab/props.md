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
[**Demo %s**](#demo-svd-props)

We verify some of the fundamental SVD properties using the built-in `svd` function.

```{index} ! MATLAB; svd
```

```{code-cell}
A = vander(1:5);
A = A(:, 1:4)
```

```{code-cell}
[U, S, V] = svd(A);
disp(sprintf("U is %d by %d. S is %d by %d. V is %d by %d.\n", size(U), size(S), size(V)))
```

We verify the orthogonality of the singular vectors as follows:

```{code-cell}
norm(U' * U - eye(5))
norm(V' * V - eye(4))
```

Here is verification of the connections between the singular values, norm, and condition number.

```{code-cell}
s = diag(S);
norm_A = norm(A)
sigma_max = s(1)
```

```{code-cell}
cond_A = cond(A)
sigma_ratio = s(1) / s(end)
```
