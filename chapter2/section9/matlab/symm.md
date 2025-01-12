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
[**Demo %s**](#demo-structure-symm)


We begin with a symmetric $\mathbf{A}$.

```{code-cell}
A_1 = [ 2     4     4     2
        4     5     8    -5
        4     8     6     2
        2    -5     2   -26 ];
```

We won't use pivoting, so the pivot element is at position (1,1). This will become the first element on the diagonal of $\mathbf{D}$. Then we divide by that pivot to get the first column of $\mathbf{L}$.

```{code-cell}
L = eye(4);
d = zeros(4, 1);
d(1) = A_1(1, 1);
L(:, 1) = A_1(:, 1) / d(1);
A_2 = A_1 - d(1) * L(:, 1) * L(:, 1)'
```

We are now set up the same way for the submatrix in rows and columns 2–4.

```{code-cell}
d(2) = A_2(2, 2);
L(:, 2) = A_2(:, 2) / d(2);
A_3 = A_2 - d(2) * L(:, 2) * L(:, 2)'
```

We continue working our way down the diagonal.

```{code-cell}
d(3) = A_3(3, 3);
L(:, 3) = A_3(:, 3) / d(3);
A_4 = A_3 - d(3) * L(:, 3) * L(:, 3)'
d(4) = A_4(4, 4);
d
L
```

We have arrived at the desired factorization, which we can validate:

```{code-cell}
norm(A_1 - (L * diag(d) * L'))
```
