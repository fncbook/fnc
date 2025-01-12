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
[**Demo %s**](#demo-pivoting-fix)

Here is the trouble-making matrix from {numref}`Demo {number} <demo-pivoting-fail>`.

```{code-cell}
A_1 = [2 0 4 3; -2 0 2 -13; 1 15 2 -4.5; -4 5 -7 -10]
```

```{index} ! MATLAB; max, ! MATLAB; \~
```

We now find the largest candidate pivot in the first column. We don't care about sign, so we take absolute values before finding the max.
```{tip}
:class: dropdown
The second output of `max` returns the location of the largest element of a vector. The `~` symbol is used to ignore the value of the first output.
```


```{code-cell}
[~, i] = max( abs(A_1(:, 1)) ) 
```

This is the row of the matrix that we extract to put into $\mathbf{U}$. That guarantees that the division used to find $\boldsymbol{\ell}_1$ will be valid.

```{code-cell}
L = zeros(4, 4);
U = zeros(4, 4);
U(1, :) = A_1(i, :);
L(:, 1) = A_1(:, 1) / U(1, 1);
A_2 = A_1 - L(:, 1) * U(1, :)
```

Observe that $\mathbf{A}_2$ has a new zero row and zero column, but the zero row is the fourth rather than the first. However, we forge on by using the largest possible pivot in column 2 for the next outer product.

```{code-cell}
[~, i] = max( abs(A_2(:, 2)) )
U(2, :) = A_2(i, :);
L(:, 2) = A_2(:, 2) / U(2, 2);
A_3 = A_2 - L(:, 2) * U(2, :)
```

Now we have zeroed out the third row as well as the second column. We can finish out the procedure.

```{code-cell}
[~, i] = max( abs(A_3(:, 3)) ) 
U(3, :) = A_3(i, :);
L(:, 3) = A_3(:, 3) / U(3, 3);
A_4 = A_3 - L(:, 3) * U(3, :)
```

```{code-cell}
[~, i] = max( abs(A_4(:, 4)) ) 
U(4, :) = A_4(i, :);
L(:, 4) = A_4(:, 4) / U(4, 4);
```

We do have a factorization of the original matrix:

```{code-cell}
A_1 - L * U
```

And $\mathbf{U}$ has the required structure:

```{code-cell}
U
```

However, the triangularity of $\mathbf{L}$ has been broken.

```{code-cell}
L
```
