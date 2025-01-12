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
[**Demo %s**](#demo-lu-derive)

For illustration, we work on a $4 \times 4$ matrix. We name it with a subscript in preparation for what comes.

```{code-cell}
A_1 = [
     2    0    4     3 
    -4    5   -7   -10 
     1   15    2   -4.5
    -2    0    2   -13
    ];
L = eye(4);
U = zeros(4, 4);
```

Now we appeal to {eq}`outer-row1`. Since $L_{11}=1$, we see that the first row of $\mathbf{U}$ is just the first row of $\mathbf{A}_1$.

```{code-cell}
U(1, :) = A_1(1, :)
```

From {eq}`outer-col1`, we see that we can find the first column of $\mathbf{L}$ from the first column of $\mathbf{A}_1$. 

```{code-cell}
L(:, 1) = A_1(:, 1) / U(1, 1)
```

 We have obtained the first term in the sum {eq}`matrixouter` for $\mathbf{L}\mathbf{U}$, and we subtract it away from $\mathbf{A}_1$.

```{code-cell}
A_2 = A_1 - L(:, 1) * U(1, :)
```

Now $\mathbf{A}_2 = \boldsymbol{\ell}_2\mathbf{u}_2^T + \boldsymbol{\ell}_3\mathbf{u}_3^T + \boldsymbol{\ell}_4\mathbf{u}_4^T.$ If we ignore the first row and first column of the matrices in this equation, then in what remains we are in the same situation as at the start. Specifically, only $\boldsymbol{\ell}_2\mathbf{u}_2^T$ has any effect on the second row and column, so we can deduce them now.

```{code-cell}
U(2, :) = A_2(2, :)
L(:, 2) = A_2(:, 2) / U(2, 2)
```

If we subtract off the latest outer product, we have a matrix that is zero in the first *two* rows and columns. 

```{code-cell}
A_3 = A_2 - L(:, 2) * U(2, :)
```

Now we can deal with the lower right $2\times 2$ submatrix of the remainder in a similar fashion.

```{code-cell}
U(3, :) = A_3(3, :);
L(:, 3) = A_3(:, 3) / U(3, 3);
A_4 = A_3 - L(:, 3) * U(3, :)
```

Finally, we pick up the last unknown in the factors.

```{code-cell}
U(4, 4) = A_4(4, 4);
```

We now have all of $\mathbf{L}$,

```{code-cell} 
L
```

and all of $\mathbf{U}$,

```{code-cell}
U
```

We can verify that we have a correct factorization of the original matrix by computing the backward error:

```{code-cell} 
A_1 - L * U
```

In floating point, we cannot expect the difference to be exactly zero as we found in this toy example. Instead, we would be satisfied to see that each element of the difference above is comparable in size to machine precision.

