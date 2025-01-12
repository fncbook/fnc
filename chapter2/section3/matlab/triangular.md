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
[**Demo %s**](#demo-systems-triangular)

```{index} ! MATLAB; tril, ! MATLAB; triu
```

It's easy to get just the lower triangular part of any matrix using the `tril` function.
```{tip}
:class: dropdown
Use `tril` to return a matrix that zeros out everything above the main diagonal. The `triu` function zeros out below the diagonal.
```

```{code-cell}
A = randi(9, 5, 5);
L = tril(A)
```

We'll set up and solve a linear system with this matrix.

```{code-cell}
b = ones(5);
x = forwardsub(L, b)
```

```{index} residual
```

It's not clear how accurate this answer is. However, the residual should be zero or comparable to $\macheps$.

```{code-cell}
b - L * x
```

```{index} ! MATLAB; diag, ! MATLAB; eye
```

Next, we'll engineer a problem to which we know the exact answer. 
```{tip}
:class: dropdown
The `eye` function creates an identity matrix. The `diag` function uses 0 as the main diagonal, positive integers as superdiagonals, and negative integers as subdiagonals.
```

```{code-cell}
alpha = 0.3;
beta = 2.2;
U = eye(5) + diag([-1 -1 -1 -1], 1);
U(1, [4, 5]) = [alpha - beta, beta]
```

```{code-cell}
x_exact = ones(5);
b = [alpha; 0; 0; 0; 1];
```

Now we use backward substitution to solve for $\mathbf{x}$, and compare to the exact solution we know already.

```{code-cell}
x = backsub(U, b);
err = x - x_exact
```

Everything seems OK here. But another example, with a different value for $\beta$, is more troubling.

```{code-cell}
alpha = 0.3;
beta = 1e12;
U = eye(5) + diag([-1 -1 -1 -1], 1);
U(1, [4, 5]) = [alpha - beta, beta];
b = [alpha; 0; 0; 0; 1];

x = backsub(U, b);
err = x - x_exact
```

It's not so good to get 4 digits of accuracy after starting with sixteen! The source of the error is not hard to track down. Solving for $x_1$ performs $(\alpha-\beta)+\beta$ in the first row. Since $|\alpha|$ is so much smaller than $|\beta|$, this a recipe for losing digits to subtractive cancellation.
