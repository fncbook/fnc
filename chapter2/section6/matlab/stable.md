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
[**Demo %s**](#demo-pivoting-stable)

We construct a linear system for this matrix with $\epsilon=10^{-12}$ and exact solution $[1, 1]$:

```{code-cell}
ep = 1e-12
A = [-ep 1; 1 -1];
b = A * [1; 1];
```

We can factor the matrix without pivoting and solve for $\mathbf{x}$.

```{code-cell}
[L, U] = lufact(A);
x = backsub( U, forwardsub(L, b) )
```

Note that we have obtained only about 5 accurate digits for $x_1$. We could make the result even more inaccurate by making $\epsilon$ even smaller:

```{code-cell}
ep = 1e-20; A = [-ep 1; 1 -1];
b = A * [1; 1];
[L, U] = lufact(A);
x = backsub( U, forwardsub(L, b) )
```

This effect is not due to ill conditioning of the problem—a solution with PLU factorization works perfectly:

```{code-cell}
A \ b
```
