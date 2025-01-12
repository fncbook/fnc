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
[**Demo %s**](#demo-pivoting-usage)

The third output of `plufact` is the permutation vector we need to apply to $\mathbf{A}$.

```{code-cell}
A = randi(20, 4, 4);
[L, U, p] = plufact(A);
A(p, :) - L * U    % should be ≈ 0
```

Given a vector $\mathbf{b}$, we solve $\mathbf{A}\mathbf{x}=\mathbf{b}$ by first permuting the entries of $\mathbf{b}$ and then proceeding as before.

```{code-cell}
b = rand(4, 1);
z = forwardsub(L, b(p));
x = backsub(U, z)
```

A residual check is successful:

```{code-cell}
b - A*x
```
