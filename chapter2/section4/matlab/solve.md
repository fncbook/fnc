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
[**Demo %s**](#demo-lu-solve)

Here are the data for a linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$. 

```{code-cell}
A = [2 0 4 3; -4 5 -7 -10; 1 15 2 -4.5; -2 0 2 -13];
b = [4; 9; 9; 4];
```

We apply {numref}`Function {number} <function-lufact>` and then do two triangular solves.

```{code-cell}
[L, U] = lufact(A)
z = forwardsub(L, b);
x = backsub(U, z);
```

A check on the residual assures us that we found the solution.

```{code-cell}
b - A * x
```
