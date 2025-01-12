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
[**Demo %s**](#demo-structure-banded)

```{index} ! MATLAB; fill, MATLAB; diagm, ! MATLAB; diag
```

Here is a small tridiagonal matrix. Note that there is one fewer element on the super- and subdiagonals than on the main diagonal.

```{code-cell}
A = [ 2 -1  0  0  0  0
      4  2 -1  0  0  0
      0  3  0 -1  0  0
      0  0  2  2 -1  0
      0  0  0  1  1 -1
      0  0  0  0  0  2 ];
```

```{index} ! MATLAB; diag
```

We can extract the elements on any diagonal using the `diag` function. The main or central diagonal is numbered zero, above and to the right of that is positive, and below and to the left is negative.
```{tip}
:class: dropdown
The `diag` function extracts the elements from a specified diagonal of a matrix.
```

```{code-cell}
diag_main = diag(A, 0)'
diag_plusone = diag(A, 1)'
diag_minusone = diag(A,-1)'
```

% We can also put whatever numbers we like onto any diagonal with `diag`.
    
```{code-cell}
A = A + diag([5 8 6 7], 2)
```

The lower and upper bandwidths of $\mathbf{A}$ are repeated in the factors from the unpivoted LU factorization. 

```{code-cell}
[L, U] = lufact(A)
```
