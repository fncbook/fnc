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
[**Demo %s**](#demo-lu-outertri)


```{index} MATLAB; tril, MATLAB; triu
```
We explore the outer product formula for two random triangular matrices.

```{code-cell}
L = tril( randi(9, 3, 3) )
```

```{code-cell}
U = triu( randi(9, 3, 3) )
```

Here are the three outer products in the sum in {eq}`matrixouter`:

```{code-cell}
L(:, 1) * U(1, :)
```

```{code-cell}
L(:, 2) * U(2, :)
```

```{code-cell}
L(:, 3) * U(3, :)
```

Simply because of the triangular zero structures, only the first outer product contributes to the first row and first column of the entire product. 
