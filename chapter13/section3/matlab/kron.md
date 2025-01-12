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
[**Demo %s**](#demo-laplace-kron)


```{code-cell}
A = [1, 2; -2, 0];
B = [1, 10, 100; -5, 5, 3];
disp("A:")
disp(A)
disp("B:")
disp(B)
```

Applying the definition manually, we get

```{code-cell}
A_kron_B = [
    A(1,1)*B  A(1,2)*B;
    A(2,1)*B  A(2,2)*B
    ]
```

```{index} ! MATLAB; kron
```

But it makes more sense to use `kron`.

```{code-cell}
kron(A, B)
```
