---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
exec(open("../../../python/FNC_init.py").read())
```
[**Demo %s**](#demo-laplace-kron)


```{code-cell}
A = array([[1, 2], [-2, 0]])
B = array([[1, 10, 100], [-5, 5, 3]])
print("A:")
print(A)
print("B:")
print(B)
```

Applying the definition manually, we get

```{code-cell}
A_kron_B = vstack([ hstack([A[0, 0] * B, A[0, 1] * B]), hstack([A[1, 0] * B, A[1, 1] * B]) ])
print(A_kron_B)
```

```{index} ! Python; kron
```

But it makes more sense to use `kron` from NumPy, or the `scipy.sparse` version when sparsity is to be preserved.

```{code-cell}
kron(A, B)
```
