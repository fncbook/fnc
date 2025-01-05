---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---

```{code-cell}
A = [1 2; -2 0]
```

```{code-cell}
B = [1 10 100; -5 5 3]
```

Applying the definition manually, we get

```{code-cell}
A_kron_B = [
    A[1, 1]*B A[1, 2]*B;
    A[2, 1]*B A[2, 2]*B
    ]
```

```{index} ! Julia; kron
```

That result should be the same as the following.

```{code-cell}
kron(A, B)
```
