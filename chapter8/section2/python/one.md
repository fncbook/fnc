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
[**Demo %s**](#demo-power-one)

Here we choose a random 5×5 matrix and a random 5-vector.

```{code-cell}
A = random.choice(range(10), (5, 5))
A = A / sum(A, 0)
x = random.randn(5)
print(x)
```

Applying matrix-vector multiplication once doesn't do anything recognizable.

```{code-cell}
y = A @ x
print(y)
```

Repeating the multiplication still doesn't do anything obvious.

```{code-cell}
z = A @ y
print(z)
```

But if we keep repeating the matrix-vector multiplication, something remarkable happens: $\mathbf{A} \mathbf{x} \approx \mathbf{x}$.

```{code-cell}
x = random.randn(5)
for j in range(6):
    x = A @ x
print(x)
print(A @ x)
```

This phenomenon is unlikely to be a coincidence!

