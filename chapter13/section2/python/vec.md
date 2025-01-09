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
[**Demo %s**](#demo-diffadv-vec)


```{code-cell}
m, n = 4, 3
x = linspace(0, 2, m+1)
y = linspace(-3, 0, n+1)

f = lambda x, y: cos(0.75 * pi * x * y - 0.5 * pi * y)
mtx, X, Y, vec, unvec, _ = FNC.tensorgrid(x, y)
F = mtx(f)
print(f"function on a {m}x{n} grid:")
with printoptions(precision=4, suppress=True):
    print(F)

print("vec(F):")
with printoptions(precision=4, suppress=True):
    print(vec(F))
```

The `unvec` operation is the inverse of vec.

```{code-cell}
print("unvec(vec(F)):")
with printoptions(precision=4, suppress=True):
    print(unvec(vec(F)))
```
