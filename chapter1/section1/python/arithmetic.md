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
[**Demo %s**](#demo-float-arithmetic)


There is no double precision number between $1$ and $1+\varepsilon_\text{mach}$. Thus, the following difference is zero despite its appearance.

```{code-cell} ipython3
eps = finfo(float).eps
e = eps/2
print((1.0 + e) - 1.0)
```

However, $1-\varepsilon_\text{mach}/2$ is a double precision number, so it and its negative are represented exactly:

```{code-cell} ipython3
print(1.0 + (e - 1.0))
```

This is now the "correct" result. But we have found a rather shocking breakdown of the associative law of addition!

