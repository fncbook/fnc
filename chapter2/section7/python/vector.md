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
[**Demo %s**](#demo-norms-vector)

```{index} ! Python; norm
```

The `norm` function from `numpy.linalg` computes vector norms.

```{code-cell} 
from numpy.linalg import norm
x = array([2, -3, 1, -1])
print(norm(x))       # 2-norm by default
```

```{code-cell} 
print(norm(x, inf))
```

```{code-cell} 
print(norm(x, 1))
```
