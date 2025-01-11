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
[**Demo %s**](#demo-pwlin-hat)

Let's define a set of four nodes (i.e., $n=3$ in our formulas).

```{index} ! Julia; annotate!
```

```{code-cell}
t = array([0, 0.075, 0.25, 0.55, 0.7, 1])
```

We plot the hat functions $H_0,\ldots,H_3$.

```{code-cell}
x = linspace(0, 1, 300)
for k in range(6):
    plot(x, FNC.hatfun(t, k)(x))
xlabel("$x$"),  ylabel("$H_k(x)$")
title("Hat functions");
```
