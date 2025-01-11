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
[**Demo %s**](#demo-galerkin-fem)

Here are the coefficient function definitions. Even though $s$ is a constant, it has to be defined as a function for {numref}`Function {number} <function-fem>` to use it.

```{code-cell}
c = lambda x: x**2
q = lambda x: 4 * ones(len(x))
f = lambda x: sin(pi * x)
```

```{code-cell}
x, u = FNC.fem(c, q, f, 0, 1, 50)
plot(x, u)
xlabel("$x$"),  ylabel("$u$")
title("Solution by finite elements");
```
