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
[**Demo %s**](#demo-nonlinear-mems)

Here is the problem definition. We use a truncated domain to avoid division by zero at $r=0$.

```{code-cell}
lamb = 0.5
phi = lambda r, w, dwdr: lamb / w**2 - dwdr / r
a, b = finfo(float).eps, 1
ga = lambda w, dw: dw
gb = lambda w, dw: w - 1
```

First we try a constant function as the initialization.

```{code-cell}
init = ones(201)
r, w1 = FNC.bvp(phi, [a, b], ga, gb, init)
plot(r, w1)
fig, ax = gcf(), gca()
xlabel("$r$"),  ylabel("$w(r)$")
title("Solution of the MEMS problem");
```

It's not necessary that the initialization satisfy the boundary conditions. In fact, by choosing a different constant function as the initial guess, we arrive at another valid solution.

```{code-cell}
r, w2 = FNC.bvp(phi, [a, b], ga, gb, 0.5 * init)
ax.plot(r, w2)
ax.set_title("Multiple solutions of the MEMS problem");
fig
```
