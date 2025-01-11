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
[**Demo %s**](#demo-polynomial-lagrange)

Here is a vector of nodes.

```{code-cell}
t = array([1, 1.5, 2, 2.25, 2.75, 3])
n = 5
```

Let's apply the definition of the cardinal Lagrange polynomial for $k=2$. First we define a polynomial $q$ that is zero at all the nodes except $i=k$. Then $\ell_2$ is found by normalizing $q$ by $q(t_k)$.

```{code-cell}
k = 2
q = lambda x: prod([x - t[i] for i in range(n + 1) if i != k])
ell_k = lambda x: q(x) / q(t[k])
```

A plot confirms the cardinal property of the result.

```{code-cell}
x = linspace(1, 3, 500)
plot(x, [ell_k(xx) for xx in x])
y = zeros(n+1)
y[k] = 1
plot(t, y, "ko")
xlabel("$x$"),  ylabel("$\\ell_2(x)$")
title(("Lagrange cardinal function"));
```

Observe that $\ell_k$ is _not_ between zero and one everywhere, unlike a hat function.
