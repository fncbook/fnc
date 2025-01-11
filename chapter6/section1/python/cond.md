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
[**Demo %s**](#demo-basics-cond)

Consider the ODEs $u'=u$ and $u'=-u$. In each case we compute $\partial f/\partial u = \pm 1$, so the condition number bound from {numref}`Theorem %s <theorem-depIC>` is $e^{b-a}$ in both problems. However, they behave quite differently. In the case of exponential growth, $u'=u$, the bound is the actual condition number.

```{code-cell}
t = linspace(0, 3, 200)
u = array([exp(t) * u0 for u0 in [0.7, 1, 1.3]])
plot(t, u.T)
xlabel("$t$")
ylabel("$u(t)$")
title(("Exponential divergence of solutions"));
```

But with $u'=-u$, solutions actually get closer together with time.

```{code-cell}
t = linspace(0, 3, 200)
u = array([exp(-t) * u0 for u0 in [0.7, 1, 1.3]])
plot(t, u.T)
xlabel("$t$")
ylabel("$u(t)$")
title(("Exponential convergence of solutions"));
```

In this case the actual condition number is one, because the initial difference between solutions is the largest over all time. Hence, the exponentially growing upper bound $e^{b-a}$ is a gross overestimate.
