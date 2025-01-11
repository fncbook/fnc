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
[**Demo %s**](#demo-stability-errfun)

We plot $|\Phi(x)|$ over the interval $[-1,1]$ with equispaced nodes for different values of $n$. 

```{code-cell} 
:tags: [hide-input]
x = linspace(-1, 1, 1601)
for n in range(10, 60, 10):
    t = linspace(-1, 1, n + 1)
    Phi = array([prod(xk - t) for xk in x])
    semilogy(x, abs(Phi), ".", markersize=2)
xlabel("$x$")
ylabel("$|\Phi(x)|$")
ylim([1e-25, 1])
title(("Effect of equispaced nodes"));
```

Each time $\Phi$ passes through zero at an interpolation node, the value on the log scale should go to $-\infty$, which explains the numerous cusps on the curves.
