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
[**Demo %s**](#demo-integration-ellipse)

```{code-cell}
f = lambda t: pi * sqrt(cos(pi * t) ** 2 + sin(pi * t) ** 2 / 4)
N = arange(4, 48, 6)
perim = zeros(N.size)
for k in range(N.size):
    h = 2 / N[k]
    t = h * arange(N[k]) - 1
    perim[k] = h * sum(f(t))
err = abs(perim - perim[-1])    # use the last value as reference
results = PrettyTable()
results.add_column("N", N)
results.add_column("perimeter", perim)
results.add_column("error", err)
results
```
The approximations gain about one digit of accuracy for each constant increment of $n$, which is consistent with spectral convergence.
