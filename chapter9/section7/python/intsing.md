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
[**Demo %s**](#demo-improper-intsing)


```{code-cell}
:tags: [hide-input]
f = lambda x: 1 / (10 * sqrt(x))
exact = 0.2
tol = array([1 / 10**d for d in arange(5, 14, 0.5)])
err = zeros((tol.size, 2))
length = zeros((tol.size, 2))
for k in range(tol.size):
    I1, x1 = FNC.intadapt(f, (tol[k]/20)**2, 1, tol[k])
    I2, x2 = FNC.intsing(f, tol[k])
    err[k] = abs(exact - array([I1, I2]))
    length[k] = [x1.size, x2.size]
loglog(length, err, "-o")
# plot(len,err,m=:o,label=["direct" "double exponential"])
n = array([100, 10000])
loglog(n, 30 / n**4, 'k--')
xlabel("number of nodes"),  ylabel("error")
title("Comparison of integration methods")
legend(["direct", "double exponential", "4th-order"], loc="lower left");
```

As in {numref}`Demo {number} <demo-improper-intinf>`, the double exponential method is more accurate than direct integration by a few orders of magnitude. Equivalently, the same accuracy can be reached with many fewer nodes.
