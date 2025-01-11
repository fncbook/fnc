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
[**Demo %s**](#demo-improper-intinf)

```{code-cell}
:tags: [hide-input]
f = lambda x: 1 / (1 + x**2)
exact = pi
tol = array([1 / 10**d for d in arange(5, 14, 0.5)])
err = zeros((tol.size, 2))
length = zeros((tol.size, 2))
for k in range(tol.size):
    I1, x1 = FNC.intadapt(f, -2/tol[k], 2/tol[k], tol[k])
    I2, x2 = FNC.intinf(f, tol[k])
    err[k] = abs(exact - array([I1, I2]))
    length[k] = [x1.size, x2.size]
loglog(length, err, "-o")
# plot(len,err,m=:o,label=["direct" "double exponential"])
n = array([100, 10000])
loglog(n, 1000 / n**4, 'k--')
xlabel("number of nodes"),  ylabel("error")
title("Comparison of integration methods")
legend(["direct", "double exponential", "4th-order"], loc="lower left");
```

Both methods are roughly fourth-order due to Simpson's formula in the underlying adaptive integration method. At equal numbers of evaluation nodes, however, the double exponential method is consistently 2--3 orders of magnitude more accurate.
