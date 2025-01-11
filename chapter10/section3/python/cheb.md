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
[**Demo %s**](#demo-diffmats-cheb)

Here is a $4\times 4$ Chebyshev differentiation matrix.

```{code-cell}
t, Dx, Dxx = FNC.diffcheb(3, [-1, 1])
print(Dx)
```

We again test the convergence rate.

```{code-cell}
f = lambda x: x + exp(sin(4 * x))
df_dx = lambda x: 1 + 4 * exp(sin(4 * x)) * cos(4 * x)
d2f_dx2 = lambda x: 4 * exp(sin(4 * x)) * (4 * cos(4 * x) ** 2 - 4 * sin(4 * x))

N = range(5, 75, 5)
err1 = zeros(len(N))
err2 = zeros(len(N))
err = zeros((len(N), 2))
for k, n in enumerate(N):
    t, Dx, Dxx = FNC.diffcheb(n, [-1, 1])
    y = f(t)
    err[k, 0] = norm(df_dx(t) - Dx @ y, inf)
    err[k, 1] = norm(d2f_dx2(t) - Dxx @ y, inf)
```

Since we expect a spectral convergence rate, we use a semi-log plot for the error.

```{code-cell}
semilogy(N, err, "-o")
xlabel("$n$"), ylabel("max error")
legend(["$f'$", "$f''$"], loc="lower left")
title("Convergence of Chebyshev derivatives");
```
