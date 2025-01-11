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
[**Demo %s**](#demo-adapt-usage)

We'll integrate the function from {numref}`Demo %s <demo-adapt-motive>`.

```{code-cell}
from scipy.integrate import quad
f = lambda x: (x + 1) ** 2 * cos((2 * x + 1) / (x - 4.3))
I, errest = quad(f, 0, 4, epsabs=1e-12, epsrel=1e-12)
print("integral:", I)    # 'exact' value
```

We perform the integration and show the nodes selected underneath the curve.

```{code-cell}
Q, t = FNC.intadapt(f, 0, 4, 0.001)
print("number of nodes:", t.size)

x = linspace(0, 4, 600)
plot(x, f(x), "k")
stem(t, f(t))
xlabel("$x$"); ylabel("$f(x)$");
```

The error turns out to be a bit more than we requested. It's only an estimate, not a guarantee.

```{code-cell}
print("error:", I - Q)
```

Let's see how the number of integrand evaluations and the error vary with the requested tolerance.

```{code-cell}
tol_ = 10.0 ** arange(-4, -12, -1)
err_ = zeros(tol_.size)
num_ = zeros(tol_.size, dtype=int)
print("    tol         error     # f-evals")
for i, tol in enumerate(tol_):
    Q, t = FNC.intadapt(f, 0, 4, tol)
    err_[i] = I - Q
    num_[i] = t.size
    print(f"  {tol:6.1e}    {err_[i]:10.3e}    {num_[i]:6d}")
```

As you can see, even though the errors are not smaller than the estimates, the two columns decrease in tandem. If we consider now the convergence not in $h$, which is poorly defined now, but in the number of nodes actually chosen, we come close to the fourth-order accuracy of the underlying Simpson scheme.

```{code-cell}
loglog(num_, abs(err_), "-o", label="results")
order4 = 0.01 * (num_ / num_[0]) ** (-4)
loglog(num_, order4, "--", label="$O(n^{-4})$")
xlabel("number of nodes"), ylabel("error")
legend()
title("Convergence of adaptive quadrature");
```
