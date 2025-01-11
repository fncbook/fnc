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
[**Demo %s**](#demo-linear-converge)

```{code-cell}
lamb = 10
exact = lambda x: sinh(lamb * x) / sinh(lamb) - 1
```

The following functions define the ODE.

```{code-cell}
p = lambda x: zeros(size(x))
q = lambda x: -(lamb**2) * ones(len(x))
r = lambda x: lamb**2 * ones(len(x))
```

We compare the computed solution to the exact one for increasing $n$.

```{code-cell}
N = array([int(2 * 10**d) for d in arange(1, 3.1, 0.25)])
err = zeros(len(N))
results = PrettyTable(["n", "error"])
for k, n in enumerate(N):
    x, u = FNC.bvplin(p, q, r, [0, 1], -1, 0, n)
    err[k] = norm(exact(x) - u, inf)
    results.add_row([n, err[k]])
print(results)
```

Each factor of 10 in $n$ reduces error by a factor of 100, which is indicative of second-order convergence.

```{code-cell}
loglog(N, err, "-o", label="observed")
loglog(N, 1 / N**2, "--", label="2nd order")
xlabel("$n$"),  ylabel("max error")
legend(),  title("Convergence of finite differences");
```
