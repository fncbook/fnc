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
[**Demo %s**](#demo-adapt-basic)

Let's run adaptive RK on  $u'=e^{t-u\sin u}$.

```{code-cell}
f = lambda t, u: exp(t - u * sin(u))
t, u = FNC.rk23(f, [0.0, 5.0], [0.0], 1e-5)
scatter(t, u[0, :])
xlabel("$t$"), ylabel("$u(t)$")
title(("Adaptive IVP solution"));
```

The solution makes a very abrupt change near $t=2.4$. The resulting time steps vary over three orders of magnitude.

```{code-cell}
dt = [t[i + 1] - t[i] for i in range(t.size - 1)]
semilogy(t[:-1], dt)
xlabel("$t$"), ylabel("time step")
title(("Adaptive step sizes"));
```

If we had to run with a uniform step size to get this accuracy, it would be

```{code-cell}
print(f"min step size was {min(dt):.2e}")
```

On the other hand, the average step size that was actually taken was

```{code-cell}
print(f"mean step size was {mean(dt):.2e}")
```

We took fewer steps by a factor of 1000! Even accounting for the extra stage per step and the occasional rejected step, the savings are clear.

