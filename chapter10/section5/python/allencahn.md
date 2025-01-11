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
[**Demo %s**](#demo-nonlinear-allencahn)


```{code-cell}
phi = lambda x, u, dudx: (u**3 - u) / epsilon
ga = lambda u, du: du
gb = lambda u, du: u - 1
```

Finding a solution is easy at larger values of $\epsilon$.

```{code-cell}
epsilon = 0.05
init = linspace(-1, 1, 141)
x, u1 = FNC.bvp(phi, [0, 1], ga, gb, init)

plot(x, u1, label="$\\epsilon = 0.05$")
fig, ax = gcf(), gca()
xlabel("$x$"),  ylabel("$u(x)$")
legend(),  title("Allen-Cahn solution");
```

Finding a good initialization is not trivial for smaller values of $\epsilon$. But the iteration succeeds if we use the first solution as the initialization at the smaller $\epsilon$.


```{code-cell}
epsilon = 0.002
x, u2 = FNC.bvp(phi, [0, 1], ga, gb, u1)
ax.plot(x, u2, label="$\\epsilon = 0.002$")
ax.legend()
fig
```

In this case we can continue further.

```{code-cell}
ϵ = 0.0005
x, u3 = FNC.bvp(phi, [0, 1], ga, gb, u2)
ax.plot(x, u3, label="$\\epsilon = 0.005$")
ax.legend()
fig
```
