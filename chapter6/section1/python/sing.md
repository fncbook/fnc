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
[**Demo %s**](#demo-basics-sing)


The equation $u'=(u+t)^2$ gives us some trouble.
```{tip}
:class: dropdown
It's a good idea to check `sol.success` after calling `solve_ivp`. If it's `False`, the solution may not be reliable. 
```

```{code-cell}
f = lambda t, u: (t + u) ** 2
sol = solve_ivp(f, [0.0, 1.0], [1.0])
if not sol.success:
    print(sol.message)
```

The warning message we received can mean that there is a bug in the formulation of the problem. But if everything has been done correctly, it suggests that the solution may not exist past the indicated time. This is a possibility in nonlinear ODEs.

```{code-cell}
semilogy(sol.t, sol.y[0, :])
xlabel("$t$")
ylabel("$u(t)$")
title(("Blowup in finite time"));
```
