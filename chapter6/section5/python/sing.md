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
[**Demo %s**](#demo-adapt-sing)

In {numref}`Demo %s <demo-basics-sing>` we saw an IVP that appears to blow up in a finite amount of time. Because the solution increases so rapidly as it approaches the blowup, adaptive stepping is required even to get close.

```{code-cell}
du_dt = lambda t, u: (t + u)**2
tspan = (0.0, 2.0)
u0 = [1.0]
t, u = FNC.rk23(du_dt, tspan, u0, 1e-5)
```

In fact, the failure of the adaptivity gives a decent idea of when the singularity occurs.

```{code-cell}
semilogy(t, u[0, :])
xlabel("$t$"),  ylabel("$u(t)$")
title("Finite-time blowup")

tf = t[-1]
axvline(x=tf, color='k', linestyle='--', label=f"t = {tf:.6f}")
legend();
```
