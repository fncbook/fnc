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
[**Demo %s**](#demo-shooting-unstable)


```{code-cell}
ga = lambda u, du : u + 1    # u=-1 at left
gb = lambda u, du : u        # u= 0 at right
init = array([-1, 0])
for lamb in range(6, 22, 4):
    phi = lambda x, u, du_dx: lamb**2 * u + lamb**2
    x, u, du_dx = FNC.shoot(phi, 0.0, 1.0, ga, gb, init)
    plot(x, u, label=f"$\\lambda$ = {lamb:.1f}")

xlabel("$x$"),  ylabel("$u(x)$"),  ylim(-1.0, 0.25)
grid(True),  legend(loc="upper left")
title("Shooting instability");
```

The numerical solutions evidently don't satisfy the right boundary condition as $\lambda$ increases, which makes them invalid. 

