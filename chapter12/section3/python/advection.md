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
[**Demo %s**](#demo-absstab-advection)

For $c=1$ we get purely imaginary eigenvalues.

```{code-cell}
:tags: [hide-input]
from scipy.linalg import eigvals
x, Dx, Dxx = FNC.diffper(40, [0, 1])
lamb = eigvals(Dx)

plot(real(lamb), imag(lamb), "o")
axis([-40, 40, -40, 40])
axis("equal")
title("Eigenvalues for pure advection");
```

Let's choose a time step of $\tau=0.1$ and compare to the stability regions of the Euler and backward Euler time steppers (shown as shaded regions):

```{code-cell}
:tags: [hide-input]
zc = exp(2j * pi * arange(361) / 360)
# points on |z|=1

z = zc - 1    # shift circle left by 1
fill(real(z), imag(z), color=(0.8, 0.8, 1))
plot(real(0.1 * lamb), imag(0.1 * lamb), "o")
axis([-5, 5, -5, 5]),  axis("equal")
title("Euler");
```

In the Euler case it's clear that *no* real value of $\tau>0$ is going to make $\zeta$ values fit within the stability region. Any method whose stability region includes none of the imaginary axis is an unsuitable choice for advection.

```{code-cell}
:tags: [hide-input]
z = zc + 1    # shift circle right by 1
fill([-6, 6, 6, -6], [-6, -6, 6, 6], color=(0.8, 0.8, 1))
fill(real(z), imag(z), color="w")
plot(real(0.1 * lamb), imag(0.1 * lamb), "o")
axis([-5, 5, -5, 5])
axis("equal")
title("Backward Euler");
```

The A-stable backward Euler time stepping tells the exact opposite story: it will be absolutely stable for any choice of the time step $\tau$.
