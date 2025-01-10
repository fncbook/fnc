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
[**Demo %s**](#demo-absstab-advdiff)

The eigenvalues of advection-diffusion are near-imaginary for $\epsilon\approx 0$ and get closer to the negative real axis as $\epsilon$ increases.

```{code-cell}
:tags: [hide-input]
from scipy.linalg import eigvals
x, Dx, Dxx = FNC.diffper(40, [0, 1])
tau = 0.1
for ep in [0.001, 0.01, 0.05]:
    lamb = eigvals(-Dx + ep * Dxx)
    plot(real(tau * lamb), imag(tau * lamb), "o", label=f"epsilon={ep:.1g}")
axis("equal")
legend()
title("Eigenvalues for advection-diffusion")
```
