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
[**Demo %s**](#demo-absstab-inflow)

Deleting the last row and column places all the eigenvalues of the discretization into the left half of the complex plane. 

```{code-cell}
from scipy.linalg import eigvals
x, Dx, Dxx = FNC.diffcheb(40, [0, 1])
A = -Dx[1:, 1:]  # leave out first row and column

lamb = eigvals(A)
plot(real(lamb), imag(lamb), "o")
xlim(-300, 100),  axis("equal"),  grid(True)
title("Eigenvalues of advection with zero inflow");
```

Note that the rightmost eigenvalues have real part at most

```{code-cell}
print(f"rightmost extent of eigenvalues: {max(real(lamb)):.3g}")
```

Consequently all solutions decay exponentially to zero as $t\to\infty$. This matches our observation of the solution: eventually, everything flows out of the domain.

