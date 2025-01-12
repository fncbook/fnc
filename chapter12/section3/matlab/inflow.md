---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
cd  /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init
```
[**Demo %s**](#demo-absstab-inflow)

Deleting the last row and column places all the eigenvalues of the discretization into the left half of the complex plane. 

```{code-cell}
[x, Dx, Dxx] = diffcheb(40, [0, 1]);
A = -Dx(2:end, 2:end);    % leave out first row and column
lambda = eig(A);
```

```{code-cell}
:tags: [hide-input]
clf
scatter(real(lambda), imag(lambda))
axis equal,  grid on 
title('Eigenvalues of advection with zero inflow')
```

Note that the rightmost eigenvalues have real part at most

```{code-cell}
max(real(lambda))
```

Consequently, all solutions decay exponentially to zero as $t\to\infty$. This matches our observation of the solution: eventually, everything flows out of the domain.

