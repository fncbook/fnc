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
[**Demo %s**](#demo-absstab-advection)

For $c=1$ we get purely imaginary eigenvalues.

```{code-cell}
:tags: [hide-input]
[x, Dx] = diffper(40, [0, 1]);
lambda = eig(Dx);
clf
scatter(real(lambda), imag(lambda))
axis equal,  grid on
title('Eigenvalues for pure advection')
```

Let's choose a time step of $\tau=0.1$ and compare to the stability regions of the Euler and backward Euler time steppers (shown as shaded regions):

```{code-cell}
:tags: [hide-input]
zc = exp( 1i * linspace(0,2*pi,361)' );    % points on |z|=1
z = zc - 1;    % shift circle left by 1
clf,  fill(real(z), imag(z), [.8, .8, 1])
hold on,  scatter(real(0.1*lambda), imag(0.1*lambda))
axis equal,  axis square
axis([-5, 5, -5, 5]),  grid on
title('Euler')
```

In the Euler case it's clear that *no* real value of $\tau>0$ is going to make $\zeta$ values fit within the stability region. Any method whose stability region includes none of the imaginary axis is an unsuitable choice for advection.

```{code-cell}
:tags: [hide-input]
clf,  fill([-6, 6, 6, -6],[-6, -6, 6, 6], [.8, .8, 1])
z = zc + 1;   % shift circle right by 1
hold on,  scatter(real(0.1*lambda), imag(0.1*lambda))
fill(real(z), imag(z), 'w')
axis equal,  axis square
axis([-5 5 -5 5]),  grid on
title('backward Euler')
```

The A-stable backward Euler time stepping tells the exact opposite story: it will be absolutely stable for any choice of the time step $\tau$.
