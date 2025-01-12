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
[**Demo %s**](#demo-absstab-advdiff)

The eigenvalues of advection-diffusion are near-imaginary for $\epsilon\approx 0$ and get closer to the negative real axis as $\epsilon$ increases.

```{code-cell}
:tags: [hide-input]
[x, Dx, Dxx] = diffper(40, [0, 1]);
tau = 0.1;
clf
for ep = [0.001 0.01 0.05]
  lambda = eig(-Dx + ep*Dxx);
  str = sprintf("\\epsilon = %.3f", ep);
  scatter(real(tau*lambda), imag(tau*lambda), displayname=str)
  hold on
end
axis equal,  grid on
legend(location='northwest')
title('Eigenvalues for advection-diffusion')
```
