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
[**Demo %s**](#demo-absstab-regions)

Euler and Backward Euler time-stepping methods were used to solve $\mathbf{u}'=\mathbf{D}_{xx}\mathbf{u}$.

```{code-cell}
m = 40;  
[x, Dx, Dxx] = diffper(m, [0, 1]);
```

The eigenvalues of this matrix are real and negative:

```{code-cell}
lambda = eig(Dxx);
clf
plot(real(lambda), imag(lambda), 'o')
axis equal, grid on
xlabel('Re \lambda'),  ylabel('Im \lambda')
title('Eigenvalues')
```

The Euler method is absolutely stable in the region $|\zeta+1| \le 1$ in the complex plane:

```{code-cell}
:tags: [hide-input]
phi = linspace(0, 2*pi, 361);
z = exp(1i*phi) - 1;   % unit circle shifted to the left by 1
fill(real(z), imag(z), [.8, .8, 1])
axis equal,  grid on
xlabel('Re \lambda'),  ylabel('Im \lambda')
title('Stability region')
```

In order to get inside this region, we have to find $\tau$ such that $\lambda \tau > -2$ for all eigenvalues $\lambda$. This is an upper bound on $\tau$. 

```{code-cell}
lambda_min = min(lambda)
max_tau = -2 / lambda_min
```

Here we plot the resulting values of $\zeta=\lambda \tau$. 

```{code-cell}
zeta = lambda * max_tau;
hold on
plot(real(zeta), imag(zeta), 'o')
title('Stability region and \zeta values')
```

In backward Euler, the region is $|\zeta-1|\ge 1$. Because they are all on the negative real axis, all of the $\zeta$ values will fit no matter what $\tau$ is chosen.

```{code-cell}
:tags: [hide-input]
clf,  fill([-6, 6, 6, -6],[-6, -6, 6, 6],[.8, .8, 1])
hold on
z = exp(1i*phi) + 1;   % unit circle shifted right by 1
fill(real(z), imag(z), 'w')
plot(real(zeta), imag(zeta), 'o')
axis equal,  axis([-4, 2, -3, 3]), grid on
xlabel('Re \lambda'),  ylabel('Im \lambda')
title('Stability region and \zeta values')
```
