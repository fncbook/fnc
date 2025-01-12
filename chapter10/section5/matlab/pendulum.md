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
[**Demo %s**](#demo-nonlinear-pendulum)

Suppose a damped pendulum satisfies the nonlinear equation $\theta'' + 0.05\theta'+\sin \theta =0$. We want to start the pendulum at $\theta=2.5$ and give it the right initial velocity so that it reaches $\theta=-2$ at exactly $t=5$. This is a boundary-value problem with Dirichlet conditions $\theta(0)=2.5$ and $\theta(5)=-2$.

The first step is to define the function $\phi$ that equals $\theta''$.

```{code-cell}
phi = @(t,theta,omega) -0.05 * omega - sin(theta);
```

Next, we define the boundary conditions.

```{code-cell}
ga = @(u, du) u - 2.5;
gb = @(u, du) u + 2;
```

The last ingredient is an initial estimate of the solution. Here we choose $n=100$ and a linear function between the endpoint values. 

```{code-cell}
init = linspace(2.5, -2, 101)';
```

We find a solution with negative initial slope, i.e., the pendulum is initially pushed back toward equilibrium.

```{code-cell}
[t, theta] = bvp(phi, 0, 5, ga, gb, init);
clf,  plot(t, theta)
xlabel('t'),  ylabel('\theta(t)')
title('Pendulum over [0,5]')
```

If we extend the time interval longer for the same boundary values, then the initial slope must adjust.

```{code-cell}
[t, theta] = bvp(phi, 0, 8, ga, gb, init);
plot(t,theta)
xlabel('t'),  ylabel('\theta(t)')
title('Pendulum over [0,8]')
```

This time, the pendulum is initially pushed toward the unstable equilibrium in the upright vertical position before gravity pulls it back down.
