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
[**Demo %s**](#demo-nonlinear-pendulum)

The first step is to define the function $\phi$ that equals $\theta''$.

```{code-cell}
phi = lambda t, theta, omega: -0.05 * omega - sin(theta)
```

Next, we define the boundary conditions.

```{code-cell}
ga = lambda u, du: u - 2.5
gb = lambda u, du: u + 2
```

The last ingredient is an initial estimate of the solution. Here we choose $n=100$ and a linear function between the endpoint values. 

```{code-cell}
init = linspace(2.5, -2, 101)
```

We find a solution with negative initial slope, i.e., the pendulum is initially pushed back toward equilibrium.

```{code-cell}
t, theta = FNC.bvp(phi, [0, 5], ga, gb, init)
plot(t, theta)
xlabel("$t$")
ylabel("$\theta(t)$")
title("Pendulum over [0,5]");
```

If we extend the time interval longer for the same boundary values, then the initial slope must adjust.

```{code-cell}

t, theta = FNC.bvp(phi, [0, 8], ga, gb, init)
plot(t, theta)
xlabel("$t$")
ylabel("$\theta(t)$")
title("Pendulum over [0,8]");
```

This time, the pendulum is initially pushed toward the unstable equilibrium in the upright vertical position before gravity pulls it back down.
