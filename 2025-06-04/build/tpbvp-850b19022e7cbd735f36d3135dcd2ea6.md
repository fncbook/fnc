---
numbering:
  enumerator: 10.1.%s
---
(section-bvp-tpbvp)=

# Two-point BVP

The initial-value problems of [Chapter 6](@chapter-ivp) are characterized by an ordinary differential equation plus a value of the solution's state at one value of the independent variable (typically, time).

In a {term}`boundary-value problem`, the state is not entirely given at any point. Instead, partial information is given at multiple values of the independent variable. We will focus on the most common type.

## Definition

::::{prf:definition} Two-point boundary-value problem: TPBVP
:label: definition-tpbvp

```{math}
:label: tpbvp
\begin{split}
u''(x) &= \phi(x,u,u'), \qquad a \le x \le b,\\
g_1(u(a),u'(a)) &= 0,\\
g_2(u(b),u'(b)) &= 0.
\end{split}
```

The functions $g_1$ and $g_2$ are called **boundary conditions**.

This TPBVP is said to be **linear** if the dependence of $\phi$, $g_1$, and $g_2$ on the solution $u(x)$ is linear. Specifically, this means that

```{math}
\phi(x, u, u') = p(x)u' + q(x)u + r(x)
```

for some coefficient functions $p$, $q$, and $r$, and

```{math}
g_i(u, u') = \alpha_i u + \beta_i u' - \gamma_i, \quad i=1,2,
```

for some constants $\alpha_i$, $\beta_i$, and $\gamma_i$.
::::

:::{note}
Often, the domain of $x$ in {eq}`tpbvp` is not explicitly stated but just implied by the appearances of $a$ and $b$ in $g_1$ and $g_2$.
:::

An IVP for the same ODE as in {eq}`tpbvp` would specify values for both $u(a)$ and $u'(a)$. Although this may look like a minor difference, the characters of IVPs and BVPs are substantially different. In a typical IVP, in which the independent variable is often time, the initial value determines everything about the future course of the solution. In a TPBVP, however, complete information is not given at any one location in the domain, and there may be more than one way (or no way) to satisfy the boundary conditions.

Certain special cases of the boundary conditions have their own nomenclature.

```{index} ! Dirichlet boundary condition, ! Neumann boundary condition, ! Robin boundary condition, homogeneous boundary condition
```

:::{prf:definition} Dirichlet, Neumann, and Robin conditions
:label: definition-bctypes
Let $g_i$ be a boundary condition in {eq}`tpbvp`, and let $\alpha,\beta$ be constants.  

{term}`Dirichlet condition`: $g_i(u,u') = u - \gamma$

{term}`Neumann condition`: $g_i(u,u') = u' - \gamma$

**Robin condition**: $g_i(u,u') = u + \beta u' - \gamma$

When $\gamma=0$ in the above, it is a {term}`homogeneous boundary condition`.
:::

::::{prf:example} Multiple solutions in a BVP
:label: example-tpbvp-pendulum
An ideal pendulum of length $L$ satisfies the ODE $\theta'' + \frac{g}{L}\sin \theta = 0$, where $\theta(t)$ is the angle of the pendulum's rod from straight downward, and $g$ is gravitational acceleration.

If we pull the pendulum bob 1 radian and release it from rest, then we have an IVP with the initial condition $\theta(0)=1$, $\theta'(0)=0$. Everything about the future trajectory of the pendulum is completely determined by that condition. But if instead we want to know how far to pull up the pendulum bob initially so that it is in the downward position 2 seconds later, then we have the Neumann boundary condition $\theta'(0)=0$ and the Dirichlet boundary condition $\theta(2)=0$.

Clearly, one solution is $\theta(0)=0$, and the pendulum never moves! But there is also an initial deflection such that $\theta(2)=0$ momentarily as the pendulum finishes its first downward swing. There could be still more solutions in which the pendulum completes one or more half-periods before passing again through the downward position at just the right moment.
::::

While time can be the independent variable in a TPBVP, as in {numref}`Example {number} <example-tpbvp-pendulum>`, it is often space, which has no intrinsic direction of information flow.

::::{prf:example} MEMS actuator
:label: example-tpbvp-mems
A micro-electromechanical system (MEMS) can have an actuator consisting of two flat, parallel, disk-shaped surfaces, one at $z=0$ and the other at $z=1$. The surface at $z=0$ is a rigid metal plate. The surface at $z=1$ is an elastic membrane fixed only at its boundary. When the surfaces are given different electric potentials, the membrane deflects in response to the induced electric field, and the field is also a position of the deflection. Assuming circular symmetry, one may derive the ordinary differential equation {cite}`peleskoEffectSmallaspectratio2006`

```{math}
:label: mems
\frac{d^2w}{d r^2} + \frac{1}{r}\, \frac{d w}{d r} = \frac{\lambda}{w^2},
```

where $w(r)$ is the vertical position of the membrane at radius $r$, and $\lambda$ is proportional to the applied electric potential. We will use $0\le r\le 1$. There are also the supplemental physical conditions

```{math}
:label: memsbc
w'(0)=0, \qquad w(1)=1,
```

derived from the circular symmetry and fixing the edge of the membrane, respectively. This is a nonlinear boundary-value problem. The solution should respect $0<w(r)\le 1$ in order to be physically meaningful.
::::

## Numerical solution

We can solve the TPBVP {eq}`tpbvp` by recasting the problem as a first-order system, as introduced in @transformation-of-high-order-systems.

::::{prf:example} Solving a TPBVP numerically
:label: demo-tpbvp-mems
As a system, the MEMS problem from {numref}`Example {number} <example-tpbvp-mems>` uses $y_1=w$, $y_2=w'$ to obtain

```{math}
:label: memssys
\begin{split}
y_1' &= y_2, \\
y_2' &= \frac{\lambda}{y_1^2} - \frac{y_2}{r}.
\end{split}
```

This has the form $\mathbf{y}' = \mathbf{f}(r,\mathbf{y})$. In the original variable, the boundary conditions are $w'(0)=0$ and $w(1)=1$. We have to encode these in the new variable as $y_2(0)=0$ and $y_1(1)-1 = 0$.

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-tpbvp-mems-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-tpbvp-mems-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-tpbvp-mems-python
:::
````
`````

::::

Characterizing the conditioning of a TPBVP theoretically is difficult. There are some numerical tools going by the name *sensitivity analysis*, but the details are too lengthy for us to explore here.

## Exercises

```{index} Bessel equation, Carrier equation
```

``````{exercise}
:label: problem-tpbvp-classify
✍ In each case, explain whether the TPBVP is linear or nonlinear.  

**(a)** $x^2 u'' +xu' + (x^2 - 1) u = 0, \quad u(0) =1,\; u(4) =0 \qquad$ [Bessel equation](wiki:Bessel_function)

**(b)** $u'' - u u' = 1, \quad u(0) = u'(1) = 1$

**(c)** $u'' + u = 1, \quad u(0) = u(1) u'(1) = 1$

**(d)** $\epsilon u'' +2(1-x^2) u +u^2  = 1, \quad u(-1) = u(1) = 0 \qquad$ *(Carrier equation {cite}`carrierSingularPerturbation1970`)*

**(e)** $u u'' = 3(u')^2, \quad u(-1) = 1, \; u(2) = \frac{1}{2}$
``````

```{index} Hermite equation, Laguerre equation, Chebyshev equation, Legendre equation
```

``````{exercise}
:label: problem-tpbvp-verify
✍ For each BVP, verify that the given solution is valid, i.e., check that the differential equation and the boundary conditions are satisfied.

**(a)** $u'' - 2xu' + 8 u = 0, \quad u(0) = 1, \; u(1) = -\frac{5}{3}\qquad$ [Hermite equation](wiki:Hermite_polynomials) 

Solution: $u(x) = \frac{4}{3}x^4-4x^2+1$

**(b)** $xu'' + (1-x) u' + 3 u = 0, \quad u(0) =1, \; u'(1) =-\frac{1}{2} \qquad$ [Laguerre equation](wiki:Laguerre_polynomials)

Solution: $u(x) = \frac{1}{6}(6-18x+9x^2-x^3)$

**(c)** $(1-x^2)u'' - xu' + 25u = 0, \quad u'(0) = 5, \; u\left(\frac{1}{2}\right) = \frac{1}{2} \qquad$ [Chebyshev equation](wiki:Chebyshev_polynomials)

Solution: $u(x) = T_5(x) = 16x^5-20x^3+5x$ 

**(d)** $(1-x^2) u'' -2xu' + 12 u = 0, \quad u(0)+2u'(0) = -3, \; u(1) = 1 \qquad$ [Legendre equation](wiki:Legendre_polynomials)

Solution: $u(x) = P_3(x) = \frac{1}{2}(5x^3-3x)$ 

**(e)** $u u'' = 3(u')^2, \quad u(-1) = 1, \; u(2) = \frac{1}{2}$ 

Solution: $u(x) =  ( x+2 )^{-1/2}$
``````

``````{exercise}
:label: problem-tpbvp-usage
⌨ For each TPBVP in @problem-tpbvp-verify, use `solve`/`bvp4c`/`solve_bvp` to find the solution, following the method in @demo-tpbvp-mems. Plot the solution and separately plot the error as a function of $x$. (Important: In some cases, you will need to truncate the domain slightly to avoid division by zero.)
``````

``````{exercise}
:label: problem-tpbvp-pendulum
⌨ Consider the pendulum from {numref}`Example {number} <example-tpbvp-pendulum>` with $g=L=1$. Suppose we want to release the pendulum from rest such that $\theta(5)=\pi/2$. Find one solution that passes through $\theta=0$ and another solution that does not. Plot $\theta(t)$ for both cases together. (Hint: Vary the initial estimate for the solution.)
``````

```{index} Allen–Cahn equation
```

``````{exercise}
:label: problem-tpbvp-allencahn
⌨  The stationary [Allen–Cahn equation](wiki:Allen-Cahn_equation) is a model of phase changes, such as the change from liquid to solid. In one spatial dimension it can be written as

$$
\epsilon u'' = u^3-u, \qquad 0 \le x \le 1, \qquad u(0)=-1, \quad u(1)=1.
$$

As $\epsilon\rightarrow 0$, the solution tends toward a step function transition between $-1$ and $1$. By symmetry, $u'(x)=-u'(1-x)$.

**(a)** Use {numref}`Function {number} <function-shoot>` with initial solution estimates $(u=-1,u'=0)$ to solve the equation for $\epsilon=0.2$. Plot the solution.

**(b)** Repeat part (a) for $\epsilon=0.02$.

**(c)** Repeat part (a) for $\epsilon=0.002$. (This is a difficult problem for the default method.) Try a few different initializations for $u$, and plot all the results. Do any seem to be valid solutions of the BVP?
``````
