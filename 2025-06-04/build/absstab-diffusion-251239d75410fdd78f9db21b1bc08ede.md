---
numbering:
  enumerator: 11.3.%s
---
(section-diffusion-absstab)=
# Absolute stability

```{index} stability; of IVP solvers
```

In {numref}`section-diffusion-methodlines` we applied several different time stepping methods to a linear, constant coefficient problem in the form 

$$
\mathbf{u}'(t)=\mathbf{A}\mathbf{u}(t).
$$

All of these methods are zero-stable in the sense of {numref}`section-ivp-zerostability`, in the limit as the time step size $\tau \to 0$.[^h2tau] Yet for some experiments with *fixed* $\tau$, as in @demo-methodlines-heatFE, we have observed exponential growth in the different limit $n\to \infty$.

[^h2tau]: In Chapter 6 we used $h$ rather than $\tau$ to denote the time step size, but now we  reserve $h$ for spacing in the $x$ direction.

```{index} eigenvalue decomposition
```
Observe that if $\mathbf{A}$ has the eigenvalue decomposition $\mathbf{A}=\mathbf{V}\mathbf{D}\mathbf{V}^{-1}$, then

\begin{align*}
  \mathbf{u}'&=(\mathbf{V}\mathbf{D}\mathbf{V}^{-1})\mathbf{u},\\
  (\mathbf{V}^{-1} \mathbf{u}') &= \mathbf{D} (\mathbf{V}^{-1} \mathbf{u}), \\
  \mathbf{y}' &= \mathbf{D} \mathbf{y},
\end{align*}

where $\mathbf{y}(t)=\mathbf{V}^{-1}\mathbf{u}(t)$. Because $\mathbf{D}$ is diagonal, the dynamics of the components of $\mathbf{y}$ are completely decoupled: each row is a self-contained equation of the form $y_j'=\lambda_j y_j$, where $\lambda_j$ is an eigenvalue of $\mathbf{A}$. 

The diagonalization argument suggests that we can look at the scalar problems 

:::{math}
:label: absstabmodel
y' = \lambda y, \quad y(0)=1,
:::

arising from the eigenvalues. These eigenvalues may not be real numbers, so in this section, $i$ stands for the imaginary unit, not an integer index. If we write $\lambda$ in real and imaginary parts as $\lambda=\alpha + i\beta$, then by Euler's identity, the exact solution of {eq}`absstabmodel` has magnitude

$$
   \bigl |e^{(\alpha+i\beta)t} \bigr| = \bigl |e^{\alpha t} \bigr| \cdot \bigl |e^{i \beta t} \bigr| = e^{\alpha t}.
$$

::::{prf:observation}
Solutions of {eq}`absstabmodel` are bounded as $t\to\infty$ if and only if $\alpha = \operatorname{Re} \lambda \le 0$. 
::::

We now consider the counterpart of this observation for the solution produced by a numerical IVP solver.


```{index} ! absolute stability
```

::::{prf:definition} Absolute stability
:label: definition-absolutestability
Let $\lambda$ be a complex number, and let $y_0,y_1,y_2,\ldots,y_n$ be the numerical solution at times $0,\tau,2\tau,\ldots,n\tau$ of {eq}`absstabmodel` using a Runge–Kutta or multistep method with fixed stepsize $\tau$. Then the method is said to have {term}`absolute stability` at $\zeta = \tau\lambda$ if $|y_n|$ is bounded above as $n\to\infty$. 
::::

The fact that absolute stability depends only on the product $\zeta = \tau\lambda$, and not independently on the individual factors, is a result of how the IVP solvers are defined, as we will see below. Since $\lambda$ has units of inverse time according to {eq}`absstabmodel`, $\zeta$ is dimensionless.
## Stability regions

Each numerical IVP solver has its own collection of $\zeta$ values for which it is absolutely stable.

```{index} ! stability region
```

::::{prf:definition} Stability region
The **stability region** of an IVP solver is the collection of all $\zeta\in\complex$ for which the method is absolutely stable.
::::

::::{prf:example}
:label: example-absstab-euler
Consider an Euler discretization of $y'=\lambda y$:
  
$$
  y_{k+1} = y_k + \tau( \lambda y_k) =   (1+ \zeta ) y_k.
$$

Given that $y_0=1$ by {eq}`absstabmodel`, we easily deduce that $y_k = (1+\zeta)^k$ for all $k$, and therefore

$$
|y_k| = |1+\zeta|^k.
$$

Hence $|y_k|$ remains bounded above as $k\to \infty$ if and only if $|1+\zeta| \le 1$. Because $\zeta$ is a complex number, it's easiest to interpret this condition geometrically:

$$
  |\zeta + 1 | = |\zeta - (-1) | \le 1.
$$

That is, the distance in the plane from $\zeta$ to the point $-1$ is less than or equal to 1. This description defines a closed disk of radius 1 centered at $(-1,0)$.
::::

::::{prf:example}
:label: example-absstab-AM1
The backward Euler method discretizes {eq}`absstabmodel` as

$$
y_{k+1} = y_k + \tau( \lambda y_{k+1}) \quad \Rightarrow \quad y_{k+1} =  \frac{1}{1-\zeta} y_k.
$$

Therefore, $y_k=(1-\zeta)^{-k}$ for all $k$, and absolute stability requires $|1-\zeta|^{-1} \le 1$, or 

$$
|\zeta-1|\ge 1.
$$

This inequality describes the region *outside* of the open disk of radius 1 centered at $1$ on the real axis of the complex plane.
::::

::::{prf:example} 
:label: example-absstab-IE2
The improved Euler method IE2 defined in {eq}`IE` discretizes {eq}`absstabmodel` as 

```{math}
{y}_{i+1} = y_i +  \zeta \left( y_i + \tfrac{1}{2}\zeta y_i \right) = (1 + \zeta + \tfrac{1}{2}\zeta^2) y_i.
```

The stability region consists of all $\zeta$ such that $| 1 + \zeta + \tfrac{1}{2}\zeta^2 | \le 1$. Although it is not elementary to describe this region geometrically, its boundary points satisfy

$$
1 - e^{i\theta} + \zeta + \tfrac{1}{2}\zeta^2 = 0
$$

for some real $\theta$, and thus we can use the quadratic formula to find all the boundary points. 
::::

Stability regions for the most common IVP integrators are given in {numref}`figure-stabreg_ab_am` and {numref}`figure-stabreg_bd_rk`.  Note that those for the implicit Adams-Moulton methods are larger than those for the explicit Adams-Bashforth methods of the same order.  For the implicit backward differentiation methods, the exteriors of the curves provide large regions of stability, but significant portions of the imaginary axis may be excluded.  Finally, while the single-step Runge-Kutta methods have smaller regions of stability, those of orders 3 and 4 do include significant portions of the imaginary axis.

```{figure} figures/stabreg_ab_am.svg
:name: figure-stabreg_ab_am
Stability regions for Adams–Bashforth methods of order 1–4 (left) and Adams–Moulton methods of order 2–5 (right). The plots are in the complex $\zeta$-plane.
```

```{figure} figures/stabreg_bd_rk.svg
:name: figure-stabreg_bd_rk
Stability regions for backward differentiation methods of order 1–4 (left, exteriors of curves) and Runge–Kutta methods of order 1–4 (right). The plots are in the complex $\zeta$-plane.
```

For any particular method and value of $\lambda$ in {eq}`absstabmodel`, we can use the stability region to deduce which, if any, values of the time step $\tau$ will give bounded solutions. Both the magnitude and the argument (angle) of $\lambda$ play a role in determining such constraints.

::::{prf:example}
:label: example-absstab-FEBE
Suppose $\lambda=-4$ and Euler's method is applied. Since the time step is always positive, $\zeta=-4\tau$ is always on the negative real axis. The only part of that line that lies within the stability region of Euler as derived in {numref}`Example {number} <example-absstab-euler>` is the real interval $[-2,0]$. Hence we require $\zeta\ge -2$, or $\tau \le 1/2$. By contrast, the stability region of backward Euler includes the entire negative real axis, so absolute stability is unconditional, i.e., assured regardless of $\tau$.

Now suppose instead that $\lambda=i$, so that $\zeta=i\tau$. Clearly $\zeta$ is always on the positive imaginary axis. But no part of this axis, aside from the origin, lies in the stability region of Euler's method, so it is unconditionally *unstable* in this circumstance. The conclusion for backward Euler is the opposite; any value of $\tau$ will do, because the entire imaginary axis is within the stability region.
::::

{numref}`Example %s <example-absstab-FEBE>` does not contradict our earlier statements about the zero stability and convergence of Euler's method in general, even for the case $\lambda=i$. But those statements are based on the limit $\tau\to 0$ for $t$ in a finite interval $[a,b]$. Both this limit and the limit $t\to \infty$ imply the number of steps $n$ goes to infinity, but the limits behave differently.

The fact that implicit methods have larger stability regions than their explicit counterparts is the primary justification for using them. While they have larger work requirements per step, they sometimes can take steps that are orders of magnitude larger than explicit methods and still remain stable.

When adaptive time stepping methods are used, as in most software for IVPs, the automatically determined time step is chosen to satisfy absolute stability requirements (otherwise errors grow exponentially). This phenomenon was manifested in @demo-methodlines-auto: in the explicit IVP method `rk23`, error control forced tiny step sizes compared to those used by `Rodas4P`, which is based on implicit methods.

## Heat equation

```{index} heat equation, method of lines
```

Now we return to the semidiscretization {eq}`heatMOL` of the heat equation, which was solved by Euler in @demo-methodlines-heatFE and backward Euler in @demo-methodlines-heatBE.

::::{prf:example} Stability regions and the heat equation
:label: demo-absstab-regions

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-absstab-regions-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-absstab-regions-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-absstab-regions-python
:::
````
`````
::::

The matrix $\mathbf{D}_{xx}$ occurring in {eq}`heatMOL` for semidiscretization of the periodic heat equation has eigenvalues that can be found explicitly. Assuming that $x\in[0,1)$ (with periodic boundary conditions), for which $h=1/m$, then the eigenvalues are (see @problem-absstab-d2eigs)

:::{math}
:label: D2eigs
\lambda_j =  -4m^2 \sin^2 \left( \frac{j\pi}{m} \right), \qquad j = 0,\ldots,m-1.   
:::

This result agrees with the observation in @demo-absstab-regions that the eigenvalues are real and negative. Furthermore, they lie within the interval $[-4m^2,0]$. In Euler time integration, this implies that $-4\tau m^2\ge -2$, or $\tau\ge 1/(2m^2)=O(m^{-2})$. For backward Euler, there is no time step restriction, and we say that backward Euler is unconditionally stable for this problem.

In summary, three things happen as $h\to 0$: 

1. The spatial discretization becomes more accurate like $O(h^2)$.
2. The size of the matrix increases like $O(h^{-1})$.
3. If we use an explicit time stepping method, then absolute stability requires $O(h^{-2})$ steps. 

The last restriction becomes rather burdensome as $h\to 0$, i.e., as we improve the spatial discretization, which is why implicit methods are preferred for diffusion. While any convergent IVP solver will get the right solution as $\tau\to 0$, the results are exponentially large nonsense until $\tau$ is small enough to satisfy absolute stability.

## Exercises

``````{exercise}
:label: problem-absdiff-diagonalize

✍ Use an eigenvalue decomposition to write the system

$$
\mathbf{u}'(t) =
\begin{bmatrix}
0 & 4 \\
-4 & 0
\end{bmatrix} \mathbf{u}(t)
$$

as an equivalent diagonal system. 
``````

``````{exercise}
:label: problem-absdiff-bounded

✍ For each system, state whether its solutions are bounded as $t\to \infty$.

**(a)** $\mathbf{u}'(t) =
\displaystyle \begin{bmatrix}
1 & 3 \\
3 & 1
\end{bmatrix} \mathbf{u}(t)$

**(b)** $\mathbf{u}'(t) =
\displaystyle \begin{bmatrix}
-1 & 3 \\
-3 & -1
\end{bmatrix} \mathbf{u}(t)$

**(c)** $\mathbf{u}'(t) =
\displaystyle \begin{bmatrix}
0 & 4 \\
-4 & 0
\end{bmatrix} \mathbf{u}(t)$
``````

``````{exercise}
:label: problem-absdiff-system1

✍ Using {numref}`figure-stabreg_ab_am` and {numref}`figure-stabreg_bd_rk`, estimate the time step restriction (if any) for the system

$$
\mathbf{u}'(t) =
\begin{bmatrix}
-4 & 0 & 0 \\
0 & -2 & 0 \\
0 & 0 & -0.5
\end{bmatrix} \mathbf{u}(t)
$$

for the following IVP methods:

**(a)** RK4 $\qquad$
**(b)** AM4 $\qquad$
**(c)** AB2
``````

``````{exercise}
:label: problem-absdiff-system2

✍ Using {numref}`figure-stabreg_ab_am` and {numref}`figure-stabreg_bd_rk`, find the time step restriction (if any) for the system

$$
\mathbf{u}'(t) =
\begin{bmatrix}
-1 & 0 & 0 \\
0 & 0 & 4 \\
0 & -4 & 0
\end{bmatrix} \mathbf{u}(t)
$$

for the following IVP methods:

**(a)** RK4 $\qquad$
**(b)** AM4 $\qquad$
**(c)** AB3
``````

``````{exercise}
:label: problem-absdiff-imaginary

✍ Of the following methods, which would be unsuitable for a problem having eigenvalues on the imaginary axis?  Justify your answer(s).

**(a)** AM2 $\qquad$
**(b)** AB2 $\qquad$
**(c)** RK2 $\qquad$
**(d)** RK3
``````

``````{exercise}
:label: problem-absdiff-negativereal

✍ Of the following methods, which would have a time step restriction for a problem with real, negative eigenvalues?  Justify your answer(s).

**(a)** AM2 $\qquad$
**(b)** AM4 $\qquad$
**(c)** BD4 $\qquad$
**(d)** RK4

``````

``````{exercise}
:label: problem-absstab-d2eigs
✍ Let $\mathbf{D}_{xx}$ be $m\times m$ and given by {eq}`heatFD22` for periodic end conditions. For any integer $k \in \{0,\ldots,m-1\}$, define $\omega = \exp(2ik\pi/m)$, and let $\mathbf{v}$ be the vector whose components are $v_j = \omega^j$ for $j=0,\ldots,m-1$.

**(a)** Show that $\omega^m = 1$. 

**(b)** Let $\mathbf{v}' = \mathbf{D}_x \mathbf{v}$. Show that for $j=1,\ldots,m-2$,

$$
v_j' = \frac{1}{h^2} \omega^{j} \left( \omega - 2 + \omega^{-1} \right).
$$

**(c)** Show that the result of part (b) holds for $j=0$ and $j=m-1$ as well.

**(d)** Explain why the above results prove that $\mathbf{v}$ is an eigenvector of $\mathbf{D}_x$ with associated eigenvalue

```{math}
:label: eq-d2eigs
\lambda =  -4 m^2 \sin^2\left( \frac{k\pi}{m} \right).
```
``````

``````{exercise}
✍ **(a)** Derive an algebraic inequality equivalent to absolute stability for the AM2 (trapezoid) formula.

✍ **(b)** Argue that the inequality in part (a) is equivalent to the restriction $\operatorname{Re} \zeta\le 0$. (Hint: Complex magnitude is equivalent to distance in the plane.)

``````
