---
numbering:
  enumerator: 4.1.%s
---
(section-nonlineqn-rootproblem)=
# The rootfinding problem

For the time being we will focus on the **rootfinding problem** for single functions of one variable.

```{index} ! rootfinding problem
```

```{prf:definition} Rootfinding problem
Given a continuous scalar function $f$ of a scalar variable, find a real number $r$, called a **root**, such that $f(r)=0$.
```

The formulation $f(x)=0$ is general enough to solve any equation. If we are trying to solve an equation $g(x)=h(x)$, we can define $f=g-h$ and find a root of $f$.

Unlike the linear problems of the earlier chapters, the usual situation here is that the root cannot be produced in a finite number of operations, even in exact arithmetic. Instead, we seek a sequence of approximations that formally converge to the root, stopping when some member of the sequence seems to be good enough, in a sense we will clarify later.

```{index} Bessel function
```

::::{prf:example} The rootfinding problem for Bessel functions
:label: demo-rootproblem-bessel
In the theory of vibrations of a circular drum, the displacement of the drumhead can be expressed in terms of pure harmonic modes, 

$$J_m(\omega_{k,m} r) \cos(m\theta) \cos(c \omega_{k,m} t),$$

where $(r,\theta)$ are polar coordinates, $0\le r\le 1$, $t$ is time, $m$ is a positive integer, $c$ is a material parameter, and $J_m$ is a _Bessel function of the first kind_. The quantity $\omega_{k,m}$ is a resonant frequency and is a positive root of the equation  

$$J_m(\omega_{k,m}) = 0,$$ 

which states that the drumhead is clamped around the rim. Bessel functions often appear in physical problems featuring radial symmetry, and tabulating approximations to the zeros of Bessel functions occupied numerous mathematician-hours before computers were on the scene.

`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #demo-rootproblem-bessel-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-rootproblem-bessel-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-rootproblem-bessel-python
:::
```` 
`````
::::

## Conditioning, error, and residual

In the rootfinding problem, the data is a continuous function $f$ and the result is a root. (This overrides our Chapter 1 notation of $f$ as the map from data to result.) How does the result change in response to perturbations in $f$? We will compute an absolute condition number rather than a relative one.

You might wonder about the relevance of perturbing a function as data to a problem. If nothing else, the values of $f$ will be represented in floating point and thus subject to rounding error. Furthermore, in many applications, $f$ might not be a simple formula but the result of a computation that uses an inexact algorithm. While there are infinitely many possible perturbations to a function, a constant perturbation is enough to get the main idea.

We assume $f$ has at least one continuous derivative near a particular root $r$. Suppose that $f$ is perturbed to $\tilde{f}(x) = f(x) + \epsilon$. As a result, the root (if it still exists) will be perturbed to $\tilde{r} = r + \delta$ such that $\tilde{f}(\tilde{r})=0$. We now compute an absolute condition number $\kappa_r$, which is the ratio $\left | \frac{\delta}{\epsilon} \right|$ as $\epsilon\to 0$.

Using Taylor's theorem,

```{math}
  0 = f(r+\delta) + \epsilon \approx f(r) + f'(r) \delta + \epsilon.
```

Since $r$ is a root, we have $f(r)=0$. This lets us relate $\delta$ to $\epsilon$, and their ratio is the condition number.

```{index} condition number; of rootfinding
```

````{prf:theorem} Condition number of rootfinding
If $f$ is differentiable at a root $r$, then the absolute condition number of $r$ with respect to constant changes in $f$ is

```{math}
:label: rootcondnum
  \kappa_r = \bigl| f'(r) \bigr|^{-1}.
```

We say $\kappa_r = \infty$ if $f'(r)=0$.
````

Equivalently, {eq}`rootcondnum` is just the magnitude of the derivative of the inverse function $f^{-1}$ at zero. 


::::{prf:example} Condition number of a rootfinding problem
:label: demo-roots-cond
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #demo-roots-cond-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-roots-cond-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-roots-cond-python
:::
```` 
`````
::::

```{index} ! residual; of rootfinding
```

We must accept that when $|f'|$ is small at the root, it may not be possible to get a small error in a computed root estimate. As always, the error is not a quantity we can compute without knowing the exact answer. There is something else we can measure, though.

::::{prf:definition} Rootfinding residual
If $\tilde{r}$ approximates a root $r$ of function $f$, then the **residual** at $\tilde{r}$ is $f(\tilde{r})$.
::::

It stands to reason that a small residual might be associated with a small error. To quantify the relationship, let $\tilde{r}$ approximate root $r$, and define the new function $g(x)=f(x)-f(\tilde{r})$. Trivially, $g(\tilde{r})=0$, meaning that $\tilde{r}$ is a true root of $g$. Since the difference $g(x)-f(x)$ is the residual value $f(\tilde{r})$, the residual is the distance to an exactly solved rootfinding problem. 

```{index} backward error
```

```{prf:observation}
The backward error in a root estimate is equal to the residual.
```

In general, it is not realistic to expect a small error in a root approximation if the condition number {eq}`rootcondnum` is large. However, we can gauge the backward error from the residual.

## Multiple roots

```{index} ! roots; multiplicity of
```

The condition number {eq}`rootcondnum` naturally leads to the question of what happens if $f'(r)=0$ at a root $r$. The following definition agrees with and extends the notion of algebraic multiplicity in a polynomial to roots of more general differentiable functions. 

::::{prf:definition} Multiplicity of a root
:label: definition-rootproblem-simple
If $f(r)=f'(r)=\cdots=f^{(m-1)}(r)=0$, but $f^{(m)}(r)\neq 0$, then we say $f$ has a root of **multiplicity** $m$ at $r$. In particular, if $f(r)=0$ and $f'(r)\neq 0$, then $m=1$ and we call $r$ a **simple root**.
::::

Another useful characterization of multiplicity $m$ is that $f(x)=(x-r)^m q(x)$ for a differentiable $q$ with $q(r)\neq 0$. 

When $r$ is a nonsimple root, the condition number {eq}`rootcondnum` is effectively infinite.[^infcond] However, even if $r$ is simple, we should expect difficulty in rootfinding if the condition number is very large. This occurs when $|f'(r)|$ is very small, which means that the quotient $q$ satisfies $q(r)\approx 0$ and another root of $f$ is very close to $r$. We made the same observation about polynomial roots all the way back in @demo-stability-roots. 

[^infcond]: Based on our definitions, this means that the relative change to the root when $f$ is changed by a perturbation of size $\epsilon$ is not $O(\epsilon)$ as $\epsilon\to 0$.

## Exercises

⌨  For each equation and given interval, do the following steps.
  
**(a)** Rewrite the equation into the standard form for rootfinding, $f(x) = 0$. Make a plot of $f$ over the given interval and determine how many roots lie in the interval. 
  
**(b)**  Use `nlsolve` to find each root, as shown in @demo-rootproblem-bessel.

**(c)** Compute the condition number of each root found in part (b). 

``````{exercise}
:label: problem-rootproblem-basic1
$x^2=e^{-x}$, over $[-2,2]$
``````

``````{exercise}
:label: problem-rootproblem-basic2
$2x = \tan x$, over $[-0.2,1.4]$
``````

``````{exercise}
:label: problem-rootproblem-basic3
$e^{x+1}=2+x$, over $[-2,2]$

``````

---

``````{exercise}
:label: problem-rootproblem-annuity
⌨ A basic safe type of investment is an annuity: one makes monthly deposits of size $P$ for $n$ months at a fixed annual interest rate $r$, and at maturity collects the amount

$$
\frac{12 P}{r} \left( \Bigl(1+\frac{r}{12}\Bigr)^n - 1\right).
$$

Say you want to create an annuity for a term of 300 months and final value of \$1,000,000. Using `nlsolve`, make a table of the interest rate you will need to get for each of the different contribution values $P=500,550,\ldots,1000$.  
``````

``````{exercise}
:label: problem-rootproblem-kepler
⌨ The most easily observed properties of the orbit of a celestial body around the sun are the period $\tau$ and the elliptical eccentricity $\epsilon$. (A circle has $\epsilon=0$.) From these, it is possible to find at any time $t$ the angle $\theta(t)$ made between the body's position and the major axis of the ellipse. This is done through

```{math}
:label: kepler1
\tan \frac{\theta}{2} = \sqrt{\frac{1+\epsilon}{1-\epsilon}}\,
\tan \frac{\psi}{2},
```

where the eccentric anomaly $\psi(t)$ satisfies Kepler's equation:

```{math}
:label: kepler2
\psi - \epsilon \sin \psi - \frac{2\pi t}{\tau} = 0.
```

Equation {eq}`kepler2` must be solved numerically to find $\psi(t)$, and then {eq}`kepler1` can be solved analytically to find $\theta(t)$.

The asteroid Eros has $\tau=1.7610$ years and $\epsilon=0.2230$. Using `nlsolve` for {eq}`kepler2`, make a plot of $\theta(t)$ for 100 values of $t$ between $0$ and $\tau$, which is one full orbit. (Note: Use `mod(θ,2π)` to put the angle between 0 and $2\pi$ if you want the result to be a continuous function.)
``````

``````{exercise}
:label: problem-rootproblem-lambertW
⌨  Lambert's $W$ function is defined as the inverse of $x e^x$. That is, $y=W(x)$ if and only if $x=ye^y$. Write a function `lambertW` that computes $W$ using `nlsolve`. Make a plot of $W(x)$ for $0\le x \le 4$.  
``````

``````{exercise}
:label: problem-rootproblem-multiplicity
✍ For each function, find the multiplicity of the given root. If it is a simple root, find its absolute condition number.

**(a)** $f(x) = x^3-2x^2+x-2$, root $r=2$

**(b)** $f(x) = (\cos x  + 1)^2$, root $r=\pi$

**(c)** $f(x) = \frac{\sin^2 x}{x}$, root $r=0$ (define $f(0) =0$)

**(d)** $f(x) =(x-1)\log(x)$, root $r=1$
``````
