---
numbering:
  enumerator: 6.5.%s
---
(section-ivp-adaptive)=
# Adaptive Runge–Kutta

```{index} Runge–Kutta method
```

The derivation and analysis of methods for initial-value problems usually assumes a fixed step size $h$. While the error behavior $O(h^p)$ is guaranteed by {numref}`Theorem %s <theorem-euler-onestepGTE>` as $h\rightarrow 0$, this bound comes with an unknowable constant, and it is not very useful as a guide to the numerical value of the error at any particular value of $h$. Furthermore, as we saw in {numref}`section-localapprox-adaptive` for numerical integration, in many problems a fixed step size is far from the most efficient strategy.

In response we will employ the basic strategy of {numref}`section-localapprox-adaptive`: estimate the error and adapt the step size in order to reach an accuracy goal. Unlike the integration problem, though, the "integrand" of an IVP is dependent on the solution itself, so the details differ greatly.
## Step size prediction

Suppose that, starting from a given value $u_i$ and using a step size $h$, we run one step of two different RK methods simultaneously: one method with order $p$, producing $u_{i+1}$, and the other method with order $p+1$, producing $\tilde{u}_{i+1}$. In most circumstances, we can expect that $\tilde{\mathbf{u}}_{i+1}$ is a much better approximation to the solution than $\mathbf{u}_{i+1}$ is. So it seems reasonable to use 

$$E_i(h)=|\tilde{\mathbf{u}}_{i+1} - \mathbf{u}_{i+1}|$$

as an estimate of the actual local error made by the $p$th-order method. For a vector IVP, we would use a norm rather than an absolute value.

<!-- If the goal is to keep global error less than some predetermined value, we could decide to accept the new solution value if $E_i$ small enough, and otherwise reject it.[^extrap]

[^extrap]: Even though the estimate $E_i$ is meant to go with the *less* accurate proposed value $\mathbf{u}_{i+1}$, it's hard to resist the temptation to keep the more accurate value instead, and this is common in practice. -->

Now we ask: looking back, what step size *should* we have taken to meet an error target of size $\epsilon$? Let's speculate, given the behavior of local truncation error as $h\rightarrow 0$, that $E_i(h)\approx C h^{p+1}$ for an unknown constant $C$. If we had used a step size $q h$ for some $q>0$, then trivially, we would expect 

$$E_i(qh)\approx C q^{p+1}h^{p+1}.$$ 

Our best guess for $q$ would therefore be to set $E_i(qh)\approx \epsilon$, or

```{math}
:label: adaptRKlocal
  q \approx \left(\frac{\epsilon}{E_i}\right)^{1/(p+1)}.
```

Perhaps, though, we should aim to control the contribution to *global* error, which is closer to $E_i(qh)/(q h)$. Then we end up with

```{math}
:label: adaptRKglobal
  q \le \left(\frac{\epsilon}{E_i}\right)^{1/p}.
```

Experts have different recommendations about whether to use {eq}`adaptRKlocal` or {eq}`adaptRKglobal`. Even though {eq}`adaptRKglobal` appears to be more in keeping with our assumptions about global errors, modern practice seems to favor {eq}`adaptRKlocal`.

```{index} adaptivity; in IVP solver
```

We now have an outline of an algorithm.

(algorithm-adaptive-adapt)=
::::{prf:algorithm} Adaptive step size for an IVP
Given a solution estimate $u_i$ at $t=t_i$, and a step size $h$, do the following:
1. Produce estimates ${u}_{i+1}$ and $\tilde{u}_{i+1}$, and estimate the error.
2. If the error is small enough, adopt $\tilde{u}_{i+1}$ as the solution value at $t=t_i+h$, then increment $i$.
3. Replace $h$ by $q h$, with $q$ given by {eq}`adaptRKlocal` or {eq}`adaptRKglobal`.
4. Repeat until $t=b$.
::::

Many details remain unspecified at this point, but we first address step 1.
## Embedded formulas

Suppose, for example, we choose to use a  pair of second- and third-order RK methods to get the $\mathbf{u}_{i+1}$ and $\tilde{\mathbf{u}}_{i+1}$ needed in {numref}`Algorithm {number} <algorithm-adaptive-adapt>`. Then we seem to need at least $2+3=5$ evaluations of $f(t,y)$ for each attempted time step. This is more than double the computational work needed by the second-order method without adaptivity. 

Fortunately, the marginal cost of adaptivity can be substantially reduced by using **embedded Runge–Kutta** formulas. Embedded RK formulas are a pair of RK methods whose stages share the same internal $f$ evaluations, combining them differently in order to get estimates of two different orders of accuracy.

A good example of an embedded method is the **Bogacki–Shampine** (BS23) formula, given by the table

```{math}
:label: bs23
\begin{array}{r|cccc}
0                  & \rule{0pt}{2.75ex} &                    &                    &                    \\
\frac{1}{2}        & \frac{1}{2}        & \rule{0pt}{2.75ex} &                    &                    \\
\frac{3}{4}        & 0                  & \frac{3}{4}        & \rule{0pt}{2.75ex} &                    \\
1                 & \frac{2}{9}        & \frac{1}{3}        & \frac{4}{9}        & \rule{0pt}{2.75ex} \\[2pt] \hline
\rule{0pt}{2.75ex} & \frac{2}{9}        & \frac{1}{3}        & \frac{4}{9}        & 0                  \\[2pt] \hline
\rule{0pt}{2.75ex} & \frac{7}{24}       & \frac{1}{4}        & \frac{1}{3}        & \frac{1}{8}
\end{array}
```

The top part of the table describes four stages in the usual RK fashion. The last two rows describe how to construct a third-order estimate $\tilde{\mathbf{u}}_{i+1}$ and a second-order estimate $\mathbf{u}_{i+1}$ by taking different combinations of those stages.

## Implementation

Our implementation of an embedded second/third-order (RK23) code is given in {numref}`Function {number} <function-rk23>`. 

(function-rk23)=
``````{prf:algorithm} rk23
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #function-rk23-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-rk23-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #function-rk23-python
:::
````
`````
``````

(demo-adapt-basic)=
::::{prf:example}
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #demo-adapt-basic-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-adapt-basic-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-adapt-basic-python
:::
```` 
`````
::::

(demo-adapt-sing)=
::::{prf:example}
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #demo-adapt-sing`-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-adapt-sing`-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-adapt-sing`-python
:::
```` 
`````
::::



```{index} stiff differential equation
```

Often the adaptively chosen steps clearly correspond to identifiable features of the solution. However, there are so-called *stiff problems* in which the time steps seem unreasonably small in relation to the observable behavior of the solution. These problems benefit from a particular type of solver that is considered in {numref}`section-ivp-implicit`.

## Exercises

1. ⌨ Using {numref}`Function {number} <function-rk23>` with an error tolerance of $10^{-8}$, solve $y'' +(1+y')^3 y = 0$ over $ 0 \le t \le 4 \pi$ with the indicated initial conditions. Plot $y(t)$ and $y'(t)$ as functions of $t$ and separately plot the time step size as a function of $t$.

    **(a)** $y(0) = 0.1, \quad y'(0) = 0$

    **(b)** $y(0) = 0.5, \quad y'(0) = 0$

    **(c)** $y(0) = 0.75, \quad y'(0) = 0$

    **(d)** $y(0) = 0.95, \quad y'(0) = 0$

2. ⌨ Solve the FitzHugh–Nagumo system from [Exercise 4.3.6](problem-systems-fitznag) for $I=0.05740$ using {numref}`Function {number} <function-rk23>` with error tolerance $10^{-2}$, $10^{-3}$, and $10^{-4}$. (This illustrates that the error tolerance is a target, not a guarantee!)

3. ✍ Derive Equation {eq}`adaptRKglobal` using the stated assumption about controlling global rather than local error.

4. ⌨ Solve the problem $u'=100u^2-u^3$, $u(0)=0.0002$, $0\le t \le 100$, and make plots that show both the solution and the time steps taken. The solution makes a quick transition between two nearly constant states. Does the step size selection behave the same in both states?

