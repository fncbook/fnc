---
numbering:
  enumerator: 6.5.%s
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
---
```{code-cell}
:tags: [remove-cell]
clear all
format short
set(0, 'defaultaxesfontsize', 12)
set(0, 'defaultlinelinewidth', 1.5)
set(0, 'defaultFunctionLinelinewidth', 1.5)
set(0, 'defaultscattermarkerfacecolor', 'flat')
gcf;
set(gcf, 'Position', [0 0 600 350])
addpath FNC-matlab
```

(section-ivp-adaptive)=

# Adaptive Runge–Kutta

```{index} Runge–Kutta method
```

The derivation and analysis of methods for initial-value problems usually assumes a fixed step size $h$. While the error behavior $O(h^p)$ is guaranteed by @theorem-euler-onestepGTE as $h\rightarrow 0$, this bound comes with an unknowable constant, and it is not very useful as a guide to the numerical value of the error at any particular value of $h$. Furthermore, as we saw in {numref}`section-localapprox-adaptive` for numerical integration, in many problems a fixed step size is far from the most efficient strategy.

In response, we will employ the basic strategy of {numref}`section-localapprox-adaptive`: estimate the error and adapt the step size in order to reach an accuracy goal. Unlike the integration problem, though, the "integrand" of an IVP is dependent on the solution itself, so the details differ greatly.

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
  q \approx \left(\frac{\epsilon}{E_i}\right)^{1/p}.
```

Experts have different recommendations about whether to use {eq}`adaptRKlocal` or {eq}`adaptRKglobal`. Even though {eq}`adaptRKglobal` appears to be more in keeping with our assumptions about global errors, modern practice seems to favor {eq}`adaptRKlocal`.

```{index} adaptivity; in IVP solver
```

We now have an outline of an algorithm.

::::{prf:algorithm} Adaptive step size for an IVP
:label: algorithm-adaptive-adapt
Given a solution estimate $u_i$ at $t=t_i$, and a step size $h$, do the following:

1. Produce estimates ${u}_{i+1}$ and $\tilde{u}_{i+1}$, and estimate the error.
2. If the error is small enough, adopt $\tilde{u}_{i+1}$ as the solution value at $t=t_i+h$, then increment $i$.
3. Replace $h$ by $q h$, with $q$ given by {eq}`adaptRKlocal` or {eq}`adaptRKglobal`.
4. Repeat until $t=b$.
::::

Many details remain unspecified at this point, but we first address step 1.

## Embedded formulas

Suppose, for example, we choose to use a pair of second- and third-order RK methods to get the $\mathbf{u}_{i+1}$ and $\tilde{\mathbf{u}}_{i+1}$ needed in {numref}`Algorithm {number} <algorithm-adaptive-adapt>`. Then we seem to need at least $2+3=5$ evaluations of $f(t,y)$ for each attempted time step. This is more than double the computational work needed by the second-order method without adaptivity.

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

``````{prf:algorithm} rk23
:label: function-rk23
```{literalinclude} ../FNC_matlab/rk23.m
:language: matlab
:linenos: true
```
::::{admonition} About the code
:class: dropdown
The check `t(i) + h == t(i)`on line 24 is to detect when $h$ has become so small that it no longer changes the floating-point value of $t_i$. This may be a sign that the underlying exact solution has a singularity near $t=t_i$, but in any case, the solver must halt by using a `break` statement to exit the loop.

On line 36, we use a combination of absolute and relative tolerances to judge the acceptability of a solution value, as in {eq}`absreltolerance`. In lines 47--49 we underestimate the step factor $q$ a bit and prevent a huge increase in the step size, since a rejected step is expensive, and then we make sure that our final step doesn't take us past the end of the domain.

Finally, line 43 exploits a subtle property of the BS23 formula called *first same as last* (FSAL).
While {eq}`bs23` calls for four stages to find the paired second- and third-order estimates, the final stage computed in stepping from $t_i$ to $t_{i+1}$ is identical to the first stage needed to step from $t_{i+1}$ to $t_{i+2}$. By repurposing `s4` as `s1` for the next pass, one of the stage evaluations comes for free, and only three evaluations of $f$ are needed per successful step.
::::

``````

::::{prf:example} Adaptive step size
:label: demo-adapt-basic

Let's run adaptive RK on  $u'=e^{t-u\sin u}$.

```{code-cell}
f = @(t, u) exp(t - u * sin(u));
u0 = 0;
a = 0;  b = 5;

[t, u] = rk23(f, [a, b], u0, 1e-5);
clf, plot(t, u)
xlabel("t");  ylabel("u(t)")
title(("Adaptive IVP solution"));
```

The solution makes a very abrupt change near $t=2.4$. The resulting time steps vary over three orders of magnitude.

```{code-cell}
Delta_t = diff(t);
semilogy(t(1:end-1), Delta_t) 
xlabel("t");  ylabel("step size")
title(("Adaptive step sizes"));
```

If we had to run with a uniform step size to get this accuracy, it would be

```{code-cell}
fprintf("minimum step size = %.2e", min(Delta_t))
```

On the other hand, the average step size that was actually taken was

```{code-cell}
fprintf("average step size = %.2e", mean(Delta_t))
```

We took fewer steps by a factor of almost 1000! Even accounting for the extra stage per step and the occasional rejected step, the savings are clear.


::::

::::{prf:example} Adaptive step size near a singularity
:label: demo-adapt-sing

In @demo-basics-sing we saw an IVP that appears to blow up in a finite amount of time. Because the solution increases so rapidly as it approaches the blowup, adaptive stepping is required even to get close.

```{code-cell}
f = @(t, u) (t + u)^2;
u0 = 1;
[t, u] = rk23(f, [0, 1], u0, 1e-5);
```

In fact, the failure of the adaptivity gives a decent idea of when the singularity occurs.

```{code-cell}
clf, semilogy(t, u)
xlabel("t");  ylabel("u(t)")
title("Adaptive solution near a singularity")

tf = t(end);
xline(tf, "linestyle", "--")
text(tf, 1e5, sprintf(" t = %.6f ", tf))
```

::::

```{index} stiff differential equation
```

Often the adaptively chosen steps clearly correspond to identifiable features of the solution. However, there are so-called *stiff problems* in which the time steps seem unreasonably small in relation to the observable behavior of the solution. These problems benefit from a particular type of solver that is considered in {numref}`section-ivp-implicit`.

## Exercises

``````{exercise}
:label: problem-adaptiverk-rk23
⌨ Using {numref}`Function {number} <function-rk23>` with an error tolerance of $10^{-8}$, solve $y'' +(1+y')^3 y = 0$ over $ 0 \le t \le 4 \pi$ with the indicated initial conditions. Plot $y(t)$ and $y'(t)$ as functions of $t$ and separately plot the time step size as a function of $t$.

**(a)** $y(0) = 0.1, \quad y'(0) = 0$

**(b)** $y(0) = 0.5, \quad y'(0) = 0$

**(c)** $y(0) = 0.75, \quad y'(0) = 0$

**(d)** $y(0) = 0.95, \quad y'(0) = 0$
``````

``````{exercise}
:label: problem-adaptiverk-fh
⌨ Solve the FitzHugh–Nagumo system from @problem-systems-fitznag for $I=0.05740$ using {numref}`Function {number} <function-rk23>` with error tolerance $10^{-2}$, $10^{-3}$, and $10^{-4}$. (This illustrates that the error tolerance is a target, not a guarantee!)
``````

``````{exercise}
:label: problem-adaptiverk-global
✍ Derive Equation {eq}`adaptRKglobal` using the stated assumption about controlling global rather than local error.
``````

``````{exercise}
:label: problem-adaptiverk-flame
⌨ Solve the problem $u'=100u^2-u^3$, $u(0)=0.0002$ for $0\le t \le 100$ using {numref}`Function {number} <function-rk23>` with error tolerance $10^{-6}$. In a 2-by-1 plot layout, plot the solution above and (on a semilog scale) the time step size below. The solution makes a quick transition between two nearly constant states. Does the step size behave the same during both of those states?
``````
