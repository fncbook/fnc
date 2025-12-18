---
numbering:
  enumerator: 6.1.%s
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---
```{code-cell}
:tags: [remove-cell]
from numpy import *
from scipy import linalg
from scipy.linalg import norm
from matplotlib.pyplot import *
from prettytable import PrettyTable
from timeit import default_timer as timer
import sys
sys.path.append('fncbook/')
import fncbook as FNC

# This (optional) block is for improving the display of plots.
# from IPython.display import set_matplotlib_formats
# set_matplotlib_formats("svg","pdf")
# %config InlineBackend.figure_format = 'svg'
rcParams["figure.figsize"] = [7, 4]
rcParams["lines.linewidth"] = 2
rcParams["lines.markersize"] = 4
rcParams['animation.html'] = "jshtml"  # or try "html5"
```

(section-ivp-basics)=

# Basics of IVPs

```{index} ! initial-value problem
```

::::{prf:definition} Initial-value problem: scalar
:label: definition-scalarivp
A scalar first-order {term}`initial-value problem` (IVP) is
  
```{math}
:label: IVP
\begin{split}
   u'(t) &= f(t,u(t)), \qquad a \le t \le b,  \\
  u(a) &=u_0.
\end{split}
```

We call $t$ the **independent variable** and $u$ the **dependent variable**. If $u'=f(t,u)=g(t)+u h(t)$, the differential equation is **linear**; otherwise, it is **nonlinear**.

A **solution** of an initial-value problem is a function $u(t)$ that makes both $u'(t)=f\bigl(t,u(t)\bigr)$ and $u(a)=u_0$ true equations.
::::

When $t$ is meant to be time, sometimes we write $\dot{u}$ (read "u-dot") instead of $u'$.

::::{prf:example}
  Suppose $u(t)$ is the size of a population at time $t$. We idealize by allowing $u$ to take any real (not just integer) value. If we assume a constant per capita birth rate (births per unit population per unit time), then
  
```{math}
\frac{d u}{d t} = k u, \qquad u(0)=u_0
```

for some $k>0$. The solution of this linear equation is $u(t)=e^{kt}u_0$, which is exponential growth.

A more realistic model would cap the growth due to finite resources. Suppose the death rate is proportional to the size of the population, indicating competition. Then

```{index} ! logistic equation
```

```{math}
:label: logistic
  \frac{d u}{d t} = ku - ru^2, \qquad u(0)=u_0.
```

This is the **logistic equation**. Although crude, it is still useful in population models.  The solution relevant for population models has the form
  
```{math}
  u(t) = \frac{k/r}{ 1 + \left( \frac{k}{r u_0} - 1 \right) e^{-k t} }.
```

For $k,r,u_0>0$, the solution smoothly varies from the initial population $u_0$ to a finite population, equal to $k/r$, that has been limited by competition.
::::

Linear problems can be solved in terms of integrals. Defining the *integrating factor* $\rho(t) = \exp\bigl[\int -h(t)\, dt \bigr]$, the solution is derived from

```{math}
  \rho(t) u(t) = u_0 + \int_a^t \rho(s) g(s) \, ds.
```

In many cases, however, the necessary integrals cannot be done in closed form. Some nonlinear ODEs, such as separable equations, may also be solvable with a short formula, perhaps with difficult integrations. Most often, though, there is no analytic formula available for the solution.

An ODE may have higher derivatives of the unknown solution present. For example, a **second-order ordinary differential equation** is often given in the form $u''(t)=f\bigl(t,u,u'\bigr)$. A second-order IVP requires two conditions at the initial time in order to specify a solution completely. As we will see in {numref}`section-ivp-systems`, we are always able to reformulate higher-order IVPs in a first-order form, so we will deal with first-order problems exclusively.

## Numerical solutions

::::{prf:example} Solving an IVP
:label: demo-basics-first

Let's use `solve_ivp` from `scipy.integrate` to define and solve the initial-value problem 

```{math}
:numbered: false
u'=\sin[(u+t)^2], \quad t \in [0,4], \quad u(0)=-1.
```

To create an initial-value problem for $u(t)$, we must supply a function that computes $u'$, an initial value for $u$, and the endpoints of the interval for $t$.

```{index} ! Python; solve_ivp
```

```{code-cell}
f = lambda t, u: sin((t + u) ** 2)
tspan = [0.0, 4.0]
u0 = array([-1.0])
```

Note above that even though this is a problem for a scalar function $u(t)$, we had to set the initial condition as a one-component vector.

```{code-cell}
from scipy.integrate import solve_ivp
sol = solve_ivp(f, tspan, u0)
```

The resulting solution object has fields `t` and `y` that contain the values of the independent and dependent variables, respectively; those field names are the same regardless of what we use in our own codes.

```{code-cell}
print("t shape:", sol.t.shape)
print("u shape:", sol.y.shape)
plot(sol.t, sol.y[0, :], "-o")
xlabel("$t$"), ylabel("$u(t)$")
title(("Solution of $u' = sin((t+u)^2)$"));
```

You can see above that the solution was not computed at enough points to make a smooth graph. There is a way to request output at times of your choosing.

```{code-cell}
sol = solve_ivp(f, tspan, u0, t_eval=linspace(0, 4, 200))
plot(sol.t, sol.y[0, :], "-")
xlabel("$t$"), ylabel("$u(t)$")
title(("Solution of $u' = sin((t+u)^2)$"));
```

Another option is to enable interpolation to evaluate the solution anywhere after the fact:

```{code-cell}
sol = solve_ivp(f, tspan, u0, dense_output=True)
for t in linspace(0, 4, 6):
    print(f"u({t:.2f}) = {sol.sol(t)[0]:.4f}")
```

::::

## Existence and uniqueness

There are simple IVPs that do not have solutions at all possible times.

::::{prf:example} Finite-time singularity
:label: demo-basics-sing


The equation $u'=(u+t)^2$ gives us some trouble.
```{tip}
:class: dropdown
It's a good idea to check `sol.success` after calling `solve_ivp`. If it's `False`, the solution may not be reliable. 
```

```{code-cell}
f = lambda t, u: (t + u) ** 2
sol = solve_ivp(f, [0.0, 1.0], [1.0])
if not sol.success:
    print(sol.message)
```

The warning message we received can mean that there is a bug in the formulation of the problem. But if everything has been done correctly, it suggests that the solution may not exist past the indicated time. This is a possibility in nonlinear ODEs.

```{code-cell}
semilogy(sol.t, sol.y[0, :])
xlabel("$t$")
ylabel("$u(t)$")
title(("Blowup in finite time"));
```

::::

We can also produce an IVP that has more than one solution.

````{prf:example}
  The functions $u(t)=u^2$ and $u(t)\equiv 0$ both satisfy the differential equation $u'=2\sqrt{u}$ and the initial condition $u(0)=0$. Thus the corresponding IVP has more than one solution.
````

The following standard theorem gives us a condition that is easy to check and guarantees that a unique solution exists. But it is not the most general possible such condition, so there are problems with a unique solution that it cannot detect. We state the theorem without proof.

````{prf:theorem} Existence and uniqueness
:label: theorem-existunique

If the derivative $\frac{\partial f}{\partial u}$ exists and $\left|\frac{\partial f}{\partial u}\right|$ is bounded by a constant $L$ for all $a\le t \le b$ and all $u$, then the initial-value problem {eq}`IVP` has a unique solution for $t\in [a,b]$.
````

## Conditioning of first-order IVPs

```{index} condition number; of initial-value problems
```

In a numerical context we have to be concerned about the conditioning of the IVP. There are two key items in {eq}`IVP` that we might consider to be the data of the initial-value ODE problem: the function $f(t,u)$, and the initial value $u_0$. It's easier to discuss perturbations to numbers than to functions, so we will focus on the effect of $u_0$ on the solution, using the following theorem that we give without proof. Happily, its conditions are identical to those in @theorem-existunique.

````{prf:theorem} Dependence on initial value
:label: theorem-depIC
If the derivative $\frac{\partial f}{\partial u}$ exists and $\left|\frac{\partial f}{\partial u}\right|$ is bounded by a constant $L$ for all $a\le t \le b$ and all $u$, then the solution $u(t;u_0+\delta)$ of $u'=f(t,u)$ with initial condition $u(0)=u_0+\delta$ satisfies
  
```{math}
:label: depIC
\left\|u(t;u_0+\delta)-u(t;u_0)\right\|_\infty \le |\delta| e^{L(b-a)}
```

for all sufficiently small $|\delta|$.
````

Numerical solutions of IVPs have errors, and those errors can be seen as perturbations to the solution. @theorem-depIC gives an upper bound of $e^{L(b-a)}$ on the infinity norm (i.e., pointwise) absolute condition number of the solution with respect to perturbations at an initial time. However, the upper bound may be a terrible overestimate of the actual sensitivity for a particular problem.

::::{prf:example} Conditioning of an IVP
:label: demo-basics-cond

Consider the ODEs $u'=u$ and $u'=-u$. In each case we compute $\partial f/\partial u = \pm 1$, so the condition number bound from @theorem-depIC is $e^{b-a}$ in both problems. However, they behave quite differently. In the case of exponential growth, $u'=u$, the bound is the actual condition number.

```{code-cell}
t = linspace(0, 3, 200)
u = array([exp(t) * u0 for u0 in [0.7, 1, 1.3]])
plot(t, u.T)
xlabel("$t$")
ylabel("$u(t)$")
title(("Exponential divergence of solutions"));
```

But with $u'=-u$, solutions actually get closer together with time.

```{code-cell}
t = linspace(0, 3, 200)
u = array([exp(-t) * u0 for u0 in [0.7, 1, 1.3]])
plot(t, u.T)
xlabel("$t$")
ylabel("$u(t)$")
title(("Exponential convergence of solutions"));
```

In this case the actual condition number is one, because the initial difference between solutions is the largest over all time. Hence, the exponentially growing upper bound $e^{b-a}$ is a gross overestimate.

::::

In general, solutions can diverge from, converge to, or oscillate around the original trajectory in response to perturbations. We won't fully consider these behaviors and their implications for numerical methods again until a later chapter.

## Exercises

``````{exercise}
:label: problem-basics-lipschitz
✍ For each IVP, determine whether the problem satisfies the conditions of @theorem-depIC. If so, determine the smallest possible value for $L$.

**(a)** $f(t,u) = 3 u,\; 0 \le t \le 1$

**(b)** $f(t,u) = -t \sin(u),\; 0 \le t \le 5$

**(c)** $f(t,u) = -(1+t^2) u^2,\; 1 \le t \le 3$

**(d)** $f(t,u) = \sqrt{u},\; 0 \le t \le 1$
``````

``````{exercise}
:label: problem-basics-usage
⌨ For each ODE in the preceding problem, assume that $u$ is initially equal to $1$ on the given interval. Solve the resulting IVP as in @demo-basics-first, and make a plot of the solution.
``````

``````{exercise}
:label: problem-basics-intfactor
✍ Use an integrating factor to find the solution of each problem in analytic form.

**(a)** $u' = -t u,\ 0 \le t \le 5,\ u(0) = 2$

**(b)** $u' - 3 u = e^{-2t},\ 0 \le t \le 1,\  u(0) = 5$
``````

``````{exercise}
:label: problem-basics-existence
✍ Consider the IVP $u'=u^2$, $u(0)=\alpha$.

**(a)** Does @theorem-existunique apply to this problem?

**(b)** Show that $u(t) = \dfrac{\alpha}{1-\alpha t}$ is a solution of the IVP.

**(c)** Does this solution necessarily exist for all $t\in[0,1]$?

``````

```{index} logistic equation
```

``````{exercise}
:label: problem-basics-logistic
⌨ Using the method in @demo-basics-first, compute solutions $x(t)$ to the [logistic equation](wiki:Logistic_function#Logistic_differential_equation) with harvesting,

```{math}
:numbered: false
x' = k (S-x)(x-M), \qquad 0\le t \le 10,
```

with $k=S=1$ and $M=0.25$, for the initial conditions $x(0)=0.9M$, $1.1M$, $1.5M$, $0.9S$, $1.1S$, $3S$. Show all the solutions together on one plot with $0\le x \le 3$. (Note: One of the solutions will throw a warning and fail to reach $t=10$, but you can plot it anyway.)
``````

``````{exercise}
:label: problem-basics-attractors
⌨ **(a)** Using the method in @demo-basics-first, solve the IVP $u'=u\cos(u) + \cos(4t)$, $0\le t \le 10$, $u(0) = u_0$ for $u_0 = -2,-1.5,-1,\ldots,1.5,2$. Plot all the solutions on a single graph. 

**(b)** All of the solutions in part (a) eventually settle into one of two periodic oscillations. To two digits of accuracy, find the value of $u_0$ in $(-1,1)$ at which the selected long-term solution changes. (This will take repeated trials, narrowing down the range for $u_0$ each time.)

``````

``````{exercise}
:label: problem-basics-caffeine
⌨ Experimental evidence (see {cite}`newtonPlasmaSalivary1981`) shows that a 300-mg oral dose of caffeine, such as might be found in a large mug of drip-brewed coffee, creates a concentration of about 8 $\mu{\rm g}$/mL in blood plasma. This boost is followed by first-order kinetics with a half-life of about 6 hours (although this rate can vary a great deal from person to person). We can model the caffeine concentration due to one drink taken over half an hour via

```{math}
:numbered: false
x'(t) = -kx + C(t),\quad x(0)=0,
```

where $k=\log(2)/6$ and

```{math}
:numbered: false
C(t) =
\begin{cases}
16, & 0\le t \le 0.5, \\
0, & t > 0.5.
\end{cases}
```

Solve the IVP and make a plot of the caffeine concentration for 12 hours. Then change $k=\log(2)/8$ (half-life of 8 hours) and plot the solution again.
``````

``````{exercise}
:label: problem-basics-skydiver
⌨ A reasonable model of the velocity $v(t)$ of a skydiver is

```{math}
:numbered: false
\frac{dv}{dt} = -g + \frac{k(t)}{m}v^2,  \quad v(0)=0,
```

where $g=9.8 \text{ m/sec}^2$ is gravitational acceleration, $m$ is the mass of the skydiver with parachute, and $k$ quantifies the effect of air resistance. At the US Air Force Academy, a training jump starts at about 1200 m and has $k=0.4875$ for $t<13$ and $k=29.16$ for $t\ge 13$. (This is an oversimplification; see {cite}`meadeDifferentialEquations1999`.) 

**(a)** Solve the IVP for $v$ for an 80-kg cadet for $t\in [0,200]$, and plot the solution.

**(b)** The total distance fallen up to time $t$ is $\displaystyle\int_0^t v(s)\, ds$. Use {numref}`Function {number} <function-intadapt>` to calculate and plot the altitude of the cadet as a function of time.

**(c)** In part (b), you should have found that the altitude becomes negative. Use {numref}`Function {number} <function-secant>` to determine accurately when the cadet reaches the ground.
``````
