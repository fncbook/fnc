---
numbering:
  enumerator: 4.7.%s
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

(section-nonlineqn-nlsq)=

# Nonlinear least squares

After the solution of square linear systems, we generalized to the case of having more constraints to satisfy than available variables. Our next step is to do the same for nonlinear equations, thus filling out this table:

|                    | **linear** | **nonlinear** |
|:---------------:|:---------:|:--------------:|
| **square**         | $\mathbf{A}\mathbf{x}=\mathbf{b}$ |$\mathbf{f}(\mathbf{x})=\boldsymbol{0}$ |
| **overdetermined** | $\min\, \bigl\|\mathbf{A}\mathbf{x} - \mathbf{b}\bigr\|_2$ | $\min\, \bigl\|\mathbf{f}(\mathbf{x}) \bigr\|_2$ |

```{index} ! nonlinear least squares
```

```{prf:definition} Nonlinear least-squares problem
:label: definition-nlsq
Given a function $\mathbf{f}(\mathbf{x})$ mapping from $\real^n$ to $\real^m$, the {term}`nonlinear least-squares problem` is to find $\mathbf{x}\in\real^n$ such that $\bigl\|\mathbf{f}(\mathbf{x})\bigr\|_2$ is minimized.
```

As in the linear case, we consider only overdetermined problems, where $m>n$. Minimizing a positive quantity is equivalent to minimizing its square, so we could also define the result as minimizing $\phi(\mathbf{x})=\mathbf{f}(\mathbf{x})^T\mathbf{f}(\mathbf{x})$.

## Gauss–Newton method

You should not be surprised to learn that we can formulate an algorithm by substituting a linear model function for $\mathbf{f}$. At a current estimate $\mathbf{x}_k$ we define

```{math}
  \mathbf{q}(\mathbf{x})  = \mathbf{f}(\mathbf{x}_k) + \mathbf{A}_k(\mathbf{x}-\mathbf{x}_k),
```

where $\mathbf{A}_k$ is the exact $m\times n$ Jacobian matrix, $\mathbf{J}(\mathbf{x}_k)$, or an approximation of it as described in {numref}`section-nonlineqn-quasinewton`.

In the square case, we solved $\mathbf{q}=\boldsymbol{0}$ to define the new value for $\mathbf{x}$, leading to the condition $\mathbf{A}_k\mathbf{s}_k=-\mathbf{f}_k$, where $\mathbf{s}_k=\mathbf{x}_{k+1}-\mathbf{x}_k$. Now, with $m>n$, we cannot expect to solve $\mathbf{q}=\boldsymbol{0}$, so instead we define $\mathbf{x}_{k+1}$ as the value that minimizes $\| \mathbf{q} \|_2$.

```{index} ! Gauss–Newton method
```

::::{prf:definition} Gauss–Newton method
:label: definition-gaussnewton
Given $\mathbf{f}$ and a starting value $\mathbf{x}_1$, for each $k=1,2,3,\ldots$

1. Compute $\mathbf{y}_k = \mathbf{f}(\mathbf{x}_k)$ and $\mathbf{A}_k$, the exact or approximate Jacobian matrix at $\mathbf{x}_k$.
2. Solve the linear least squares problem $\argmin \| \mathbf{A}_k\mathbf{s}_k  + \mathbf{y}_k\|_2$ for $\mathbf{s}_k$.
3. Let $\mathbf{x}_{k+1} = \mathbf{x}_k + \mathbf{s}_k$.

Under appropriate conditions, the sequence $\mathbf{x}_{1}, \mathbf{x}_{2},\dots$ converges to a local minimizer of $\|\mathbf{f}(\mathbf{x})\|_2$.
::::

In brief, Gauss–Newton solves a series of linear least-squares problems in order to solve a nonlinear least-squares problem.

Surprisingly, {numref}`Function {number} <function-newtonsys>` and {numref}`Function {number} <function-levenberg>`, which were introduced for the case of $m=n$ nonlinear equations, work *without modification* as the Gauss–Newton method for the overdetermined case! The reason is that the backslash operator applies equally well to the linear system and linear least-squares problems, and nothing else in those functions was written with explicit reference to $n$.

## Convergence

In the multidimensional Newton method for a nonlinear system, we expect quadratic convergence to a solution in the typical case. For the Gauss–Newton method, the picture is more complicated.

As always in least-squares problems, the residual $\mathbf{f}(\mathbf{x})$ will not necessarily be zero when $\|\mathbf{f}\|$ is minimized. Suppose that the minimum value of $\|\mathbf{f}\|$ is $R>0$. In general, we might observe quadratic-like convergence until the iterate $\|\mathbf{x}_k\|$ is within distance $R$ of a true minimizer, and linear convergence thereafter. When $R$ is not sufficiently small, the convergence can be quite slow.

::::{prf:example} Convergence of nonlinear least squares
:label: demo-nlsq-converge

We will observe the convergence of {numref}`Function {number} <function-levenberg>` for different levels of the minimum least-squares residual. We start with a function mapping from $\real^2$ into $\real^3$, and a point that will be near the optimum.

```{code-cell}
g = lambda x: array([sin(x[0] + x[1]), cos(x[0] - x[1]), exp(x[0] - x[1])])
p = array([1, 1])
```

The function $\mathbf{g}(\mathbf{x}) - \mathbf{g}(\mathbf{p})$ obviously has a zero residual at $\mathbf{p}$. We'll make different perturbations of that function in order to create nonzero residuals.

```{code-cell}
for R in [1e-3, 1e-2, 1e-1]:
    # Define the perturbed function.
    f = lambda x: g(x) - g(p) + R * array([-1, 1, -1]) / sqrt(3)
    x = FNC.levenberg(f, [0, 0])
    r = x[-1]
    err = [norm(x[j] - r) for j in range(len(x) - 1)]
    normres = norm(f(r))
    semilogy(err, label=f"R={normres:.2g}")
title("Convergence of Gauss–Newton")
xlabel("iteration"), ylabel("error")
legend();
```

In the least perturbed case, where the minimized residual is less than $10^{-3}$, the convergence is plausibly quadratic. At the next level up, the convergence starts similarly but suddenly stagnates for a long time. In the most perturbed case, the quadratic phase is nearly gone and the overall shape looks linear.

::::

## Nonlinear data fitting

```{index} data fitting; nonlinear
```

In {numref}`section-leastsq-fitting` we saw how to fit functions to data values, provided that the set of candidate fitting functions depends linearly on the undetermined coefficients. We now have a tool to generalize that process to fitting functions that depend nonlinearly on unknown parameters.

Suppose that $(t_i,y_i)$ for $i=1,\ldots,m$ are given points. We wish to model the data by a function $g(t,\mathbf{x})$ that depends on unknown parameters $x_1,\ldots,x_n$ in an arbitrary way. A standard approach is to minimize the discrepancy between the model and the observations, in a least-squares sense. Define

```{math}
:label: nlsq-misfit
\mathbf{f}(\mathbf{x}) = \left[\, g(t_i,\mathbf{x})-y_i  \, \right]_{\,i=1,\ldots,m}.
```

We call $\mathbf{f}$ a **misfit** function. By minimizing $\bigl\| \mathbf{f}(\mathbf{c}) \bigr\|^2$, we get the best possible fit to the data. If an explicit Jacobian matrix is desired for the minimization, we can compute

```{math}
\mathbf{f}{\,}'(\mathbf{x}) = \left[ \frac{\partial}{\partial x_j} g(t_i,\mathbf{x}) \right]_{\,i=1,\ldots,m;\,j=1,\ldots,n.}
```

The form of $g$ is up to the modeler. There may be compelling theoretical choices, or you may just be looking for enough algebraic power to express the data well. Naturally, in the special case where the dependence on $\mathbf{c}$ is linear, i.e.,

```{math}
  g(t,\mathbf{c}) = c_1 g_1(t) + c_2 g_2(t) + \cdots + c_m g_m(t),
```

then the misfit function is also linear in $\mathbf{c}$ and the fitting problem reduces to linear least squares.

::::{prf:example} Nonlinear data fitting
:label: demo-nlsq-MM
Inhibited enzyme reactions often follow what are known as *Michaelis–Menten_ kinetics, in which a reaction rate $w$ follows a law of the form

$$w(s) = \frac{V s}{K_m + s},$$

where $s$ is the concentration of a substrate. The real values $V$ and $K_m$ are parameters that are free to fit to data. For this example, we cook up some artificial data with $V=2$ and $K_m=1/2$.

```{code-cell}
m = 25
V, Km = 2, 0.5
s = linspace(0.05, 6, m)
model = lambda x: V * x / (Km + x)
w = model(s) + 0.15 * cos(2 * exp(s / 16) * s)    # noise added

fig, ax = subplots()
ax.scatter(s, w, label="data")
ax.plot(s, model(s), 'k--', label="unperturbed model")
xlabel("s"), ylabel("w")
legend();
```

```{index} ! Python; destructuring
```

The idea is to pretend that we know nothing of the origins of this data and use nonlinear least squares to find the parameters in the theoretical model function $v(s)$. In {eq}`nlsq-misfit`, the $s$ variable plays the role of $t$, and $v$ plays the role of $g$.
```{tip}
:class: dropdown
Putting comma-separated values on the left of an assignment will **destructure** the right-hand side, drawing individual assignments from entries of a vector, for example.
```

```{code-cell}
def misfit(c):
    V, Km = c  # rename components for clarity
    f = V * s / (Km + s) - w
    return f
```

In the Jacobian the derivatives are with respect to the parameters in $\mathbf{x}$.

```{code-cell}
def misfitjac(x):
    V, Km = x   # rename components for clarity
    J = zeros([m, 2])
    J[:, 0] = s / (Km + s)          # d/d(V)
    J[:, 1] = -V * s / (Km + s)**2  # d/d(Km)
    return J
```

```{code-cell}
x1 = [1, 0.75]
x = FNC.newtonsys(misfit, misfitjac, x1)
V, Km = x[-1]  # final values
print(f"estimates are V = {V:.3f}, Km = {Km:.3f}")
```

The final values are reasonably close to the values $V=2$, $K_m=0.5$ that we used to generate the noise-free data. Graphically, the model looks close to the original data.

```{code-cell}
# since V and Km have been updated, model() is too
ax.plot(s, model(s), label="nonlinear fit")
```

For this particular model, we also have the option of linearizing the fit process. Rewrite the model as 

```{math}
:enumerated: false
\frac{1}{w} = \frac{\alpha}{s} + \beta = \alpha \cdot s^{-1} + \beta
```

for the new fitting parameters $\alpha=K_m/V$ and $\beta=1/V$. This corresponds to the misfit function whose entries are

$$f_i([\alpha,\beta]) = \left(\alpha \cdot \frac{1}{s_i} + \beta\right) - \frac{1}{w_i}$$
for $i=1,\ldots,m$. Although this misfit is nonlinear in $s$ and $w$, it's linear in the unknown parameters $\alpha$ and $\beta$. This lets us pose and solve it as a linear least-squares problem.

```{code-cell}
from numpy.linalg import lstsq
A = array( [[1 / s[i], 1.0] for i in range(len(s))] )
z = lstsq(A, 1 / w, rcond=None)[0]
alpha, beta = z
print("alpha:", alpha, "beta:", beta)
```

The two fits are different; they do not optimize the same quantities.

```{code-cell}
linmodel = lambda x: 1 / (beta + alpha / x)
ax.plot(s, linmodel(s), label="linear fit")
ax.legend()
fig
```

The truly nonlinear fit is clearly better in this case. It optimizes a residual for the original measured quantity rather than a transformed one we picked for algorithmic convenience.

::::

## Exercises

``````{exercise}
:label: problem-nlsq-onevar
✍ Define $\mathbf{f}(x)=[ x-8, \; x^2-4 ]$.

**(a)** Write out the linear model of $\mathbf{f}$ at $x=2$.

**(b)** Find the estimate produced by one step of the Gauss–Newton method, starting at $x=2$.
``````

``````{exercise}
:label: problem-nlsq-quadratic
✍ (Continuation of @problem-nlsq-onevar.) The Gauss–Newton method replaces $\mathbf{f}(\mathbf{x})$ by a linear model and minimizes the norm of its residual. An alternative is to replace $\| \mathbf{f}(\mathbf{x}) \|_2^2$ by a scalar quadratic model $q(\mathbf{x})$, and minimize $q$.

**(a)** Using $\mathbf{f}(x) = [ x-8, \; x^2-4 ]$, let $q(x)$ be defined by the first three terms in the Taylor series for $\phi(x) = \| \mathbf{f}(x) \|_2^2$ at $x=2$.

**(b)** Find the unique $x$ that minimizes $q(x)$. Is the result the same as the estimate produced by Gauss–Newton?
``````

``````{exercise}
:label: problem-nlsq-epidemic
⌨  A famous result by Kermack and McKendrick in 1927 {cite}`kermackContributionMathematical1927` suggests that in epidemics that kill only a small fraction of a susceptible population, the death rate as a function of time is well modeled by

```{math}
w'(t) = A \operatorname{sech}^2[B(t-C)]
```

for constant values of the parameters $A,B,C$. Since the maximum of sech is $\operatorname{sech}(0)=1$, $A$ is the maximum death rate and $C$ is the time of peak deaths. You will use this model to fit the deaths per week from plague recorded in Mumbai 
during 1906:

``` julia
5, 10, 17, 22, 30, 50, 51, 90, 120, 180, 292, 395, 445, 775, 780,
700, 698, 880, 925, 800, 578, 400, 350, 202, 105, 65, 55, 40, 30, 20
```

**(a)** Use {numref}`Function {number} <function-levenberg>` to find the best least-squares fit to the data using the $\operatorname{sech}^2$ model. Make a plot of the model fit superimposed on the data. 

**(b)** Repeat part (a) using only the first 15 data values. (The fit will not be as good.)
``````

``````{exercise}
:label: problem-nlsq-circlefit
⌨  (Variation on @problem-newtonsys-circlefit.) Suppose the points $(x_i,y_i)$ for $i=1,\ldots,m$ are given, and the goal is to find the circle with center $(a,b)$ and radius $r$ that best fits the points. Define 

$$f_i(a,b,r) = (a-x_i)^2 + (b-y_i)^2 - r^2, \qquad i=1,\ldots,m.$$ 

Then we can define the best circle as the one that minimizes $\|\mathbf{f}\|$.  

Define data points as follows:

``` python
import numpy as np
m = 30;
t = 2*np.pi * np.random.rand(m)
x = -2 + 5 * np.cos(t) + 0.2 * np.random.randn(m)
y =  1 + 5 * np.sin(t) + 0.2 * np.random.randn(m)
```

Use {numref}`Function {number} <function-levenberg>` to find the best-fit circle, and make a plot of the circle superimposed on the points.
``````

``````{exercise}
:label: problem-nlsq-blink
⌨ The position of the upper lid during an eye blink can be measured from high-speed video {cite}`wuEffectsMild2014`, and it may be possible to classify blinks based in part on fits to the lid position {cite}`broschBlinkCharacterization2017`. The lid position functions proposed to fit blinks is a product of a monomial or polynomial multiplying a decaying exponential {cite}`berkeKineticsLid1998`.  In this problem, you will generate representative data, add a small amount of noise to it, and then perform nonlinear least-squares fits to the data.

**(a)** Consider the fitting function $y(\mathbf{a}) = a_1 t^2 \exp \left( -a_2 t^{a_3} \right)$, using the vector of coefficients $\mathbf{a} = [a_1,a_2,a_3]$, and create synthetic eyelid position data as follows:

``` python
import numpy as np
N = 20                                  # number of t values
t = (np.arange(1, N+1)) / N;            # equally spaced to t=1
a = [10, 10, 2]                         # baseline values
y = a[0] * t**2 * np.exp(-a[1]*t**a[2]) # ideal data
noise = 0.03                            # amplitude of noise
ym = y + noise * np.random.rand(N)      # add noise
```

**(b)** Using the data `(t,ym)`, find the nonlinear least-squares fit using {numref}`Function {number} <function-levenberg>`.

**(c)** Plot the fit at 300 points in $t\in[0,1]$, together with the measured data points `ym`.

**(d)** Increase the noise to 5% and 10%. You may have to increase the number of measured points `N` and/or the maximum number of iterations.  How close are the coefficients?  Plot the data and the resulting fit for each case.
``````

``````{exercise}
:label: problem-nlsq-blink2
⌨ Repeat the previous problem using the fitting function $y(\mathbf{a}) = (a_1+a_2 t + a_3 t^2) t^2 \exp \left( -a_4 t^{a_5} \right)$, using the vector of coefficients $\mathbf{a} = [a_1,\ldots,a_5]$. (This was the choice used in Brosch et al {cite}`broschBlinkCharacterization2017`.)  Use `a = [20, -10, -8, 7, 2]` to create the data and as an initial guess for the coefficients for the fit to the noisy data.
``````
