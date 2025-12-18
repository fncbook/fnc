---
numbering:
  enumerator: 6.7.%s
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

(section-ivp-implicit)=

# Implementation of multistep methods

```{index} multistep method; implementation of
```

We now consider some of the practical issues that arise when multistep formulas are used to solve IVPs. In this section we emphasize the vector IVP, $\mathbf{u}'=\mathbf{f}(t,\mathbf{u})$, and use boldface in the difference formula {eq}`multistep` as well.

## Explicit methods

As a concrete example, the AB4 method is defined by the formula

```{math}
:label: ab4
\mathbf{u}_{i+1} = \mathbf{u}_i + h\, ( \tfrac{55}{24}\mathbf{f}_i - \tfrac{59}{24} \mathbf{f}_{i-1} + \tfrac{37}{24}\mathbf{f}_{i-2} - \tfrac{9}{24}\mathbf{f}_{i-3}), \quad i=3,\ldots,n-1.
```

{numref}`Function {number} <function-ab4>` shows a basic implementation of AB4.

Observe that {numref}`Function {number} <function-rk4>` is used to find the starting values $\mathbf{u}_1,\mathbf{u}_2,\mathbf{u}_3$ that are needed before the iteration formula takes over. As far as RK4 is concerned, it needs to solve  (the same step size as in the AB4 iteration). These results are then used to find $\mathbf{f}_0,\ldots,\mathbf{f}_3$ and get the main iteration started.

``````{prf:algorithm} ab4
:label: function-ab4
```{literalinclude} chapter06.py
:filename: ab4.py
:start-at: def ab4
:end-at: return t
:language: python
:linenos: true
```
::::{admonition} About the code
:class: dropdown
Line 15 sets `σ` to be the coefficients of the generating polynomial $\sigma(z)$ of AB4. Lines 19--21 set up the IVP over the time interval $a \le t \le a+3 h$, call `rk4` to solve it using the step size $h$, and use the result to fill the first four values of the solution. Then line 24 computes the vector $[f_2,f_1,f_0]$.

Line 28 computes $f_i$, based on the most recent solution value and time. That goes into the first spot of `f`, followed by the three values that were previously most recent. These are the four values that appear in {eq}`ab4`. Each particular $f_i$ value starts at the front of `f`, moves through each position in the vector over three iterations, and then is forgotten.
::::
``````

::::{prf:example} Convergence of Adams–Bashforth
:label: demo-implicit-ab4

We study the convergence of AB4 using the IVP $u'=\sin[(u+t)^2]$ over $0\le t \le 4$, with $u(0)=-1$. As usual, `solve_ivp` is called to give an accurate reference solution.

```{code-cell}
from scipy.integrate import solve_ivp
du_dt = lambda t, u: sin((t + u)**2)
tspan = (0.0, 4.0)
u0 = [-1.0]
u_ref = solve_ivp(du_dt, tspan, u0, dense_output=True, rtol=1e-13, atol=1e-13).sol
```

Now we perform a convergence study of the AB4 code.

```{code-cell}
n = array([int(4 * 10**k) for k in linspace(0, 3, 7)])
err = []
results = PrettyTable(["n", "AB4 error"])
for i in range(len(n)):
    t, u = FNC.ab4(du_dt, tspan, u0, n[i])
    err.append( abs(u_ref(4)[0] - u[0][-1]) )
    results.add_row([n[i], err[-1]])

print(results)
```

The method should converge as $O(h^4)$, so a log-log scale is appropriate for the errors.

```{code-cell}
loglog(n, err, "-o", label="AB4")
loglog(n, 0.5 * err[-1] * (n / n[-1])**(-4), "--", label="4th order")

xlabel("$n$"),  ylabel("final error")
legend(), title("Convergence of AB4");
```

::::

## Implicit methods

```{index} implicit IVP solver
```

The implementation of an implicit multistep method is more difficult. Consider the second-order implicit formula AM2, also known as the trapezoid method. To advance from step $i$ to $i+1$, we need to solve

```{math}
:label: AM2solve
  \mathbf{z} - \tfrac{1}{2} h f(t_{i+1},\mathbf{z})  = \mathbf{u}_i + \tfrac{1}{2} h \mathbf{f}(t_i,\mathbf{u}_i)
```

for $\mathbf{z}$. This equation can be written as $\mathbf{g}(\mathbf{z})=\boldsymbol{0}$, so the rootfinding methods of Chapter 4 can be used. The new value $\mathbf{u}_{i+1}$ is equal to the root of this equation.  

An implementation of AM2 using {numref}`Function {number} <function-levenberg>` from {numref}`section-nonlineqn-quasinewton` is shown in {numref}`Function {number} <function-am2>`.

``````{prf:algorithm} am2
:label: function-am2
```{literalinclude} chapter06.py
:filename: am2.py
:start-at: def am2
:end-at: return t
:language: python
:linenos: true
```
::::{admonition} About the code
:class: dropdown
Lines 22-23 define the function $\mathbf{g}$ and call `levenberg` to find the new solution value, using an Euler half-step as its starting value. A robust code would have to intercept the case where `levenberg` fails to converge, but we have ignored this issue for the sake of brevity.
::::
``````

## Stiff problems

At each time step in {numref}`Function {number} <function-am2>`, or any implicit IVP solver, a rootfinding iteration of uncertain expense is needed, requiring multiple calls to evaluate the function $\mathbf{f}$. This fact makes the cost of an implicit method much greater on a per-step basis than for an explicit one. Given this drawback, you are justified to wonder whether implicit methods are ever competitive! The answer is emphatically yes, as @demo-implicit-stiff demonstrates.

```{index} stiff differential equation
```

::::{prf:example} Stiffness
:label: demo-implicit-stiff

The following simple ODE uncovers a surprise.

```{code-cell}
f = lambda t, u: u**2 - u**3
u0 = array([0.005])
tspan = [0, 400]
```

We will solve the problem first with the implicit AM2 method using $n=200$ steps.

```{code-cell}
tI, uI = FNC.am2(f, [0.0, 400.0], u0, 200)
fig, ax = subplots()
ax.plot(tI, uI[0], label="AM2")
xlabel("$t$"), ylabel("$y(t)$");
```

So far, so good. Now we repeat the process using the explicit AB4 method.

```{code-cell}
tE, uE = FNC.ab4(f, [0.0, 400.0], u0, 200)
ax.scatter(tE, uE[0], label="AB4")
ax.set_ylim([-4, 2]), ax.legend()
fig
```

Once the solution starts to take off, the AB4 result goes catastrophically wrong.

```{code-cell}
uE[0, 104:111]
```

We hope that AB4 will converge in the limit $h\to 0$, so let's try using more steps.

```{code-cell}
plot(tI, uI[0], color="k", label="AM2")
tE, uE = FNC.ab4(f, [0, 400], u0, 1000)
plot(tE, uE[0], ".-", label="AM4, n=1000")
tE, uE = FNC.ab4(f, [0, 400], u0, 1600)
plot(tE, uE[0], ".-", label="AM4, n=1600")
xlabel("$t$"),  ylabel("$u(t)$")
legend();
```

So AB4, which is supposed to be _more_ accurate than AM2, actually needs something like 8 times as many steps to get a reasonable-looking answer!

::::

Although the result of @demo-implicit-stiff may seem counter-intuitive, there is no contradiction. A fourth-order explicit formula is more accurate than a second-order implicit one, in the limit $h\to 0$. But there is another limit to consider, $t\to \infty$ with $h$ fixed, and in this one the implicit method wins.

Problems for which implicit methods are much more efficient than explicit counterparts are called **stiff**. A complete mathematical description will wait for Chapter 11, but a sure sign of stiffness is the presence of phenomena on widely different time scales. In @demo-implicit-stiff, for instance, there are two slow periods during which the solution changes very little, interrupted by a very fast transition in the state. An explicit method "thinks" that the step size must always be dictated by the timescale of the fast transition, whereas an implicit method can take large steps during the slow periods.

```{index} adaptivity; in IVP solver
```

## Adaptivity

As with RK methods, we can run two time stepping methods simultaneously in order to estimate the error and adjust the step size. For example, we could pair AB3 with AB4 as practically no cost because the methods differ only in how they include known information from the recent past. The more accurate AB4 value should allow an accurate estimate of the local error in the AB3 value, and so on.

Because multistep methods rely on the solution history, though, changing the step size is more algebraically complicated than for RK methods. If $h$ is changed, then the historical values $\mathbf{u}_{i-1},\mathbf{u}_{i-2}\ldots$ and $\mathbf{f}_{i-1},\mathbf{f}_{i-2}\ldots$ are no longer given at the right moments in time to apply the iteration formula. A typical remedy is to use interpolation to re-evaluate the historical values at the appropriate times. The details are important but not especially illuminating, and we do not give them here.

## Exercises

% Must stay as #1

``````{exercise}
:label: problem-implicit-ab4tests
⌨ For each IVP, solve the problem using {numref}`Function {number} <function-ab4>` with $n=100$, and plot the solution and the error $u-\hat{u}$ on separate plots.

**(a)** $u' = -2t u, \ 0 \le t \le 2, \ u(0) = 2;\  \hat{u}(t) = 2e^{-t^2}$

**(b)** $u' = u + t, \ 0 \le t \le 1, \ u(0) = 2;\  \hat{u}(t) = -1-t+e^{3t}$

**(c)** $u' = x^2/[u(1+x^3)],\ 0 \le x \le 3, \ u(0) =1;\ \hat{u}(x) =[1+(2/3)\ln (1+x^3)]^{1/2}$

**(d)** $u''+ 9u = 9t, \: 0< t< 2\pi, \: u(0) =1,\: u'(0) = 1; \: \hat{u}(t) = t+\cos (3t)$

**(e)** $u''+ 9u = \sin(2t), \: 0< t< 2\pi, \: u(0) =2,\: u'(0) = 1$;
$\quad \hat{u}(t) = (1/5) \sin(3t) + 2 \cos (3t)+ (1/5) \sin (2t)$

**(f)** $u''- 9u = 9t \: 0< t< 1, \: u(0) =2,\: u'(0) = -1; \: \hat{u}(t) = e^{3t} + e^{-3t}-t$

**(g)** $u''+ 4u'+ 4u = t, \: 0< t< 4, \: u(0) =1,\: u'(0) = 3/4; \: \hat{u}(t) = (3t+5/4)e^{-2t} + (t-1)/4$

**(h)** $x^2 u'' +5xu' + 4u = 0,\: 1<x<e^2, \: u(1) =1, \: u'(1) = -1; \: \hat{u}(x) = x^{-2}( 1 + \ln x)$

**(i)** $2 x^2 u'' +3xu' - u = 0,\: 1<x<16, \: u(1) =4, \: u'(1) = -1$;
$\quad \hat{u}(x) = 2(x^{1/2} + x^{-1})$

``````

``````{exercise}
:label: problem-implicit-ab4converge
⌨ For each IVP in @problem-implicit-ab4tests, use {numref}`Function {number} <function-ab4>` for $n=10\cdot2^d$ and $d=1,\ldots,10$. Make a log-log convergence plot for the final time error $|u_n-\hat{u}(t_n)|$ versus $n$, and add a straight line indicating fourth-order convergence.
``````

``````{exercise}
:label: problem-implicit-am2
⌨ Repeat @problem-implicit-ab4tests using {numref}`Function {number} <function-am2>`.
``````

``````{exercise}
:label: problem-implicit-am2converge
⌨  Repeat @problem-implicit-ab4converge using {numref}`Function {number} <function-am2>` and comparing to second-order rather than fourth-order convergence.
``````

``````{exercise}
:label: problem-implicit-bd2
⌨ Using {numref}`Function {number} <function-am2>` as a model, write a function `bd2` that applies the BD2 method to solve an IVP. Use one step of backward Euler to get the needed starting value $u_1$. Test the convergence of your function on the IVP in part (c) of @problem-implicit-ab4tests.
``````

``````{exercise}
:label: problem-implicit-stiff
⌨ Within machine epsilon, the exact solution of the IVP in @demo-implicit-stiff satisfies $\hat{u}(400)\approx 1$. For $n=600,800,1000,\ldots,2000,$ solve that IVP using {numref}`Function {number} <function-ab4>` and {numref}`Function {number} <function-am2>`, keeping track of the error $|u_n-1|$ in each case. Make a log-log convergence plot of the error versus $n$ for both methods on the same graph. (They are dramatcally different, for reasons discussed in Chapter 11.)
``````

``````{exercise}
:label: problem-implicit-ivpimag
Consider the IVP

```{math}
\mathbf{u}'(t) = \mathbf{A} \mathbf{u}(t), \quad \mathbf{A}=
\begin{bmatrix}
0&-4\\4&0
\end{bmatrix}, \quad \mathbf{u}(0) =
\begin{bmatrix}
1\\0
\end{bmatrix}.
```

**(a)** ✍ Define $E(t) = \bigl\|\mathbf{u}(t)\bigr\|_2^2$. Show that $E(t)$ is constant. (Hint: differentiate $\mathbf{u}^T\mathbf{u}$ with respect to time and simplify it.)

**(b)** ⌨ Use {numref}`Function {number} <function-ab4>` to solve the IVP for $t\in[0,20]$ with $n=100$ and $n=150$. Plot $|E(t)-E(0)|$ versus time for both solutions on a single log-linear graph. You should see exponential growth in time. (In this regime, AB4 is acting unstably in a sense discussed in {numref}`section-diffusion-absstab`.)

**(c)** ⌨ Repeat part (b) with $n=400$ and $n=600$, but on a linear-linear plot. Now you should see only linear growth of $|E(t)-E(0)|$. (In this regime, AB4 is fully stable.)

**(d)** ⌨ Repeat part (b) with AM2 instead of AB4, on a linear-linear plot. You will find that AM2 conserves energy, just like the exact solution. 
``````

``````{exercise}
:label: problem-implicit-ab2
⌨ **(a)** Modify {numref}`Function {number} <function-ab4>` to implement the AB2 method.

**(b)** Repeat part (b) of @problem-implicit-ivpimag, using AB2 in place of AB4.

**(c)** Repeat part (c) of @problem-implicit-ivpimag, using AB2 in place of AB4.
``````

``````{exercise}
:label: problem-implicit-am1
⌨ **(a)** Modify {numref}`Function {number} <function-am2>` to implement the backward Euler (AM1) method.

**(b)** Repeat part (d) of @problem-implicit-ivpimag, using AM1 in place of AM2 and $n=400,800$. Does the AM1 method conserve energy?
``````
