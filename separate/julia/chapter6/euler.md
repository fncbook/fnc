---
numbering:
  enumerator: 6.2.%s
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.12
---
```{code-cell}
:tags: [remove-cell]
import Pkg; Pkg.activate("/Users/driscoll/Documents/GitHub/fnc")

using FNCFunctions

using Plots
default(
    titlefont=(11,"Helvetica"),
    guidefont=(11,"Helvetica"),
    linewidth = 2,
    markersize = 3,
    msa = 0,
    size=(500,320),
    label="",
    html_output_format = "svg"
)

using PrettyTables, LaTeXStrings, Printf
using LinearAlgebra

@ptconf backend = Val(:html) tf = tf_html_simple
```

(section-ivp-euler)=

# Euler's method

Let a first-order initial-value problem be given in the form

```{math}
:label: euler-ivp
\begin{split}
  u'(t) &= f\bigl(t,u(t)\bigr), \qquad a \le t \le b,\\
  u(a)& =u_0.
\end{split}
```

We represent a numerical solution of an IVP by its values at a finite collection of nodes, which for now we require to be equally spaced:

```{math}
:label: nodes-euler
t_i = a + ih, \qquad h=\frac{b-a}{n}, \qquad i=0,\ldots,n.
```

The number $h$ is called the **step size**.

Because we don't get exactly correct values of the solution at the nodes, we need to take some care with the notation. From now on we let $\hat{u}(t)$ denote the exact solution of the IVP. The approximate value at $t_i$ computed at the nodes by our numerical methods will be denoted by $u_i\approx \hat{u}(t_i)$. Because we are given the initial value $u(a)=u_0$ exactly, there is no need to distinguish whether we mean $u_0$ as the exact or the numerical solution.

Consider a piecewise linear interpolant to the (as yet unknown) values $u_0,u_1,\ldots$, $u_n$. For $t_i < t < t_{i+1}$, its slope is

```{math}
\frac{u_{i+1} - u_{i}}{t_{i+1}-t_i} = \frac{u_{i+1}-u_i}{h}.
```

We can connect this derivative to the differential equation by following the model of $u'=f(t,u)$:

```{math}
\frac{u_{i+1}-u_i}{h} = f(t_i,u_i), \qquad i=0,\ldots,n-1.
```

```{index} ! Euler's method
```

The left-hand side is a forward-difference approximation to $u'(t)$ at $t=t_i$. We can rearrange the equation to get **Euler's method**, our first method for IVPs.

::::{prf:definition} Euler's method for an IVP
:label: definition-eulerivp
Given the IVP $u'=f(t,u)$, $u(a)=u_0$, and the nodes {eq}`nodes-euler`, iteratively compute the sequence

```{math}
:label: euler1
  u_{i+1}=u_i + h f(t_i,u_i), \qquad i=0,\ldots,n-1.
```

Then $u_i$ is approximately the value of the solution at $t=t_i$.
::::

Euler's method marches ahead in $t$, obtaining the solution at a new time level explicitly in terms of the latest value.

A basic implementation of Euler's method is shown in {numref}`Function {number} <function-euler>`.

``````{prf:algorithm} euler
:label: function-euler

```{literalinclude} chapter06.jl
:filename: euler.jl
:start-after: # begin euler
:end-before: # end euler
:language: julia
:linenos: true
```
::::{admonition} About the code
:class: dropdown
The `ivp` input argument is an `ODEProblem`, like in @demo-basics-first. It has fields `ivp.f`, `ivp.tspan`, `ivp.u0`, and `ivp.p` that fully define the problem. The outputs are vectors of the nodes and approximate solution values at those nodes.
::::
``````

## Local truncation error

Let $\hat{u}(t)$ be the exact solution of the IVP {eq}`euler-ivp`, and suppose that somehow we have access to it at $t=t_i$, so that $u_i=\hat{u}(t_i)$. How good is $u_{i+1}$ as an approximation to $\hat{u}(t_{i+1})$? The answer is revealed through a Taylor series:

```{math}
:label: eulerLTE
\begin{split}
  \hat{u}(t_{i+1}) - \bigl[ u_i + hf(t_i,u_i) \bigr]
 &=  \hat{u}(t_{i+1}) - \bigl[ \hat{u}(t_i) + hf\bigl(t_i,\hat{u}(t_i)\bigr) \bigr] \\
 &= \bigl[ \hat{u}(t_i) + h \hat{u}'(t_i) + \tfrac{1}{2}h^2 \hat{u}''(t_i) + O(h^3) \bigr] - \bigl[ \hat{u}(t_i) + h\hat{u}'(t_i) \bigr] \notag \\
  &= \tfrac{1}{2}h^2 \hat{u}''(t_i) + O(h^3),
\end{split}
```

where we used the fact that $\hat{u}$ satisfies the differential equation.

We now introduce some formalities.

```{index} ! initial-value problem; one-step method for
```

::::{prf:definition} One-step IVP method
:label: definition-ivponestep
A **one-step method** for the IVP {eq}`euler-ivp` is a formula of the form

```{math}
:label: onestepODE
{u}_{i+1} = u_i + h\phi(t_i,u_i,h), \qquad i=0,\ldots,n-1.
```

::::

Euler's method is the particular case of {eq}`onestepODE` with $\phi(t,u,h) = f(t,u)$, but we will see other one-step methods in later sections.

```{index} ! truncation error; of a one-step IVP solver
```

In close analogy with {numref}`section-localapprox-fd-converge`, we define truncation error as the residual of {eq}`onestepODE` when the exact solution is inserted.

::::{prf:definition} Truncation error of a one-step IVP method
:label: definition-eulerlte
The {term}`local truncation error` (LTE) of the one-step method {eq}`onestepODE` is

```{math}
:label: onestepLTE
  \tau_{i+1}(h) := \frac{\hat{u}(t_{i+1})-\hat{u}(t_i)}{h} - \phi\bigl(t_i,\hat{u}(t_i),h\bigr).
```

The method is called **consistent** if $\tau_{i+1}(h)\to 0$ as $h\to 0$.
::::

One result follows immediately from the definitions:

:::{prf:lemma}
If $\phi(t,u,0)=f(t,u)$ for any function $u$, then the method {eq}`onestepODE` is consistent.
::::

## Convergence

While the local truncation error is straightforward to calculate from its definition, it is not the quantity we want to know about and control.

```{index} ! global error
```

::::{prf:definition} Global error of an IVP solution
:label: definition-globalerror
Given an IVP whose exact solution is $\hat{u}(t)$, the **global error** of approximate solution values $u_0,u_1,\ldots,u_n$ at times $t_i$ in {eq}`nodes-euler` is the vector $[ \hat{u}(t_i) - u_i ]_{\,i=0,\ldots,n}$.
::::

At times the term *global error* may be interpreted as the max-norm of the global error vector, or as its final value.

By our definitions, the local error in stepping from $t_i$ to $t_{i+1}$ is $h\tau_{i+1}(h)$. To reach the time $t=b$ from $t=a$ with step size $h$, we need to take $n=(b-a)/h$ steps. If we want to reach, say, $t=(a+b)/2$, then we would have to take $n/2$ steps, and so on. In fact, to reach any fixed time in the interval, we need to take $O(n)=O(h^{-1})$ steps. By expressing the local error with a factor of $h$ taken out, the LTE $\tau$ itself is accounting for the simple accumulation of error caused by taking $O(n)$ steps.[^altLTE]

[^altLTE]: Another point of view is that we can of course make local errors smaller by chopping $h$ in half, but then we have to take twice as many steps. The important quantity, then, is local error *per unit step length*, which is how $\tau$ is defined.

However, global error is not as simple as a sum of local errors. As explained in @theorem-depIC and illustrated in @demo-basics-cond, each step causes a perturbation of the solution that can grow as $t$ advances. Thus, we have to account for the flow evolution of individual step truncation errors as well as their mere accumulation. That is the subject of the following theorem.

````{prf:theorem}
:label: theorem-euler-onestepGTE
Suppose that the unit local truncation error of the one-step method {eq}`onestepODE` satisfies
  
```{math}
:label: ULTEbound
  |\tau_{i+1}(h)| \le C h^p,
```

and that

```{math}
:label: GTELip
\left| \frac{\partial \phi}{\partial u} \right| \le L
```

for all $t\in[a,b]$, all $u$, and all $h>0$. Then the global error satisfies

```{math}
:label: GTEbound
|\hat{u}(t_i) - u_i| \le \frac{Ch^p}{L} \left[ e^{L(t_i-a)} - 1
\right] = O(h^p),
```

as $h\rightarrow 0$.
````

````{prf:proof}
:enumerated: false

Define the global error sequence $ϵ_i=\hat{u}(t_i)-u_i$. Using {eq}`onestepODE`, we obtain
  
```{math}
ϵ_{i+1} - ϵ_i = \hat{u}(t_{i+1}) - \hat{u}(t_i) - ( {u}_{i+1} - u_i ) =
\hat{u}(t_{i+1}) - \hat{u}(t_i) - h\phi(t_i,u_i,h),
```

or

```{math}
  ϵ_{i+1} = ϵ_i + [\hat{u}(t_{i+1}) - \hat{u}(t_i) - h\phi(t_i,\hat{u}(t_i),h)] +
  h[\phi(t_i,\hat{u}(t_i),h)- \phi(t_i,u_i,h)].
```

We apply the triangle inequality,  {eq}`onestepLTE`, and {eq}`ULTEbound` to find

```{math}
  |ϵ_{i+1}| \le |ϵ_i| + Ch^{p+1} + h \left| \phi(t_i,\hat{u}(t_i),h)- \phi(t_i,u_i,h)\right|.
```

The Fundamental Theorem of Calculus implies that

```{math}
\begin{split}
  \left| \phi(t_i,\hat{u}(t_i),h)- \phi(t_i,u_i,h)\right|
      & = \left|  \int_{u_i}^{\hat{u}(t_i)} \frac{\partial \phi}{\partial u} \,du  \right|\\
    & \le  \int_{u_i}^{\hat{u}(t_i)} \left|\frac{\partial \phi}{\partial u}\right| \,du \\[1mm]
    & \le L | \hat{u}(t_i)-u_i| = L\, |ϵ_i|.
\end{split}
```

Thus

```{math}
\begin{split}
  |ϵ_{i+1}| &\le Ch^{p+1} + (1 + hL) |ϵ_i| \\
  &\le Ch^{p+1} + (1 + hL) \bigl[ Ch^{p+1} + (1 + hL) |ϵ_{i-1}|
  \bigr]\\
  &\;\vdots \\
  &\le Ch^{p+1} \left[ 1 + (1+hL) + (1+hL)^2 + \cdots + (1+hL)^i
  \right].
\end{split}
```

To get the last line we applied the inequality recursively until reaching $ϵ_0$, which is zero. Replacing $i+1$ by $i$ and simplifying the geometric sum, we get

```{math}
  |ϵ_i| \le Ch^{p+1}\frac{(1+hL)^i - 1}{(1+hL)-1} = \frac{Ch^p}{L}
  \left[ (1+hL)^i - 1 \right].
```

We observe that $1+x \le e^x$ for $x\ge 0$ (see @problem-euler-inequality). Hence $(1+hL)^i \le e^{i h L}$, which completes the proof.
````

```{index} ! order of accuracy; of a one-step IVP method
```

The theorem justifies one more general definition.

::::{prf:definition} Order of accuracy of a one-step IVP method
:label: definition-orderonestep
If the local truncation error of the one-step method {eq}`onestepODE` satisfies $\tau_{i+1}(h)=O(h^p)$ for a positive integer $p$, then $p$ is the **order of accuracy** of the formula.
::::

We could restate @theorem-euler-onestepGTE as saying that the global error has the same order of accuracy as the LTE. Note, however, that the $O(h^p)$ convergence hides a leading constant that grows exponentially in time. When the time interval is bounded as $h\to 0$, this does not interfere with the conclusion, but the behavior as $t\to\infty$ contains no such guarantee.

::::{prf:example} Convergence of Euler's method
:label: demo-euler-converge

We consider the IVP

```{math}
:numbered: false
u'=\sin[(u+t)^2], \quad t \in [0,4], \quad u(0)=-1.
```

```{code-cell}
using OrdinaryDiffEq
f(u, p, t) = sin((t + u)^2);
tspan = (0.0, 4.0);
u0 = -1.0;
ivp = ODEProblem(f, u0, tspan)
```

Here is the call to {numref}`Function {number} <function-euler>`.

```{code-cell}
using Plots
t, u = FNC.euler(ivp, 20)
plot(t, u;
    m=2,  label="n=20", 
    xlabel=L"t",  ylabel=L"u(t)",
    title="Solution by Euler's method")
```

We could define a different interpolant to get a smoother picture above, but the derivation of Euler's method assumed a piecewise linear interpolant, and that sets the limit of its accuracy. We can instead request more steps to make the interpolant look smoother.

```{code-cell}
t, u = FNC.euler(ivp, 50)
plot!(t, u, m=2, label="n=50")
```

Increasing $n$ changed the solution noticeably. Since we know that interpolants and finite differences become more accurate as $h\to 0$, we should anticipate the same behavior from Euler's method. We don't have an exact solution to compare to, so we will use a `DifferentialEquations` solver to construct an accurate reference solution.

```{code-cell}
u_exact = solve(ivp, Tsit5(), reltol=1e-14, abstol=1e-14)
plot!(u_exact, l=(2, :black), label="reference")
```

Now we can perform a convergence study.

```{code-cell}
n = [round(Int, 5 * 10^k) for k in 0:0.5:3]
err = []
for n in n
    t, u = FNC.euler(ivp, n)
    push!(err, norm(u_exact.(t) - u, Inf))
end
@pt :header=["n", "inf-norm error"] [n err]
```

The error is approximately cut by a factor of 10 for each increase in $n$ by the same factor. A log-log plot also confirms first-order convergence. Keep in mind that since $h=(b-a)/n$, it follows that $O(h)=O(n^{-1})$.

```{code-cell}
plot(n, err;
    m=:o, label="results",
    xaxis=(:log10, L"n"),  yaxis=(:log10, "inf-norm global error"),
    title="Convergence of Euler's method")

# Add line for perfect 1st order.
plot!(n, 0.5 * err[end] * (n / n[end]) .^ (-1), l=:dash, label=L"O(n^{-1})")
```

::::

Euler's method is the ancestor of the two major families of IVP methods presented in this chapter. Before we describe them, though, we generalize the initial-value problem itself in a crucial way.

## Exercises

``````{exercise}
:label: problem-euler-byhand
✍ Do two steps of Euler's method for the following problems using the given step size $h$. Then, compute the error using the given exact solution.

**(a)** $u' = -2t u, \ u(0) = 2;\ h=0.1;\ \hat{u}(t) = 2e^{-t^2}$

**(b)** $u' = u + t, \ u(0) = 2;\ h=0.2;\ \hat{u}(t) = -1-t+3e^t$

**(c)** $t u' + u = 1, \ u(1) = 6, \ h = 0.25;\ \hat{u}(t) = 1+5/t$

**(d)** $u' - 2u(1-u) = 0, \ u(0) = 1/2, \ h = 0.25; \ \hat{u}(t) = 1/(1 + e^{-2t})$
``````

``````{exercise}
:label: problem-euler-usage
⌨ For each IVP, solve the problem using {numref}`Function {number} <function-euler>`. (i) Plot the solution for $n=320$. (ii) For $n=10\cdot2^k$, $k=2,3,\ldots,10$, compute the error at the final time and make a log-log convergence plot, including a reference line for first-order convergence.

**(a)** $u' = -2t u, \ 0 \le t \le 2, \ u(0) = 2;\  \hat{u}(t) = 2e^{-t^2}$

**(b)** $u' = u + t, \ 0 \le t \le 1, \ u(0) = 2;\  \hat{u}(t) = -1-t+3e^t$

**(c)** $(1+t^3)uu' = t^2,\ 0 \le t \le 3, \ u(0) =1;\ \hat{u}(t) = [1+(2/3)\ln (1+t^3)]^{1/2}$

**(d)** $u' - 2u(1-u) = 0, \ 0 \le t \le 2, \ u(0) = 1/2; \ \hat{u}(t) = 1/(1 + e^{-2t})$

**(e)** $v' - (1+x^2) v = 0, \ 1 \le x \le 3, \ v(1) = 1, \ \hat{v}(x) = e^{(x^3+3x-4)/3}$

**(f)** $v' + (1+x^2) v^2 = 0, \ 0 \le x \le 2, \ v(0) = 2, \ \hat{v}(x) = 6/(2x^3+6x+3)$

**(g)** $u' = 2(1+t)(1+u^2), \ 0 \le t \le 0.5, \ u(0) = 0,  \ \hat{u}(t) = \tan(2t + t^2)$
``````

``````{exercise}
:label: problem-euler-improved
 ✍ Here is an alternative to Euler's method:

```{math}
\begin{split}
v_{i+1} &= u_i + h f(t_i,u_i),\\
u_{i+1} &= u_i + hf(t_{i}+h,v_{i+1}).
\end{split}
```

**(a)** Write out the method explicitly in the general one-step form {eq}`onestepODE` (i.e., clarify what $\phi$ is for this method).

**(b)** Show that the method is consistent.
``````

``````{exercise}
:label: problem-euler-bounded
✍ Consider the problem $u'=ku$, $u(0)=1$ for a real constant $k$ and $t>0$.

**(a)** Find an explicit formula in terms of $h$, $k$, and $i$ for the Euler solution $u_i$ at $t=ih$.

**(b)** Find values of $k$ and $h$ such that $|u_i|\to\infty$ as $i\to\infty$ while the exact solution $\hat{u}(t)$ is bounded as $t\to\infty$.

``````

``````{exercise}
:label: problem-euler-inequality
✍ Prove the fact, used in the proof of @theorem-euler-onestepGTE, that $1+x\le e^x$ for all $x\ge 0$.
``````

``````{exercise}
:label: problem-euler-roundoff
✍ Suppose that the error in making a step is also subject to roundoff error, so that the total local error per unit step is $Ch^p + \delta_{i+1} h^{-1}$. Assume that $|\delta_{i+1}| \le \macheps$ for all $i$ and that the initial condition is known exactly. Generalize @theorem-euler-onestepGTE for this case.
``````
