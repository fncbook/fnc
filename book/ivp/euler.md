# Euler's method

```{index} Euler's method
```

Let a first-order initial-value problem be given in the form

```{math}
\begin{split}
  u'(t) &= f\bigl(t,u(t)\bigr), \qquad a \le t \le b,\\
  u(a)& =u_0.
\end{split}
```

```{margin}
A numerical method represents the solution of an IVP by its values at a finite collection of times.
```

We represent a numerical solution of an IVP by its values at a finite collection of nodes, which for now we require to be equally spaced:

```{math}
t_i = a + ih, \qquad h=\frac{b-a}{n}, \qquad i=0,\ldots,n.
```

The number $h$ is called the {term}`step size`.

Because we don't get exactly correct values of the solution at the nodes, we need to take some care with the notation. From now on we let $\hat{u}(t)$ denote the exact solution of the IVP. The approximate value at $t_i$ computed at the nodes by our numerical methods will be denoted by $u_i\approx \hat{u}(t_i)$. Because we are given the initial value $u(a)=u_0$ exactly, there is no need to distinguish whether we mean $u_0$ as the exact or the numerical solution.

Consider a piecewise linear interpolant to the (as yet unknown) values $u_0,u_1,\ldots$, $u_n$. Its derivative is piecewise constant with values

```{math}
\frac{u_{i+1} - u_{i}}{t_{i+1}-t_i} = \frac{u_{i+1}-u_i}{h}, \qquad t_i < t < t_{i+1},
```

where $i=0,\ldots,n-1$. We can connect this derivative to the differential equation by following the model of $u'=f(t,u)$:

```{math}
\frac{u_{i+1}-u_i}{h} = f(t_i,u_i), \qquad i=0,\ldots,n-1.
```

We could also view the left-hand side as a forward-difference approximation to $u'(t)$ at $t=t_i$. Either way, we can rearrange to get

```{math}
  :label: euler1
  u_{i+1}=u_i + h f(t_i,u_i), \qquad i=0,\ldots,n-1.
```

Together with the starting value $u_0$ from the initial condition, this formula defines an iteration known as {term}`Euler's method`.  It is an explicit method, meaning that the answer at the new time level is explicitly given in terms of the older time level(s).

A basic implementation of Euler's method is shown in {numref}`Function {number}<function-euler>`. It expects the problem to be specified in the form of a function $f$ of two arguments, an interval defining the time domain, and an initial condition. It also requires the number of intervals $n$ defined by the nodes (or equivalently, the number of steps in the iteration). The output of {numref}`Function {number}<function-euler>` is a vector of the nodes and a vector of approximate solution values at those nodes.

(function-euler)=

````{proof:function} euler
**Euler's method for an initial-value problem**

```{code-block} julia
:lineno-start: 1
"""
euler(ivp,n)

Apply Euler's method to solve the given IVP using `n` time steps.
Returns a vector of times and a vector of solution values.
"""
function euler(ivp,n)
    # Time discretization.
    a,b = ivp.tspan
    h = (b-a)/n
    t = [ a + i*h for i in 0:n ]

    # Initial condition and output setup.
    u = fill(float(ivp.u0),n+1)

    # The time stepping iteration.
    for i in 1:n
        u[i+1] = u[i] + h*ivp.f(u[i],ivp.p,t[i])
    end
    return t,u
end
```
````

## Local truncation error

It should be clear that Euler's method cannot get the exact solution unless it happens to be piecewise linear. More generally, suppose that an oracle has granted us the exact solution at $t=t_i$, so that $u_i=\hat{u}(t_i)$. How would we then fare in obtaining $u_{i+1}$ as an approximation to $\hat{u}(t_{i+1})$? The answer is revealed through a Taylor series:

```{math}
  :label: eulerLTE
\begin{split}
  \hat{u}(t_{i+1}) - \bigl[ u_i + hf(t_i,u_i) \bigr]
 &=  \hat{u}(t_{i+1}) - \bigl[ \hat{u}(t_i) + hf\bigl(t_i,\hat{u}(t_i)\bigr) \bigr] \\
 &= \bigl[ \hat{u}(t_i) + h \hat{u}'(t_i) + \tfrac{1}{2}h^2 \hat{u}''(t_i) + \cdots \bigr] \\
 &\qquad\quad  - \bigl[ \hat{u}(t_i) + h\hat{u}'(t_i) \bigr] \notag \\
  &= \tfrac{1}{2}h^2 \hat{u}''(t_i) + O(h^3),
\end{split}
```

where used the fact that $\hat{u}$ satisfies the differential equation.

We formalize this calculation as follows. Euler's method may be written in the abstract form

```{math}
:label: onestepODE
\mathbf{u}_{i+1} = u_i + h\phi(t_i,u_i,h), \qquad i=0,\ldots,n-1,
```

```{index} initial-value problem; one-step method for solving
```

```{index} truncation error; of a one-step IVP solver
```

which we call a general **one-step method**. Euler's method is the particular case $\phi(t,u,h) = f(t,u)$, but we will see other one-step methods in future sections. When we substitute the exact solution at $t=t_i$, calculate the resulting error at $t_{i+1}$, and divide by $h$, we get a quantity called the {term}`local truncation error` (LTE) of the one-step formula. In the general one-step formula this is

```{math}
  :label: onestepLTE
  \tau_{i+1}(h) := \frac{\hat{u}(t_{i+1})-\hat{u}(t_i)}{h} - \phi\bigl(t_i,\hat{u}(t_i),h\bigr).
```

Compared to {eq}`eulerLTE` there is an extra division by $h$, which we explain below. First, though, note that in the limit $h\rightarrow 0$ in {eq}`onestepLTE`, we obtain

```{math}
  \hat{u}'(t_i) - \phi(t_i,\hat{u}(t_i),0),
```

which through the ODE is the same as

```{math}
  f(t_i,\hat{u}(t_i)) - \phi(t_i,\hat{u}(t_i),0).
```

It seems reasonable to expect the LTE to vanish as the step size goes to zero for any ODE, which implies that $\phi(t,u,0)=f(t,u)$ for any function $u$ whatsoever. This condition on the one-step formula is called **consistency**. It is trivially true for Euler's method.

## Convergence

The local truncation error measures the effect of a single step of the numerical method. It's straightforward to calculate from the formula, but the practical quantity of interest is the {term}`global error`, $\hat{u}(t_i) - u_i$, over the entire time interval. By {eq}`onestepLTE`, $h\tau_{i+1}(h)$ describes how much error is made by taking a single step, starting from the exact value. If there were no other sources of or effects on the error, we would add up all of those local errors to get the global error.

To reach the time $t=b$ from $t=a$ with step size $h$, we need to take $n=(b-a)/h$ steps. If we want to reach, say, $t=(a+b)/2$, then we would have to take $n/2$ steps, and so on. The point is that to reach any fixed time in the interval, we need to take $O(n)=O(h^{-1})$ steps. That is why we express the error made in one step as $h\tau_{i+1}(h)$, with that extra factor of $h$ taken out. By this reasoning, for instance, the LTE of Euler computed in {eq}`eulerLTE` implies a global error that is $O(h)$.

However, global error is not just a simple sum of local errors. As each step causes a perturbation of the solution, we jump from one solution curve to a new one. The new curve will have its own trajectory, i.e., the error will propagate through the ODE (see {prf:ref}`demos-basics-cond`). This phenomenon is precisely the subject of {pref:ref}`theorem-depIC`: jumping to a different solution curve incurs a condition number at time $t>t_i$ of $e^{L(t-t_i)}$, which is constant at fixed time as $h\to 0$.

The following theorem puts our observations above on a rigorous footing.

````{prf:theorem}
:label: theorem-onestepGTE
Suppose that the unit local truncation error of the one-step method {eq}`onestepODE` satisfies
  
```{math}
  :label: ULTEbound
  |\tau_{i+1}(h)| \le C h^p
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
Define the global error sequence $E_i=\hat{u}(t_i)-u_i$. Using {eq}`onestepODE`, we obtain
  
```{math}
  E_{i+1} - E_i = \hat{u}(t_{i+1}) - \hat{u}(t_i) - ( \mathbf{u}_{i+1} - u_i ) =
  \hat{u}(t_{i+1}) - \hat{u}(t_i) - h\phi(t_i,u_i,h),
```

or

```{math}
  E_{i+1} = E_i + \hat{u}(t_{i+1}) - \hat{u}(t_i) - h\phi(t_i,\hat{u}(t_i),h) +
  h[\phi(t_i,\hat{u}(t_i),h)- \phi(t_i,u_i,h)].
```

We apply the triangle inequality,  {eq}`onestepLTE`, and {eq}`ULTEbound` to find

```{math}
  |E_{i+1}| \le |E_i| + Ch^{p+1} + h \left| \phi(t_i,\hat{u}(t_i),h)- \phi(t_i,u_i,h)\right|.
```

The Fundamental Theorem of Calculus implies that

```{math}
\begin{split}
  \left| \phi(t_i,\hat{u}(t_i),h)- \phi(t_i,u_i,h)\right|
      & = \left|  \int_{u_i}^{\hat{u}(t_i)} \frac{\partial \phi}{\partial u} \,du  \right|\\
    & \le  \int_{u_i}^{\hat{u}(t_i)} \left|\frac{\partial \phi}{\partial u}\right| \,du \\
    &= \le L | \hat{u}(t_i)-u_i| = L\, |E_i|.
\end{split}
```

Thus

```{math}
\begin{split}
  |E_{i+1}| &\le Ch^{p+1} + (1 + hL) |E_i| \\
  &\le Ch^{p+1} + (1 + hL) \bigl[ Ch^{p+1} + (1 + hL) |E_{i-1}|
  \bigr]\\
  &\;\vdots \\
  &\le Ch^{p+1} \left[ 1 + (1+hL) + (1+hL)^2 + \cdots + (1+hL)^i
  \right].
\end{split}
```

To get the last line we applied the inequality recursively until reaching $E_0$, which is zero. Replacing $i+1$ by $i$ and simplifying the geometric sum, we get

```{math}
  |E_i| \le Ch^{p+1}\frac{(1+hL)^i - 1}{(1+hL)-1} = \frac{Ch^p}{L}
  \left[ (1+hL)^i - 1 \right].
```

We observe that $1+x \le e^x$ for $x\ge 0$ (see [this exercise](problem-expdominate)). Hence $(1+hL)^i \le e^{i h L}$, which completes the proof.
````

```{margin}
The local truncation error of a one-step method has the same order of accuracy as the global error.
```

```{index} order of accuracy; of a one-step IVP solver
```

```{prf:example} Julia demo
:class: demo
:label: demos-euler-converge
{doc}`demos/euler-converge`
```

The theorem justifies a general definition of {term}`order of accuracy` as the leading exponent of $h$ in $\tau_{i+i}(h)$: the local truncation error of a one-step method has the same order of accuracy as the global error. This agrees with the first-order convergence we observed experimentally for Euler in {prf:ref}`demos-euler-converge`. Note, however, that the $O(h^p)$ convergence hides a leading constant that grows exponentially in time. When the time interval is bounded as $h\to 0$, this does not interfere with the conclusion, but the behavior as $t\to\infty$ contains no such guarantee.

## Exercises

1. ✍ Do two steps of Euler's method for the following problems using the given step size $h$. Then, compute the error using the given exact solution.

    **(a)** $u' = -2t u, \ u(0) = 2;\ h=0.1;\ \hat{u}(t) = 2e^{-t^2}$

    **(b)** $u' = u + t, \ u(0) = 2;\ h=0.2;\ \hat{u}(t) = -1-t+3e^t$

    **(c)** $t u' + u = 1, \ u(1) = 6, \ h = 0.25;\ \hat{u}(t) = 1+5/t$

    **(d)** $u' - 2u(1-u) = 0, \ u(0) = 1/2, \ h = 0.25; \ \hat{u}(t) = 1/(1 + e^{-2t})$

    ````{only} solutions
    ````

2. ⌨ For each IVP, solve the problem using {numref}`Function {number}<function-euler>`. (i) Plot the solution for $n=320$. (ii) For $n=10\cdot2^k$, $k=2,3,\ldots,10$, compute the error at the final time and make a log--log convergence plot, including a reference line for first-order convergence.

    **(a)** $u' = -2t u, \ 0 \le t \le 2, \ u(0) = 2;\  \hat{u}(t) = 2e^{-t^2}$

    **(b)** $u' = u + t, \ 0 \le t \le 1, \ u(0) = 2;\  \hat{u}(t) = -1-t+3e^t$

    **(c)** $(1+t^3)uu' = t^2,\ 0 \le xt \le 3, \ u(0) =1;\ \hat{u}(t) = [1+(2/3)\ln (1+xt^3)]^{1/2}$

    **(d)** $u' - 2u(1-u) = 0, \ 0 \le t \le 2, \ u(0) = 1/2; \ \hat{u}(t) = 1/(1 + e^{-2t})$

    **(e)** $v' - (1+x^2) v = 0, \ 1 \le x \le 3, \ v(1) = 1, \ \hat{v}(x) = e^{-(\pi/4)+\arctan(x)}$

    **(f)** $v' + (1+x^2) v^2 = 0, \ 0 \le x \le 2, \ v(0) = 2, \ \hat{v}(x) = 1/(0.5+\arctan(x))$

    **(g)** $u' = 2(1+t)(1+u^2), \ 0 \le t \le 0.5, \ u(0) = 0,  \ \hat{u}(t) = \tan(2t + t^2)$

    ````{only} solutions
    ````

3. (deleted)

4. ✍ Here is an alternative to Euler's method:

    ```{math}
    \begin{split}
      v_{i+1} &= u_i + h f(t_i,u_i)\\
      u_{i+1} &= u_i + hf(t_{i}+h,v_{i+1}).
    \end{split}
    ```

    **(a)** Write out the method explicitly in the general one-step form {eq}`onestepODE` (i.e., clarify what $\phi$ is for this method).

    **(b)** Show that the method is consistent.
  
    ````{only} solutions
    ````

5. ✍ Consider the problem $u'=ku$, $u(0)=1$ for constant $k$ and $t>0$.

      **(a)** Find an explicit formula in terms of $h$, $k$, and $i$ for the Euler solution $u_i$ at $t=ih$.

      **(b)** Find values of $k$ and $h$ such that $|u_i|\to\infty$ as $i\to\infty$ while the exact solution $\hat{u}(t)$ is bounded as $t\to\infty$.
  
    ````{only} solutions
    ````

    (problem-expdominate)=
6. ✍ Prove the fact, used in the proof of {prf:ref}`theorem-onestepGTE`, that $1+x\le e^x$ for all $x\ge 0$.

    ````{only} solutions
    ````

7. ✍ Suppose that the error in making a step is also subject to roundoff error $\epsilon_{i+1}$, so that $\tau_{i+1}(h) = Ch^p+\epsilon_{i+1} h^{-1}$; assume that $|\epsilon_{i+1}| \le \epsilon$ is the largest roundoff error in the computation and that the initial condition is known exactly. Generalize {prf:ref}`theorem-onestepGTE` for this case.

    ````{only} solutions
    ````
