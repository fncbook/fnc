# The rootfinding problem

For the time being we will focus on single equations in one variable.

```{prf:definition} Rootfinding problem
Given a continuous scalar function $f$ of a scalar variable, find a real number $r$ such that $f(r)=0$.
```

```{index} roots
```

We call $r$ a **root** of the function $f$. The formulation $f(x)=0$ is general enough to solve any equation, for if we are given an equation $g(x)=h(x)$, we can define $f=g-h$ and find a root of $f$.

````{prf:example} Julia demo
:class: demo
:label: demos-roots-bessel
{doc}`demos/roots-bessel`
````

Unlike the linear problems of the earlier chapters, the usual situation here is that the root cannot be produced in a finite number of operations, even in exact arithmetic. Instead, we seek a sequence of approximations that formally converge to the root, stopping when some member of the sequence seems to be "good enough" (more on that later). The `NLsolve` package for Julia has a function `nlsolve` for general-purpose rootfinding.

## Conditioning, error, and residual

In the rootfinding problem the data is a continuous function $f$ whose root we seek. Let's assume $f$ has at least one continuous derivative near a particular root $r$. Say that $f$ is perturbed to $\tilde{f}(x) = f(x) + \epsilon$. As a result, the root will be perturbed to $\tilde{r} = r + \delta$ satisfying, by definition, $\tilde{f}(\tilde{r})=0$. We now compute an absolute condition number $\kappa_r$, which is the ratio $\left | \frac{\delta}{\epsilon} \right|$ as $\epsilon\to 0$.

Using Taylor series expansions,

```{math}
  0 = f(r+\delta) + \epsilon = f(r) + f'(r) \delta + \epsilon + O(\delta^2).
```

Given that $f(r)=0$, this implies

```{math}
  :label: rootcondnum
  \kappa_r = \bigl| f'(r) \bigr|^{-1}.
```

```{margin}
The condition number of the rootfinding problem is equivalent to that for evaluating the inverse function.
```

```{index} condition number; of rootfinding
```

````{prf:example} Julia demo
:class: demo
:label: demos-roots-cond
{doc}`demos/roots-cond`
````

Recall that the absolute condition number for the evaluation of $f$ at $x$ is simply $|f'(x)|$. The condition number of the rootfinding problem is equivalent to that of evaluating the inverse function. The influence of $|f'(r)|$ on rootfinding is easily illustrated.

```{index} residual
```

We must accept that when $|f'|$ is small at the root, it may not be possible to get a small error in a computed root estimate. As always, the error is not a quantity we can compute without knowing the exact answer. But we can compute the {term}`residual` of a root estimate, which is the value of $f$ there. Since $f(r)=0$, it stands to reason that a small residual might be associated with a small error.

Suppose that we find an approximation $\tilde{r}$ to the actual root $r$. Define the new function $g(x)=f(x)-f(\tilde{r})$. Trivially, $g(\tilde{r})=0$, meaning that the root estimate is a true root of $g$. Since the difference between $g$ and the original $f$ is the residual value $f(\tilde{r})$, the residual is the distance to a rootfinding problem that our root estimate solves exactly. That is, the residual is the backward error of the estimate. 

```{margin}
It is not always possible to get a small error in a root approximation. But the residual equals the backward error.
```

```{index} backward error
```

To summarize: In general, it is not always realistic to expect a small error in a root approximation. However, the backward error is the same as the residual of the estimate.

## Multiple roots

```{index} roots; nonsimple
```

The condition number {eq}`rootcondnum` naturally leads to the question of what happens if $f'(r)=0$ at a root $r$. Suppose first that $f$ is a polynomial of degree $n>0$, so that

```{math}
  f(x) = (x-r)q(x),
```

for a polynomial $q$ of degree $n-1$. If $r$ is a simple root of $f$—that is, it appears just once in the list of the $n$ roots—then it follows that $q(r)\neq 0$. Conversely, if $q(r)=0$, then $r$ appears among the roots of $q$ and is a multiple root of $f$. However, we don't need to know the quotient polynomial $q$ explicitly in order to make the determination. Consider that

```{math}
  f'(x) = (x-r)q'(x) + q(x),
```

so that $f'(r) = q(r)$. Hence $r$ is simple if and only if $f'(r)\neq 0$.

```{margin}
$r$ is a simple root of $f$ if and only if $f'(r)\neq 0$.
```

This conclusion extends to non-polynomial differentiable functions $f$. If $r$ is a root of $f$, define $q(x)=f(x)/(x-r)$. By L'Hôpital's rule, $g$ is well defined at $x=r$ as long as $f'$ is. Now we can again write $f(x)=(x-r)q(x)$ for $x\neq r$, and by continuity it works at $x=r$ as well. So  the reasoning we applied to polynomials can be repeated: $r$ is a simple root of $f$ if and only if $f'(r)\neq 0$.

Now suppose that $f'(r)=q(r)=0$, so that $r$ is not simple. If $q$ is differentiable, we may apply the same logic to it that we did to $f$. Hence $r$ is not simple for $q$ if and only if $q'(r)=0$. Now observe that

```{math}
  f''(x) = (x-r)q''(x) + 2q'(x),
```

and thus $f''(r)=2q'(r)$, so $f''(r)=0$ if and only if $r$ is a multiple root of $q$. In general we define $r$ as a **root of multiplicity $m$** if $f(r)=f'(r)=\cdots=f^{(m-1)}(r)=0$, but $f^{(m)}(r)\neq 0$. If $m=1$, we say $r$ is a {term}`simple root`.

It's useful to think through the consequences of these definitions for the Taylor series at the point $r$,

```{math}
  f(x) = a_0 + a_1(x-r) + a_2(x-r)^2 + \cdots,
```

where $a_n=f^{(n)}(r)/n!$. The fact that $r$ is a root implies $f(r)=a_0=0$. If $r$ is a simple root, then $a_1\neq 0$, and conversely. If $r$ is a double root, then $a_2\neq 0$, and so on. Simply put, if $r$ is a root of order $m$, then the series expansion begins with $(x-r)^m$.

When $r$ is a multiple root, the condition number {eq}`rootcondnum` is apparently infinite.[^infcond] However, even if $r$ is technically simple, we should expect difficulty if the condition number is very large. This occurs when $|f'(r)|$ is very small, which means that quotient $q$ satisfies $q(r)\approx 0$ and another root of $f$ is very close to $r$. The situation is reminiscent of the linear system problem: the degenerate case (singular matrix/multiple root) is mathematically isolated when considered exactly, but the effect on fixed precision computation is just as drastic in a neighborhood of the singularity.

[^infcond]: Based on our definitions, this means that the relative change to the root when $f$ is changed by a perturbation of size $\epsilon$ is not $O(\epsilon)$ as $\epsilon\to 0$.

## Exercises

1. ✍ For each function, find the multiplicity of the given root. If it is a simple root, find its absolute condition number.
  
    **(a)** $f(x) = x^3-2x^2+x-2$, root $r=2$

    **(b)** $f(x) = (\cos x  + 1)^2$, root $r=\pi$

    **(c)** $f(x) = \frac{\sin^2 x}{x}$, root $r=0$ (define $f(0) =0$)

    **(d)** $f(x) =(x-1)\log(x)$, root $r=1$
  
    %(a) We need to calculate the derivative $f'(x)=3x^2-4x+1$ at $x=r=2$; we find $f'(2)=5$.
    %Since $f'(r)\neq 0$, it is a simple root with multiplicity 1.  For condition number, we also need $||f||_\infty$ on an interval around the root.  Using $x\in[3/2,5/2]$, and either (i) a lot of points in Matlab or (ii) analyzing the function in this interval, one can find that $||f||_\infty=|f(5/2)|=3.625$.  Then
    %\[ \kappa (r) = \frac{||f||_\infty}{|r||f'(r)|} = \frac{3.625}{2\times 5} = 0.3625. \]
    %Thus the conditioning for this problem is good.\\
    %
    %(b) For this case, we calculate that $f'(i)=-2-4i\neq 0$, so that this is also a simple root with multiplicity 1. Using an interval of unit width centered around $r=i$ in either the imaginary or real direction
    %results in a maximum magnitude $\max |f|=3.625$ (still).  So
    %\[ \kappa (r) = \frac{||f||_\infty}{|r||f'(r)|} = \frac{3.625}{1\times \sqrt{20}} \approx 0.81. \]
    %The conditioning is also good for this root.

    (problem-conddoubleroot)=
2. For any $\epsilon>0$, let $f_\epsilon(x) = \sin[(x-1+\epsilon)(x-1)]$. This function has roots at $1-\epsilon$ and $1$.
  
    **(a)** ✍ Find $|f_\epsilon'(1)|$. According to  {eq}`rootcondnum`, the condition number of the root $r=1$ is inversely proportional to this quantity.

    **(b)** ⌨ Define a perturbation function by $g(x) = \cos[10x+\sin(20x)]$. Verify that $f_\epsilon(x)+10^{-10}g(x)$ has a root in the interval $[1-\epsilon,1+10\epsilon]$ for  $\epsilon=10^{-3}$.

    **(c)** ⌨ For $\epsilon=10^{-3},10^{-4},10^{-5},10^{-6}$, use `nlsolve` to find the root in the interval given in part (b). Make a table of $\epsilon$, $1/|f_\epsilon'(1)|$, $|r-1|$, and $|(r-1)f_\epsilon'(1)|$. The last value should be approximately constant.
  
    % partial answer
    <!-- ep = 1e-6;
    f = @(x) sin( (x-1+ep).*(x-1) );
    g = @(x) cos(10*x+sin(20*x));
    r = fzero( @(x) f(x)+1e-10*g(x), [1-ep,1+10*ep] );
    [ep, abs(r-1), ep*abs(r-1) ]
     -->

3. (continuation) The condition number theory {eq}`rootcondnum` suggests that if $f$ has simple roots at $r_1$ and $r_2$ that are close to each other, then the condition number is large. But then the function $\tilde{f}(x) = f(x)/(x-r_2)$ no longer has a root at $r_2$ and should have a much better condition number if there are no other nearby roots. This trick (called *deflation*) can work even if the division factor does not use the exact root $r_2$. For the steps below, use the same definitions as in the preceding problem.
  
    **(a)** ✍ Define $\tilde{f}_\epsilon(x) = f_\epsilon(x)/(x+1.01)$. Find $|\tilde{f}_\epsilon'(1)|$. (You may want to use computer algebra for this.)

    **(b)** ⌨ Repeat part (c) of the preceding problem, but using the interval $[1-\epsilon/200,1+\epsilon/100]$ each time. By what factor are
    the errors improved over using $f_\epsilon$? Do you still find that $|(x-1)\tilde{f}_\epsilon'(1)|$ is roughly constant?
  
    % incomplete
    <!-- ep = 1e-6;
    ft = @(x) f(x)/(x-.99);
    r = fzero(@(x) ft(x)+1e-10*g(x),[1-ep/200,1+ep/100]);
    err = 1 - r
    ep*err -->
