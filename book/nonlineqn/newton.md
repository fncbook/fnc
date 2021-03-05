# Newton's method in one variable

```{index} Newton's method
```

````{prf:example} Julia demo
:class: demo
{doc}`demos/newton-line`
````

Newton's method is the cornerstone of rootfinding. We introduce the key idea with an example in {doc}`demos/newton-line`.

Using general notation, if we have a root approximation $x_k$, we can construct a **linear model** of $f(x)$ using the classic formula for the tangent line of a differentiable function,

```{math}
  :label: tangentline
  q(x) = f(x_k) + f'(x_k)(x-x_k).
```

Finding the root of $q(x)=0$ is trivial. We define the next approximation by the condition $q(x_{k+1})=0$, or

```{math}
  :label: newton
  x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)}.
```

Starting with an initial estimate $x_1$, this formula defines a sequence of estimates $x_2,x_3,\ldots$. The iteration so defined is what we call {term}`Newton's method`.

## Quadratic convergence

The graphs of {doc}`demos/newton-line` suggest why the Newton iteration may converge to a root: Any differentiable function looks more and more like its tangent line as we zoom in to the point of tangency. Yet it is far from clear that it *must* converge, or at what rate it will do so. The matter of the convergence rate is fairly straightforward to resolve. Define the error sequence

```{math}
:label: errorseq
\epsilon_k = r - x_k, \quad k=1,2,\ldots,
```

where $r$ is the limit of the sequence and $f(r)=0$. Exchanging $x$-values for $\epsilon$-values in {eq}`newton` gives

```{math}
  \epsilon_{k+1} = \epsilon_k + \frac{f(r-\epsilon_k)}{f'(r-\epsilon_k)}.
```

We assume that $|\epsilon_k|\to 0$; eventually, the errors remain as small as we please forever. Then a Taylor expansion of $f$ about $x=r$ gives

```{math}
  \epsilon_{k+1} = \epsilon_k + \frac{ f(r) - \epsilon_kf'(r) + \frac{1}{2}\epsilon_k^2f''(r) +
    O(\epsilon_k^3)}{ f'(r) - \epsilon_kf''(r) + O(\epsilon_k^2)}.
```

We use the fact that $f(r)=0$ and additionally assume now that $f'(r)\ne 0$. Then

```{math}
\epsilon_{k+1} = \epsilon_k - \epsilon_k \left[ 1 - \dfrac{1}{2}\dfrac{f''(r)}{f'(r)} \epsilon_k
+ O(\epsilon_k^2)\right] \, \left[ 1 -  \dfrac{f''(r)}{f'(r)}\epsilon_k + O(\epsilon_k^2)\right]^{-1}.
```

The series in the denominator is of the form $(1+z)^{-1}$. Provided $|z|<1$, this is the limit of the geometric series $1-z+z^2-z^3 + \cdots$. Keeping only the lowest-order terms, we derive

```{math}
:label: newtonerr
\begin{split}
\epsilon_{k+1} &= \epsilon_k - \epsilon_k \left[ 1 - \dfrac{1}{2}\dfrac{f''(r)}{f'(r)} \epsilon_k + O(\epsilon_k^2) \right] \, \left[ 1 + \dfrac{f''(r)}{f'(r)}
\epsilon_k + O(\epsilon_k^2) \right]\\
&= -\frac{1}{2}\, \frac{f''(r)}{f'(r)} \epsilon_k^2 + O(\epsilon_k^3).
\end{split}
```

```{margin}
Asymptotically, each iteration of Newton's method roughly squares the error.
```

```{index} quadratic convergence
```

```{index} convergence rate; quadratic
```

Equation {eq}`newtonerr` suggests that eventually, each iteration of Newton's method roughly squares the error. This behavior is called {term}`quadratic convergence`. The formal definition of quadratic convergence is that there exists a number $\alpha>0$ such that

```{math}
  :label: quadratic-convergence
  \lim_{k\to\infty} \frac{|x_{k+1}-r|}{|x_k-r|^2} = \alpha.
```

Recall that linear convergence is identifiable by trending toward a straight line on a log--linear plot of the error. When the convergence is quadratic, no such straight line exists—the convergence keeps getting steeper. Alternatively, note that (neglecting high-order terms)

```{math}
  \log(|\epsilon_{k+1}|) \approx 2 \log(|\epsilon_{k}|) + \text{constant},
```

````{prf:example} Julia demo
:class: demo
{doc}`demos/newton-converge`
````

which is equivalent to saying that the number of accurate digits approximately doubles at each iteration, once the errors become small enough.

```{index} roots; nonsimple
```

Let's revisit the assumptions made to derive quadratic convergence as given by {eq}`newtonerr`:

1. The residual function $f$ has to have enough continuous derivatives to make the Taylor series expansion valid. Often this is stated as $f$ being "smooth enough." This is usually not a problem, but see [this exercise](problem-newtonalternate).
2. We required $f'(r)\neq 0$—that is, $r$ must be a *simple* root. See [this exercise](problem-newtonmultiple) to investigate what happens at a multiple root.
3. We assumed that the sequence converged, which is not easy to guarantee in any particular case. In fact,
finding a starting guess from which the Newton iteration converges is
often the most challenging part of a rootfinding problem. We will try to deal with this issue in \secref{quasinewton}.

## Implementation

(function-newton)=

````{proof:function} newton
**Newton's method for a scalar rootfinding problem.**

```{code-block} julia
:lineno-start: 1
"""
newton(f,dfdx,x1)

Use Newton's method to find a root of `f` starting from `x1`, where
`dfdx` is the derivative of `f`. Returns a vector of root estimates.
"""
function newton(f,dfdx,x1)
    # Operating parameters.
    funtol = 100*eps();  xtol = 100*eps();  maxiter = 40;

    x = [x1]
    y = f(x1)
    dx = Inf   # for initial pass below
    k = 1

    while (abs(dx) > xtol) && (abs(y) > funtol) && (k < maxiter)
        dydx = dfdx(x[k])
        dx = -y/dydx            # Newton step
        push!(x,x[k]+dx)        # append new estimate

        k = k+1
        y = f(x[k])
    end

    if k==maxiter
        @warn "Maximum number of iterations reached."
    end

    return x
end
```
````

````{prf:example} Julia demo
:class: demo
{doc}`demos/newton-usage`
````

Our implementation of Newton's iteration is given in {ref}`function-newton`. It accepts mathematical functions $f$ and $f'$ and the starting guess $x_1$ as input arguments. Beginning programmers are tempted to embed the functions directly into the code, but there are two good reasons not to do so. First, you would need a new copy of the whole code for each new instance of the problem, even though very little code may need to change. Second, you may want to try more than one rootfinding implementation for a particular problem, and keeping the definition of the problem separate from the algorithm for its solution makes this task much easier. As a practical issue in MATLAB, you can define short functions inline as in {doc}`demos/newton-converge`. For longer functions, you can write separate function files and pass them as arguments to {ref}`function-newton` by affixing an at-sign \verb+@+ to the front of the name of the function file.

```{index} backward error
```

```{index} residual
```

The function {ref}`function-newton` also deals with a thorny practical issue: how to stop the iteration. It adopts a three-part criterion. First, it monitors the difference between successive root estimates, $|x_k-x_{k-1}|$, which is used as a stand-in for the unknown error $|r-x_k|$. In addition, it monitors the residual $|f(x_k)|$, which is equivalent to the backward error and more realistic to control in badly conditioned problems (see {doc}`rootproblem`). If either of these quantities is considered to be sufficiently small, the iteration ends. Finally, we need to protect against the possibility of a nonconvergent iteration, so the procedure terminates with a warning if a maximum number of iterations is exceeded.[^inputs]

[^inputs]: In more practical codes, the thresholds used to make these decisions are controllable through additional user inputs to the procedure.

## Exercises

For each of problems 1--3, do the following steps.
  
**(a)** ✍ Rewrite the equation into the standard form for rootfinding, $f(x) = 0$, and compute $f'(x)$. **(b)** ⌨  Make a plot of $f$ over the given interval and determine how many roots lie in the interval. **(c)** ⌨ Use `nlsolve` to find an "exact" value for each root. **(d)** ⌨ Use {ref}`function-newton` to find each root. **(e)** ⌨ For one of the roots, define `e` as a vector of the errors in the Newton sequence. Determine numerically whether the convergence is roughly quadratic.

1. $x^2=e^{-x}$, over $[-2,2]$

    ````{only} solutions

    ``` matlab
    %% Part (a)
    %
    % To write the function in standard form, put everything on one side of the
    % equation to define
    %
    % $$ f(x) = -x^2+e^{-x}. $$
    %
    % An anonymous function definition follows, as well as a plot.

    f = @(x) -x.^2+exp(-x);
    xx = linspace(-2,2,101);
    plot(xx,f(xx),'b-',[-2,2],[0,0],'g--','LineWidth',2)
    xlabel('x'), ylabel('y')

    %% Part (b)
    %
    % Use fzero to find the roots in $$ x \in [-2,2]$.

    [r,fatr,iflag] = fzero(f,0.13)

    %% Part (c)
    %
    % Solve for the root

    fp = @(x) -2*x-exp(-x);
    rn = newton(f,fp,0.25)

    %% Part (d)
    %
    % Is the convergence quadratic?
    err = abs(rn-r)
    disp('If the exponents double, the convergence is quadratic')

    disp(' ')
    disp('Use ratios; tending to a constant shows quadratic convergence.')
    err(2:end)./(err(1:end-1).^2)

    % Make a figure with a model sequence for comparison
    figure
    a = 2.^(-2.^(1:5));
    loglog(err(1:end-1),err(2:end),'-o',a(1:end-1),a(2:end),'g--','LineWidth',2)
    xlabel('e_n'); ylabel('e_{n+1}');
    title('Slope 2 is quadratic!','FontSize',14);
    legend('Error','Model Seq','Location','NorthWest');

    disp(' ')
    disp('The convergence appears to be quadratic until the finite precision')
    disp('interferes.')
    ```
    ````

2. $2x = \tan x$, over $[-0.2,1.4]$

3. $e^{x+1}=2+x$, over $[-2,2]$

    ````{only} solutions

    ``` matlab
    %% Problem 4.1.3
    %

    %% Setup
    clear all, close all, clc, format compact, format short e

    %% Part (a)
    %
    % To write the function in standard form, put everything on one side of the
    % equation to define
    %
    % $$ f(x) = x + 2 - e^{x+1}. $$
    %
    % An anonymous function definition follows, as well as a plot.

    f = @(x) x+2-exp(x+1);
    xx = linspace(-2,2,101);
    plot(xx,f(xx),'b-','LineWidth',2)
    xlabel('x'), ylabel('y')

    %% Part (b)
    %
    % Use fzero to find the roots in $$ x \in [-2,2]$.

    [r,fatr,iflag] = fzero(f,-0.13)

    % This approach fails because it is a double root at x=-1, and fzero can't
    % solve it because it doesn't change sign across the root.

    %% Part (c)
    %
    % The condition number may be computed as follows.  We want to know about
    % the function on the interval $$ x \in [-2,2] $$.  For the infinity norm
    % of f, use the following.  Note that the root is r=-1 and f'(r)=0;
    % will cause an infinite condition number.

    f_inf_norm = norm(f(xx),inf);
    fp = @(x) 1-exp(x+1);
    r = -1;
    cond = f_inf_norm/abs(r)/abs(fp(r))
    ```

4. ⌨  Consider the equation $f(x)=x^{-2} - \sin x=0$ on the interval $x \in [0.1,4\pi]$.  Use a plot to approximately locate the roots of $f$. To which roots do the following initial guesses converge when using {ref}`function-newton`?  Is the root obtained the one that is closest to that guess?

    **(a)** $x_0 = 1.5,\quad$
    **(b)** $x_0 = 2,\quad$
    **(c)** $x_0 = 3.2,\quad$
    **(d)** $x_0 = 4,\quad$
    **(e)** $x_0 = 5,\quad$
    **(f)** $x_0 = 2\pi$.

    ````{only} solutions

    ``` matlab
    fplot(@(x) 1./x.^2,[.1 4*pi])
    hold on
    fplot(@sin,[.1 4*pi])
    ylim([-3 3])

    %%
    % There are four roots in the interval. Though the function is close
    % to another at the right endpoint, there is no root nearby.

    f = @(x) 1./x.^2 - sin(x);
    dfdx = @(x) -2./x.^3 - cos(x);
    all_roots = [];
    init_ = [1.5 2 3.2 4 5 2*pi];
    for init = init_
        x = newton(f,dfdx,init);
        all_roots = [all_roots; x(end)];
    end

    %%
    format long
    all_roots

    %%
    % Find the four unique ones.
    r = sort(all_roots);
    r(diff(r) < 1e-11) = []

    %%
    % Distance from each initial guess to each root
    fprintf('   Init ')
    fprintf(' %9.3f',r)
    fprintf('     Found\n-----------------------------------------------------------\n')
    for k = 1:6
        init = init_(k);
        fprintf('   %.2f ',init)
        fprintf('%9.1f ',abs(init-r))
        fprintf('  %7.2f\n',all_roots(k))
    end

    %%
    % By inspection, it follows that cases (b)
    % and (e) found a root that was not closest to the initial guess.
    ```

5. ✍ Show that if $f(x)=x^{-1}-b$ for nonzero $b$, then Newton's iteration converging to the root $r=1/b$ can be implemented without performing any divisions. 

    ````{only} solutions

    ``` matlab
    %% (a)
    % $f(x) = bx - 1$ will not work! Instead use
    %
    % $$ f(x) = b - (1/x).$$

    %% (b)
    % $$x_{n+1} = x_n - \displaystyle\frac{b-1/x_n}{1/x_n^2} = x_n - \left( bx_n^2 -x_n \right) = x_n(2-bx_n).$$
    ```

    (problem-newtonalternate)=
6. ✍ Discuss what happens when Newton's method is applied to find a root of $f(x) = \operatorname{sign}(x) \sqrt{|x|}$, starting at $x_0\ne 0$.

    (problem-newtonmultiple)=
7. ✍ In the case of a multiple root, where $f(r)=f'(r)=0$, the derivation of the quadratic error convergence in {eq}`newtonerr` is invalid. Redo the derivation to show that in this circumstance and with $f''(r)\neq 0$, the error converges only linearly.

8. ✍ In {ref}`function-newton` and elsewhere, the actual error is not available, so we use $|x_k-x_{k-1}|$ as an approximate indicator of error to determine when to stop the iteration. Find an example that foils this indicator; that is, a sequence $\{x_k\}$ such that
  
    ```{math}
    \lim_{k\rightarrow \infty} (x_k-x_{k-1}) = 0,
    ```

    but $\{x_k\}$ diverges. (Hint: You have seen such sequences in calculus.) Hence the need for residual tolerances and escape hatches in the code!
