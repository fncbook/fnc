# Fixed point iteration

```{index} fixed point problem
```

````{prf:example} Julia demo
:class: demo
{doc}`demos/fp-spiral`
````

The rootfinding problem $f(x)=0$ can always be transformed into another form, $g(x)=x$, known as the {term}`fixed point problem`. Given $f$, one such transformation is to define $g(x)=x-f(x)$. Then the fixed point equation is true at, and only at, a root of $f$. {doc}`demos/fp-spiral` shows that evaluations of the function $g$ can be used to try to locate a fixed point.

This is our first example of an **iterative algortihm**. The idea is to generate not a single answer but a sequence of values that one hopes will converge to the correct result. Often the iteration is constructed by defining a formula to map one member of the sequence to the next one. In this case we have

```{math}
  :label: fixedpointiter
  x_{k+1} = g(x_k), \qquad k=1,2,\ldots,
```

```{index} fixed point iteration
```

which is known as the {term}`fixed point iteration`. In order to fully define the process, we must also provide a starting value $x_1$.  Then {eq}`fixedpointiter` defines the rest of the sequence $x_2,x_3,\ldots$ for as many terms as we care to generate.

## Series analysis

In {doc}`demos/fp-spiral`, the two computed iterations differ only in the choice of $x_1$. In the first case we evidently generated a sequence that converged to one of the fixed points. In the second case, however, the generated sequence diverged.\footnote{We can only ever generate a finite sample from an infinite sequence, which in principle does not guarantee anything whatsoever about the limit or divergence of that sequence. However, in practical computing one usually assumes that well-established trends in the sequence will continue, and we complement the observed experiences with rigorous theory where possible.} The easiest way to uncover the essential difference between the two cases is to use Taylor series expansions.

Suppose a fixed point $r$ is the desired limit of an iteration $x_1,x_2,\ldots$. It's often easier to express quantities in terms of the error sequence $\epsilon_1,\epsilon_2,\ldots,$ where $\epsilon_k=x_k-r$. Starting from {eq}`fixedpointiter`, we have

```{math}
\begin{split}
  \epsilon_{k+1}+r = g( \epsilon_{k}+r ) = g(r) + g'(r) \epsilon_k + \frac{1}{2}g''(r) \epsilon_k^2 + \cdots,
\end{split}
```

assuming that $g$ has at least two continuous derivatives. But by definition, $g(r)=r$, so

```{math}
  :label: fpconverge
  \epsilon_{k+1} = g'(r) \epsilon_k + O(\epsilon_k^2).
```

If the iteration is to converge to $r$, the errors must approach zero. In this case we can neglect the second-order term and conclude that $\epsilon_{k+1} \approx g'(r) \epsilon_k$. This is satisfied if $\epsilon_k\approx C \bigl[g'(r)\bigr]^k$ for a constant $C$ and all sufficiently large $k$. 

Hence if $|g'(r)|>1$, we are led to the contradictory conclusion that the errors must *grow*, not vanish. More precisely, if a term in the sequence really does get close to the fixed point $r$, the iteration will start producing new values that are farther away from it. However, if $|g'(r)|<1$ we do see errors that converge to zero.

(example-fprate)=

````{prf:example}
  The role of $g'(r)$ is clear in {doc}`demos/fp-spiral`. We have $g(x) = -x^2+5x-3.5$ and $g'(x)=-2x+5$. For the first fixed point, near $2.71$, we get $g'(r)\approx-0.42$, indicating convergence. For the second fixed point, near 1.29, we get $g'(r)\approx 2.42$, which is consistent with the observed divergence.
````

## Linear convergence

```{index} linear convergence
```

In numerical computation we want to know not just whether an iteration converges but also the *rate* at which convergence occurs, i.e. how quickly the errors approach zero. Other things being equal, faster convergence is preferred to slower convergence, as it usually implies that the computation will take less time to achieve a desired accuracy.

The prediction of the series analysis above is that if the fixed point iteration converges, the errors satisfy

```{math}
  :label: fplinear
  |\epsilon_k| = |x_k-r| \approx C \sigma^k, \quad \sigma = |g'(r)| < 1.
```

Taking logs, we get

```{math}
  \log |\epsilon_k| \approx k(\log \sigma) + (\log C).
```

```{index} convergence rate; linear
```

This is in the form $\log |\epsilon_k| \approx \alpha k + \beta$, which is a linear relationship, and we refer to this situation as {term}`linear convergence`. The formal definition of linear convergence for the sequence $x_1,x_2,\ldots$ to the number $r$ is that there exists a number $0<\sigma<1$ such that

```{math}
  :label: linearconvergence
  \lim_{k\to\infty} \frac{|x_{k+1}-r|}{|x_k-r|} = \sigma.
```

Practically speaking, linear convergence is identified by two different observations about the errors $\epsilon_k$:

1. The errors lie on a straight line on a log-linear graph, as implied by $\log |\epsilon_k| \approx \alpha k + \beta$.
2. The error is reduced by a constant factor in each step, as implied by $|\epsilon_k| \approx C \sigma^k$.

```{prf:example} Julia demo
:class: demo
{doc}`demos/fp-converge`
```

Both statements are approximate and only apply for sufficiently large values of $k$, so a certain amount of judgment has to be applied.

## Contraction maps

The convergence condition $\sigma=|g'(r)|<1$ derived by series expansion is a special case of a more general condition.

```{index} Lipschitz condition
```

````{prf:definition}
A function $g$ is said to satisfy a **Lipschitz condition** with constant $L$ on the interval $S\subset\mathbb{R}$ if
  
```{math}
    :label: lipschitz
    \bigl| g(s)-g(t) \bigr| \le L \bigl| s-t \bigr|, \; \text{for all } s,t\in S.
```
````

```{index} contraction mapping
```

It can be shown that a function satisfying {eq}`lipschitz` is continuous in $S$. If $L<1$ we call $g$ a **contraction mapping**, because distances between points decrease after an application of $g$.

From the Fundamental Theorem of Calculus, which asserts that $g(s)-g(t)=\int_s^t g'(x)\, dx$, it's not difficult to conclude that an upper bound of $|g'(x)|\le L$ for all $x$ results in {eq}`lipschitz`. But the weaker Lipschitz condition is enough to guarantee the success of fixed point iteration.

(theorem-contraction)=

````{proof:theorem} Contraction mapping
Suppose that $g$ satisfies {eq}`lipschitz` with $L<1$ on an interval $S$. Then $S$ contains exactly one fixed point $r$ of $g$. If $x_1,x_2,\ldots$ are generated by the fixed point iteration {eq}`fixedpointiter`, and $x_1,x_2,\ldots$ all lie in $S$, then $|x_k-r|\le L^{k-1} |x_1-r|$ for all $k>1$.
````

````{proof:proof}
(partial proof)  First we show there is at most one fixed point in $S$. Suppose $f(r)=r$ and $f(s)=s$ in $S$. Then by {eq}`lipschitz`, $|r-s|=|g(r)-g(s)|\le L|r-s|$, which for $L<1$ is possible only if $|r-s|=0$, so $r=s$.

Now suppose that for some $r\in S$, $g(r)=r$. By the definition of the fixed point iteration and the Lipschitz condition,
  
```{math}
|x_{k+1} - r | = |g(x_k) - g(r)| \le L |x_k-r|,
```

which shows that $x_k\to r$ as $k\to \infty$. To show that $r$ must exist and complete the proof, one needs to apply the Cauchy theory of convergence of a sequence, which is beyond the scope of this book.
````

There are stronger and more general statements of [the contraction mapping theorem](theorem-contraction). For instance, it's possible to show that all initial $x_1$ that are sufficiently close to the fixed point will lead to convergence of the iteration. Algorithmically the main virtue of the fixed point iteration is that it is incredibly easy to apply. However, as we are about to discover, it's not the fastest option.

## Exercises

(problem-fixedptsimple)=

1. ✍ In each case, show that the given $g(x)$ has a fixed point at the given $r$ and use {eq}`fpconverge` to show that fixed point iteration can converge to it.
  
    **(a)** $g(x) = \frac{1}{2}\bigl(x + \frac{9}{x}\bigr)$, $r=3$.

    **(b)** $g(x) = \pi + \frac{1}{4}\sin(x)$, $r=\pi$.

    **(c)** $g(x) = x+1-\tan(x/4)$, $r=\pi$

2. ⌨ For each case in the preceding problem, apply fixed point iteration in Julia and use a log--linear graph of the error to verify linear convergence. Then use numerical values of the error to determine an approximate value for $\sigma$ in {eq}`fplinear`.

3. ✍  In each case, show that the given $g(x)$ has a fixed point of the given $r$. Then determine analytically whether the fixed point iteration could converge to that point given a close enough starting value.
  
    **(a)** $g(x) = 3+x-x^2$, $r=\sqrt{3}$

    **(b)** $g(x) = \sqrt{1+x}$, $r=(1+\sqrt{5})/2$

    **(c)** $g(x) = -\sqrt{1+x}$, $r=(1-\sqrt{5})/2$

    **(d)** $g(x) = x+1-\tan(\pi x)$, $r=1/4$
  
4. In {doc}`demos/fp-spiral` we used $g(x)=x-f(x)$ to find a fixed point of the polynomial $f(x)=x^2 - 4x + 3.5$.
  
    **(a)** ✍ Why does the iteration "spiral in" to the fixed point? (Refer to the series analysis.)

    **(b)** ✍ Show that if $\hat{g}(x) = (x^2+3.5)/4$, then any fixed point of $g$ is a root of $f$.

    **(c)** ⌨ Use fixed point iteration on $\hat{g}$ to try to find both roots of $f$, and note which case(s), if either, converge.

    **(d)** ✍ Use {eq}`fpconverge` to explain the observed behavior in part~(c).
  
5. ✍ The $m$th root of a positive real number $a$ is a fixed point of the function
  
    ```{math}
    g(x) = \frac{a}{x^{m-1}}.
    ```

    For what positive integer values of $m$ will the fixed point iteration for $g$ converge (for close enough initial guesses)?

6. **(a)** ✍ Show that $r=1/3$ is a fixed point of $g(x) = 2x-3x^2$.

    **(b)** ✍ Find $g'(1/3)$. How does this affect {eq}`fpconverge`?

    **(c)** ⌨ Do an experiment with fixed point iteration on $g$ to converge to $r=1/3$. Is the convergence linear?
  
    ```{only} solutions
    %% (b)
    % g'(1/3)=0, so the linear term in the analysis drops out.
    %% (c)
    x = .5;
    for i = 1:6
        x = 2*x-3*x^2;
        x - 1/3
    end
    %% The convergence is quadratic.
    ```

7. ✍  Consider the iteration
  
    ```{math}
    x_{k+1} = x_k - \frac{f(x_k)}{c}, \qquad k=0,1,\ldots.
    ```

    Suppose $f(r)=0$ and $f'(r)$ exist. Find one or more conditions on $c$ such that the iteration converges to $r$.
