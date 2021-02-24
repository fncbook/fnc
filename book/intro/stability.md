# Stability

If we solve a problem using a computer algorithm and see a large error in the result, we might suspect poor conditioning in the original mathematical problem. But algorithms can also be sources of errors. When error in the result of an algorithm exceeds what conditioning can explain, we call the algorithm {term}`unstable`.

## Case study

In {ref}`sec-conditioning-multidim` we showed that finding the roots of a quadratic polynomial $ax^2 + b x+c$ is poorly conditioned if and only if the roots are close to each other relative to their size. Thus, for the polynomial

```{math}
:label: quadunstable
p(x) = (x-10^6)(x-10^{-6}) = x^2 - (10^6+10^{-6})x + 1,
```

finding roots is a very well conditioned problem. An obvious algorithm for finding those roots is to directly apply the quadratic formula:

```{math}
:label: quadform
x_1 = \frac{-b + \sqrt{b^2-4ac}}{2a}, \qquad
x_2 = \frac{-b - \sqrt{b^2-4ac}}{2a}.
```

```{sidebar} Demo
:class: demo
{doc}`demos/stability-quadbad`
```

```{index} subtractive cancellation
```

{doc}`demos/stability-quadbad` suggests that the venerable quadratic formula is an *unstable* means of computing roots in finite precision. The roots themselves were not sensitive to the data or arithmetic—it's the specific computational path we chose that caused the huge growth in errors. 
<!-- ```{sidebar} demo
:class: demo
{doc}`demos/stability-quadgood`
``` 
-->

```{sidebar} Demo
:class: demo
{doc}`demos/stability-quadgood`
```

We can confirm this conclusion by finding a different path that avoids subtractive cancellation. A little algebra using {eq}`quadform` confirms the additional formula $x_1x_2=c/a$.  So given one root $r$, we compute the other root using $c/ar$, which has only multiplication and division and therefore creates no numerical trouble.

```{margin}
The sensitivity of the true problem is governed only by its condition number, but the sensitivity of an algorithm depends on the condition numbers of all of its steps.
```

Both algorithms for calculating roots are equivalent when using real numbers and exact arithmetic. When perturbations are introduced into each intermediate value, though, their effects may depend dramatically on the specific sequence of steps taken. The sensitivity of the target problem $f(x)$ is governed only by $\kappa_f$, but the sensitivity of an algorithm depends on the condition numbers of all of its steps. In this example, the direct application of the quadratic formula included a poorly conditioned step—subtractive cancellation—that can be avoided.

This situation may sound hopelessly complicated. But the elementary operations we take for granted are quite well conditioned in most circumstances. The glaring exceptions occur when $|f(x)|$ is much smaller than $|x|$, which is not inherently hard to detect; however, not all such situations create sensitivity. A practical characterization of instability is that results are much less accurate than the conditioning of the problem suggests. Typically one should apply an algorithm to test problems whose answers are well known, or for which other programs are known to work well, in order to spot likely instabilities. In the rest of this book we will see some specific ways in which instability is manifested for different types of problems.

## Backward error

In the presence of poor conditioning, even a good algorithm $\tilde{f}$ for a problem $f$ may not have a small error $|\tilde{f}(x)-f(x)|/|f(x)|$. Just the act of rounding the data to floating point may introduce a large change in the result. There is another way to characterize the error that can be a useful alternative measurement, as illustrated in {numref}`fig-backwarderror`.

```{figure} figures/backwarderror.svg
:name: fig-backwarderror
Illustration of backward error.
```

```{margin}
Backward error measures what change to the original data would reproduce the result found by an algorithm.
```

Let $\tilde{y} = \tilde{f}(x)$ be a computed result for the original data $x$. If there is a value $\tilde{x}$ such that

```{math}
:label: backwarderror
f(\tilde{x}) = \tilde{y} = \tilde{f}(x),
```

```{sidebar} Demo
:class: demo
{doc}`demos/stability-roots`
```

then we call $|\tilde{x}-x|/|x|$ the (relative) {term}`backward error` of the result. Instead of asking, "How close is to the true answer is your answer?", backward error asks, "How close to the true question is the question you answered?"

```{margin}
In an ill-conditioned problem, we can only hope for small backward error, not small error.
```

Small backward error is the best we can hope for in a poorly conditioned problem. Without getting into the formal details, know that if an algorithm always produces small backward errors, then it is stable. But the converse is not always true: some stable algorithms may produce a large backward error.

(example-not-backward-stable)=

````{proof:example}
  One stable algorithm that is not backward stable is floating point evaluation for our old standby, $f(x)=x+1$. If $|x|<\epsilon_\text{mach}/2$, then the computed result is $\tilde{f}(x)=1$, since there are no floating point numbers between $1$ and $1+\epsilon_\text{mach}$. Hence the only possible choice for a real number $\tilde{x}$ satisfying {eq}`backwarderror` is $\tilde{x}=0$. But then $|\tilde{x}-x|/|x|=1$, which indicates 100\% backward error!
````

## Exercises

1. The formulas

    ```{math}
    f(x)=\frac{1-\cos x}{\sin x},\qquad g(x) = \frac{2\sin^2(x/2)}{\sin(x)},
    ```

    are mathematically equivalent in exact arithmetic but suggest algorithms that can behave quite differently in floating point.

    **(a)** ✍ Show that the relative condition numbers for these expressions are the same. (You might want to use computer algebra to help; try simplifying their ratio to 1.)

    **(b)** ✍ Estimate the number of digits lost to subtractive cancellation in the numerator of $f$ when $x=10^{-6}$.

    **(c)** ✍ The formula for $g$ uses only multiplication, division, and evaluation of the sine function; therefore, any potential loss of accuracy would be due to evaluating sin. Find

    ```{math}
    \lim_{x\to 0} \kappa_\text{sin}(x).
    ```

    **(d)** ⌨ Evaluate $f$ and $g$ at $x=10^{-6}$ and compute the relative difference.
  
    ```{only} solutions
    **(a)** Both are equivalent to $x\csc(x)$.

    **(b)** Since $\cos(x)\approx 1-x^2/2$ near zero, the conditioning of the subtraction step is very close to $2/x^2=2\times10^{12}$. That is, 12 digits are lost.

    **(c)** $\lim_{x\to0} |x\cot(x)| = 1. (Hence there is no loss of precision.)

    **(d)** The results are `5.000444502912538e-07` and `5.000000000000416e-07`, for a relative difference of about $8.89\times 10^{-5}$.
    ```

2. Let $f(x) = (e^x-1)/x$.
  
    **(a)** ✍ Find the condition number $\kappa_f(x)$. What is the maximum of $\kappa_f(x)$ over $-1\le x \le 1$?
  
    **(b)** ⌨  Use the "obvious" algorithm

    ``` julia
    y = (exp(x)-1) / x
    ```

    to compute $f(x)$ at $x=10^{-2},10^{-3},10^{-4},\ldots,10^{-11}$.  

    **(c)** ⌨ Use the first 8 terms of the Taylor series

    ```{math}
    f(x) = 1 + \frac{1}{2!}x + \frac{1}{3!}x^2 + \frac{1}{4!}x^3 + \cdots
    ```

    to create a second algorithm, and evaluate it at the same set of points.
  
    **(d)** ⌨  Make a table of the relative difference between the two algorithms as a function of $x$. Which algorithm is more accurate, and why?
  
3. ⌨ The function
  
    ```{math}
    x = \cosh(y) = \frac{e^t + e^{-t}}{2}
    ```

    can be inverted to yield a formula for $\operatorname{acosh}(x)$:
  
    ```{math}
    :label: acosh
    y = \log\bigl(x-\sqrt{x^2-1}\bigr).
    ```

    For the steps below, define $y_i=10^{-4i}$ and $x_i=\cosh(y_i)$ for $i=1,\dots,4$. Hence $y_i=\operatorname{acosh}(x_i)$.

    **(a)** Find the relative condition number of evaluating $f(x) = \operatorname{acosh}(x)$. (You can use {eq}`acosh`) or look up a formula for $f'$ in a calculus book.)  Evaluate $\kappa_f$ at all the $x_i$. You will find that the problem is well-conditioned at these inputs.

    **(b)** Use {eq}`acosh` to approximate $f(x_i)$ for all $i$. Compute the relative accuracy of the results. What is the source of the observed inaccuracy?

    **(c)** An alternative fornula is

    ```{math}
    :label: acosh2
    t = -2\log\left(\sqrt{\frac{x+1}{2}} + \sqrt{\frac{x-1}{2}}\right).
    ```

    Apply {eq}`acosh2` to approximate $f(x_i)$ for all $i$, again computing the relative accuracy of the results.

4. ⌨ (Continuation of [a previous problem](problem-samplevar). Adapted from {cite}`highamAccuracyStability2002`.) One problem with the formula {eq}`samplevar` for sample variance is that one computes a sum for $\overline{x}$, then another sum to find $s^2$. Some statistics textbooks quote a one-pass formula
  
    ```{math}
    \begin{split}
    s^2 &= \frac{1}{n-1} \left( u - \tfrac{1}{n}v^2 \right),\\
    u & = \sum_{i=1}^n x_i^2, \\
    v &= \sum_{i=1}^n x_i.
    \end{split}
    ```

    "One-pass" means that both $u$ and $v$ can be computed in a single loop. Try this formula for the two datasets

    ``` julia
    x = [ 1e6, 1+1e6, 2+1e6]
    x = [ 1e9, 1+1e9, 2+1e9]
    ```

    computing the relative difference from the output of your earlier `samplevar` function. Explain the results.
