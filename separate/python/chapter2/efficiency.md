---
numbering:
  enumerator: 2.5.%s
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

(section-linsys-efficiency)=

# Efficiency of matrix computations

Predicting how long an algorithm will take to solve a particular problem, on a particular computer, as written in a particular way in a particular programming language, is an enormously difficult undertaking. It's more practical to predict how the required time will scale as a function of the size of the problem. In the case of a linear system of equations, the problem size is $n$, the number of equations and variables in the system.  Because expressions of computational time are necessarily approximate, it's customary to suppress all but the term that is dominant as $n\to\infty$. We first need to build some terminology for these expressions.

## Asymptotic analysis

```{index} ! asymptotic notation
```

::::{prf:definition} Asymptotic notation
:label: definition-asymptotic-notation
Let $f(n)$ and $g(n)$ be positive-valued functions. We say $f(n)=O(g(n))$ (read "$f$ is {term}`big-O` of $g$") as $n\rightarrow \infty$ if $f(n)/g(n)$ is bounded above as $n\to\infty$.

We say $f(n)\sim g(n)$ (read "$f$ is {term}`asymptotic` to $g$") as $n\rightarrow \infty$ if $f(n)/g(n)\rightarrow 1$ as $n\rightarrow\infty$.
::::

One immediate consequence is that $f\sim g$ implies $f=O(g)$.[^sets]

[^sets]: More precisely, $O(g)$ and $\sim g$ are *sets* of functions, and $\sim g$ is a subset of $O(g)$. That we write $f=O(g)$ rather than $f\in O(g)$ is a quirk of convention.

````{prf:example}
Consider the functions $f(n) = a_1 n^3 + b_1 n^2 + c_1 n$ and $g(n) = a_2 n^3$ in the limit $n\to \infty$.  Then
  
```{math}
    \lim_{n \to \infty} \frac{f(n)}{g(n)}
    = \lim_{n \to \infty} \frac{a_1 + b_1n^{-1} + c_1n^{-2}}{a_2} =
    \frac{a_1}{a_2} .
```

Since $a_1/a_2$ is a constant, $f(n) = O(g(n))$; if $a_1=a_2$, then $f \sim g$.
````

````{prf:example}
  Consider $f(n) = \sin (1/n)$, $g(n)=1/n$ and $h(n) = 1/n^2$. For large $n$, Taylor's theorem with remainder implies that
  
```{math}
f(n) = \frac{1}{n} - \cos(1/\xi)\frac{1}{6 n^3},
```

where $n<\xi<\infty$.  But
  
```{math}
\lim_{n\to \infty} \frac{f}{g} = \lim_{n\to \infty} 1-\cos(1/\xi)\frac{1}{6 n^2} = 1,
```

and so $f \sim g$.  On the other hand, comparing $f$ and $h$, we find

```{math}
\lim_{n\to \infty} \frac{f}{h} = \lim_{n\to \infty}  n-\cos(1/\xi)\frac{1}{6 n} = \infty,
```

so we cannot say that $f = O(h)$. A consideration of $h/f$ will show that $h = O(f)$, however.
````

It's conventional to use asymptotic notation that is as specific as possible. For instance, while it is true that $n^2+n=O(n^{10})$, it's more informative, and usually expected, to say $n^2+n=O(n^2)$. There are additional notations that enforce this requirement strictly, but we will just stick to the informal understanding.

There is a memorable way to use asymptotic notation to simplify sums:

```{math}
:label: sumflops
\begin{split}
  \sum_{k=1}^n k&\sim \frac{n^2}{2} = O(n^2), \text{ as $n\to\infty$}, \\
  \sum_{k=1}^n k^2 &\sim \frac{n^3}{3} = O(n^3), \text{ as $n\to\infty$}, \\
  &\vdots \\
  \sum_{k=1}^n k^p &\sim \frac{n^{p+1}}{p+1} = O(n^{p+1}), \text{ as $n\to\infty$}.
\end{split}
```

These formulas greatly resemble the definite integral of $x^p$.

::::{prf:example}
:label: example-efficiency-sums

```{math}
:numbered: false
\begin{align*}
\sum_{k=1}^{n-1} 4k^2 + 3 & = 4 \left( \sum_{k=1}^{n-1} k^2\right)  + 3 \sum_{k=1}^{n-1} 1\\
&\sim 4 \left( \frac{1}{3} (n-1)^3 \right) + 3(n-1) \\
& = \frac{4}{3} (n^3 - 3n^2 + 3n - 1)  + 3n - 3 \\
&\sim \frac{4}{3} n^3.
\end{align*}
```

::::

## Flop counting

```{index} ! flops
```

Traditionally, in numerical linear algebra, we count arithmetic operations to measure an algorithm's running time.

::::{prf:definition} Floating-point operation
:label: definition-flop
A **floating-point operation** or {term}`flop` is a single addition, subtraction, multiplication, division, or square root of two floating-point numbers.
::::

:::{important}
Flop counts are, to put it mildly, a very rough measure of the time a computer program takes to run. Many other factors, such as memory access and copying, have significant effects as well. Still, flop counts can be a useful way to compare different algorithms and to estimate how the runtime requirements scale with the problem size.
:::

::::{prf:example} Floating-point operations in matrix-vector multiplication
:label: demo-flops-mvmult


Here is a straightforward implementation of matrix-vector multiplication.

```{code-cell} 
n = 6
A = random.rand(n, n)
x = ones(n)
y = zeros(n)
for i in range(n):
    for j in range(n):
        y[i] += A[i, j] * x[j]   # 2 flops
```

Each of the loops implies a summation of flops. The total flop count for this algorithm is

$$
\sum_{i=1}^n \sum_{j=1}^n 2 = \sum_{i=1}^n 2n = 2n^2.
$$

Since the matrix $\mathbf{A}$ has $n^2$ elements, all of which have to be involved in the product, it seems unlikely that we could get a flop count that is smaller than $O(n^2)$ in general.

Let's run an experiment with the built-in matrix-vector multiplication. We assume that flops dominate the computation time and thus measure elapsed time. 

```{code-cell} 
N = 400 * arange(1, 11)
t = []
print("  n           t")
for i, n in enumerate(N):
    A = random.randn(n, n)  
    x = random.randn(n)
    start = timer()
    for j in range(50): A @ x
    t.append(timer() - start)
    print(f"{n:5}   {t[-1]:10.3e}")
```

The reason for doing multiple repetitions at each value of $n$ above is to avoid having times so short that the resolution of the timer is a factor.

Looking at the timings just for $n=2000$ and $n=4000$, they have ratio:

```{code-cell} 
print(t[9] / t[4])
```

If the run time is dominated by flops, then we expect this ratio to be 

$$
\frac{2(4000)^2}{2(2000)^2}=4.
$$

::::

Suppose that the running time $t$ of an algorithm obeys a function that is $O(n^p)$. For sufficiently large $n$, $t\approx Cn^p$ for a constant $C$ should be a good approximation. Hence

```{math}
:label: loglogfit
t \approx Cn^p \qquad \Longleftrightarrow \qquad \log t \approx p(\log n) + \log C.
```

So, we expect that a graph of $\log t$ as a function of $\log n$ will be a straight line of slope $p$.

::::{prf:example} Asymptotics in log-log plots
:label: demo-flops-loglog

Let's repeat the experiment of the previous example for more, and larger, values of $n$.

```{code-cell} 
N = arange(400, 6200, 200)
t = zeros(len(N))
for i, n in enumerate(N):
    A = random.randn(n,n)  
    x = random.randn(n)
    start = timer()
    for j in range(20): A@x
    t[i] = timer() - start
```

Plotting the time as a function of $n$ on log-log scales is equivalent to plotting the logs of the variables, but is formatted more neatly. 

```{code-cell} 
fig, ax = subplots()
ax.loglog(N, t, "-o", label="observed")
ylabel("elapsed time (sec)");
xlabel("$n$");
title("Timing of matrix-vector multiplications");
```

You can see that while the full story is complicated, the graph is trending to a straight line of positive slope. For comparison, we can plot a line that represents $O(n^2)$ growth exactly. (All such lines have slope equal to 2.)

```{code-cell} 
ax.loglog(N, t[-1] * (N/N[-1])**2, "--", label="$O(n^2)$")
ax.legend();  fig
```

::::

## Solution of linear systems

Recall the steps of {numref}`Algorithm {number} <algorithm-lu-solve>` for the system $\mathbf{A}\mathbf{x}=\mathbf{b}$:

1. Factor $\mathbf{L}\mathbf{U}=\mathbf{A}$ using Gaussian elimination.
2. Solve $\mathbf{L}\mathbf{z}=\mathbf{b}$ for $\mathbf{z}$ using forward substitution.
3. Solve $\mathbf{U}\mathbf{x}=\mathbf{z}$ for $\mathbf{x}$ using backward substitution.

The second and third steps are solved by {numref}`Function {number} <function-forwardsub>` and {numref}`Function {number} <function-backsub>`.

Take `forwardsub`, for instance. Line 11 computes
```python
s = L[i, :i] @ x[:i]
```

This is the inner product between the first $i$ elements of the $i$th row of `L` and the first $i$ elements of `x`, requiring $i$ multiplications and $(i-1)$ additions. Line 212 adds two more flops, for a total of $2i+1$ in each pass through the loop. In this loop, $i$ ranges from $0$ to $n-1,$ so the total count is

```{math}
:numbered: false
\sum_{i=0}^{n-1} (2i+1) = \sum_{i=1}^n (2i-1) = -n + 2 \sum_{i=1}^n i.
```

It is not hard to find an exact formula for the sum above, but we use the asymptotic expression {eq}`sumflops` to simplify it to $\sim n^2$. After all, since flop counting is only an approximation of true running time, why bother with the more complicated exact expression? An analysis of backward substitution yields the same result.

```{prf:lemma}
Solving a triangular $n\times n$ system by forward or backward substitution takes $\sim n^2$ flops asymptotically.
```

Before counting flops for the LU factorization, we have to admit that {numref}`Function {number} <function-lufact>` is not written as economically as it could be. Recall from our motivating example in @demo-lu-derive that we zero out the first row and column of $\mathbf{A}$ with the first outer product, the second row and column with the second outer product, and so on. There is no good reason to do multiplications and additions with values known to be zero.

Suppose we replace lines 14–17 of `lufact` in @function-lufact with
```{code-block} python
:lineno-start: 14
for k in range(n-1):
    U[k, k:n] = A_k[k, k:n]
    L[k:n, k] = A_k[k:n, k] / U[k, k]
    A_k[k:n, k:n] -= np.outer(L[k:n, k], U[k, k:n])
```

We will use the following handy fact.
:::{prf:observation}
The Numpy slice `k:n`, where $k\le n$, has $n-k$ elements.
:::

Line 16 above divides each element of the vector `A_k[k:n, k] / U[k, k]` by a scalar. Hence, the number of flops equals the length of the vector, which is $n-k$.

Line 17 has an outer product followed by a matrix subtraction. The definition @definition-outerprod of the outer product makes it clear that that computation takes one flop (multiplication) per element of the result, which here results in $(n-k)^2$ flops. The number of subtractions is identical.

Altogether, the factorization takes

```{math}
:numbered: false
\sum_{k=0}^{n-2} n-k + 2(n-k)^2 = \sum_{k=1}^{n-1} n-k + 1 + 2(n-k+1)^2
```

flops.


There are different ways to simplify the total count,

```{math}
:label: gecount1
\sum_{k=1}^{n-1} n-k + 1 + 2(n-k+1)^2.
```

We will make a change of summation index using $j=n-k$. The endpoints of the sum are $j=n-1$ when $k=1$ and $j=1$ when $k=n-1$. Since the order of terms in a sum doesn't matter, we get

```{math}
\begin{align*}
\sum_{j=1}^{n-1} 1+j+2(j+1)^2 &=  \sum_{j=1}^{n-1} 3 + 5j + 2j^2 \\
  & \sim  3(n-1) + \frac{5}{2}(n-1)^2 + \frac{2}{3}(n-1)^3 \\
  & \sim \frac{2}{3}n^3.
\end{align*}
```

We have proved the following.

```{prf:theorem} Efficiency of LU factorization
The LU factorization of an $n\times n$ matrix takes $\sim\frac{2}{3}n^3$ flops as $n\to \infty$. This dominates the flops for solving an $n\times n$ linear system.
```

::::{prf:example} Floating-point operations in LU factorization
:label: demo-flops-lufact


We'll test the conclusion of $O(n^3)$ flops experimentally using the `lu` function imported from `scipi.linalg`.

```{code-cell}
from scipy.linalg import lu
N = arange(200, 2600, 200)
t = zeros(len(N))
for i, n in enumerate(N):
    A = random.randn(n,n)  
    start = timer()
    for j in range(5): lu(A)
    t[i] = timer() - start
```

We plot the timings on a log-log graph and compare it to $O(n^3)$. The result could vary significantly from machine to machine, but in theory the data should start to parallel the line as $n\to\infty$.

```{code-cell}
loglog(N, t, "-o", label="obseved")
loglog(N, t[-1] * (N / N[-1])**3, "--", label="$O(n^3)$")
legend();
xlabel("$n$");
ylabel("elapsed time (sec)");
title("Timing of LU factorizations");
```

::::

In practice, flops are not the only aspect of an implementation that occupies significant time. Our position is that counting flops as a measure of performance is a useful oversimplification. We will assume that LU factorization (and as a result, the solution of a linear system of $n$ equations) requires a real-world time that is roughly $O(n^3)$. This growth rate is a great deal more tolerable than, say, $O(2^n)$, but it does mean that for (at this writing) $n$ greater than 10,000 or so, something other than general LU factorization will have to be used.

## Exercises

``````{exercise}
:label: problem-efficiency-asmptotic1
✍ The following are asymptotic assertions about the limit $n\rightarrow\infty$. In each case, prove the statement true or false.

**(a)** $n^2 = O(\log n),\quad$ 
**(b)** $n^{a} = O(n^b)$ if $a\le b,\quad$
**(c)** $e^n \sim e^{2n},\quad$
**(d)** $n+\sqrt{n}\sim n+2\sqrt{n}$.
``````

``````{exercise}
:label: problem-efficiency-asmptotic2
✍ The following are asymptotic assertions about the limit $h\to 0$. In each case, prove the statement true or false.

**(a)** $h^2\log(h) = O(h^3),\quad$
**(b)** $h^{a} = O(h^b)$ if $a < b,\quad$
**(c)** $\sin(h) \sim h,\quad$
**(d)** $(e^{2h}-1)\sim h$.
``````

``````{exercise}
:label: problem-efficiency-innerproduct
✍ Show that the inner product of two $n$-vectors takes exactly $2n-1$ flops.
``````

``````{exercise}
:label: problem-efficiency-matrixmult
✍ Show that the multiplication of two $n\times n$ matrices takes $\sim 2n^3$ flops.
``````

``````{exercise}
:label: problem-efficiency-polynomial
✍ This problem is about evaluation of a polynomial $c_1 + c_2 x + \cdots + c_{n}x^{n-1}$.

**(a)** Here is a little code to do the evaluation.


``` python
y = c[0]
xpow = 1
for i in range(1, n):
    xpow *= x
    y += c[i]*xpow
```


Assuming that `x` is a scalar, how many flops does this function take, as a function of $n$?

**(b)** Compare the count from (a) to the flop count for Horner's algorithm, {numref}`Function {number} <function-horner>`.
``````

``````{exercise}
:label: problem-efficiency-exact

The exact sums for $p=1,2$ in @sumflops are as follows:

```{math}
:numbered: false
\sum_{k=1}^{n} k = \frac{n(n+1)}{2}, \qquad 
\sum_{k=1}^{n} k^2 = \frac{n(n+1)(2n+1)}{6}.
```

**(a)** ✍  Use these to find the exact result for @gecount1.

**(b)** ⌨ Plot the ratio of your result from (a) and the asymptotic result $2n^3/3$ for all $n=10^{1+0.03i}$, $i=0,\dots,100$, using a log scale for $n$ and a linear scale for the ratio. (The curve should approach 1 asymptotically.)
``````

``````{exercise}
:label: problem-efficiency-power
✍ Show that for any nonnegative constant integer $m$,

```{math}
:numbered: false
\sum_{k=0}^{n-m} k^p \sim \frac{n^{p+1}}{p+1}.
```
``````

``````{exercise}
:label: problem-efficiency-triangular
⌨ *(Julia or MATLAB only)* Because triangular systems are so much faster to solve than general systems, it's worthwhile for the backslash operator in Julia and MATLAB to detect triangular structure and use forward or backward substitution when possible. 

Define a random $1000\times 1000$ matrix $\mathbf{A}$, and a length-1000 random column vector $\mathbf{b}$. Then define `T = tril(A)` to get a lower-triangular matrix. Compare how long it takes to solve the systems $\mathbf{A}\mathbf{x}=\mathbf{b}$ and $\mathbf{T}\mathbf{x}=\mathbf{b}$ using the backslash operator. (You should solve 100 or more systems each time to get reliable timings.)
``````
